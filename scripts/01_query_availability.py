#!/usr/bin/env python3
"""
01_query_availability.py
Survey DENV complete genome availability on NCBI and Nextstrain.
Counts only — no FASTA download. Run 02_download_genomes.py to fetch sequences.

Requires:
  module load ncbi-datasets/v2
  conda activate ./env   (pandas, requests, zstandard)

Output (data/raw/):
  ncbi_denv_metadata.tsv          -- per-accession metadata from NCBI summary
  nextstrain_metadata_<sero>.tsv  -- Nextstrain metadata per serotype
  availability_report.tsv         -- counts by source / serotype / host

Usage:
  python scripts/01_query_availability.py [--outdir data/raw]
                                          [--ncbi-only | --nextstrain-only]
"""

import argparse
import json
import logging
import os
import shutil
import subprocess
import sys
import time
from io import BytesIO
from pathlib import Path

import pandas as pd
import requests

# ── Logging ───────────────────────────────────────────────────────────────────
logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(message)s",
    handlers=[logging.StreamHandler(sys.stdout)],
)
log = logging.getLogger(__name__)

# ── Constants ─────────────────────────────────────────────────────────────────
# DENV species taxid — covers all 4 serotypes; serotype assigned post-hoc
DENV_TAXON = "dengue virus"

# Serotype name patterns (from virus.organism_name field)
# Handles both Arabic (type 1) and Roman numerals (type I) — both appear in NCBI
SEROTYPE_MAP = {
    "dengue virus type 1": "DENV1", "dengue virus type i": "DENV1",
    "dengue virus 1": "DENV1",
    "dengue virus type 2": "DENV2", "dengue virus type ii": "DENV2",
    "dengue virus 2": "DENV2",
    "dengue virus type 3": "DENV3", "dengue virus type iii": "DENV3",
    "dengue virus 3": "DENV3",
    "dengue virus type 4": "DENV4", "dengue virus type iv": "DENV4",
    "dengue virus 4": "DENV4",
}

# Nextstrain pre-curated metadata — zstd-compressed TSVs
NEXTSTRAIN_BASE = "https://data.nextstrain.org/files/workflows/dengue"
NEXTSTRAIN_META = {
    "DENV1": f"{NEXTSTRAIN_BASE}/metadata_denv1.tsv.zst",
    "DENV2": f"{NEXTSTRAIN_BASE}/metadata_denv2.tsv.zst",
    "DENV3": f"{NEXTSTRAIN_BASE}/metadata_denv3.tsv.zst",
    "DENV4": f"{NEXTSTRAIN_BASE}/metadata_denv4.tsv.zst",
}


# ── Helpers ───────────────────────────────────────────────────────────────────

def check_dependencies():
    if not shutil.which("datasets"):
        log.error(
            "'datasets' not found on PATH. Load the module first:\n"
            "  module load ncbi-datasets/v2"
        )
        sys.exit(1)
    log.info("datasets CLI: %s", shutil.which("datasets"))


def assign_serotype(organism_name: str) -> str:
    """Assign DENV serotype label from organism_name string."""
    name = organism_name.lower()
    for pattern, label in SEROTYPE_MAP.items():
        if pattern in name:
            return label
    # Fallback: check for bare digits
    for digit, label in [("1", "DENV1"), ("2", "DENV2"), ("3", "DENV3"), ("4", "DENV4")]:
        if f"type {digit}" in name or f"denv{digit}" in name:
            return label
    return "unknown"


def assign_host_category(host_name: str) -> str:
    """Categorise host organism name into broad groups."""
    name = host_name.lower()
    if "homo sapiens" in name or "human" in name:
        return "human"
    if "aedes" in name or "culex" in name or "mosquito" in name or "stegomyia" in name:
        return "mosquito"
    if name:
        return "other"
    return "unknown"


def flatten_record(rec: dict) -> dict:
    """Flatten a datasets summary JSON record to a flat dict."""
    host_info = rec.get("host", {})
    isolate   = rec.get("isolate", {})
    location  = rec.get("location", {})
    virus     = rec.get("virus", {})
    organism  = virus.get("organism_name", "")
    host_name = host_info.get("organism_name", "")
    return {
        "accession":       rec.get("accession", ""),
        "serotype":        assign_serotype(organism),
        "host_organism":   host_name,
        "host_category":   assign_host_category(host_name),
        "length":          rec.get("length", ""),
        "completeness":    rec.get("completeness", ""),
        "collection_date": isolate.get("collection_date", ""),
        "isolate_name":    isolate.get("name", ""),
        "isolate_source":  isolate.get("source", ""),
        "country":         location.get("geographic_location", ""),
        "region":          location.get("geographic_region", ""),
        "virus_name":      organism,
        "source_db":       rec.get("source_database", ""),
        "release_date":    rec.get("release_date", ""),
    }


# ── NCBI survey ───────────────────────────────────────────────────────────────

def query_ncbi(outdir: Path) -> pd.DataFrame:
    """
    Single bulk query for all DENV complete genomes.
    Filters by host category AFTER retrieving all records (no --host flag).
    """
    api_key = os.environ.get("NCBI_API_KEY", "")
    cmd = [
        "datasets", "summary", "virus", "genome",
        "taxon", DENV_TAXON,
        # NOTE: --complete-only is NOT used here — "COMPLETE" is submitter-assigned
        # and inconsistent. Nextstrain data shows ~8k DENV2 sequences ≥10kb, but
        # --complete-only returns only ~1.7k. We query all and filter by length downstream.
        "--as-json-lines",
    ]
    if api_key:
        cmd += ["--api-key", api_key]
        log.info("Using NCBI API key (10 req/s limit)")
    else:
        log.warning("NCBI_API_KEY not set — limited to 3 req/s. Set it to speed up.")

    log.info("Running: %s", " ".join(cmd))
    result = subprocess.run(cmd, capture_output=True, text=True)

    if result.returncode != 0:
        log.error("datasets CLI error:\n%s", result.stderr.strip())
        return pd.DataFrame()

    records = []
    for line in result.stdout.strip().splitlines():
        if line:
            try:
                records.append(flatten_record(json.loads(line)))
            except (json.JSONDecodeError, Exception) as e:
                log.warning("Skipping malformed line: %s", e)

    if not records:
        log.warning("No NCBI records returned.")
        return pd.DataFrame()

    meta_df = pd.DataFrame(records)
    out_path = outdir / "ncbi_denv_metadata.tsv"
    meta_df.to_csv(out_path, sep="\t", index=False)
    log.info("NCBI: %d records → %s", len(meta_df), out_path)

    # Summary counts
    counts = (
        meta_df.groupby(["serotype", "host_category"])
        .size()
        .reset_index(name="n_genomes")
    )
    counts.insert(0, "source", "NCBI_GenBank")
    log.info("\n%s", counts.to_string(index=False))
    return counts


# ── Nextstrain survey ─────────────────────────────────────────────────────────

def fetch_nextstrain(serotype: str, url: str, outdir: Path) -> pd.DataFrame | None:
    """
    Streaming zstd decompression — avoids loading the whole file into memory.
    Uses curl --continue-at - for resume on partial downloads.
    """
    try:
        import zstandard as zstd
    except ImportError:
        log.error("zstandard not installed. Run: conda activate ./env")
        return None

    out_path = outdir / f"nextstrain_metadata_{serotype.lower()}.tsv"

    # Use requests with streaming
    log.info("Fetching Nextstrain %s …", serotype)
    for attempt in range(1, 4):
        try:
            with requests.get(url, stream=True, timeout=120) as r:
                r.raise_for_status()
                dctx = zstd.ZstdDecompressor()
                # Stream-decompress directly into pandas
                with dctx.stream_reader(r.raw) as reader:
                    df = pd.read_csv(reader, sep="\t", low_memory=False)
            df.to_csv(out_path, sep="\t", index=False)
            log.info("  %s: %d sequences → %s", serotype, len(df), out_path)
            return df
        except Exception as e:
            log.warning("Attempt %d/3 failed for %s: %s", attempt, serotype, e)
            time.sleep(2 ** attempt)

    log.error("All retries failed for %s", serotype)
    return None


def query_nextstrain(outdir: Path) -> pd.DataFrame:
    """Fetch all Nextstrain DENV metadata and return counts table."""
    counts = []
    for sero, url in NEXTSTRAIN_META.items():
        df = fetch_nextstrain(sero, url, outdir)
        if df is None:
            continue
        # Identify host column (column names vary across Nextstrain builds)
        host_col = next(
            (c for c in ["host", "Host", "host_species", "host_name"] if c in df.columns),
            None
        )
        if host_col:
            for host_val, grp in df.groupby(host_col):
                counts.append({
                    "source":         "Nextstrain",
                    "serotype":       sero,
                    "host_category":  assign_host_category(str(host_val)),
                    "host_organism":  str(host_val),
                    "n_genomes":      len(grp),
                })
        else:
            log.warning("  No host column found in Nextstrain %s (columns: %s)", sero, list(df.columns[:10]))
            counts.append({
                "source":        "Nextstrain",
                "serotype":      sero,
                "host_category": "not_annotated",
                "host_organism": "",
                "n_genomes":     len(df),
            })
    return pd.DataFrame(counts) if counts else pd.DataFrame()


# ── Main ──────────────────────────────────────────────────────────────────────

def main():
    parser = argparse.ArgumentParser(description="DENV genome availability survey (counts only)")
    parser.add_argument("--outdir", default="data/raw", help="Output directory (default: data/raw)")
    parser.add_argument("--ncbi-only",       action="store_true")
    parser.add_argument("--nextstrain-only", action="store_true")
    args = parser.parse_args()

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    check_dependencies()
    log.info("=== DENV availability survey | output → %s ===", outdir)

    dfs = []

    if not args.nextstrain_only:
        ncbi_counts = query_ncbi(outdir)
        if not ncbi_counts.empty:
            dfs.append(ncbi_counts)

    if not args.ncbi_only:
        ns_counts = query_nextstrain(outdir)
        if not ns_counts.empty:
            dfs.append(ns_counts)

    if dfs:
        report = pd.concat(dfs, ignore_index=True)
        report_path = outdir / "availability_report.tsv"
        report.to_csv(report_path, sep="\t", index=False)
        log.info("=== Report → %s ===", report_path)
        print("\n" + report.to_string(index=False))
    else:
        log.warning("No data collected.")


if __name__ == "__main__":
    main()
