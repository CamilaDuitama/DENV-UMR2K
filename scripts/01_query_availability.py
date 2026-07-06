#!/usr/bin/env python3
"""
01_query_availability.py
Survey DENV complete genome availability on NCBI (via ncbi-datasets CLI)
and Nextstrain for all four serotypes, split by host (human / mosquito).

Requires:
  module load ncbi-datasets/v2   (datasets CLI, no Python package needed)
  conda env denv-umr2k            (pandas, requests, zstandard)

Output (data/raw/):
  ncbi_denv_counts.tsv            -- counts by source/serotype/host
  ncbi_denv_metadata.tsv          -- per-accession metadata from NCBI
  nextstrain_metadata_<sero>.tsv  -- Nextstrain metadata per serotype
  availability_report.tsv         -- merged comparison table

Usage:
  python scripts/01_query_availability.py [--ncbi-only] [--nextstrain-only]
"""

import json
import logging
import subprocess
import sys
import time
import argparse
from io import BytesIO
from pathlib import Path

import pandas as pd
import requests

# ── Config ───────────────────────────────────────────────────────────────────
OUTDIR = Path("data/raw")
OUTDIR.mkdir(parents=True, exist_ok=True)

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(message)s",
    handlers=[logging.StreamHandler(sys.stdout)],
)
log = logging.getLogger(__name__)

# NCBI taxon names accepted by datasets CLI
SEROTYPES = {
    "DENV1": "dengue virus type 1",
    "DENV2": "dengue virus type 2",
    "DENV3": "dengue virus type 3",
    "DENV4": "dengue virus type 4",
}

# Host queries: datasets CLI accepts scientific names
# "mosquito" is not a recognized NCBI host term; use species names
HOSTS = {
    "human":            "Homo sapiens",
    "Aedes_aegypti":    "Aedes aegypti",
    "Aedes_albopictus": "Aedes albopictus",
}

# Nextstrain pre-curated metadata (correct URLs — .tsv.zst compressed)
NEXTSTRAIN_BASE = "https://data.nextstrain.org/files/workflows/dengue"
NEXTSTRAIN_META = {
    "DENV1": f"{NEXTSTRAIN_BASE}/metadata_denv1.tsv.zst",
    "DENV2": f"{NEXTSTRAIN_BASE}/metadata_denv2.tsv.zst",
    "DENV3": f"{NEXTSTRAIN_BASE}/metadata_denv3.tsv.zst",
    "DENV4": f"{NEXTSTRAIN_BASE}/metadata_denv4.tsv.zst",
}


# ── NCBI datasets CLI ─────────────────────────────────────────────────────────

def run_datasets(taxon: str, host: str) -> list[dict]:
    """
    Call `datasets summary virus genome taxon <taxon> --host <host>
    --complete-only --as-json-lines` and return parsed records.
    Requires `datasets` on PATH (module load ncbi-datasets/v2).
    """
    cmd = [
        "datasets", "summary", "virus", "genome", "taxon", taxon,
        "--host", host,
        "--complete-only",
        "--as-json-lines",
    ]
    log.info("Running: %s", " ".join(cmd))
    result = subprocess.run(cmd, capture_output=True, text=True)
    if result.returncode != 0:
        log.error("datasets CLI error: %s", result.stderr.strip())
        return []
    records = []
    for line in result.stdout.strip().splitlines():
        if line:
            try:
                records.append(json.loads(line))
            except json.JSONDecodeError as e:
                log.warning("JSON parse error: %s", e)
    return records


def flatten_record(rec: dict, serotype: str, host_label: str) -> dict:
    """Extract flat metadata fields from a datasets JSON record."""
    host_info = rec.get("host", {})
    isolate   = rec.get("isolate", {})
    location  = rec.get("location", {})
    virus     = rec.get("virus", {})
    return {
        "accession":       rec.get("accession", ""),
        "serotype":        serotype,
        "host_query":      host_label,
        "host_organism":   host_info.get("organism_name", ""),
        "length":          rec.get("length", ""),
        "completeness":    rec.get("completeness", ""),
        "collection_date": isolate.get("collection_date", ""),
        "isolate_name":    isolate.get("name", ""),
        "isolate_source":  isolate.get("source", ""),
        "country":         location.get("geographic_location", ""),
        "region":          location.get("geographic_region", ""),
        "virus_name":      virus.get("organism_name", ""),
        "source_db":       rec.get("source_database", ""),
        "release_date":    rec.get("release_date", ""),
    }


def query_ncbi() -> tuple[pd.DataFrame, pd.DataFrame]:
    """Query NCBI for all serotype x host combinations. Returns (counts, metadata)."""
    counts = []
    all_records = []

    for sero_label, taxon in SEROTYPES.items():
        for host_label, host_term in HOSTS.items():
            records = run_datasets(taxon, host_term)
            n = len(records)
            log.info("  %s / %s: %d complete genomes", sero_label, host_label, n)
            counts.append({
                "source":    "NCBI_GenBank",
                "serotype":  sero_label,
                "host":      host_label,
                "n_genomes": n,
            })
            for r in records:
                all_records.append(flatten_record(r, sero_label, host_label))
            time.sleep(0.5)  # polite pause between CLI calls

    counts_df = pd.DataFrame(counts)
    meta_df   = pd.DataFrame(all_records)

    counts_df.to_csv(OUTDIR / "ncbi_denv_counts.tsv", sep="\t", index=False)
    meta_df.to_csv(OUTDIR / "ncbi_denv_metadata.tsv", sep="\t", index=False)
    log.info("NCBI metadata → %s (%d rows)", OUTDIR / "ncbi_denv_metadata.tsv", len(meta_df))

    return counts_df, meta_df


# ── Nextstrain ────────────────────────────────────────────────────────────────

def fetch_nextstrain(serotype: str, url: str) -> pd.DataFrame | None:
    """Download a Nextstrain metadata .tsv.zst file and return as DataFrame."""
    try:
        import zstandard as zstd
    except ImportError:
        log.error("zstandard not installed. Activate conda env: conda activate denv-umr2k")
        return None

    log.info("Fetching Nextstrain %s: %s", serotype, url)
    try:
        r = requests.get(url, timeout=120, stream=True)
        r.raise_for_status()
        raw = b"".join(r.iter_content(chunk_size=65536))
        dctx = zstd.ZstdDecompressor()
        tsv_bytes = dctx.decompress(raw)
        df = pd.read_csv(BytesIO(tsv_bytes), sep="\t", low_memory=False)
        out_path = OUTDIR / f"nextstrain_metadata_{serotype.lower()}.tsv"
        df.to_csv(out_path, sep="\t", index=False)
        log.info("  %s: %d sequences → %s", serotype, len(df), out_path)
        return df
    except Exception as e:
        log.warning("Nextstrain %s failed: %s", serotype, e)
        return None


def query_nextstrain() -> pd.DataFrame:
    """Fetch all Nextstrain DENV metadata and return a counts table."""
    counts = []
    for sero, url in NEXTSTRAIN_META.items():
        df = fetch_nextstrain(sero, url)
        if df is not None:
            host_col = next(
                (c for c in ["host", "Host", "host_species"] if c in df.columns),
                None
            )
            if host_col:
                for host_val, grp in df.groupby(host_col):
                    counts.append({
                        "source":    "Nextstrain",
                        "serotype":  sero,
                        "host":      str(host_val),
                        "n_genomes": len(grp),
                    })
            else:
                counts.append({
                    "source":    "Nextstrain",
                    "serotype":  sero,
                    "host":      "not_annotated",
                    "n_genomes": len(df),
                })
    return pd.DataFrame(counts)


# ── Main ──────────────────────────────────────────────────────────────────────

def main():
    parser = argparse.ArgumentParser(description="DENV genome availability survey")
    parser.add_argument("--ncbi-only",      action="store_true", help="Only query NCBI")
    parser.add_argument("--nextstrain-only", action="store_true", help="Only query Nextstrain")
    args = parser.parse_args()

    log.info("=== DENV sequence availability survey ===")

    dfs = []

    if not args.nextstrain_only:
        ncbi_counts, _ = query_ncbi()
        dfs.append(ncbi_counts)

    if not args.ncbi_only:
        ns_counts = query_nextstrain()
        if not ns_counts.empty:
            dfs.append(ns_counts)

    if dfs:
        report = pd.concat(dfs, ignore_index=True)
        report_path = OUTDIR / "availability_report.tsv"
        report.to_csv(report_path, sep="\t", index=False)
        log.info("=== Report → %s ===", report_path)
        print("\n" + report.to_string(index=False))
    else:
        log.warning("No data collected.")


if __name__ == "__main__":
    main()
