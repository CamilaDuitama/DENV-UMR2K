#!/usr/bin/env python3
"""
02_download_genomes.py
Download DENV full genomes prioritising Nextstrain (curated, clade-annotated),
then supplement with NCBI GenBank sequences not already captured.

Strategy
--------
1. Nextstrain  — download sequences_denv{1-4}.fasta.zst + use metadata TSVs
                  from 01_query_availability.py (already downloaded to data/raw/)
2. NCBI        — download ALL dengue (taxid 12637, no --complete-only),
                  subtract Nextstrain accessions, keep remaining ≥10 kb
3. Merge       — cat Nextstrain + NCBI complement → length filter → dedup

Requires:
  module load ncbi-datasets/v2
  module load SeqKit/2.8.2
  conda activate ./env   (pandas, requests, zstandard)

Output (data/processed/):
  nextstrain_seqs_full.fasta          -- Nextstrain sequences ≥10 kb
  nextstrain_metadata_full.tsv        -- corresponding metadata
  ncbi_complement.fasta               -- NCBI sequences ≥10 kb not in Nextstrain
  ncbi_complement_metadata.tsv        -- corresponding metadata
  denv_final.fasta                    -- merged, deduped full-genome dataset
  denv_final_metadata.tsv             -- merged metadata
"""

import argparse
import io
import json
import logging
import os
import shutil
import subprocess
import sys
import time
import zipfile
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
DENV_TAXID = "12637"
MIN_LEN = 10_000
MAX_LEN = 12_000

# Vaccine / passage / chimeric strains not caught by is_lab_host=True.
# Identified by scanning strain names in Nextstrain metadata (2026-07-10).
# PDK = Primary Dog Kidney passage (attenuated); rDEN = recombinant chimeric.
VACCINE_ACCESSIONS = {
    "MW945952",  # RDENV1-WP-1A     — recombinant DENV1
    "AF180818",  # 16007 (PDK-13)   — DENV1, PDK-13 passage attenuated
    "U87412",    # PDK-53           — DENV2, PDK-53 attenuated vaccine
    "KU725664",  # PDK53            — DENV2, PDK-53 variant
    "M84728",    # 16681-PDK53      — DENV2, PDK-53 variant
    "KJ160505",  # rDENV3-4         — DENV3/4 chimeric recombinant
    "MW793459",  # PDK48            — DENV4, PDK-48 attenuated vaccine
    "KJ160504",  # rDENV4           — recombinant DENV4
}

NEXTSTRAIN_BASE = "https://data.nextstrain.org/files/workflows/dengue"
NEXTSTRAIN_SEQS = {
    "DENV1": f"{NEXTSTRAIN_BASE}/sequences_denv1.fasta.zst",
    "DENV2": f"{NEXTSTRAIN_BASE}/sequences_denv2.fasta.zst",
    "DENV3": f"{NEXTSTRAIN_BASE}/sequences_denv3.fasta.zst",
    "DENV4": f"{NEXTSTRAIN_BASE}/sequences_denv4.fasta.zst",
}

SEROTYPE_MAP = {
    "dengue virus type 1": "DENV1", "dengue virus type i":   "DENV1", "dengue virus 1": "DENV1",
    "dengue virus type 2": "DENV2", "dengue virus type ii":  "DENV2", "dengue virus 2": "DENV2",
    "dengue virus type 3": "DENV3", "dengue virus type iii": "DENV3", "dengue virus 3": "DENV3",
    "dengue virus type 4": "DENV4", "dengue virus type iv":  "DENV4", "dengue virus 4": "DENV4",
}


# ── Helpers ───────────────────────────────────────────────────────────────────

def check_dependencies():
    missing = []
    if not shutil.which("datasets"):
        missing.append("datasets  →  module load ncbi-datasets/v2")
    if not shutil.which("seqkit"):
        missing.append("seqkit    →  module load SeqKit/2.8.2")
    if missing:
        log.error("Missing:\n  " + "\n  ".join(missing))
        sys.exit(1)


def run(cmd: list[str], desc: str) -> subprocess.CompletedProcess:
    log.info("[%s] %s", desc, " ".join(cmd))
    r = subprocess.run(cmd, capture_output=True, text=True)
    if r.returncode != 0:
        log.error("[%s] FAILED:\n%s", desc, r.stderr.strip())
        sys.exit(1)
    return r


def count_seqs(fasta: Path) -> int:
    r = subprocess.run(["seqkit", "stats", "-T", str(fasta)], capture_output=True, text=True)
    lines = r.stdout.strip().splitlines()
    if len(lines) > 1:
        return int(lines[1].split("\t")[3])
    return 0


def assign_serotype(name: str) -> str:
    n = name.lower()
    for pat, lab in SEROTYPE_MAP.items():
        if pat in n:
            return lab
    return "unknown"


def assign_host_category(host: str) -> str:
    h = host.lower()
    if "homo sapiens" in h or "human" in h:
        return "human"
    if "aedes" in h or "culex" in h or "mosquito" in h or "stegomyia" in h:
        return "mosquito"
    return "other" if h else "unknown"


def download_zst_fasta(url: str, out_fasta: Path, label: str) -> bool:
    """Stream-decompress a .fasta.zst URL directly to a FASTA file."""
    try:
        import zstandard as zstd
    except ImportError:
        log.error("zstandard not installed. conda activate ./env")
        return False

    log.info("Downloading %s …", label)
    for attempt in range(1, 4):
        try:
            with requests.get(url, stream=True, timeout=300) as r:
                r.raise_for_status()
                dctx = zstd.ZstdDecompressor()
                with dctx.stream_reader(r.raw) as reader, open(out_fasta, "wb") as fh:
                    shutil.copyfileobj(reader, fh)
            log.info("  %s → %s", label, out_fasta)
            return True
        except Exception as e:
            log.warning("  Attempt %d/3 failed: %s", attempt, e)
            time.sleep(2 ** attempt)
    log.error("All retries failed for %s", label)
    return False


# ── Phase 1: Nextstrain sequences ─────────────────────────────────────────────

def download_nextstrain_seqs(raw_dir: Path, proc_dir: Path) -> tuple[Path, pd.DataFrame]:
    """
    Download Nextstrain FASTA for all serotypes, cat, length-filter, collect metadata.
    Returns (filtered_fasta_path, metadata_df).
    """
    log.info("=== Phase 1: Nextstrain sequences ===")

    # Download per-serotype FASTAs
    per_sero_fastas = []
    for sero, url in NEXTSTRAIN_SEQS.items():
        out = raw_dir / f"nextstrain_seqs_{sero.lower()}.fasta"
        if download_zst_fasta(url, out, f"Nextstrain {sero}"):
            per_sero_fastas.append((sero, out))

    if not per_sero_fastas:
        log.error("No Nextstrain sequences downloaded.")
        sys.exit(1)

    # Concatenate
    all_ns_fasta = raw_dir / "nextstrain_seqs_all.fasta"
    with open(all_ns_fasta, "wb") as out_fh:
        for _, fa in per_sero_fastas:
            with open(fa, "rb") as fh:
                shutil.copyfileobj(fh, out_fh)
    log.info("Concatenated %d serotype FASTAs → %s (%d seqs)",
             len(per_sero_fastas), all_ns_fasta, count_seqs(all_ns_fasta))

    # Length filter
    ns_filtered = proc_dir / "nextstrain_seqs_full.fasta"
    run(["seqkit", "seq", "-m", str(MIN_LEN), "-M", str(MAX_LEN),
         str(all_ns_fasta), "-o", str(ns_filtered)], "seqkit length filter [Nextstrain]")
    n_ns = count_seqs(ns_filtered)
    log.info("Nextstrain full genomes (≥%d bp): %d", MIN_LEN, n_ns)

    # Load metadata TSVs (produced by 01_query_availability.py)
    meta_dfs = []
    for sero in ["DENV1", "DENV2", "DENV3", "DENV4"]:
        meta_path = raw_dir / f"nextstrain_metadata_{sero.lower()}.tsv"
        if meta_path.exists():
            df = pd.read_csv(meta_path, sep="\t", low_memory=False)
            df["serotype"] = sero
            meta_dfs.append(df)
        else:
            log.warning("Nextstrain metadata not found: %s — run 01_query_availability.py first", meta_path)

    if not meta_dfs:
        log.error("No Nextstrain metadata found. Run 01_query_availability.py first.")
        sys.exit(1)

    ns_meta = pd.concat(meta_dfs, ignore_index=True)

    # Remove vaccine / chimeric strains not caught by is_lab_host flag
    vax_mask = ns_meta["accession"].isin(VACCINE_ACCESSIONS)
    if vax_mask.sum():
        log.info("Removing %d vaccine/chimeric accessions from Nextstrain metadata", vax_mask.sum())
        ns_meta = ns_meta[~vax_mask].copy()

    # Keep only accessions retained after length filter
    retained_ids = set(subprocess.run(
        ["seqkit", "seq", "--name", "--only-id", str(ns_filtered)],
        capture_output=True, text=True
    ).stdout.strip().splitlines())

    # Nextstrain uses 'accession' or 'strain' as ID; try both
    id_col = "accession" if "accession" in ns_meta.columns else "strain"
    ns_meta_full = ns_meta[ns_meta[id_col].isin(retained_ids)].copy()

    # Standardise host columns
    host_col = next((c for c in ["host", "Host", "host_species"] if c in ns_meta_full.columns), None)
    if host_col:
        ns_meta_full["host_category"] = ns_meta_full[host_col].fillna("").apply(assign_host_category)
    else:
        ns_meta_full["host_category"] = "unknown"

    out_meta = proc_dir / "nextstrain_metadata_full.tsv"
    ns_meta_full.to_csv(out_meta, sep="\t", index=False)
    log.info("Nextstrain metadata → %s (%d rows)", out_meta, len(ns_meta_full))

    return ns_filtered, ns_meta_full


# ── Phase 2: NCBI complement ──────────────────────────────────────────────────

def download_ncbi_complement(
    raw_dir: Path, proc_dir: Path, ns_accessions: set[str]
) -> tuple[Path, pd.DataFrame]:
    """
    Download all DENV from NCBI (no --complete-only), subtract Nextstrain accessions,
    length-filter, return complement FASTA + metadata.
    """
    log.info("=== Phase 2: NCBI complement ===")

    zip_path = raw_dir / "ncbi_dengue_all.zip"
    unzip_dir = raw_dir / "ncbi_dataset_download"

    api_key = os.environ.get("NCBI_API_KEY", "")
    cmd = [
        "datasets", "download", "virus", "genome",
        "taxon", DENV_TAXID,
        # No --complete-only: length filter applied by seqkit below
        # "info" is not a valid --include value for virus downloads;
        # data_report.jsonl is included automatically with every download
        "--include", "genome",
        "--filename", str(zip_path),
    ]
    if api_key:
        cmd += ["--api-key", api_key]
    else:
        log.warning("NCBI_API_KEY not set — downloads may be throttled")

    run(cmd, "datasets download NCBI")
    log.info("Downloaded → %s (%.1f MB)", zip_path, zip_path.stat().st_size / 1e6)

    if unzip_dir.exists():
        shutil.rmtree(unzip_dir)
    with zipfile.ZipFile(zip_path, "r") as zf:
        zf.extractall(unzip_dir)
    zip_path.unlink()

    # Locate FASTA and metadata
    fna_candidates = list(unzip_dir.rglob("genomic.fna"))
    jsonl_candidates = list(unzip_dir.rglob("data_report.jsonl"))
    if not fna_candidates or not jsonl_candidates:
        log.error("NCBI package missing genomic.fna or data_report.jsonl")
        sys.exit(1)

    raw_fasta  = raw_dir / "ncbi_all_raw.fasta"
    shutil.copy2(fna_candidates[0], raw_fasta)
    log.info("NCBI raw FASTA: %d seqs", count_seqs(raw_fasta))

    # Parse metadata
    records = []
    with open(jsonl_candidates[0]) as fh:
        for line in fh:
            line = line.strip()
            if not line:
                continue
            try:
                rec = json.loads(line)
            except json.JSONDecodeError:
                continue
            host_name = rec.get("host", {}).get("organism_name", "")
            organism  = rec.get("virus", {}).get("organism_name", "")
            isolate   = rec.get("isolate", {})
            location  = rec.get("location", {})
            records.append({
                "accession":       rec.get("accession", ""),
                "serotype":        assign_serotype(organism),
                "host_organism":   host_name,
                "host_category":   assign_host_category(host_name),
                "length":          rec.get("length", ""),
                "completeness":    rec.get("completeness", ""),
                "collection_date": isolate.get("collection_date", ""),
                "isolate_name":    isolate.get("name", ""),
                "country":         location.get("geographic_location", ""),
                "region":          location.get("geographic_region", ""),
                "virus_name":      organism,
                "source":          "NCBI_GenBank",
            })
    ncbi_meta = pd.DataFrame(records)
    log.info("NCBI metadata: %d records", len(ncbi_meta))

    # Subtract Nextstrain accessions
    complement_acc = set(ncbi_meta["accession"]) - ns_accessions
    complement_meta = ncbi_meta[ncbi_meta["accession"].isin(complement_acc)].copy()
    log.info("NCBI complement (not in Nextstrain): %d sequences", len(complement_meta))

    # Write accession allowlist for seqkit grep
    acc_file = proc_dir / "ncbi_complement_accessions.txt"
    complement_meta["accession"].to_csv(acc_file, index=False, header=False)

    # grep + length filter
    ncbi_grep   = proc_dir / "ncbi_complement_grep.fasta"
    ncbi_filtered = proc_dir / "ncbi_complement.fasta"
    run(["seqkit", "grep", "-f", str(acc_file), str(raw_fasta), "-o", str(ncbi_grep)],
        "seqkit grep [NCBI complement]")
    run(["seqkit", "seq", "-m", str(MIN_LEN), "-M", str(MAX_LEN),
         str(ncbi_grep), "-o", str(ncbi_filtered)], "seqkit length filter [NCBI]")
    ncbi_grep.unlink(missing_ok=True)
    acc_file.unlink(missing_ok=True)

    n_ncbi = count_seqs(ncbi_filtered)
    log.info("NCBI complement full genomes (≥%d bp): %d", MIN_LEN, n_ncbi)

    # Keep only retained accessions in metadata
    retained = set(subprocess.run(
        ["seqkit", "seq", "--name", "--only-id", str(ncbi_filtered)],
        capture_output=True, text=True
    ).stdout.strip().splitlines())
    complement_meta = complement_meta[complement_meta["accession"].isin(retained)].copy()

    out_meta = proc_dir / "ncbi_complement_metadata.tsv"
    complement_meta.to_csv(out_meta, sep="\t", index=False)
    log.info("NCBI complement metadata → %s (%d rows)", out_meta, len(complement_meta))

    return ncbi_filtered, complement_meta


# ── Phase 3: Merge + validate ─────────────────────────────────────────────────

def merge_and_validate(
    ns_fasta: Path, ns_meta: pd.DataFrame,
    ncbi_fasta: Path, ncbi_meta: pd.DataFrame,
    proc_dir: Path
):
    log.info("=== Phase 3: Merge and validate ===")

    # Cat FASTAs
    merged_fasta = proc_dir / "denv_final.fasta"
    with open(merged_fasta, "wb") as out:
        for fa in [ns_fasta, ncbi_fasta]:
            if fa.exists() and fa.stat().st_size > 0:
                with open(fa, "rb") as fh:
                    shutil.copyfileobj(fh, out)

    # Final dedup by ID
    dedup_fasta = proc_dir / "denv_final_dedup.fasta"
    run(["seqkit", "rmdup", "-n", str(merged_fasta), "-o", str(dedup_fasta)], "seqkit rmdup")
    merged_fasta.unlink()
    dedup_fasta.rename(merged_fasta)

    n_final = count_seqs(merged_fasta)
    log.info("Final dataset: %d sequences → %s", n_final, merged_fasta)

    # Merge metadata — align columns as best we can
    ns_meta["source"] = "Nextstrain"
    shared_cols = ["accession", "serotype", "host_category", "source"]
    for col in shared_cols:
        if col not in ns_meta.columns:
            ns_meta[col] = ""
        if col not in ncbi_meta.columns:
            ncbi_meta[col] = ""

    merged_meta = pd.concat([ns_meta, ncbi_meta], ignore_index=True)

    # Cross-check: keep only sequences present in final FASTA
    retained = set(subprocess.run(
        ["seqkit", "seq", "--name", "--only-id", str(merged_fasta)],
        capture_output=True, text=True
    ).stdout.strip().splitlines())

    id_col = next((c for c in ["accession", "strain"] if c in merged_meta.columns), None)
    if id_col:
        merged_meta = merged_meta[merged_meta[id_col].isin(retained)]

    meta_out = proc_dir / "denv_final_metadata.tsv"
    merged_meta.to_csv(meta_out, sep="\t", index=False)
    log.info("Final metadata → %s (%d rows)", meta_out, len(merged_meta))

    # Final summary
    summary = merged_meta.groupby(["serotype", "host_category", "source"]).size().reset_index(name="n")
    log.info("\n=== FINAL DATASET SUMMARY ===\n%s", summary.to_string(index=False))
    log.info("=== Done. Final FASTA → %s ===", merged_fasta)


# ── Main ──────────────────────────────────────────────────────────────────────

def main():
    parser = argparse.ArgumentParser(description="Download DENV full genomes (Nextstrain-first)")
    parser.add_argument("--outdir",   default="data", help="Base output dir (default: data)")
    parser.add_argument("--skip-ncbi", action="store_true", help="Only download Nextstrain sequences")
    args = parser.parse_args()

    raw_dir  = Path(args.outdir) / "raw"
    proc_dir = Path(args.outdir) / "processed"
    raw_dir.mkdir(parents=True, exist_ok=True)
    proc_dir.mkdir(parents=True, exist_ok=True)

    check_dependencies()
    log.info("=== DENV genome download | Nextstrain-first strategy ===")
    log.info("Length filter: %d–%d bp | Output: %s", MIN_LEN, MAX_LEN, args.outdir)

    # Phase 1: Nextstrain
    ns_fasta, ns_meta = download_nextstrain_seqs(raw_dir, proc_dir)
    ns_accessions = set(ns_meta.get("accession", ns_meta.get("strain", pd.Series())).dropna())

    if args.skip_ncbi:
        log.info("--skip-ncbi: skipping NCBI download")
        # Still save as final
        ns_fasta.rename(proc_dir / "denv_final.fasta")
        ns_meta.to_csv(proc_dir / "denv_final_metadata.tsv", sep="\t", index=False)
        log.info("=== Done (Nextstrain only) ===")
        return

    # Phase 2: NCBI complement
    ncbi_fasta, ncbi_meta = download_ncbi_complement(raw_dir, proc_dir, ns_accessions)

    # Phase 3: Merge
    merge_and_validate(ns_fasta, ns_meta, ncbi_fasta, ncbi_meta, proc_dir)


if __name__ == "__main__":
    main()
