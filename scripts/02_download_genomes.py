#!/usr/bin/env python3
"""
02_download_genomes.py
Download ALL DENV complete genomes from NCBI, filter by host and length,
deduplicate, and produce a clean FASTA + metadata TSV ready for analysis.

Requires:
  module load ncbi-datasets/v2
  conda activate ./env   (pandas, biopython, zstandard)
  seqkit on PATH         (from conda env)

Output (data/processed/):
  denv_all_raw.fasta              -- all downloaded sequences (unfiltered)
  denv_all_metadata.tsv           -- parsed metadata from data_report.jsonl
  denv_filtered.fasta             -- length + host filtered + deduped
  denv_filtered_metadata.tsv      -- metadata for retained sequences

Usage:
  python scripts/02_download_genomes.py [--outdir data/processed]
                                        [--min-len 10000] [--max-len 12000]
                                        [--hosts human mosquito]
                                        [--keep-zip]
"""

import argparse
import json
import logging
import os
import shutil
import subprocess
import sys
import zipfile
from pathlib import Path

import pandas as pd

# ── Logging ───────────────────────────────────────────────────────────────────
logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(message)s",
    handlers=[logging.StreamHandler(sys.stdout)],
)
log = logging.getLogger(__name__)

# ── Constants ─────────────────────────────────────────────────────────────────
DENV_TAXID   = "12637"   # Dengue virus species (all 4 serotypes)
MIN_LEN_DEFAULT = 10_000
MAX_LEN_DEFAULT = 12_000

SEROTYPE_MAP = {
    "dengue virus type 1": "DENV1", "dengue virus 1": "DENV1",
    "dengue virus type 2": "DENV2", "dengue virus 2": "DENV2",
    "dengue virus type 3": "DENV3", "dengue virus 3": "DENV3",
    "dengue virus type 4": "DENV4", "dengue virus 4": "DENV4",
}


# ── Helpers ───────────────────────────────────────────────────────────────────

def check_dependencies():
    missing = []
    if not shutil.which("datasets"):
        missing.append("datasets  →  module load ncbi-datasets/v2")
    if not shutil.which("seqkit"):
        missing.append("seqkit    →  conda activate ./env")
    if missing:
        log.error("Missing dependencies:\n  " + "\n  ".join(missing))
        sys.exit(1)
    log.info("datasets: %s", shutil.which("datasets"))
    log.info("seqkit:   %s", shutil.which("seqkit"))


def assign_serotype(organism_name: str) -> str:
    name = organism_name.lower()
    for pattern, label in SEROTYPE_MAP.items():
        if pattern in name:
            return label
    for digit, label in [("1","DENV1"),("2","DENV2"),("3","DENV3"),("4","DENV4")]:
        if f"type {digit}" in name or f"denv{digit}" in name:
            return label
    return "unknown"


def assign_host_category(host_name: str) -> str:
    name = host_name.lower()
    if "homo sapiens" in name or "human" in name:
        return "human"
    if "aedes" in name or "culex" in name or "mosquito" in name or "stegomyia" in name:
        return "mosquito"
    return "other" if name else "unknown"


def run(cmd: list[str], desc: str) -> subprocess.CompletedProcess:
    log.info("[%s] %s", desc, " ".join(cmd))
    result = subprocess.run(cmd, capture_output=True, text=True)
    if result.returncode != 0:
        log.error("[%s] FAILED:\n%s", desc, result.stderr.strip())
        sys.exit(1)
    return result


# ── Phase 1: Download ─────────────────────────────────────────────────────────

def download_ncbi(outdir: Path, keep_zip: bool) -> Path:
    """
    Download all complete DENV genomes using datasets download.
    Returns path to unzipped ncbi_dataset directory.
    """
    zip_path = outdir / "ncbi_dengue_all.zip"
    unzip_dir = outdir / "ncbi_dataset_download"

    api_key = os.environ.get("NCBI_API_KEY", "")
    cmd = [
        "datasets", "download", "virus", "genome",
        "taxon", DENV_TAXID,
        "--complete-only",
        "--include", "genome,cds,gff3,info",
        "--filename", str(zip_path),
    ]
    if api_key:
        cmd += ["--api-key", api_key]
        log.info("Using NCBI API key")
    else:
        log.warning("NCBI_API_KEY not set — downloads may be slower")

    log.info("=== Phase 1: Downloading NCBI genomes ===")
    run(cmd, "datasets download")
    log.info("Downloaded → %s (%.1f MB)", zip_path, zip_path.stat().st_size / 1e6)

    # Unzip
    log.info("Unzipping …")
    if unzip_dir.exists():
        shutil.rmtree(unzip_dir)
    with zipfile.ZipFile(zip_path, "r") as zf:
        zf.extractall(unzip_dir)
    log.info("Unzipped → %s", unzip_dir)

    if not keep_zip:
        zip_path.unlink()
        log.info("Removed zip archive")

    return unzip_dir


# ── Phase 2: Parse metadata ───────────────────────────────────────────────────

def parse_metadata(unzip_dir: Path, raw_outdir: Path) -> pd.DataFrame:
    """
    Parse data_report.jsonl from the downloaded package.
    Fields: accession, serotype, host_organism, host_category,
            country, region, collection_date, length, source_db, etc.
    """
    log.info("=== Phase 2: Parsing metadata ===")

    # data_report.jsonl can live in different locations depending on datasets version
    candidates = list(unzip_dir.rglob("data_report.jsonl"))
    if not candidates:
        log.error("data_report.jsonl not found under %s", unzip_dir)
        sys.exit(1)
    report_path = candidates[0]
    log.info("Metadata file: %s", report_path)

    records = []
    with open(report_path) as fh:
        for line in fh:
            line = line.strip()
            if not line:
                continue
            try:
                rec = json.loads(line)
            except json.JSONDecodeError:
                continue

            host_info = rec.get("host", {})
            isolate   = rec.get("isolate", {})
            location  = rec.get("location", {})
            virus     = rec.get("virus", {})
            organism  = virus.get("organism_name", "")
            host_name = host_info.get("organism_name", "")

            records.append({
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
            })

    meta_df = pd.DataFrame(records)
    out_path = raw_outdir / "denv_all_metadata.tsv"
    meta_df.to_csv(out_path, sep="\t", index=False)
    log.info("Parsed %d records → %s", len(meta_df), out_path)

    # Summary by serotype × host
    summary = meta_df.groupby(["serotype", "host_category"]).size().reset_index(name="n")
    log.info("\nRaw counts:\n%s", summary.to_string(index=False))

    return meta_df


# ── Phase 3: Copy + locate FASTA ─────────────────────────────────────────────

def locate_fasta(unzip_dir: Path, raw_outdir: Path) -> Path:
    """Find and copy the merged genomic FASTA."""
    log.info("=== Phase 3: Locating FASTA ===")
    candidates = list(unzip_dir.rglob("genomic.fna"))
    if not candidates:
        log.error("genomic.fna not found under %s", unzip_dir)
        sys.exit(1)
    src = candidates[0]
    dst = raw_outdir / "denv_all_raw.fasta"
    shutil.copy2(src, dst)
    log.info("FASTA → %s", dst)
    return dst


# ── Phase 4: Filter and deduplicate ──────────────────────────────────────────

def filter_and_dedup(
    raw_fasta: Path,
    meta_df: pd.DataFrame,
    processed_dir: Path,
    host_categories: list[str],
    min_len: int,
    max_len: int,
) -> None:
    """
    1. Build accession allowlist from metadata (host filter).
    2. seqkit grep to keep only allowed accessions.
    3. seqkit seq to length-filter.
    4. seqkit rmdup to remove exact duplicates by sequence ID.
    5. Write filtered metadata TSV.
    """
    log.info("=== Phase 4: Filtering and deduplication ===")
    processed_dir.mkdir(parents=True, exist_ok=True)

    # Host filter
    keep_df = meta_df[meta_df["host_category"].isin(host_categories)].copy()
    log.info("Host filter (%s): %d / %d sequences retained",
             "+".join(host_categories), len(keep_df), len(meta_df))

    if keep_df.empty:
        log.warning("No sequences pass host filter — check host_categories argument.")
        return

    # Write accession list for seqkit grep
    acc_list = processed_dir / "keep_accessions.txt"
    keep_df["accession"].to_csv(acc_list, index=False, header=False)

    step1 = processed_dir / "host_filtered.fasta"
    step2 = processed_dir / "length_filtered.fasta"
    step3 = processed_dir / "denv_filtered.fasta"

    # Step 1 — seqkit grep by accession
    run(
        ["seqkit", "grep", "-f", str(acc_list), str(raw_fasta), "-o", str(step1)],
        "seqkit grep"
    )
    log.info("After host filter: %s", _count_seqs(step1))

    # Step 2 — seqkit length filter
    run(
        ["seqkit", "seq", "-m", str(min_len), "-M", str(max_len),
         str(step1), "-o", str(step2)],
        "seqkit length filter"
    )
    log.info("After length filter [%d–%d]: %s", min_len, max_len, _count_seqs(step2))

    # Step 3 — seqkit rmdup (by sequence ID, not sequence content)
    run(
        ["seqkit", "rmdup", "-n", str(step2), "-o", str(step3)],
        "seqkit rmdup"
    )
    log.info("After deduplication: %s", _count_seqs(step3))

    # Clean up intermediates
    step1.unlink(missing_ok=True)
    step2.unlink(missing_ok=True)
    acc_list.unlink(missing_ok=True)

    # ── Phase 5: Validation & filtered metadata ───────────────────────────────
    log.info("=== Phase 5: Validation ===")

    # Extract retained accessions from filtered FASTA
    result = subprocess.run(
        ["seqkit", "seq", "--name", "--only-id", str(step3)],
        capture_output=True, text=True
    )
    retained_ids = set(result.stdout.strip().splitlines())
    log.info("Retained unique accessions: %d", len(retained_ids))

    # Cross-check against metadata
    filtered_meta = keep_df[keep_df["accession"].isin(retained_ids)].copy()
    meta_out = processed_dir / "denv_filtered_metadata.tsv"
    filtered_meta.to_csv(meta_out, sep="\t", index=False)
    log.info("Filtered metadata → %s (%d rows)", meta_out, len(filtered_meta))

    # Check for duplicates
    dup_count = len(retained_ids) - filtered_meta["accession"].nunique()
    if dup_count != 0:
        log.warning("Accession count mismatch — check for duplicate IDs in FASTA")
    else:
        log.info("No duplicate accessions in final dataset")

    # Final summary
    summary = filtered_meta.groupby(["serotype", "host_category"]).size().reset_index(name="n")
    log.info("\nFinal counts:\n%s", summary.to_string(index=False))
    log.info("=== Done. Final FASTA → %s ===", step3)


def _count_seqs(fasta: Path) -> str:
    r = subprocess.run(["seqkit", "stats", "-T", str(fasta)], capture_output=True, text=True)
    lines = r.stdout.strip().splitlines()
    return lines[1] if len(lines) > 1 else "?"


# ── Main ──────────────────────────────────────────────────────────────────────

def main():
    parser = argparse.ArgumentParser(description="Download + filter DENV complete genomes")
    parser.add_argument("--outdir",    default="data",     help="Base output dir (default: data)")
    parser.add_argument("--min-len",   type=int, default=MIN_LEN_DEFAULT)
    parser.add_argument("--max-len",   type=int, default=MAX_LEN_DEFAULT)
    parser.add_argument("--hosts",     nargs="+", default=["human", "mosquito"],
                        choices=["human", "mosquito", "other", "unknown"],
                        help="Host categories to retain (default: human mosquito)")
    parser.add_argument("--keep-zip",  action="store_true", help="Keep the downloaded zip")
    args = parser.parse_args()

    raw_dir       = Path(args.outdir) / "raw"
    processed_dir = Path(args.outdir) / "processed"
    raw_dir.mkdir(parents=True, exist_ok=True)
    processed_dir.mkdir(parents=True, exist_ok=True)

    check_dependencies()
    log.info("=== DENV genome download pipeline ===")
    log.info("Hosts: %s | Length: %d–%d bp", args.hosts, args.min_len, args.max_len)

    unzip_dir = download_ncbi(raw_dir, args.keep_zip)
    meta_df   = parse_metadata(unzip_dir, raw_dir)
    raw_fasta  = locate_fasta(unzip_dir, raw_dir)

    filter_and_dedup(
        raw_fasta, meta_df, processed_dir,
        host_categories=args.hosts,
        min_len=args.min_len,
        max_len=args.max_len,
    )


if __name__ == "__main__":
    main()
