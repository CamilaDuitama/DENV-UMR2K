#!/usr/bin/env python3
"""
01_query_availability.py
Query NCBI Virus (nuccore) and Nextstrain for DENV full genome availability
by serotype and host (human / mosquito).

Output:
    data/raw/ncbi_denv_summary.tsv   -- per-accession metadata from NCBI
    data/raw/availability_report.tsv -- counts table by source/serotype/host
    data/raw/nextstrain_denv2_metadata.tsv -- Nextstrain DENV2 metadata if reachable

Usage (on cluster after activating conda env):
    python scripts/01_query_availability.py
"""

import os
import sys
import time
import csv
import json
import logging
import requests
import pandas as pd
from pathlib import Path
from Bio import Entrez

# ── Config ──────────────────────────────────────────────────────────────────
Entrez.email = "cduitama@pasteur.fr"   # required by NCBI
OUTDIR = Path("data/raw")
OUTDIR.mkdir(parents=True, exist_ok=True)

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(message)s",
    handlers=[logging.StreamHandler(sys.stdout)],
)
log = logging.getLogger(__name__)

# NCBI search parameters
# Complete genome = RefSeq + GenBank complete genomes
# Host terms: human (Homo sapiens) or mosquito (Aedes, Culex, Stegomyia)
SEROTYPES = {
    "DENV1": "Dengue virus 1[Organism]",
    "DENV2": "Dengue virus 2[Organism]",
    "DENV3": "Dengue virus 3[Organism]",
    "DENV4": "Dengue virus 4[Organism]",
}
HOSTS = {
    "human":    '("Homo sapiens"[Host] OR "human"[Host])',
    "mosquito": '("Aedes"[Host] OR "Culex"[Host] OR "mosquito"[Host] OR "Aedes aegypti"[Host] OR "Aedes albopictus"[Host])',
}
GENOME_FILTER = '"complete genome"[Title] OR "complete sequence"[Title]'
DB = "nuccore"
RETMAX = 100_000  # upper bound; NCBI will cap at actual count

# Nextstrain DENV2 metadata (public S3 or data.nextstrain.org)
NEXTSTRAIN_DENV2_META = (
    "https://data.nextstrain.org/files/dengue/denv2/metadata.tsv.zst"
)

# ── Helpers ──────────────────────────────────────────────────────────────────

def esearch_count(query: str) -> int:
    """Return hit count for an NCBI esearch query."""
    with Entrez.esearch(db=DB, term=query, retmax=0) as h:
        rec = Entrez.read(h)
    return int(rec["Count"])


def efetch_summaries(query: str, retmax: int = RETMAX) -> list[dict]:
    """Fetch DocSummary records for a query. Returns list of summary dicts."""
    log.info("esearch: %s", query)
    with Entrez.esearch(db=DB, term=query, retmax=retmax, usehistory="y") as h:
        search_results = Entrez.read(h)
    count = int(search_results["Count"])
    webenv = search_results["WebEnv"]
    query_key = search_results["QueryKey"]
    log.info("  → %d records found", count)

    records = []
    batch = 500
    for start in range(0, min(count, retmax), batch):
        log.info("  fetching %d–%d …", start, start + batch)
        with Entrez.esummary(
            db=DB,
            webenv=webenv,
            query_key=query_key,
            retstart=start,
            retmax=batch,
            rettype="json",
        ) as h:
            data = Entrez.read(h)
        for uid in data["DocumentSummarySet"]["DocumentSummary"]:
            records.append(dict(uid))
        time.sleep(0.35)  # NCBI rate limit: ≤3 req/s without API key
    return records


def summarize_record(rec: dict, serotype: str, host_label: str) -> dict:
    return {
        "accession": rec.get("AccessionVersion", rec.get("Caption", "")),
        "title": rec.get("Title", ""),
        "length": rec.get("Slen", ""),
        "organism": rec.get("Organism", ""),
        "host": rec.get("SubType", "") + "|" + rec.get("SubName", ""),
        "country": "",   # not in DocSummary by default
        "collection_date": "",
        "serotype": serotype,
        "host_query": host_label,
    }


# ── Main ─────────────────────────────────────────────────────────────────────

def query_ncbi() -> pd.DataFrame:
    rows = []
    counts = []

    for st, org_term in SEROTYPES.items():
        for host_label, host_term in HOSTS.items():
            query = f"({org_term}) AND ({host_term}) AND ({GENOME_FILTER})"
            count = esearch_count(query)
            log.info("%s / %s: %d sequences", st, host_label, count)
            counts.append({"source": "NCBI", "serotype": st, "host": host_label, "n_sequences": count})

    # Full metadata fetch for DENV2 only (main target)
    denv2_records = []
    for host_label, host_term in HOSTS.items():
        query = f"({SEROTYPES['DENV2']}) AND ({host_term}) AND ({GENOME_FILTER})"
        recs = efetch_summaries(query)
        for r in recs:
            denv2_records.append(summarize_record(r, "DENV2", host_label))

    meta_df = pd.DataFrame(denv2_records)
    meta_path = OUTDIR / "ncbi_denv2_metadata.tsv"
    meta_df.to_csv(meta_path, sep="\t", index=False)
    log.info("Saved DENV2 metadata → %s (%d rows)", meta_path, len(meta_df))

    counts_df = pd.DataFrame(counts)
    return counts_df


def query_nextstrain() -> pd.DataFrame | None:
    """Download Nextstrain DENV2 metadata (zstd-compressed TSV)."""
    try:
        import zstandard as zstd  # noqa: F401 – optional
    except ImportError:
        log.warning("zstandard not installed; trying uncompressed fallback")

    # Try plain TSV first (some nextstrain datasets expose uncompressed)
    url_plain = "https://data.nextstrain.org/files/dengue/denv2/metadata.tsv"
    try:
        log.info("Fetching Nextstrain DENV2 metadata …")
        r = requests.get(url_plain, timeout=60, stream=True)
        if r.status_code == 200:
            meta_path = OUTDIR / "nextstrain_denv2_metadata.tsv"
            with open(meta_path, "wb") as f:
                for chunk in r.iter_content(chunk_size=65536):
                    f.write(chunk)
            df = pd.read_csv(meta_path, sep="\t", low_memory=False)
            log.info("Nextstrain DENV2: %d sequences", len(df))
            # Count by host
            if "host" in df.columns:
                host_counts = df.groupby("host").size().reset_index(name="n_sequences")
                host_counts["source"] = "Nextstrain"
                host_counts["serotype"] = "DENV2"
                return host_counts
            return pd.DataFrame([{"source": "Nextstrain", "serotype": "DENV2", "host": "unknown", "n_sequences": len(df)}])
        else:
            log.warning("Nextstrain returned HTTP %s", r.status_code)
    except Exception as e:
        log.warning("Nextstrain fetch failed: %s", e)
    return None


def main():
    log.info("=== DENV sequence availability survey ===")

    ncbi_counts = query_ncbi()

    ns_counts = query_nextstrain()
    if ns_counts is not None:
        all_counts = pd.concat([ncbi_counts, ns_counts], ignore_index=True)
    else:
        all_counts = ncbi_counts

    report_path = OUTDIR / "availability_report.tsv"
    all_counts.to_csv(report_path, sep="\t", index=False)
    log.info("=== Availability report saved → %s ===", report_path)
    print(all_counts.to_string(index=False))


if __name__ == "__main__":
    main()
