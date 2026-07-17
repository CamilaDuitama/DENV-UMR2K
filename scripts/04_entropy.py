#!/usr/bin/env python3
"""
04_entropy.py
Per-site Shannon entropy (amino acid level) for all DENV genes,
split by serotype and host (human vs mosquito).

Pipeline per serotype — follows Testa et al. 2026 (Nat Ecol Evol):
  1. Filter metadata → human + mosquito sequences only
  2. cd-hit-est 98% NT clustering to reduce redundancy
  3. Translate polyprotein (from NT pos 97) with EMBOSS transeq
  4. MAFFT --amino to align all polyprotein AA sequences
     (reference AA prepended as anchor)
  5. Map reference AA positions → alignment columns
  6. Extract per-gene AA sub-alignments and split by host
  7. Compute per-site Shannon entropy + standardised entropy (H/log2(20))
  8. Write TSV: serotype/gene/host/site/entropy/entropy_std/n_seqs/n_informative

No NT alignment needed — entropy is computed purely at the AA level.
RevTrans (codon-aware NT alignment) is only required for phydms (WP2).

Requires:
  module load blast+/2.12.0
  module load cd-hit/4.8.1
  module load EMBOSS/6.6.0
  module load mafft/7.526
  conda activate ./env   (pandas, biopython)

Usage:
  python scripts/04_entropy.py [--skip-cdhit] [--threads 8]
  python scripts/04_entropy.py --serotypes DENV2 --skip-cdhit   # quick test
"""

import argparse
import logging
import math
import shutil
import subprocess
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from Bio import SeqIO, AlignIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from Bio.Align import MultipleSeqAlignment

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(message)s",
    handlers=[logging.StreamHandler(sys.stdout)],
)
log = logging.getLogger(__name__)

# ── Gene coordinates ──────────────────────────────────────────────────────────
# NT coordinates in DENV2 reference NC_001474.2 (1-based, inclusive)
POLYPROTEIN_START_NT = 97   # C protein starts here

NT_GENE_COORDS = {
    "C":    (97,   411),
    "prM":  (412,  936),
    "E":    (937,  2421),
    "NS1":  (2422, 3477),
    "NS2A": (3478, 4131),
    "NS2B": (4132, 4521),
    "NS3":  (4522, 6375),
    "NS4A": (6376, 6825),
    "2K":   (6826, 6891),
    "NS4B": (6892, 7650),
    "NS5":  (7651, 10272),
}

# Convert to AA positions within the polyprotein (1-based)
AA_GENE_COORDS = {
    gene: (
        (nt_s - POLYPROTEIN_START_NT) // 3 + 1,
        (nt_e - POLYPROTEIN_START_NT) // 3,
    )
    for gene, (nt_s, nt_e) in NT_GENE_COORDS.items()
}

REF_DIR = Path("data/references")

VACCINE_ACCESSIONS = {
    "MW945952", "AF180818", "U87412", "KU725664",
    "M84728",   "KJ160505", "MW793459", "KJ160504",
}

MIN_SEQS = 5    # minimum sequences per group to compute entropy
LOG2_20  = math.log2(20)   # for entropy standardisation


# ── Helpers ───────────────────────────────────────────────────────────────────

def check_tools() -> None:
    missing = [t for t in ["cd-hit-est", "transeq", "mafft"]
               if not shutil.which(t)]
    if missing:
        log.error("Missing tools: %s — load required modules.", missing)
        sys.exit(1)


def run(cmd: list, desc: str) -> subprocess.CompletedProcess:
    log.info("[%s] %s", desc, " ".join(str(c) for c in cmd))
    r = subprocess.run(cmd, capture_output=True, text=True)
    if r.returncode != 0:
        log.error("[%s] FAILED:\n%s", desc, r.stderr.strip())
        sys.exit(1)
    return r


def write_fasta(fasta: Path, ids: set[str], out: Path) -> int:
    """Write subset matching ids. Returns count written."""
    recs = [r for r in SeqIO.parse(fasta, "fasta")
            if r.id.split(".")[0] in ids]
    if recs:
        SeqIO.write(recs, out, "fasta")
    return len(recs)


def cdhit(fasta: Path, out: Path, threads: int) -> None:
    run(["cd-hit-est", "-i", str(fasta), "-o", str(out),
         "-c", "0.98", "-n", "8", "-T", str(threads),
         "-M", "4000", "-d", "0"], "cd-hit-est 98%")


def translate_polyprotein(nt_fasta: Path, aa_out: Path) -> int:
    """
    For each NT sequence: extract from POLYPROTEIN_START_NT, trim to codon
    boundary, translate with EMBOSS transeq, return number translated.
    """
    poly_nt = aa_out.with_suffix(".poly_nt.fasta")
    n = 0
    with open(poly_nt, "w") as fh:
        for rec in SeqIO.parse(nt_fasta, "fasta"):
            poly = rec.seq[POLYPROTEIN_START_NT - 1:]
            trim = (len(poly) // 3) * 3
            if trim < 90:    # < 30 codons — skip
                continue
            SeqIO.write(
                SeqRecord(poly[:trim], id=rec.id, description=""),
                fh, "fasta"
            )
            n += 1
    if n == 0:
        return 0
    run(["transeq", "-sequence", str(poly_nt), "-outseq", str(aa_out),
         "-frame", "1", "-clean", "-osformat2", "fasta"],
        "transeq polyprotein")
    return n


def mafft_align_aa(aa_fasta: Path, ref_aa: Path, out_aln: Path,
                   threads: int) -> None:
    """Prepend reference AA, run MAFFT --amino."""
    combined = out_aln.parent / (out_aln.stem + "_combined.fasta")
    with open(combined, "w") as fh:
        for rec in SeqIO.parse(ref_aa, "fasta"):
            SeqIO.write(rec, fh, "fasta")
        for rec in SeqIO.parse(aa_fasta, "fasta"):
            SeqIO.write(rec, fh, "fasta")
    result = subprocess.run(
        ["mafft", "--amino", "--auto",
         "--thread", str(threads), "--quiet", str(combined)],
        capture_output=True, text=True,
    )
    if result.returncode != 0:
        log.error("MAFFT AA failed:\n%s", result.stderr)
        sys.exit(1)
    with open(out_aln, "w") as fh:
        fh.write(result.stdout)


def ref_aa_col_map(aln_fasta: Path, ref_prefix: str = "REF_") -> dict[int, int]:
    """Return map: reference AA position (1-based) → alignment column (0-based)."""
    aln = AlignIO.read(aln_fasta, "fasta")
    ref = next(r for r in aln if r.id.startswith(ref_prefix))
    pos_map: dict[int, int] = {}
    aa_pos = 0
    for col, aa in enumerate(str(ref.seq)):
        if aa != "-":
            aa_pos += 1
            pos_map[aa_pos] = col
    return pos_map


def extract_gene_aa(aln_fasta: Path, pos_map: dict[int, int],
                    gene: str, host_ids: set[str]) -> list[str] | None:
    """
    Extract AA alignment columns for gene from aln, restricted to host_ids.
    transeq appends '_1' to sequence IDs — strip it when matching.
    Returns list of AA strings, or None if too few sequences.
    """
    aa_start, aa_end = AA_GENE_COORDS[gene]
    cols = [pos_map[p] for p in range(aa_start, aa_end + 1) if p in pos_map]
    if not cols:
        return None
    aln = AlignIO.read(aln_fasta, "fasta")
    seqs = []
    for rec in aln:
        if rec.id.startswith("REF_"):
            continue
        # transeq appends "_1" to ID — strip it
        base_id = rec.id.rsplit("_", 1)[0] if rec.id.endswith("_1") else rec.id
        base_id = base_id.split(".")[0]
        if base_id not in host_ids:
            continue
        seqs.append("".join(str(rec.seq)[c] for c in cols))
    return seqs if len(seqs) >= MIN_SEQS else None


def shannon_entropy(seqs: list[str]) -> pd.DataFrame:
    """
    Per-site Shannon entropy from a list of equal-length AA strings.
    Returns: site / entropy / entropy_std / entropy_mm / n_seqs / n_informative
      entropy_std = H / log2(20)  — standardised to [0,1], Testa et al. 2026
      entropy_mm  = Miller-Madow bias-corrected entropy:
                    H_MM = H_ML + (K_obs - 1) / (2 * N)
                    where K_obs = number of distinct AAs, N = informative seqs
                    Corrects for downward bias with small sample sizes.
    """
    if not seqs:
        return pd.DataFrame()
    n_sites = len(seqs[0])
    rows = []
    for i in range(n_sites):
        col = [s[i] for s in seqs]
        informative = [aa for aa in col if aa not in "-X*x"]
        n_inf = len(informative)
        if n_inf < 2:
            h = h_std = h_mm = np.nan
        else:
            freqs: dict[str, int] = {}
            for aa in informative:
                freqs[aa] = freqs.get(aa, 0) + 1
            k_obs = len(freqs)
            h = round(-sum((c / n_inf) * math.log2(c / n_inf)
                           for c in freqs.values()), 6)
            h_std = round(h / LOG2_20, 6)
            # Miller-Madow correction: H_MM = H + (K-1)/(2N)
            h_mm = round(h + (k_obs - 1) / (2 * n_inf), 6)
        rows.append({"site": i + 1, "entropy": h, "entropy_std": h_std,
                     "entropy_mm": h_mm,
                     "n_seqs": len(col), "n_informative": n_inf})
    return pd.DataFrame(rows)


# ── Main ──────────────────────────────────────────────────────────────────────

def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--fasta",      default="data/processed/denv_final.fasta")
    parser.add_argument("--meta",       default="data/processed/denv_final_metadata.tsv")
    parser.add_argument("--outdir",     default="data/processed/entropy")
    parser.add_argument("--workdir",    default="data/entropy_work")
    parser.add_argument("--serotypes",  nargs="+",
                        default=["DENV1","DENV2","DENV3","DENV4"])
    parser.add_argument("--hosts",      nargs="+", default=["Human","Mosquito"])
    parser.add_argument("--genes",      nargs="+",
                        default=list(AA_GENE_COORDS.keys()))
    parser.add_argument("--skip-cdhit", action="store_true")
    parser.add_argument("--threads",    type=int, default=8)
    args = parser.parse_args()

    check_tools()

    fasta   = Path(args.fasta)
    out_dir = Path(args.outdir)
    work    = Path(args.workdir)
    out_dir.mkdir(parents=True, exist_ok=True)
    work.mkdir(parents=True, exist_ok=True)

    meta = pd.read_csv(args.meta, sep="\t", low_memory=False)
    meta["accession"] = meta["accession"].astype(str).str.split(".").str[0]
    meta = meta[~meta["accession"].isin(VACCINE_ACCESSIONS)]
    meta = meta[meta["is_lab_host"].fillna("").astype(str).str.lower() != "true"]

    all_results: list[pd.DataFrame] = []

    for sero in args.serotypes:
        ref_nt = REF_DIR / f"ref_{sero}.fasta"
        if not ref_nt.exists():
            log.warning("Reference missing for %s — skipping", sero)
            continue

        sero_dir = work / sero
        sero_dir.mkdir(exist_ok=True)

        # ── Translate reference polyprotein ────────────────────────────────
        ref_aa_path = sero_dir / f"ref_{sero}_aa.fasta"
        ref_poly_nt = sero_dir / f"ref_{sero}_poly_nt.fasta"
        for rec in SeqIO.parse(ref_nt, "fasta"):
            poly = rec.seq[POLYPROTEIN_START_NT - 1:]
            trim = (len(poly) // 3) * 3
            SeqIO.write(
                SeqRecord(poly[:trim], id=f"REF_{sero}", description=""),
                ref_poly_nt, "fasta"
            )
        run(["transeq", "-sequence", str(ref_poly_nt), "-outseq", str(ref_aa_path),
             "-frame", "1", "-clean", "-osformat2", "fasta"], f"transeq ref {sero}")

        # ── Collect sequences ──────────────────────────────────────────────
        sero_ids = set(
            meta[(meta["serotype"] == sero) &
                 (meta["host_type"].isin(args.hosts))]["accession"]
        )
        sero_fasta = sero_dir / "all_nt.fasta"
        n_all = write_fasta(fasta, sero_ids, sero_fasta)
        log.info("%s: %d sequences", sero, n_all)
        if n_all == 0:
            continue

        # ── cd-hit ─────────────────────────────────────────────────────────
        if args.skip_cdhit:
            clustered = sero_fasta
        else:
            clustered = sero_dir / "clustered_nt.fasta"
            cdhit(sero_fasta, clustered, args.threads)

        clustered_ids = {r.id.split(".")[0]
                         for r in SeqIO.parse(clustered, "fasta")}
        n_clustered = len(clustered_ids)
        log.info("%s: %d after cd-hit", sero, n_clustered)

        # Build host → set of retained accessions
        host_id_map: dict[str, set[str]] = {}
        for host in args.hosts:
            host_id_map[host] = set(
                meta[(meta["serotype"] == sero) &
                     (meta["host_type"] == host)]["accession"]
            ) & clustered_ids
            log.info("%s / %s: %d sequences", sero, host,
                     len(host_id_map[host]))

        # ── Translate polyprotein ──────────────────────────────────────────
        aa_raw = sero_dir / "all_aa_raw.fasta"
        n_trans = translate_polyprotein(clustered, aa_raw)
        log.info("%s: %d sequences translated", sero, n_trans)
        if n_trans == 0:
            continue

        # ── Align AA polyprotein ───────────────────────────────────────────
        aa_aln = sero_dir / "all_aa_aligned.fasta"
        mafft_align_aa(aa_raw, ref_aa_path, aa_aln, args.threads)
        log.info("%s: AA alignment complete", sero)

        # ── Reference AA position → alignment column map ───────────────────
        pos_map = ref_aa_col_map(aa_aln)

        # ── Per-gene, per-host entropy ─────────────────────────────────────
        for gene in args.genes:
            if gene not in AA_GENE_COORDS:
                continue
            for host in args.hosts:
                host_ids = host_id_map.get(host, set())
                if len(host_ids) < MIN_SEQS:
                    log.info("%s / %s / %s: only %d seqs — skipping",
                             sero, gene, host, len(host_ids))
                    continue

                seqs = extract_gene_aa(aa_aln, pos_map, gene, host_ids)
                if seqs is None:
                    log.info("%s / %s / %s: too few after extraction",
                             sero, gene, host)
                    continue

                df = shannon_entropy(seqs)
                if df.empty:
                    continue

                df.insert(0, "host",          host)
                df.insert(0, "gene",          gene)
                df.insert(0, "serotype",      sero)
                df["n_seqs_total"]   = len(seqs)
                df["n_seqs_cluster"] = n_clustered

                all_results.append(df)
                log.info("%s / %s / %s: %d sites, %d seqs, mean H=%.3f, "
                         "mean H_std=%.3f",
                         sero, gene, host, len(df), len(seqs),
                         df["entropy"].mean(), df["entropy_std"].mean())

    if not all_results:
        log.warning("No results produced.")
        return

    result_df = pd.concat(all_results, ignore_index=True)
    out_path = out_dir / "entropy_per_site.tsv"
    result_df.to_csv(out_path, sep="\t", index=False)
    log.info("=== Results → %s (%d rows) ===", out_path, len(result_df))

    summary = (result_df.groupby(["serotype", "gene", "host"])
               .agg(n_sites=("site", "count"),
                    n_seqs=("n_seqs_total", "first"),
                    mean_H=("entropy", "mean"),
                    mean_H_std=("entropy_std", "mean"),
                    max_H=("entropy", "max"))
               .round(4))
    log.info("\n%s", summary.to_string())


if __name__ == "__main__":
    main()
