# DENV NS1 Evolution — CNRS UMR2K Seed Grant 2026

**Author:** Camila Duitama (Institut Pasteur)

## Overview

Bioinformatics component of the UMR2K 2026 project integrating deep mutational scanning (DMS) data of DENV NS1 into phylogenetic models to dissect host-specific selective pressures (human vs. mosquito).

**WP1:** Compare experimental MFEs against natural DENV NS1 diversity across serotypes and host environments.  
**WP2:** Benchmark experimentally-informed phylogenetic models (phydms) against classical substitution models (phyloMAd).

## Data source

**Nextstrain DENV** (primary) — pre-curated, QC'd, clade-annotated sequences from GenBank.  
**NCBI GenBank** (complement) — sequences not in Nextstrain, downloaded via `ncbi-datasets` CLI.

Full-genome criterion: `length ≥ 10,000 bp` AND `genome_coverage ≥ 0.95` AND `is_lab_host ≠ True`, human + mosquito hosts only.

## Current status

Download complete: **21,630 full genomes** in `data/processed/denv_final.fasta` (Nextstrain-only; 220 MB).  
Interactive metadata report: **https://camiladuitama.github.io/DENV-UMR2K/**

## Repository structure

```
UMR2K/
├── data/raw/          # survey metadata (Nextstrain TSVs, NCBI summary)
├── data/processed/    # filtered FASTA + metadata (post-download)
├── scripts/
│   ├── 01_query_availability.py  # metadata survey (Nextstrain + NCBI)
│   ├── 02_download_genomes.py    # Nextstrain-first download + seqkit filter
│   ├── 03_build_site.py          # GitHub Pages report generator
│   ├── run_query.slurm           # SLURM: survey job
│   └── run_download.slurm        # SLURM: download job (64G, 8 CPU)
├── docs/index.html    # GitHub Pages site
├── environment.yaml   # conda env
└── env/               # local conda env (gitignored)
```

## Cluster usage

```bash
module load ncbi-datasets/v2
module load SeqKit/2.8.2
mamba env create -f environment.yaml --prefix ./env
conda activate ./env

sbatch scripts/run_query.slurm      # availability survey
python scripts/03_build_site.py     # rebuild site (submit node)
sbatch scripts/run_download.slurm   # full download (after approval)
```

## License

MIT
