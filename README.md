# DENV NS1 Evolution — CNRS UMR2K Seed Grant 2026

**Bioinformatician:** Camila Duitama (Institut Pasteur)

## Overview

Bioinformatics component of the UMR2K 2026 project integrating deep mutational scanning (DMS) data of dengue virus (DENV) NS1 into phylogenetic models to dissect host-specific selective pressures (human vs. mosquito).

## Objectives

- **WP1:** Compare experimental mutation fitness effects (MFEs) against natural DENV NS1 diversity across serotypes and host environments.
- **WP2:** Benchmark experimentally-informed phylogenetic models (phydms) against classical substitution models using model adequacy (phyloMAd).

## Data Collection (this repo)

Focus: DENV **full genome** sequences from **human** and **mosquito** hosts.

| Serotype | Priority |
|----------|----------|
| DENV2 | Primary |
| DENV1/3/4 | Survey (availability-dependent) |

### Sources under evaluation

| Source | Access | Notes |
|--------|--------|-------|
| [NCBI Virus / GenBank](https://www.ncbi.nlm.nih.gov/labs/virus/) | Public | Curated full genomes, rich metadata |
| [NCBI SRA](https://www.ncbi.nlm.nih.gov/sra) | Public | Raw reads, requires assembly |
| [Nextstrain / nextstrain.org](https://nextstrain.org/dengue) | Public | Pre-filtered, phylogenetically curated |
| [ViPR](https://www.viprbrc.org/) | Public | Virus Pathogen Resource, DENV-specific |

## Repository Structure

```
UMR2K/
├── data/
│   ├── raw/          # downloaded sequences and metadata
│   └── processed/    # filtered, aligned datasets
├── scripts/          # download and QC scripts
├── docs/             # GitHub Pages site source
└── README.md
```

## GitHub Pages

A metadata summary website will be auto-generated at `docs/` and published via GitHub Pages, showing sequence counts, host distribution, geography, and temporal coverage of the final dataset.

## Dependencies

- Python ≥ 3.10
- `Biopython`, `pandas`, `requests`
- NCBI `datasets` CLI / `Entrez Direct`
- Nextstrain CLI (optional)

## License

MIT
