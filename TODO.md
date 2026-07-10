# TODO — DENV NS1 UMR2K Project

## WP1 — Data collection

- [x] Repo, git, environment.yaml (conda env at `./env`)
- [x] `01_query_availability.py` — survey metadata via `ncbi-datasets` CLI + Nextstrain
- [x] `02_download_genomes.py` — Nextstrain-first download, NCBI complement, seqkit filter
- [x] `03_build_site.py` — GitHub Pages interactive report (Plotly)
- [x] Availability survey run — 58,641 NCBI sequences; ~15k Nextstrain full genomes (human+mosquito, ≥10kb, cov≥0.95)
- [x] GitHub Pages live: https://camiladuitama.github.io/DENV-UMR2K/
- [x] Beatriz review received (2026-07-10):
  - **Mosquito species**: keep all, filter later if needed; primary focus *Aedes aegypti*
  - **Gene of interest**: ⚠️ possibly **NS4A-2K-NS4B** instead of NS1 — TBC once mutagenesis system is optimised
  - **Geography / sampling**: include everything now; apply subsampling per-analysis later
  - **Vaccine-derived sequences**: check whether these exist in DENV Nextstrain data and remove if significant
  - **Initial analysis requested**: Shannon entropy per site (amino acid level) across all proteins
  - **Reference**: https://www.nature.com/articles/s41559-026-02993-8#Sec10 (methods: natural diversity analysis)
- [ ] **Check for vaccine-derived sequences** in Nextstrain DENV metadata (is_vaccine_strain field?)
- [ ] Run download job (`sbatch scripts/run_download.slurm`)
- [ ] QC final dataset: deduplication, length/coverage checks

## WP1 — Natural diversity analysis

- [ ] Multiple sequence alignment (`module load mafft`; full genome per serotype)
- [ ] Extract coding regions per gene — reference DENV2 NC_001474 coordinates:
  - C: 97–411; prM: 412–936; E: 937–2421; NS1: 2422–3477
  - NS2A: 3478–4131; NS2B: 4132–4521; NS3: 4522–6375
  - NS4A: 6376–6825; 2K: 6826–6891; NS4B: 6892–7650; NS5: 7651–10272
  - ⚠️ **NS4A-2K-NS4B** likely the target gene — confirm with Beatriz
- [ ] Translate to amino acids per gene
- [ ] **Shannon entropy per site** per protein (per host: human vs mosquito)
- [ ] Phylogeny (IQ-TREE or FastTree) for QC and context
- [ ] Compare site-level diversity with experimental MFEs once gene is confirmed

## WP2 — DMS data integration

- [ ] ⚠️ **Wait for Beatriz to confirm target gene** (NS1 vs NS4A-2K-NS4B) before sourcing DMS data
- [ ] Find published DMS dataset for confirmed target gene
  - If NS1: Dolan et al. 2021 *Cell Host Microbe*; check Bloom lab GitHub
  - If NS4A-2K-NS4B: literature search needed
  - Note: `dms_tools2` is legacy; current tools are `dms-variants` / `polyclonal`
- [ ] Obtain MFE files (amino-acid preference CSVs, one per host condition)
- [ ] Integrate MFEs into `phydms` codon model
- [ ] Benchmark phydms vs. standard models (M0, M8, GY94)

## WP2 — Phylogenetic modeling

- [ ] Install and test `phydms` (https://jbloomlab.github.io/phydms/)
- [ ] Test `phyloMAd` for model adequacy (dev — contact authors)
- [ ] Run at serotype scale (DENV2) then species scale (all 4 serotypes)

## Infrastructure

- [x] GitHub repo: https://github.com/CamilaDuitama/DENV-UMR2K
- [x] GitHub Pages live from `docs/`
- [x] GitHub Actions to auto-rebuild site on push
- [ ] Add NCBI API key as GitHub secret
- [ ] Tag data freeze version when dataset is locked
- [ ] Add `CITATION.cff`
