# TODO — DENV NS1 UMR2K Project

## WP1 — Data collection

- [x] Repo, git, environment.yaml (conda env at `./env`)
- [x] `01_query_availability.py` — survey metadata via `ncbi-datasets` CLI + Nextstrain
- [x] `02_download_genomes.py` — Nextstrain-first download, NCBI complement, seqkit filter
- [x] `03_build_site.py` — GitHub Pages interactive report (Plotly)
- [x] Availability survey: 58,641 NCBI sequences; 21,630 Nextstrain full genomes (human+mosquito, ≥10kb, cov≥0.95)
- [x] GitHub Pages live: https://camiladuitama.github.io/DENV-UMR2K/
- [x] Beatriz review received (2026-07-10):
  - **Mosquito species**: keep all for now, filter later; primary focus *Aedes aegypti*
  - **Gene of interest**: ⚠️ possibly **NS4A-2K-NS4B** instead of NS1 — TBC once mutagenesis system optimised
  - **Geography / sampling**: include everything; subsampling per-analysis later
  - **Initial analysis requested**: Shannon entropy per site (AA level) across all proteins
  - **Reference methods**: https://www.nature.com/articles/s41559-026-02993-8#Sec10
- [x] Vaccine/lab sequences checked and excluded:
  - 170 via `is_lab_host = True` (lab cell-line passages)
  - 8 PDK-passage / chimeric strains (hardcoded in `02_download_genomes.py`)
- [x] **Download complete** — `data/processed/denv_final.fasta` (21,634 seqs, 220 MB, Nextstrain-only)
  - NCBI complement pending: bug fixed (`--include genome`); resubmit `run_download.slurm`
- [x] QC added to website: seqkit stats, serotype × host breakdown, length distribution

## WP1 — Natural diversity analysis

> **Note**: pipeline below follows the reference methods (Testa et al. 2026 Nat Ecol Evol).
> Steps marked ⚙️ are confirmed; steps marked 🔄 may change depending on boss/collaborator decisions.

- ⚙️ **cd-hit clustering at 98% nt** to remove near-identical sequences  
  `module load cd-hit/4.8.1`  
  Reduces ~21k → ~8-10k non-redundant sequences; required for feasible phydms analysis
- ⚙️ **Per-gene extraction** from aligned full genomes (DENV2 NC_001474 reference coordinates):
  - C: 97–411 | prM: 412–936 | E: 937–2421 | NS1: 2422–3477
  - NS2A: 3478–4131 | NS2B: 4132–4521 | NS3: 4522–6375
  - NS4A: 6376–6825 | 2K: 6826–6891 | **NS4B: 6892–7650** | NS5: 7651–10272
  - ⚠️ **NS4A-2K-NS4B region** likely the study target — confirm with Beatriz
- ⚙️ **Codon-aware alignment**: MAFFT (`module load mafft`) per gene per serotype
- ⚙️ **Translation to AA**: EMBOSS `transeq` (`module load EMBOSS/6.6.0`)
- ⚙️ **Shannon entropy per site** per protein, split by host (human vs mosquito)  
  Formula: $H = -\sum_i p_i \log_2(p_i)$ where $p_i$ = AA frequency at site  
  Custom Python script (using numpy)
- 🔄 **Phylogeny**: IQ-TREE 3.1.0 (`module load IQ-TREE/3.1.0`) for QC tree  
  (Paper used RAxML 8.2.12 with PROTGAMMAWAG; both available)
- ⚙️ Compare site-level entropy/diversity with experimental MFEs once gene confirmed

## WP2 — DMS data integration

- [ ] ⚠️ **Wait for Beatriz to confirm target gene** (NS1 vs NS4A-2K-NS4B)
- [ ] Find published DMS dataset for target gene:
  - If NS1: Dolan et al. 2021 *Cell Host Microbe* / Bloom lab GitHub
  - If NS4A-2K-NS4B: search Bloom lab + literature
  - Note: `dms_tools2` legacy; current tools are `dms-variants` / `polyclonal`
- [ ] Obtain MFE files (AA preference CSVs, one per host condition)
- [ ] **cd-hit at 95% AA** before phydms (paper used this to stay within phydms limits)
- [ ] Integrate MFEs → site preferences → `phydms`
- [ ] Benchmark phydms vs. standard models (M0, M8, GY94)

## WP2 — Phylogenetic modeling

- [ ] Install and test `phydms` (https://jbloomlab.github.io/phydms/)
- [ ] Test `phyloMAd` for model adequacy (dev — contact authors)
- [ ] Run at serotype scale (DENV2) then species scale (all 4 serotypes)

## Infrastructure

- [x] GitHub repo: https://github.com/CamilaDuitama/DENV-UMR2K
- [x] GitHub Pages: https://camiladuitama.github.io/DENV-UMR2K/
- [x] GitHub Actions auto-deploy on push
- [ ] Add NCBI API key as GitHub secret (10× faster NCBI downloads)
- [ ] Tag data freeze version when dataset is locked
- [ ] Add `CITATION.cff`
