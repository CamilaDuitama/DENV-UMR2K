# TODO — DENV NS1 UMR2K Project

## WP1 — Data collection

- [x] Repo, git, environment.yaml (conda env at `./env`)
- [x] `01_query_availability.py` — survey metadata via `ncbi-datasets` CLI + Nextstrain
- [x] `02_download_genomes.py` — Nextstrain-first download, NCBI complement, seqkit filter
- [x] `03_build_site.py` — GitHub Pages interactive report (Plotly)
- [x] Availability survey run — 58,641 NCBI sequences; ~15k Nextstrain full genomes (human+mosquito, ≥10kb, cov≥0.95)
- [x] GitHub Pages live: https://camiladuitama.github.io/DENV-UMR2K/
- [x] Email sent to Beatriz for dataset review
- [ ] **WAITING: Beatriz reply** — confirm mosquito species, genes of interest, exclusions
- [ ] Run download job after approval (`sbatch scripts/run_download.slurm`)
- [ ] QC final dataset: deduplication, length/coverage checks

## WP1 — Natural diversity analysis

- [ ] Multiple sequence alignment (`module load mafft`)
- [ ] Extract NS1 coding region (nt 2389–4149, DENV2 reference NC_001474)
- [ ] Phylogeny (IQ-TREE or FastTree)
- [ ] Natural amino-acid frequencies per NS1 site
- [ ] Compare natural frequencies with experimental MFEs

## WP2 — DMS data integration

- [ ] Find published DENV NS1 DMS dataset
  - Dolan et al. 2021 *Cell Host Microbe* (flavivirus NS1)
  - Check Bloom lab GitHub: https://github.com/jbloomlab
  - Note: `dms_tools2` is legacy; newer tools are `dms-variants` / `polyclonal`
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
- [ ] GitHub Actions to auto-rebuild site on push
- [ ] Add NCBI API key as GitHub secret
- [ ] Tag data freeze version when dataset is locked
- [ ] Add `CITATION.cff`
