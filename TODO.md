# TODO — DENV NS1 UMR2K Project

## WP1 — Data collection

- [x] Set up repo scaffold and git
- [x] Write `environment.yaml` (plotly, zstandard, biopython, fasttree)
- [x] Rewrite `01_query_availability.py` to use `ncbi-datasets` CLI + Nextstrain
- [ ] **Run availability survey** (`sbatch scripts/run_query.slurm`) and review counts
- [ ] Decide final host filter: `Homo sapiens` + `Aedes aegypti` + `Aedes albopictus`
- [ ] Download full DENV2 genomes (FASTA + metadata) using `datasets download virus genome`
- [ ] Download other serotypes (DENV1/3/4) for context — assess counts from survey
- [ ] QC filter: completeness, length (≥ 10 kb), remove duplicates, remove low-quality

## WP1 — Natural diversity analysis

- [ ] Multiple sequence alignment (module load mafft; DENV2 full genome, then NS1 region)
- [ ] Extract NS1 coding region (nt 2389–4149 in DENV2 reference NC_001474)
- [ ] Infer phylogeny (IQ-TREE or FastTree) for QC and serotype-level analyses
- [ ] Compute natural amino-acid frequencies per NS1 site
- [ ] Compare natural frequencies with experimental MFEs (WP1 core analysis)

## WP2 — DMS data integration

- [ ] **Understand `dms_tools2` pipeline** (Bloom lab):
  - Repo: https://github.com/jbloomlab/dms_tools2
  - Docs: https://jbloomlab.github.io/dms_tools2/
  - Key tools: `dms2_bcsubamp` (read processing), `dms2_prefs` (amino-acid preferences),
    `dms2_batch_diffsel` (differential selection between conditions/hosts)
  - **Input data format**: paired-end Illumina reads from barcoded subamplicon library.
    Each library = a pool of all codon-mutants of the gene of interest, passaged in
    one condition (e.g., mosquito cells C6/36 or human cells Vero/Huh7).
    After sequencing: FASTQ → `dms2_bcsubamp` → per-codon counts CSV → `dms2_prefs`
    → amino-acid preference files (CSV) → site-preference logo plots.
  - **NOTE**: `dms_tools2` is legacy (Python 3.6+). Current Bloom lab tools are
    `dms-variants` and `polyclonal`. Check whether the DENV NS1 DMS dataset
    (Dolan et al. or equivalent) was processed with `dms_tools2` or newer tools.
  - **TODO**: Find the published DENV NS1 DMS dataset. Likely candidates:
    - Dolan et al. 2021 Cell Host Microbe (flavivirus NS1)
    - Check Bloom lab GitHub for DENV NS1 specific repo
- [ ] Obtain MFE data files (amino-acid preferences CSVs, one per host condition)
- [ ] Integrate MFEs into phydms codon model framework
- [ ] Benchmark phydms vs. standard models (M0, M8, GY94) on DENV2 NS1 tree

## WP2 — Phylogenetic modeling

- [ ] Install and test `phydms` (https://jbloomlab.github.io/phydms/)
- [ ] Test `phyloMAd` for model adequacy assessment (dev version — contact authors)
- [ ] Run phydms at serotype scale (DENV2 only) then species scale (all 4 serotypes)

## GitHub Pages website

- [ ] Design site layout: summary stats, interactive Plotly charts, metadata table
- [ ] Implement `scripts/build_site.py` → renders `docs/index.html` from Jinja2 + Plotly
- [ ] Charts to include:
  - Bar chart: n genomes by serotype × host × source (NCBI vs Nextstrain)
  - Choropleth map: geographic distribution of sequences
  - Timeline: collection date histogram
- [ ] Set up GitHub Actions to auto-rebuild site on push
- [ ] Enable GitHub Pages from `docs/` folder

## Infrastructure / reproducibility

- [ ] Register GitHub repo (camiladuitama/DENV-UMR2K or similar)
- [ ] Add NCBI API key as GitHub secret (speeds up `datasets` calls 3x → 10 req/s)
- [ ] Tag data freeze version when final dataset is locked
- [ ] Add `CITATION.cff`
