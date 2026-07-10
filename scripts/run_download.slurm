#!/bin/bash
#SBATCH --job-name=denv_download
#SBATCH --partition=edid
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --output=logs/denv_download_%j.out
#SBATCH --error=logs/denv_download_%j.err
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=cduitama@pasteur.fr

# ── Environment ───────────────────────────────────────────────────────────────
set -euo pipefail
WORKDIR="/pasteur/helix/scratch/cduitama/EDID/UMR2K"
cd "$WORKDIR"
mkdir -p logs data/raw data/processed

echo "== Job started on $(hostname) at $(date) =="

# Load modules
module load ncbi-datasets/v2
module load SeqKit/2.8.2

# Activate local conda env
source /pasteur/appa/homes/cduitama/miniconda3/etc/profile.d/conda.sh
conda activate /pasteur/helix/scratch/cduitama/EDID/UMR2K/env

echo "== Python:   $(which python) =="
echo "== datasets: $(which datasets) =="
echo "== seqkit:   $(which seqkit) =="

# Optional: set NCBI API key for faster downloads (10 req/s vs 3 req/s)
# export NCBI_API_KEY="your_key_here"

# ── Run download + filter pipeline ───────────────────────────────────────────
python scripts/02_download_genomes.py --outdir data

echo "== Job finished at $(date) =="
