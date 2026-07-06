#!/bin/bash
#SBATCH --job-name=denv_query
#SBATCH --partition=edid
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --time=04:00:00
#SBATCH --output=logs/denv_query_%j.out
#SBATCH --error=logs/denv_query_%j.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=cduitama@pasteur.fr

# ── Environment ──────────────────────────────────────────────────────────────
set -euo pipefail
WORKDIR="/pasteur/helix/scratch/cduitama/EDID/UMR2K"
cd "$WORKDIR"
mkdir -p logs

echo "== Job started on $(hostname) at $(date) =="
echo "== Working dir: $WORKDIR =="

# Load ncbi-datasets CLI (provides `datasets` command)
module load ncbi-datasets/v2

# Activate conda env
source /pasteur/appa/homes/cduitama/miniconda3/etc/profile.d/conda.sh
conda activate denv-umr2k

echo "== Python: $(which python) =="
echo "== Python version: $(python --version) =="

# ── Run ──────────────────────────────────────────────────────────────────────
python scripts/01_query_availability.py

echo "== Job finished at $(date) =="
