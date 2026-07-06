#!/bin/bash
#SBATCH --job-name=denv_query
#SBATCH --partition=edid
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
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

# Load modules
module load ncbi-datasets/v2

# Activate conda env (installed locally in project dir)
source /pasteur/appa/homes/cduitama/miniconda3/etc/profile.d/conda.sh
conda activate /pasteur/helix/scratch/cduitama/EDID/UMR2K/env

echo "== Python: $(which python) =="
echo "== Python version: $(python --version) =="

# ── Run ──────────────────────────────────────────────────────────────────────
python scripts/01_query_availability.py

echo "== Job finished at $(date) =="
