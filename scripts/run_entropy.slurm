#!/bin/bash
#SBATCH --job-name=denv_entropy
#SBATCH --partition=edid
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --output=logs/denv_entropy_%j.out
#SBATCH --error=logs/denv_entropy_%j.err
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=cduitama@pasteur.fr

set -euo pipefail
WORKDIR="/pasteur/helix/scratch/cduitama/EDID/UMR2K"
cd "$WORKDIR"
mkdir -p logs data/processed/entropy data/entropy_work

echo "== Job started on $(hostname) at $(date) =="

module load blast+/2.12.0
module load fasta
module load ruby
module load mafft/7.526
module load cd-hit/4.8.1
module load EMBOSS/6.6.0
module load SeqKit/2.8.2

source /pasteur/appa/homes/cduitama/miniconda3/etc/profile.d/conda.sh
conda activate /pasteur/helix/scratch/cduitama/EDID/UMR2K/env

echo "== Python:  $(which python) =="
echo "== mafft:   $(which mafft) =="
echo "== cd-hit:  $(which cdhit) =="
echo "== transeq: $(which transeq) =="

# Run entropy pipeline
# --skip-cdhit can be added for a quick test run
python scripts/04_entropy.py \
    --fasta   data/processed/denv_final.fasta \
    --meta    data/processed/denv_final_metadata.tsv \
    --outdir  data/processed/entropy \
    --workdir data/entropy_work \
    --threads 8

echo "== Entropy done. Building HTML page =="

python scripts/05_build_entropy_page.py

echo "== Job finished at $(date) =="
