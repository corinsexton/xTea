#!/bin/bash
#SBATCH -J mv2_e2e035
#SBATCH -A park_contrib
#SBATCH -p priopark
#SBATCH -c 4
#SBATCH --mem=24G
#SBATCH -t 1:00:00
#SBATCH -o /n/scratch/users/c/cos689/xtea2_review/mosaic_v2/e2e/logs/%x_%j.out
# Rescore the end-to-end WashU_100X run at the fixed 0.35 ceiling (and 0.30) on its own post-filter
# inputs: postF_mosaic + shared tail rerun with the counting_v2 tune code copy.
set -euo pipefail
source /n/app/conda/miniforge3/24.11.3-0/etc/profile.d/conda.sh; conda activate xtea2; export PYTHONNOUSERSITE=1
export PATH=/home/cos689/software/bin:$PATH
B=/n/scratch/users/c/cos689/xtea2_review/mosaic_v2
python $B/scripts/sweep_fixed_e2e.py --run $B/e2e/runs/WashU_100X --name WashU_100X --mode mosaic --param disc \
  --values=0.30,0.35 --out $B/e2e/fixed --tsv $B/e2e/fixed_sweep.tsv -n 4
