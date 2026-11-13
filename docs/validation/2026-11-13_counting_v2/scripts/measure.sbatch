#!/bin/bash
#SBATCH -J hclip_meas
#SBATCH -A park_contrib
#SBATCH -p priopark
#SBATCH -c 1
#SBATCH --mem=4G
#SBATCH -t 1:00:00
#SBATCH -o /n/scratch/users/c/cos689/xtea2_review/counting_v2/logs/%x_%A_%a.out
set -euo pipefail
G=/n/scratch/users/c/cos689/xtea2_testdata/gtval
read -r ID NAME BAM REF < <(awk -v i="$SLURM_ARRAY_TASK_ID" '$1==i' $G/manifest.tsv)
source /n/app/conda/miniforge3/24.11.3-0/etc/profile.d/conda.sh; conda activate xtea2
export PYTHONNOUSERSITE=1
M=/n/scratch/users/c/cos689/xtea2_review/counting_v2/measure
python $M/measure_hclip.py $G/runs_bug2_tuned "$NAME" "$BAM" "$REF" $M/$NAME
