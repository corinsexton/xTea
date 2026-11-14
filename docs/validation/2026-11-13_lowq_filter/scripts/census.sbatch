#!/bin/bash
#SBATCH -J lqcensus
#SBATCH -A park_contrib
#SBATCH -p priopark
#SBATCH -c 1
#SBATCH --mem=8G
#SBATCH -t 1:00:00
#SBATCH -o /n/scratch/users/c/cos689/xtea2_review/lowq_filter/logs/census_%A_%a.out
set -euo pipefail
G=/n/scratch/users/c/cos689/xtea2_testdata/gtval
read -r ID NAME BAM REF < <(awk -v i="$SLURM_ARRAY_TASK_ID" '$1==i' $G/manifest.tsv)
RL=$(grep -o 'cns-read-length [0-9]*' $G/runs_b15/$NAME/driver.log | head -1 | awk '{print $2}')
source /n/app/conda/miniforge3/24.11.3-0/etc/profile.d/conda.sh; conda activate xtea2
export PYTHONNOUSERSITE=1
python /n/scratch/users/c/cos689/xtea2_review/lowq_filter/scripts/lowq_reads.py $G/runs_b15/$NAME/$NAME "$BAM" "$REF" $RL /n/scratch/users/c/cos689/xtea2_review/lowq_filter/census/$NAME.tsv.gz
echo done
