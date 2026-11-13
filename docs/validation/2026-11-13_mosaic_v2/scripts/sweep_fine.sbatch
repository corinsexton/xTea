#!/bin/bash
#SBATCH -J mv2_sweep
#SBATCH -A park
#SBATCH -p short
#SBATCH -c 4
#SBATCH --mem=24G
#SBATCH -t 2:00:00
#SBATCH -o /n/scratch/users/c/cos689/xtea2_review/mosaic_v2/logs/%x_%A_%a.out
# Fine disc_upper sweep 0.20-0.70 on the B15 mosaic inputs (counting_v2 code copy), for the
# per-sample ceiling fit.
set -euo pipefail
source /n/app/conda/miniforge3/24.11.3-0/etc/profile.d/conda.sh; conda activate xtea2; export PYTHONNOUSERSITE=1
export PATH=/home/cos689/software/bin:$PATH
C=/n/scratch/users/c/cos689/xtea2_review/counting_v2/mosaic
NAMES=(x BCM_300X UW_100X WashU_100X)
N=${NAMES[$SLURM_ARRAY_TASK_ID]}
D=$C/inputs
t0=$(date +%s)
python $C/sweep.py --run /n/scratch/users/c/cos689/xtea2_testdata/vafmix/runs/$N --name $N --mode mosaic --param disc \
  --values=0.20,0.21,0.22,0.23,0.24,0.25,0.26,0.27,0.28,0.29,0.30,0.31,0.32,0.33,0.34,0.35,0.36,0.37,0.38,0.39,0.40,0.41,0.42,0.43,0.44,0.45,0.46,0.47,0.48,0.49,0.50,0.51,0.52,0.53,0.54,0.55,0.56,0.57,0.58,0.59,0.60,0.61,0.62,0.63,0.64,0.65,0.66,0.67,0.68,0.69,0.70 --mosaic-inputs L1=$D/${N}_L1.input_cv2.txt,Alu=$D/${N}_Alu.input_cv2.txt,SVA=$D/${N}_SVA.input_cv2.txt \
  --out /n/scratch/users/c/cos689/xtea2_review/mosaic_v2/sweep_runs/$N --tsv /n/scratch/users/c/cos689/xtea2_review/mosaic_v2/sweeps/${N}_fine.tsv -n 4
echo "WALL $(( $(date +%s)-t0 ))s"
