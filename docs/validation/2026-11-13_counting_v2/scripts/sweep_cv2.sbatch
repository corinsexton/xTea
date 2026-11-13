#!/bin/bash
#SBATCH -J cv2_mosaic
#SBATCH -A park_contrib
#SBATCH -p priopark
#SBATCH -c 8
#SBATCH --mem=32G
#SBATCH -t 1:00:00
#SBATCH -o /n/scratch/users/c/cos689/xtea2_review/counting_v2/mosaic/logs/%x_%A_%a.out
# Mosaic HapMap-mixture check at MOSAIC_AF_BOUNDS disc_upper 0.35, counting-v2 code copy.
# tasks 1-3: counting-v2 inputs (B14 inputs + per-record cv2-minus-0553edb recount delta)
# tasks 4-6: B14 inputs unchanged (reproduction check of the code copy vs the B14 sweep)
set -euo pipefail
source /n/app/conda/miniforge3/24.11.3-0/etc/profile.d/conda.sh; conda activate xtea2; export PYTHONNOUSERSITE=1
export PATH=/home/cos689/software/bin:$PATH
M=/n/scratch/users/c/cos689/xtea2_review/counting_v2/mosaic
T=/n/scratch/users/c/cos689/xtea2_review/bug2/validation/tune
NAMES=(x BCM_300X UW_100X WashU_100X BCM_300X UW_100X WashU_100X)
N=${NAMES[$SLURM_ARRAY_TASK_ID]}
if [ $SLURM_ARRAY_TASK_ID -le 3 ]; then K=cv2; D=$M/inputs; S=input_cv2; else K=b14; D=$T/mosaic_inputs; S=input_new; fi
t0=$(date +%s)
python $M/sweep.py --run /n/scratch/users/c/cos689/xtea2_testdata/vafmix/runs/$N --name $N --mode mosaic --param disc \
  --values=0.35 --mosaic-inputs L1=$D/${N}_L1.$S.txt,Alu=$D/${N}_Alu.$S.txt,SVA=$D/${N}_SVA.$S.txt \
  --out $M/runs/$N/$K --tsv $M/sweeps/${N}_$K.tsv -n 8
echo "WALL $(( $(date +%s)-t0 ))s"
