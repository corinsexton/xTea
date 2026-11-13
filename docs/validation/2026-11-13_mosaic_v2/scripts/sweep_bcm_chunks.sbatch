#!/bin/bash
#SBATCH -J mv2_sweepbcm
#SBATCH -A park
#SBATCH -p short
#SBATCH -c 4
#SBATCH --mem=32G
#SBATCH -t 2:00:00
#SBATCH -o /n/scratch/users/c/cos689/xtea2_review/mosaic_v2/logs/%x_%A_%a.out
# BCM_300X fine disc_upper sweep in 4 chunks (one per task).
set -euo pipefail
source /n/app/conda/miniforge3/24.11.3-0/etc/profile.d/conda.sh; conda activate xtea2; export PYTHONNOUSERSITE=1
export PATH=/home/cos689/software/bin:$PATH
C=/n/scratch/users/c/cos689/xtea2_review/counting_v2/mosaic; N=BCM_300X; D=$C/inputs
VALS=$(python -c "i=$SLURM_ARRAY_TASK_ID; print(','.join('%.2f'%(0.20+0.01*j) for j in range(51) if j%4==i-1))")
python $C/sweep.py --run /n/scratch/users/c/cos689/xtea2_testdata/vafmix/runs/$N --name $N --mode mosaic --param disc \
  --values=$VALS --mosaic-inputs L1=$D/${N}_L1.input_cv2.txt,Alu=$D/${N}_Alu.input_cv2.txt,SVA=$D/${N}_SVA.input_cv2.txt \
  --out /n/scratch/users/c/cos689/xtea2_review/mosaic_v2/sweep_runs/${N}_c$SLURM_ARRAY_TASK_ID --tsv /n/scratch/users/c/cos689/xtea2_review/mosaic_v2/sweeps/${N}_fine_c$SLURM_ARRAY_TASK_ID.tsv -n 4
