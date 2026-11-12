#!/bin/bash
#SBATCH -J tune_sweep
#SBATCH -A park_contrib
#SBATCH -p priopark
#SBATCH -c 8
#SBATCH --mem=32G
#SBATCH -t 1:00:00
#SBATCH -o /n/scratch/users/c/cos689/xtea2_review/bug2/validation/tune/logs/%x_%A_%a.out
set -euo pipefail
source /n/app/conda/miniforge3/24.11.3-0/etc/profile.d/conda.sh; conda activate xtea2; export PYTHONNOUSERSITE=1
export PATH=/home/cos689/software/bin:$PATH
T=/n/scratch/users/c/cos689/xtea2_review/bug2/validation/tune
G=/n/scratch/users/c/cos689/xtea2_testdata/gtval
IFS=$'\t' read -r ID MODE NAME PARAM VALUES INP < <(awk -F'\t' -v i="$SLURM_ARRAY_TASK_ID" '$1==i' $T/sweep_tasks.tsv)
echo "task $ID $MODE $NAME $PARAM inputs=$INP node $(hostname) started $(date '+%F %T')"
t0=$(date +%s)
if [ "$MODE" = germline ]; then
  python $T/sweep.py --run $G/runs_bug2/$NAME --name $NAME --mode germline --param $PARAM --values="$VALUES" \
    --out $T/runs/$NAME/$PARAM --tsv $T/sweeps/${NAME}_${PARAM}.tsv -n 8
else
  M=/n/scratch/users/c/cos689/xtea2_testdata/vafmix/runs/$NAME
  EXTRA=""
  [ "$INP" != old ] && EXTRA="--mosaic-inputs L1=$T/mosaic_inputs/${NAME}_L1.input_new.txt,Alu=$T/mosaic_inputs/${NAME}_Alu.input_new.txt,SVA=$T/mosaic_inputs/${NAME}_SVA.input_new.txt"
  python $T/sweep.py --run $M --name $NAME --mode mosaic --param disc --values="$VALUES" $EXTRA \
    --out $T/runs/$NAME/disc_$INP --tsv $T/sweeps/${NAME}_disc_${INP}.tsv -n 8
fi
echo "WALL $(( $(date +%s)-t0 ))s"
