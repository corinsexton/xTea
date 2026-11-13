#!/bin/bash
#SBATCH -J cv2_recount
#SBATCH -A park_contrib
#SBATCH -p priopark
#SBATCH -c 16
#SBATCH --mem=48G
#SBATCH -t 1:00:00
#SBATCH -o /n/scratch/users/c/cos689/xtea2_review/counting_v2/mosaic/logs/%x_%A_%a.out
set -euo pipefail
source /n/app/conda/miniforge3/24.11.3-0/etc/profile.d/conda.sh; conda activate xtea2; export PYTHONNOUSERSITE=1
M=/n/scratch/users/c/cos689/xtea2_review/counting_v2/mosaic
S=/n/data1/hms/dbmi/park/julia/SMaHT/benchmarking/seqdata/SMaHT/HapMapMix_Illumina/subsampled
case $SLURM_ARRAY_TASK_ID in
 1) N=BCM_300X; B=$S/BCM_Illumina_300X.bam; IS=659;;
 2) N=UW_100X; B=$S/UW_Illumina_100X.bam; IS=430;;
 3) N=WashU_100X; B=$S/WashU_Illumina_100X.bam; IS=619;;
esac
t0=$(date +%s)
python $M/recount_mosaic_cv2.py --name $N --bam $B \
  --ref /n/data1/hms/dbmi/park-smaht_dac/ref/GRCh38_no_alt/hg38_no_alt.fa --is $IS --rlen 151 -n 16 --out $M/inputs 2> >(grep -v hts_idx_load >&2)
echo "WALL $(( $(date +%s)-t0 ))s"
