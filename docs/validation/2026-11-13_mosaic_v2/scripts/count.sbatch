#!/bin/bash
#SBATCH -J mv2_count
#SBATCH -A park
#SBATCH -p park,short
#SBATCH -c 8
#SBATCH --mem=24G
#SBATCH -t 1:00:00
#SBATCH -o /n/scratch/users/c/cos689/xtea2_review/mosaic_v2/logs/count_%A_%a.out
set -euo pipefail
M=/n/scratch/users/c/cos689/xtea2_review/mosaic_v2
V=/n/scratch/users/c/cos689/xtea2_testdata/vafsite
read -r ID NAME BAM REF SITES < <(awk -v i="$SLURM_ARRAY_TASK_ID" '$1==i' $M/counts_manifest.tsv)
source /n/app/conda/miniforge3/24.11.3-0/etc/profile.d/conda.sh; conda activate xtea2; export PYTHONNOUSERSITE=1
export TMPDIR=/n/scratch/users/c/cos689/tmp; mkdir -p $TMPDIR
export XTEA_DIR=/n/scratch/users/c/cos689/xtea2_wt/vaf_site_v2/xtea
echo "$NAME $BAM $(date)"; t0=$(date +%s)
python $M/scripts/count_at_sites.py --bam "$BAM" --ref "$REF" --sites $M/sites/$SITES.tsv --out $M/counts/$NAME.tsv -n 8 --isize-from $V/counts/$NAME.tsv.isize 2> >(grep -v hts_idx_load >&2)
echo "WALL $(( $(date +%s)-t0 ))s"
