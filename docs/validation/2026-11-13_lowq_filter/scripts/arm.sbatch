#!/bin/bash
#SBATCH -J lqarm
#SBATCH -A park_contrib
#SBATCH -p priopark
#SBATCH -c 20
#SBATCH --mem=48G
#SBATCH -t 1:00:00
#SBATCH -o /n/scratch/users/c/cos689/xtea2_review/lowq_filter/logs/arm_%x_%A_%a.out
# Downstream-only rerun of a runs_b15 dataset with the cns low-MAPQ clip counts replaced by one
# variant column of census/<name>.counts.tsv (scan_pass + cns_emit adopted via --resume).
# ARM=<arm dir name>, COL=<counts column: scanner|est|exact|est_alt|exact_alt|zero>.
set -euo pipefail
G=/n/scratch/users/c/cos689/xtea2_testdata/gtval
L=/n/scratch/users/c/cos689/xtea2_review/lowq_filter
WT=/n/scratch/users/c/cos689/xtea2_review/lowq_filter/code_2af82e4
: "${ARM:?}"; : "${COL:?}"
read -r ID NAME BAM REF < <(awk -v i="$SLURM_ARRAY_TASK_ID" '$1==i' $G/manifest.tsv)
OUT=$L/arms/$ARM/$NAME
echo "arm $ARM col $COL  $NAME  code 2af82e4 frozen copy  $(date '+%F %T')"
rm -rf "$OUT"; mkdir -p "$(dirname $OUT)"
cp -r $G/runs_b15/$NAME "$OUT"
S=$OUT/$NAME
# rewrite absolute paths? scan/cns outputs are referenced by path in the run; re-adopted below
rm -rf $S/.resume $S/$NAME.vcf $S/merged
for f in L1 Alu SVA; do
  awk -v F=$f -v C=$COL 'NR==1{for(i=1;i<=NF;i++) h[$i]=i; next} $1==F{v=(C=="zero")?0:$(h[C]); print $2, $3, v}' \
    $L/census/$NAME.counts.tsv > $S/$f/tmp/cns/low_mapq_clip_cnt.txt
done
source /n/app/conda/miniforge3/24.11.3-0/etc/profile.d/conda.sh
conda activate xtea2
export PATH=/home/cos689/software/bin:$PATH
export PYTHONNOUSERSITE=1
unset XTEA_BASIC_INFO_SEED XTEA_TD_DECOY_CONCAT XTEA_TD_DECOY_STRUCTURAL || true
t0=$(date +%s)
cd "$WT"
python xtea/run_xtea.py \
  -i "$BAM" -s "$NAME" -o "$OUT" \
  --reference "$REF" \
  --rep-lib-dir /home/cos689/xTea/rep_lib_annotation \
  --gene-gff3 /home/cos689/data1/corinne/ref/gencode.v45.basic.annotation.gff3 \
  --blacklist /home/cos689/xTea/rep_lib_annotation/blacklist/hg38/centromere.bed \
  --types L1,Alu,SVA \
  --scan-bin "$WT/src/csrc/xtea_scan" \
  --bwa /home/cos689/software/bin/bwa \
  --scan-min-clip 2 \
  --mate-in-rep-families alu,sva \
  --cns-emit-scan \
  --td-gather-scan \
  --parallel-downstream \
  --keep-intermediates \
  --ignore-rg \
  --resume --adopt-existing \
  -n "${SLURM_CPUS_PER_TASK:-20}" \
  > "$OUT/driver.log" 2>&1 || { rc=$?; echo "FAIL: exit $rc"; tail -60 "$OUT/driver.log"; exit $rc; }
grep -c "RESUMED\|\[resume\]" $OUT/driver.log || true
grep "\[resume\]" $OUT/driver.log | head -5
printf '  merged %s\n' "$(grep -vc '^#' "$S/$NAME.vcf" 2>/dev/null || echo NA)"
echo "WALL $(( $(date +%s) - t0 ))s"
