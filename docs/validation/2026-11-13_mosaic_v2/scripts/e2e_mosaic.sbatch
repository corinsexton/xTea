#!/bin/bash
#SBATCH -J mv2_e2e
#SBATCH -A park
#SBATCH -p short
#SBATCH -c 20
#SBATCH --mem=64G
#SBATCH -t 10:00:00
#SBATCH -o /n/scratch/users/c/cos689/xtea2_review/mosaic_v2/e2e/logs/%x_%j.out
# Copy of xtea2_testdata/vafmix/mosaic.sbatch for the mosaic_v2 end-to-end check: xTea2 (vaf_site_v2 worktree), --mosaic, on SMaHT HapMap mixture BAMs (VAF calibration). NAME, BAM from env.
# Task id -> row of gtval/manifest.tsv (id, name, input, reference). Command line identical to
# xtea2_testdata/v2.0.0_hg002.sbatch except --ignore-rg (the GIAB 2x250 BAM has no @RG).
set -euo pipefail
G=/n/scratch/users/c/cos689/xtea2_testdata/gtval
NAME=${NAME:?}; BAM=${BAM:?}; REF=/n/data1/hms/dbmi/park-smaht_dac/ref/GRCh38_no_alt/hg38_no_alt.fa; ID=$NAME
WT=${WT:-/n/scratch/users/c/cos689/xtea2_wt/vaf_site_v2}
OUT=/n/scratch/users/c/cos689/xtea2_review/mosaic_v2/e2e/runs/$NAME
REPLIB=/home/cos689/xTea/rep_lib_annotation
SAMPLE=$NAME
SCAN_BIN="$WT/src/csrc/xtea_scan"
CORES="${SLURM_CPUS_PER_TASK:-24}"
echo "task $ID $NAME  in=$BAM  ref=$REF"
echo "node $(hostname)  code $WT ($(git -C "$WT" rev-parse --short HEAD))  started $(date '+%F %T')"
[ -x "$SCAN_BIN" ] || { echo "FATAL: $SCAN_BIN not built"; exit 1; }
[ -z "$(git -C "$WT" status --porcelain)" ] || { echo "FATAL: worktree is dirty"; exit 1; }
unset XTEA_BASIC_INFO_SEED XTEA_TD_DECOY_CONCAT XTEA_TD_DECOY_STRUCTURAL || true
source /n/app/conda/miniforge3/24.11.3-0/etc/profile.d/conda.sh
conda activate xtea2
export PATH=/home/cos689/software/bin:$PATH
[ "${PYTHONNOUSERSITE:-}" = 1 ] || { echo "FATAL: PYTHONNOUSERSITE is not 1"; exit 1; }
rm -rf "$OUT"; mkdir -p "$OUT"
t0=$(date +%s)
cd "$WT"
python xtea/run_xtea.py \
  -i "$BAM" -s "$SAMPLE" -o "$OUT" \
  --reference "$REF" \
  --rep-lib-dir "$REPLIB" \
  --gene-gff3 /home/cos689/data1/corinne/ref/gencode.v45.basic.annotation.gff3 \
  --blacklist "$REPLIB/blacklist/hg38/centromere.bed" \
  --types L1,Alu,SVA \
  --scan-bin "$SCAN_BIN" \
  --bwa /home/cos689/software/bin/bwa \
  --scan-min-clip 2 \
  --mate-in-rep-families alu,sva \
  --cns-emit-scan \
  --td-gather-scan \
  --parallel-downstream \
  --keep-intermediates \
  --ignore-rg \
  --mosaic \
  -n "$CORES" \
  > "$OUT/driver.log" 2>&1 || { rc=$?; echo "FAIL: exit $rc"; tail -80 "$OUT/driver.log"; exit $rc; }
wall=$(( $(date +%s) - t0 ))
for pair in L1:LINE1 Alu:ALU SVA:SVA; do
  t="${pair%%:*}"; alt="${pair#*:}"
  f=$(ls "$OUT/$SAMPLE/merged"/*_${alt}.vcf 2>/dev/null | head -1)
  printf '  %-4s %s\n' "$t" "$([ -n "$f" ] && grep -vc '^#' "$f" || echo NA)"
done
printf '  merged %s\n' "$(grep -vc '^#' "$OUT/$SAMPLE/$SAMPLE.vcf" 2>/dev/null || echo NA)"
printf 'WALL %02d:%02d:%02d  finished %s\n' $((wall/3600)) $((wall%3600/60)) $((wall%60)) "$(date '+%F %T')"
