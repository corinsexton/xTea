#!/bin/bash
#SBATCH -J cv2
#SBATCH -A park_contrib
#SBATCH -p priopark
#SBATCH -c 24
#SBATCH --mem=64G
#SBATCH -t 1:30:00
#SBATCH -o /n/scratch/users/c/cos689/xtea2_review/counting_v2/validation/logs/%x_%A_%a.out
# Copy of gtval/xtea_final.sbatch for the counting-v2 (right-clip map_pos + supplementary AF clips) revalidation.
# Differences: WT and RUNS come from the environment (--export), output goes to gtval/$RUNS,
# and since the worktree carries the uncommitted fix, the "clean worktree" check is replaced
# by a check that the worktree state (diff vs HEAD + untracked files) hashes to EXPECT_STATE.
set -euo pipefail
G=/n/scratch/users/c/cos689/xtea2_testdata/gtval
read -r ID NAME BAM REF < <(awk -v i="$SLURM_ARRAY_TASK_ID" '$1==i' $G/manifest.tsv)
: "${WT:?set WT}"; : "${RUNS:?set RUNS}"; : "${EXPECT_STATE:?set EXPECT_STATE}"
case "$RUNS" in runs_final|runs_v2.0.0|runs_vaf|runs_bug2*|runs_ctrl_ca5b843) echo "FATAL: refusing to write $RUNS"; exit 1;; esac
OUT=$G/$RUNS/$NAME
REPLIB=/home/cos689/xTea/rep_lib_annotation
SAMPLE=$NAME
SCAN_BIN="$WT/src/csrc/xtea_scan"
CORES="${SLURM_CPUS_PER_TASK:-24}"
state() { ( cd "$WT"; git diff HEAD; git ls-files --others --exclude-standard | sort | xargs -r sha256sum ) | sha256sum | cut -c1-16; }
echo "task $ID $NAME  in=$BAM  ref=$REF  runs=$RUNS"
echo "node $(hostname)  code $WT ($(git -C "$WT" rev-parse --short HEAD), state $(state))  started $(date '+%F %T')"
[ -x "$SCAN_BIN" ] || { echo "FATAL: $SCAN_BIN not built"; exit 1; }
[ "$(state)" = "$EXPECT_STATE" ] || { echo "FATAL: worktree state $(state) != $EXPECT_STATE"; exit 1; }
git -C "$WT" status --porcelain
unset XTEA_BASIC_INFO_SEED XTEA_TD_DECOY_CONCAT XTEA_TD_DECOY_STRUCTURAL || true
source /n/app/conda/miniforge3/24.11.3-0/etc/profile.d/conda.sh
conda activate xtea2
export PATH=/home/cos689/software/bin:$PATH
export PYTHONNOUSERSITE=1
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
  -n "$CORES" \
  > "$OUT/driver.log" 2>&1 || { rc=$?; echo "FAIL: exit $rc"; tail -80 "$OUT/driver.log"; exit $rc; }
wall=$(( $(date +%s) - t0 ))
[ "$(state)" = "$EXPECT_STATE" ] || echo "WARNING: worktree state changed during the run: $(state)"
for pair in L1:LINE1 Alu:ALU SVA:SVA; do
  t="${pair%%:*}"; alt="${pair#*:}"
  f=$(ls "$OUT/$SAMPLE/merged"/*_${alt}.vcf 2>/dev/null | head -1)
  printf '  %-4s %s\n' "$t" "$([ -n "$f" ] && grep -vc '^#' "$f" || echo NA)"
done
printf '  merged %s\n' "$(grep -vc '^#' "$OUT/$SAMPLE/$SAMPLE.vcf" 2>/dev/null || echo NA)"
printf 'WALL %02d:%02d:%02d  finished %s\n' $((wall/3600)) $((wall%3600/60)) $((wall%60)) "$(date '+%F %T')"
