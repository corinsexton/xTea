#!/bin/bash
#SBATCH -J tdev_wgs
#SBATCH -c 24
#SBATCH -t 6:00:00
#SBATCH --mem=64G
#SBATCH -p priopark
#SBATCH --account=park_contrib
#SBATCH -o /n/scratch/users/c/cos689/xtea2_testdata/gate_transduction_evidence_wgs_%j.out
#SBATCH --mail-type=NONE
# =============================================================================
# The transduction decision (xtea/transduction_evidence.py) at WGS scale. ONE arm, against
# the FROZEN BASELINE-0 -- not two arms.
#
# WHY ONE ARM IS THE RIGHT SHAPE HERE, when the barrier gates all ran three. Those gates
# compared two topologies, so they needed a same-node control to compare timings against.
# This change is a pure refactor whose *only* correctness claim is "no artifact moves", and
# the chr1 gate (`gate_transduction_evidence_chr1.sbatch`) already establishes that against
# the parent commit over ~98 byte-identical files. What chr1 cannot supply is the
# POPULATION: it yields 17 calls and no created records at all, where WGS yields
# 188/1223/82 (merged 1493) including the 52 records that exist only through site creation.
# BASELINE-0 is frozen with sha256s and is the canonical statement of those numbers, so one
# arm against it answers the population question at half the compute. A second arm would
# only re-measure what chr1 measured, an hour more expensively.
#
# WHAT MUST BE TRUE:
#   1. every frozen BASELINE-0 VCF is byte-identical. This is the whole correctness claim.
#   2. the ledger exists once per FLANKED family (L1, SVA -- Alu has no transduction flank,
#      so the step never runs for it), is non-empty, and every reason in it is in the
#      module's own vocabulary with no unknown phase. A refactor that produced an empty
#      ledger would pass check 1 in silence, because nothing reads the ledger.
#   3. the call counts are 188 / 1223 / 82 / 1493.
#   4. `downstream_parallel` is in the band the last WGS gates measured. NOT a tight
#      assertion -- job 49599838's node was ~10% slower than job 49550925's (OFF 511.75s
#      against 462.7s), so absolutes are not comparable between jobs. The bound here is
#      loose on purpose: it exists to catch the decision loop having become structurally
#      expensive (an evidence object per locus over ~40K loci per family), not to measure
#      a few seconds of jitter.
#
# Frozen worktree, not the live checkout -- run_xtea spawns a fresh x_TEA_main per step, so
# a job sharing the checkout picks up mid-run edits.
#
# USAGE  sbatch [--dependency=afterok:<chr1 job>] gate_transduction_evidence_wgs.sbatch <worktree>
# =============================================================================
set -euo pipefail

SCRATCH=/n/scratch/users/c/cos689/xtea2_testdata
WT="${1:?usage: sbatch gate_transduction_evidence_wgs.sbatch <worktree> [types]}"
TYPES="${2:-L1,Alu,SVA}"
OUT="${XTEA2_OUT:-$SCRATCH/gate_transduction_evidence_wgs_${SLURM_JOB_ID:-manual}}"
BASE0="${XTEA2_BASELINE0:-$SCRATCH/baseline0_hg002_L1-Alu-SVA}"
TESTDATA=/n/data1/hms/dbmi/park/corinne/claude_sandbox/test_data/HG002
BAM="$TESTDATA/alignment/HG002.novaseq.pcr-free.30x.dedup.grch38.bam"
SAMPLE=HG002.wgs
REF=/n/data1/hms/dbmi/park-smaht_dac/ref/GRCh38_no_alt/hg38_no_alt.fa
CORES="${SLURM_CPUS_PER_TASK:-24}"
MAX_DOWNSTREAM_S="${XTEA2_MAX_DOWNSTREAM_S:-700}"

source /n/app/conda/miniforge3/24.11.3-0/etc/profile.d/conda.sh
conda activate xtea2
export PATH=/home/cos689/software/bin:$PATH
export XTEA_BASIC_INFO_SEED=42   # pin the insert-size draw, as BASELINE-0 did
test "$PYTHONNOUSERSITE" = 1     # the env's guard against ~/.local's onnxruntime

[ -d "$WT" ]                    || { echo "FATAL: worktree $WT missing"; exit 1; }
[ -x "$WT/src/csrc/xtea_scan" ] || { echo "FATAL: $WT/src/csrc/xtea_scan not built"; exit 1; }
[ -f "$BAM" ]                   || { echo "FATAL: $BAM missing"; exit 1; }
[ -f "$BAM.bai" ]               || { echo "FATAL: $BAM.bai missing"; exit 1; }
[ -z "$(git -C "$WT" status --porcelain)" ] || { echo "FATAL: worktree is dirty"; exit 1; }
[ -d "$BASE0/BASELINE-0" ]      || { echo "FATAL: $BASE0/BASELINE-0 missing"; exit 1; }
# The code under test has to be present BEFORE an hour of compute, not after.
[ -f "$WT/xtea/transduction_evidence.py" ] || {
    echo "FATAL: worktree has no transduction_evidence.py; wrong commit"; exit 1; }
grep -q 'transduction_evidence.select_clip_sites' "$WT/xtea/x_transduction.py" || {
    echo "FATAL: worktree's step is not wired to the gates; wrong commit"; exit 1; }

mkdir -p "$OUT"
cd "$WT"

echo "=============================================================="
echo "TRANSDUCTION EVIDENCE GATE (WGS)  --  types=$TYPES"
echo "  job      : ${SLURM_JOB_ID:-<none>} on $(hostname)  cores=$CORES"
echo "  commit   : $(git rev-parse --short HEAD) $(git log -1 --format=%s | cut -c1-60)"
echo "  worktree : $WT"
echo "  baseline : $BASE0"
echo "  out      : $OUT"
echo "  started  : $(date '+%F %T')"
echo "=============================================================="

t0=$(date +%s)
python xtea/run_xtea.py \
  -i "$BAM" -s "$SAMPLE" -o "$OUT" \
  --reference "$REF" \
  --rep-lib-dir /home/cos689/xTea/rep_lib_annotation \
  --gene-gff3 /home/cos689/data1/corinne/ref/gencode.v45.basic.annotation.gff3 \
  --blacklist /home/cos689/xTea/rep_lib_annotation/blacklist/hg38/centromere.bed \
  --types "$TYPES" \
  --scan-bin "$WT/src/csrc/xtea_scan" \
  --bwa /home/cos689/software/bin/bwa \
  --scan-min-clip 2 \
  --mate-in-rep-families alu,sva \
  --cns-emit-scan \
  --td-gather-scan \
  --parallel-downstream \
  --keep-intermediates \
  -n "$CORES" \
  > "$OUT/driver.log" 2>&1 || {
      echo "FAIL: the run exited $?"; tail -80 "$OUT/driver.log"; exit 1; }
wall=$(( $(date +%s) - t0 ))
printf '  run wall %02d:%02d:%02d\n' $((wall/3600)) $((wall%3600/60)) $((wall%60))

RC=0
B="$OUT/$SAMPLE"

echo
echo "--- check 1: every frozen BASELINE-0 VCF is byte-identical ---"
n_same=0; n_diff=0
for f in "$BASE0"/BASELINE-0/*.vcf; do
    b=$(find "$B" -name "$(basename "$f")" -type f | head -1)
    if [ -z "$b" ]; then echo "  MISSING: $(basename "$f")"; n_diff=$((n_diff+1)); continue; fi
    if cmp -s "$f" "$b"; then
        n_same=$((n_same+1)); echo "  $(basename "$f"): IDENTICAL"
    else
        echo "  $(basename "$f"): DIFFERS"; n_diff=$((n_diff+1))
        diff <(grep -v '^#' "$f" | cut -f1,2) <(grep -v '^#' "$b" | cut -f1,2) | head -20 || true
    fi
done
echo "  identical: $n_same    differing: $n_diff"
[ "$n_diff" -eq 0 ] && [ "$n_same" -eq 4 ] || {
    echo "  FAIL: the refactor moved a call at WGS scale"; RC=1; }

echo
echo "--- check 2: the ledgers, at WGS population ---"
mapfile -t VOCAB < <(python - "$WT" <<'EOF'
import sys, os
sys.path.insert(0, os.path.join(sys.argv[1], "xtea"))
import transduction_evidence as te
for r in sorted(te.PHASE_OF_REASON):
    print(r)
EOF
)
echo "  vocabulary (${#VOCAB[@]} reasons)"
[ "${#VOCAB[@]}" -ge 14 ] || { echo "    FATAL: could not read the reason vocabulary"; RC=1; }
n_led=0
while IFS= read -r led; do
    n_led=$((n_led+1))
    rows=$(( $(wc -l < "$led") - 1 ))
    fam=$(basename "$(dirname "$led")")
    echo "  $fam: $rows rows, $(wc -c < "$led") bytes"
    [ "$(head -1 "$led")" = "#chrom	pos	phase	reason	detail" ] || {
        echo "    FATAL: unexpected header"; RC=1; }
    [ "$rows" -gt 0 ] || { echo "    FATAL: empty ledger"; RC=1; }
    bad=$(awk -F'\t' 'NR>1 {print $4}' "$led" | sort -u | grep -vxF -f <(printf '%s\n' "${VOCAB[@]}") || true)
    [ -z "$bad" ] || { echo "    FATAL: reasons outside the vocabulary: $bad"; RC=1; }
    n_q=$(awk -F'\t' 'NR>1 && $3=="?"' "$led" | wc -l)
    [ "$n_q" -eq 0 ] || { echo "    FATAL: $n_q rows with an unknown phase"; RC=1; }
    echo "    by phase : $(awk -F'\t' 'NR>1 {c[$3]++} END {for (p in c) printf "%s=%d ", p, c[p]}' "$led")"
    echo "    by reason: $(awk -F'\t' 'NR>1 {c[$4]++} END {for (r in c) printf "%s=%d ", r, c[r]}' "$led")"
done < <(find "$B" -name '*.decision_ledger.tsv' | sort)
echo "  ledgers found: $n_led (expect 2 -- the flanked families)"
[ "$n_led" -eq 2 ] || RC=1
grep -ah 'TD-LEDGER' "$B"/*.steps.log 2>/dev/null || echo "  (no [TD-LEDGER] stderr lines found)"

echo
echo "--- check 3: call counts must be 188 / 1223 / 82 / 1493 ---"
l1=$(grep -vc '^#' "$B/merged/HG002.novaseq.pcr-free.30x.dedup.grch38_LINE1.vcf" 2>/dev/null || echo 0)
alu=$(grep -vc '^#' "$B/merged/HG002.novaseq.pcr-free.30x.dedup.grch38_ALU.vcf" 2>/dev/null || echo 0)
sva=$(grep -vc '^#' "$B/merged/HG002.novaseq.pcr-free.30x.dedup.grch38_SVA.vcf" 2>/dev/null || echo 0)
mrg=$(grep -vc '^#' "$B/$SAMPLE.vcf" 2>/dev/null || echo 0)
echo "  L1=$l1 Alu=$alu SVA=$sva merged=$mrg"
[ "$l1" -eq 188 ] && [ "$alu" -eq 1223 ] && [ "$sva" -eq 82 ] && [ "$mrg" -eq 1493 ] || {
    echo "  FAIL: call counts moved"; RC=1; }

echo
echo "--- check 4: downstream_parallel, loosely bounded ---"
# Not a measurement of the change (one arm, one node, no control): a bound that catches the
# decision loop having become structurally expensive. The last WGS gates put this at 462.7s
# (job 49550925) and 511.75s (job 49599838, a ~10% slower node) with no barrier flags.
dp=$(grep -a 'downstream_parallel' "$OUT/driver.log" | tail -3 || true)
echo "$dp" | sed 's/^/  /'
secs=$(grep -aoE 'downstream_parallel[^0-9]*([0-9]+\.[0-9]+)s' "$OUT/driver.log" | tail -1 |
       grep -oE '[0-9]+\.[0-9]+' | tail -1 || true)
if [ -n "$secs" ]; then
    echo "  downstream_parallel = ${secs}s (bound: <= ${MAX_DOWNSTREAM_S}s)"
    awk -v s="$secs" -v m="$MAX_DOWNSTREAM_S" 'BEGIN {exit !(s <= m)}' || {
        echo "  FAIL: downstream_parallel above the loose bound"; RC=1; }
else
    echo "  (could not parse downstream_parallel from the driver log -- reporting, not failing)"
fi

echo
echo "=============================================================="
if [ "$RC" -eq 0 ]; then
    echo "TRANSDUCTION EVIDENCE GATE (WGS): PASS"
else
    echo "TRANSDUCTION EVIDENCE GATE (WGS): FAIL"
fi
echo "  finished: $(date '+%F %T')"
echo "=============================================================="
exit "$RC"
