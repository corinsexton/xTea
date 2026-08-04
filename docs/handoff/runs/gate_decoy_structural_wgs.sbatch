#!/bin/bash
#SBATCH -J decoy_wgs
#SBATCH -c 24
#SBATCH -t 6:00:00
#SBATCH --mem=64G
#SBATCH -p priopark
#SBATCH --account=park_contrib
#SBATCH -o /n/scratch/users/c/cos689/xtea2_testdata/gate_decoy_structural_wgs_%j.out
#SBATCH --mail-type=NONE
# =============================================================================
# The SVA decoy fix at WGS scale. ONE arm, WITH `--td-decoy-structural` ON, against the frozen
# BASELINE-0 -- and this is the leg the fix exists for, because chr1's gate B has 0 SVA
# survivors and therefore no population to move.
#
# WHY ONE ARM, AND WHY THE FROZEN BASELINE IS THE OTHER ONE. The containment claim ("with the
# flag off nothing moves") is `gate_decoy_structural_chr1.sbatch`'s AFTER_OFF arm, over the whole
# output tree. BASELINE-0 is frozen with sha256s and was reproduced byte-for-byte by the last WGS
# gate at the behavioural parent of this commit (job 49703322, 4/4 VCFs), so it IS the off-arm
# here at no compute. That makes check 1 a sharper statement than "some things moved":
#
#   **L1 and Alu must be byte-identical to BASELINE-0 WITH THE FIX ON, and SVA must not be.**
#
# L1 flank ids collide with no decoy tag (`TD_DECOY_LINE` is "LINE1", not "L1" -- the probe
# measures 0 collisions in both L1 libraries), and Alu never runs the transduction step. So a
# moved L1 or Alu file falsifies the classification's premise, and an unmoved SVA file means the
# flag did not reach the step. Both directions are failures, which is what makes it a gate.
#
# WHAT IS BEING MEASURED, and check 5 is the payoff:
#   1. L1 + Alu identical to BASELINE-0; SVA and the merged VCF differ. Stated as a prediction.
#   2. L1 = 188 and Alu = 1223 exactly; SVA and merged reported, not predicted.
#   3. the ledgers: L1's must be byte-identical to nothing in particular but must not change its
#      B_source population, SVA's is expected to move. Vocabulary-clean, no unknown phase.
#   4. `downstream_parallel` inside a loose bound. NOT a measurement -- one arm on one node has
#      no control. Prior no-flag WGS runs: 462.7s, 511.75s (a ~10% slower node), 463.3s, 485.8s.
#      Expect the fix to make the SVA step SLOWER, not faster: hits that used to be discarded
#      at the decoy branch now reach the source extraction and the dominance pick.
#   5. **THE PAYOFF.** `probe_td_source_classes.py`'s run half must report SVA reference-copy
#      picked sources > 0. It is 0 today and 0 on master and 0 on every genome build, out of 105
#      picked sources, against a 5,454-record library built to be exactly that source set. If
#      this is still 0, the fix does not do the one thing it was written to do.
#   6. the truth-set score, BASELINE-0 vs this arm, from `benchmark/benchmark_truth_set.py`.
#      REPORTED, NOT GATED, and read §7.1 before reading it as a verdict: the HG002 truth set
#      enumerates MEI insertion POSITIONS and its derived transduction denominator is 0 (L1) and
#      2 (SVA), so it cannot score transduction calls as transductions. What it can say is
#      whether SVA position sensitivity/precision (85.6% / 93.9% at BASELINE-0, 90 truth
#      records) got better or worse, and whether new SVA calls land on truth.
#   7. the FN attribution, rerun. Same status: information, not a gate.
#
# Frozen worktree, not the live checkout -- run_xtea spawns a fresh x_TEA_main per step, so a
# job sharing the checkout picks up mid-run edits.
#
# USAGE  sbatch [--dependency=afterok:<chr1 job>] gate_decoy_structural_wgs.sbatch <worktree>
# =============================================================================
set -euo pipefail

SCRATCH=/n/scratch/users/c/cos689/xtea2_testdata
WT="${1:?usage: sbatch gate_decoy_structural_wgs.sbatch <worktree> [types]}"
TYPES="${2:-L1,Alu,SVA}"
OUT="${XTEA2_OUT:-$SCRATCH/gate_decoy_structural_wgs_${SLURM_JOB_ID:-manual}}"
BASE0="${XTEA2_BASELINE0:-$SCRATCH/baseline0_hg002_L1-Alu-SVA}"
TESTDATA=/n/data1/hms/dbmi/park/corinne/claude_sandbox/test_data/HG002
BAM="$TESTDATA/alignment/HG002.novaseq.pcr-free.30x.dedup.grch38.bam"
SAMPLE=HG002.wgs
REF=/n/data1/hms/dbmi/park-smaht_dac/ref/GRCh38_no_alt/hg38_no_alt.fa
CORES="${SLURM_CPUS_PER_TASK:-24}"
MAX_DOWNSTREAM_S="${XTEA2_MAX_DOWNSTREAM_S:-900}"
# L1's ledger row count on job 49703322, which is the behavioural parent of this commit.
# Overridable so a future run with different inputs can restate it rather than lie.
EXPECT_LEDGER_L1="${XTEA2_EXPECT_LEDGER_L1:-21158}"
TRUTH_DIR="${XTEA2_TRUTH_DIR:-$TESTDATA/truth}"
REPLIB=/home/cos689/xTea/rep_lib_annotation
VCF_PREFIX=HG002.novaseq.pcr-free.30x.dedup.grch38

source /n/app/conda/miniforge3/24.11.3-0/etc/profile.d/conda.sh
conda activate xtea2
export PATH=/home/cos689/software/bin:$PATH
export XTEA_BASIC_INFO_SEED=42   # pin the insert-size draw, as BASELINE-0 did
test "$PYTHONNOUSERSITE" = 1     # the env's guard against ~/.local's onnxruntime
# The flag is passed on the command line below; an inherited env var would make it impossible to
# tell which channel armed it, and the driver banner is check 0's evidence either way.
unset XTEA_TD_DECOY_STRUCTURAL

[ -d "$WT" ]                    || { echo "FATAL: worktree $WT missing"; exit 1; }
[ -x "$WT/src/csrc/xtea_scan" ] || { echo "FATAL: $WT/src/csrc/xtea_scan not built"; exit 1; }
[ -f "$BAM" ]                   || { echo "FATAL: $BAM missing"; exit 1; }
[ -f "$BAM.bai" ]               || { echo "FATAL: $BAM.bai missing"; exit 1; }
[ -z "$(git -C "$WT" status --porcelain)" ] || { echo "FATAL: worktree is dirty"; exit 1; }
[ -d "$BASE0/BASELINE-0" ]      || { echo "FATAL: $BASE0/BASELINE-0 missing"; exit 1; }
# The code under test has to be present BEFORE an hour of compute, not after.
if ! grep -q 'def decoy_tags' "$WT/xtea/transduction_evidence.py"; then
    echo "FATAL: worktree has no decoy_tags; wrong commit"; exit 1
fi
if ! grep -q 'td-decoy-structural' "$WT/xtea/run_xtea.py"; then
    echo "FATAL: worktree has no --td-decoy-structural flag; wrong commit"; exit 1
fi

mkdir -p "$OUT"
cd "$WT"

echo "=============================================================="
echo "GATE: the SVA decoy fix, WGS, FLAG ON  --  types=$TYPES"
echo "  job      : ${SLURM_JOB_ID:-<none>} on $(hostname)  cores=$CORES"
echo "  commit   : $(git rev-parse --short HEAD) $(git log -1 --format=%s | cut -c1-60)"
echo "  worktree : $WT"
echo "  baseline : $BASE0  (the frozen OFF arm)"
echo "  out      : $OUT"
echo "  started  : $(date '+%F %T')"
echo "=============================================================="

t0=$(date +%s)
python xtea/run_xtea.py \
  -i "$BAM" -s "$SAMPLE" -o "$OUT" \
  --reference "$REF" \
  --rep-lib-dir "$REPLIB" \
  --gene-gff3 /home/cos689/data1/corinne/ref/gencode.v45.basic.annotation.gff3 \
  --blacklist /home/cos689/xTea/rep_lib_annotation/blacklist/hg38/centromere.bed \
  --types "$TYPES" \
  --scan-bin "$WT/src/csrc/xtea_scan" \
  --bwa /home/cos689/software/bin/bwa \
  --scan-min-clip 2 \
  --mate-in-rep-families alu,sva \
  --cns-emit-scan \
  --td-gather-scan \
  --td-decoy-structural \
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
echo "--- check 0: this arm really ran with the fix ON ---"
n=$(grep -ac '\[decoy\] the transduction decoy test is STRUCTURAL' "$OUT/driver.log" || true)
echo "  [decoy] banner: $n (expect 1)"
[ "$n" -eq 1 ] || { echo "  FATAL: the flag did not reach the driver -- this arm is a no-op"; RC=1; }

echo
echo "--- check 1: THE PREDICTION -- L1 and Alu identical to BASELINE-0, SVA not ---"
n_same=0; n_diff=0
for f in "$BASE0"/BASELINE-0/*.vcf; do
    bn=$(basename "$f")
    b=$(find "$B" -name "$bn" -type f | head -1)
    if [ -z "$b" ]; then echo "  MISSING: $bn"; n_diff=$((n_diff+1)); continue; fi
    if cmp -s "$f" "$b"; then s=IDENTICAL; n_same=$((n_same+1)); else s=DIFFERS; n_diff=$((n_diff+1)); fi
    echo "  $bn: $s"
    case "$bn" in
      *_LINE1.vcf|*_ALU.vcf)
        [ "$s" = IDENTICAL ] || {
            echo "    FAIL: the fix moved a family that collides with no decoy tag (L1) or"
            echo "    never runs the transduction step at all (Alu). Either the classification"
            echo "    is not what the probe measured, or something else in this arm moved."
            diff <(grep -v '^#' "$f" | cut -f1,2) <(grep -v '^#' "$b" | cut -f1,2) | head -20 || true
            RC=1; } ;;
      *_SVA.vcf)
        [ "$s" = DIFFERS ] || {
            echo "    FAIL: SVA is byte-identical to BASELINE-0, so the fix reached nothing."
            echo "    The whole claim is that 5,454 hg38 SVA flank records stop being read as"
            echo "    the consensus decoy; if no SVA call moves at WGS, it did not happen."
            RC=1; }
        echo "    positions gained/lost vs BASELINE-0:"
        diff <(grep -v '^#' "$f" | cut -f1,2) <(grep -v '^#' "$b" | cut -f1,2) |
            grep -E '^[<>]' | sed 's/^/      /' | head -40 || true ;;
    esac
done
echo "  identical: $n_same    differing: $n_diff  (expect 2 and 2 -- L1+Alu, SVA+merged)"

echo
echo "--- check 2: call counts. L1 = 188 and Alu = 1223 are gated; SVA is reported ---"
l1=$(grep -vc '^#' "$B/merged/${VCF_PREFIX}_LINE1.vcf" 2>/dev/null || echo 0)
alu=$(grep -vc '^#' "$B/merged/${VCF_PREFIX}_ALU.vcf" 2>/dev/null || echo 0)
sva=$(grep -vc '^#' "$B/merged/${VCF_PREFIX}_SVA.vcf" 2>/dev/null || echo 0)
mrg=$(grep -vc '^#' "$B/$SAMPLE.vcf" 2>/dev/null || echo 0)
sva0=$(grep -vc '^#' "$BASE0/BASELINE-0/${VCF_PREFIX}_SVA.vcf" 2>/dev/null || echo 0)
echo "  L1=$l1 (BASELINE-0 188)  Alu=$alu (1223)  SVA=$sva (BASELINE-0 $sva0)  merged=$mrg (1493)"
{ [ "$l1" -eq 188 ] && [ "$alu" -eq 1223 ]; } || { echo "  FAIL: L1 or Alu moved"; RC=1; }
echo "  SVA delta: $(( sva - sva0 )) calls"

echo
echo "--- check 3: the ledgers ---"
mapfile -t VOCAB < <(python - "$WT" <<'EOF'
import sys, os
sys.path.insert(0, os.path.join(sys.argv[1], "xtea"))
import transduction_evidence as te
for r in sorted(te.PHASE_OF_REASON):
    print(r)
EOF
)
echo "  vocabulary (${#VOCAB[@]} reasons) -- unchanged by this fix, which adds no reason"
[ "${#VOCAB[@]}" -ge 18 ] || { echo "    FATAL: expected >=18 reasons, got ${#VOCAB[@]}"; RC=1; }
n_led=0
while IFS= read -r led; do
    n_led=$((n_led+1))
    rows=$(( $(wc -l < "$led") - 1 ))
    fam=$(basename "$(dirname "$led")")
    echo "  $fam: $rows rows  (job 49703322: L1 21158, SVA 27787)"
    # L1's ledger is a second, finer statement of the same prediction as check 1: the fix cannot
    # touch L1, so not one of its 21,158 rows may move. It is a stronger test than the VCF,
    # because a locus can be dropped at a different gate and still not be called either way.
    if [ "$fam" = L1 ]; then
        [ "$rows" -eq "$EXPECT_LEDGER_L1" ] || {
            echo "    FAIL: L1's ledger moved ($rows vs $EXPECT_LEDGER_L1). L1 collides with no"
            echo "    decoy tag, so the fix must be invisible to every L1 gate, not just to the"
            echo "    calls it emits."; RC=1; }
    fi
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

echo
echo "--- check 4: downstream_parallel, loosely bounded ---"
# Prior no-flag WGS runs: 462.7s (49550925), 511.75s (49599838, a ~10% slower node), 463.3s
# (49634290), 485.8s (49703322). The fix does strictly MORE work in SVA's disc parse, so the
# bound is deliberately looser than those gates' 700s.
secs=$(grep -aoE 'downstream_parallel[^0-9]*([0-9]+\.[0-9]+)s' "$OUT/driver.log" | tail -1 |
       grep -oE '[0-9]+\.[0-9]+' | tail -1 || true)
grep -a 'downstream_parallel' "$OUT/driver.log" | tail -3 | sed 's/^/  /' || true
if [ -n "$secs" ]; then
    echo "  downstream_parallel = ${secs}s (bound: <= ${MAX_DOWNSTREAM_S}s)"
    awk -v s="$secs" -v m="$MAX_DOWNSTREAM_S" 'BEGIN {exit !(s <= m)}' || {
        echo "  FAIL: downstream_parallel above the loose bound"; RC=1; }
else
    echo "  (could not parse downstream_parallel -- reporting, not failing)"
fi
echo "  the SVA transduction step's own phases, for where the extra work landed:"
grep -a 'TD-PHASE' "$B"/*.steps.log 2>/dev/null | awk -F'\t' '$3==4 || $3=="4"' |
    sed 's/^/    /' | head -20 || echo "    (no [TD-PHASE] lines found)"

echo
echo "--- check 5: THE PAYOFF -- SVA reference-copy sources must no longer be zero ---"
probe_out="$OUT/probe_td_source_classes.txt"
python "$WT/docs/handoff/runs/probe_td_source_classes.py" \
    --rep-lib-dir "$REPLIB" --run-dir "$B" > "$probe_out" 2>&1 || {
        echo "  the probe exited nonzero:"; sed 's/^/  /' "$probe_out"; RC=1; }
sed 's/^/  /' "$probe_out"
# The probe prints `SVA  picked sources N  reference copies R  polymorphic P`. R is the number
# this whole change exists to move off zero.
n_ref_sva=$(awk '/^  SVA / {for (i=1;i<=NF;i++) if ($i=="copies") {print $(i+1); exit}}' "$probe_out")
echo "  SVA reference-copy picked sources: ${n_ref_sva:-<unparsed>} (BASELINE-0: 0 of 105)"
if [ -z "${n_ref_sva:-}" ]; then
    echo "    FAIL: could not parse the probe's SVA line -- the check did not run"; RC=1
elif [ "$n_ref_sva" -le 0 ]; then
    echo "    FAIL: still zero. Every hit to the 5,454-record hg38 SVA flank library is still"
    echo "    being discarded, so whatever this arm changed, it was not the decoy test."; RC=1
else
    echo "    PASS: reference SVA sources are reachable for the first time"
fi

echo
echo "--- check 6: the truth-set score, BASELINE-0 vs this arm (REPORTED, not gated) ---"
# Read §7.1 before reading this as a verdict: the derived transduction denominator on HG002 is
# 0 (L1) and 2 (SVA), so this scores POSITIONS, not transductions. SVA position sensitivity /
# precision at BASELINE-0 is 85.6% / 93.9% over 90 truth records; 5 FPs is the whole FP budget,
# so a handful of new calls can move precision several points in either direction.
if [ -d "$TRUTH_DIR" ]; then
    for arm in BASELINE-0 THIS_ARM; do
        if [ "$arm" = BASELINE-0 ]; then v="$BASE0/BASELINE-0/HG002.wgs.vcf"; else v="$B/$SAMPLE.vcf"; fi
        echo "  === $arm ($v) ==="
        if python "$WT/benchmark/benchmark_truth_set.py" score \
                --calls "$v" --truth-dir "$TRUTH_DIR" --families L1,Alu,SVA \
                > "$OUT/score_$arm.txt" 2>&1; then
            sed -n '/--- L1 ---/,$p' "$OUT/score_$arm.txt" | sed 's/^/    /'
        else
            echo "    (the scorer exited nonzero -- reporting, not failing)"
            tail -20 "$OUT/score_$arm.txt" | sed 's/^/    /'
        fi
    done
else
    echo "  (no truth dir at $TRUTH_DIR -- skipped; set XTEA2_TRUTH_DIR)"
fi

echo
echo "--- check 7: the FN attribution, rerun (REPORTED, not gated) ---"
if [ -d "$TRUTH_DIR" ]; then
    python "$WT/benchmark/ledger_attribution.py" \
        --run-dir "$B" --calls "$B/$SAMPLE.vcf" --truth-dir "$TRUTH_DIR" \
        > "$OUT/attribution.txt" 2>&1 && sed 's/^/  /' "$OUT/attribution.txt" || {
            echo "  (the attribution tool exited nonzero -- reporting, not failing)"
            sed 's/^/  /' "$OUT/attribution.txt"; }
else
    echo "  (no truth dir -- skipped)"
fi

echo
echo "--- check 8: the funnel probe still reconciles (the gates must not have shifted) ---"
funnel_out="$OUT/probe_funnel.txt"
if python "$WT/docs/handoff/runs/probe_decision_funnel.py" --run-dir "$B" > "$funnel_out" 2>&1; then
    sed 's/^/  /' "$funnel_out"
    grep -q 'reconciliation checks: FAILED' "$funnel_out" && {
        echo "    FAIL: the ledger and the site files disagree about a gate's population"; RC=1; } || true
else
    echo "  the probe exited nonzero:"; sed 's/^/  /' "$funnel_out"; RC=1
fi

echo
echo "=============================================================="
if [ "$RC" -eq 0 ]; then
    echo "GATE decoy (WGS): PASS"
else
    echo "GATE decoy (WGS): FAIL"
fi
echo "  finished: $(date '+%F %T')"
echo "=============================================================="
exit "$RC"
