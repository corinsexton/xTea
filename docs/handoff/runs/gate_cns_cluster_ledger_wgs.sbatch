#!/bin/bash
#SBATCH -J b2_wgs
#SBATCH -c 24
#SBATCH -t 6:00:00
#SBATCH --mem=64G
#SBATCH -p priopark
#SBATCH --account=park_contrib
#SBATCH -o /n/scratch/users/c/cos689/xtea2_testdata/gate_cns_cluster_ledger_wgs_%j.out
#SBATCH --mail-type=NONE
# =============================================================================
# Gate [B2] at WGS scale. ONE arm, against the FROZEN BASELINE-0 -- and this is the leg that
# actually exercises the population, because chr1's B2 domain is 5 loci for L1 and empty for SVA.
#
# WHY ONE ARM. The correctness claim is "no call moves", and the chr1 gate
# (`gate_cns_cluster_ledger_chr1.sbatch`) establishes that against the parent commit over the
# whole output tree, including the stronger statement that the ledger changes by exactly its new
# B2 rows and in nothing else. BASELINE-0 is frozen with sha256s and is the canonical statement
# of the WGS numbers, so one arm answers the population question at half the compute.
#
# WHAT MUST BE TRUE, and checks 5 and 6 are the ones this gate exists for:
#   1. every frozen BASELINE-0 VCF is byte-identical, and the calls are 188 / 1223 / 82 / 1493.
#      Instrumentation must not move a call.
#   2. both ledgers are present, non-empty, vocabulary-clean, no unknown phase.
#   3. the ledgers now carry B2 rows at WGS population.
#   4. `downstream_parallel` inside a loose bound. NOT a measurement -- one arm on one node has
#      no control, and the last three WGS runs put this at 462.7s, 511.75s (a ~10% slower node)
#      and 463.3s. The bound exists to catch the gate having become structurally expensive.
#   5. **THE PREDICTION.** The funnel probe measured `picked sources` minus `cluster-consistent`
#      at 135 (L1) and 49 (SVA) on job 49634290, from the step's own site FILES and with no
#      ledger involved. This run's B2 rows must be those same numbers. Two derivations from
#      different sources agreeing is the whole point of instrumenting the gate rather than
#      inferring its population -- and if they disagree, the inference in every prior handoff
#      about this function is what was wrong, which is worth knowing.
#   6. the probe's own reconciliation, run on this arm, must print OK rather than FAILED or
#      SKIPPED. It compares the two set-wise, not just by count.
#
# Frozen worktree, not the live checkout -- run_xtea spawns a fresh x_TEA_main per step, so a
# job sharing the checkout picks up mid-run edits.
#
# USAGE  sbatch [--dependency=afterok:<chr1 job>] gate_cns_cluster_ledger_wgs.sbatch <worktree>
# =============================================================================
set -euo pipefail

SCRATCH=/n/scratch/users/c/cos689/xtea2_testdata
WT="${1:?usage: sbatch gate_cns_cluster_ledger_wgs.sbatch <worktree> [types]}"
TYPES="${2:-L1,Alu,SVA}"
OUT="${XTEA2_OUT:-$SCRATCH/gate_cns_cluster_ledger_wgs_${SLURM_JOB_ID:-manual}}"
BASE0="${XTEA2_BASELINE0:-$SCRATCH/baseline0_hg002_L1-Alu-SVA}"
TESTDATA=/n/data1/hms/dbmi/park/corinne/claude_sandbox/test_data/HG002
BAM="$TESTDATA/alignment/HG002.novaseq.pcr-free.30x.dedup.grch38.bam"
SAMPLE=HG002.wgs
REF=/n/data1/hms/dbmi/park-smaht_dac/ref/GRCh38_no_alt/hg38_no_alt.fa
CORES="${SLURM_CPUS_PER_TASK:-24}"
MAX_DOWNSTREAM_S="${XTEA2_MAX_DOWNSTREAM_S:-700}"
# The funnel probe's own measurement on job 49634290, from two site files per family and no
# ledger. Overridable so a future run with different inputs can restate them rather than lie.
EXPECT_B2_L1="${XTEA2_EXPECT_B2_L1:-135}"
EXPECT_B2_SVA="${XTEA2_EXPECT_B2_SVA:-49}"

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
if ! grep -q 'def gate_cns_disc_cluster' "$WT/xtea/transduction_evidence.py"; then
    echo "FATAL: worktree has no gate_cns_disc_cluster; wrong commit"; exit 1
fi
if ! grep -q 'ledger\.record(' "$WT/xtea/x_clip_disc_filter.py"; then
    echo "FATAL: worktree's check_seprt_disc_consistency is not instrumented"; exit 1
fi

mkdir -p "$OUT"
cd "$WT"

echo "=============================================================="
echo "GATE B2 (cns-cluster ledger, WGS)  --  types=$TYPES"
echo "  job      : ${SLURM_JOB_ID:-<none>} on $(hostname)  cores=$CORES"
echo "  commit   : $(git rev-parse --short HEAD) $(git log -1 --format=%s | cut -c1-60)"
echo "  worktree : $WT"
echo "  baseline : $BASE0"
echo "  expect B2: L1=$EXPECT_B2_L1  SVA=$EXPECT_B2_SVA (job 49634290's site-file measurement)"
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
    echo "  FAIL: instrumenting gate B2 moved a call at WGS scale"; RC=1; }

echo
echo "--- check 2: call counts must be 188 / 1223 / 82 / 1493 ---"
l1=$(grep -vc '^#' "$B/merged/HG002.novaseq.pcr-free.30x.dedup.grch38_LINE1.vcf" 2>/dev/null || echo 0)
alu=$(grep -vc '^#' "$B/merged/HG002.novaseq.pcr-free.30x.dedup.grch38_ALU.vcf" 2>/dev/null || echo 0)
sva=$(grep -vc '^#' "$B/merged/HG002.novaseq.pcr-free.30x.dedup.grch38_SVA.vcf" 2>/dev/null || echo 0)
mrg=$(grep -vc '^#' "$B/$SAMPLE.vcf" 2>/dev/null || echo 0)
echo "  L1=$l1 Alu=$alu SVA=$sva merged=$mrg"
[ "$l1" -eq 188 ] && [ "$alu" -eq 1223 ] && [ "$sva" -eq 82 ] && [ "$mrg" -eq 1493 ] || {
    echo "  FAIL: call counts moved"; RC=1; }

echo
echo "--- check 3: the ledgers, at WGS population ---"
mapfile -t VOCAB < <(python - "$WT" <<'EOF'
import sys, os
sys.path.insert(0, os.path.join(sys.argv[1], "xtea"))
import transduction_evidence as te
for r in sorted(te.PHASE_OF_REASON):
    print(r)
EOF
)
echo "  vocabulary (${#VOCAB[@]} reasons)"
[ "${#VOCAB[@]}" -ge 18 ] || { echo "    FATAL: expected >=18 reasons, got ${#VOCAB[@]}"; RC=1; }
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
echo "--- check 4: downstream_parallel, loosely bounded ---"
# One arm, one node, no control: a bound that catches the gate having become structurally
# expensive, not a measurement. Prior no-flag WGS runs: 462.7s (49550925), 511.75s (49599838,
# a ~10% slower node), 463.3s (49634290).
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
echo "--- check 5: THE PREDICTION -- B2 rows must be $EXPECT_B2_L1 (L1) and $EXPECT_B2_SVA (SVA) ---"
# Those two numbers were measured from the step's own site files on job 49634290, with no
# ledger in existence. This is the first time the same population is counted the other way.
for fam in L1 SVA; do
    led="$B/$fam/candidate_disc_filtered_cns2.txt.all_non_sibling_td.txt.decision_ledger.tsv"
    if [ ! -f "$led" ]; then echo "  $fam: NO LEDGER at $led"; RC=1; continue; fi
    got=$(awk -F'\t' 'NR>1 && $3=="B2_cns_cluster"' "$led" | wc -l)
    eval "want=\$EXPECT_B2_$fam"
    echo "  $fam: B2 rows = $got   (site-file prediction: $want)"
    echo "    split: $(awk -F'\t' 'NR>1 && $3=="B2_cns_cluster" {c[$4]++} END {for (r in c) printf "%s=%d ", r, c[r]}' "$led")"
    [ "$got" -eq "$want" ] || {
        echo "    FAIL: the ledger and the site files disagree about gate B2's population."
        echo "    One of the two is wrong, and every prior statement about this function"
        echo "    rested on the site-file derivation."; RC=1; }
done

echo
echo "--- check 6: the funnel probe's own set-wise reconciliation ---"
probe_out="$OUT/probe_funnel.txt"
if python "$WT/docs/handoff/runs/probe_decision_funnel.py" --run-dir "$B" > "$probe_out" 2>&1; then
    sed 's/^/  /' "$probe_out"
    if grep -q 'reconciliation checks: FAILED' "$probe_out"; then
        echo "    FAIL: the probe's reconciliation failed"; RC=1
    fi
    if grep -q 'SKIPPED, not passed' "$probe_out"; then
        echo "    FAIL: a family with picked sources produced no B2 rows at all"; RC=1
    fi
else
    echo "  the probe exited nonzero:"; sed 's/^/  /' "$probe_out"; RC=1
fi

echo
echo "--- check 7: the FN attribution, rerun with gate B2 in the vocabulary ---"
# Not pass/fail: the reason this gate is worth its compute. Until now every attributed miss was
# in gate A and 306 of 322 FNs had no ledger row. Some of those 306 may now attribute to B2 --
# or may not, since B2's domain is only ~160 loci per family. Either answer is information.
TRUTH_DIR="${XTEA2_TRUTH_DIR:-/n/data1/hms/dbmi/park/corinne/claude_sandbox/test_data/HG002/truth}"
if [ -d "$TRUTH_DIR" ]; then
    python "$WT/benchmark/ledger_attribution.py" \
        --run-dir "$B" --calls "$B/$SAMPLE.vcf" --truth-dir "$TRUTH_DIR" \
        > "$OUT/attribution.txt" 2>&1 && sed 's/^/  /' "$OUT/attribution.txt" || {
            echo "  (the attribution tool exited nonzero -- reporting, not failing)"
            sed 's/^/  /' "$OUT/attribution.txt"; }
else
    echo "  (no truth dir at $TRUTH_DIR -- skipped; set XTEA2_TRUTH_DIR)"
fi

echo
echo "=============================================================="
if [ "$RC" -eq 0 ]; then
    echo "GATE B2 (WGS): PASS"
else
    echo "GATE B2 (WGS): FAIL"
fi
echo "  finished: $(date '+%F %T')"
echo "=============================================================="
exit "$RC"
