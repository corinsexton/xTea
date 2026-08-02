#!/bin/bash
#SBATCH -c 24
#SBATCH -t 10:00:00
#SBATCH --mem=64G
#SBATCH -p priopark
#SBATCH --account=park_contrib
#SBATCH -o /n/scratch/users/c/cos689/xtea2_testdata/stage2c_family_assign_%j.out
#SBATCH --mail-type=NONE

# =============================================================================
# STAGE 2c.3 GATE -- family assignment moves from the postF_hc merge up to
# cns_filter (--family-assign-early).
#
#     sbatch docs/handoff/runs/stage2c_family_assign_wgs.sbatch
#
# THIS IS THE FIRST GATE IN THE OVERHAUL THAT A BYTE-DIFF CANNOT DECIDE. Stages 0
# and 1 were mechanics: any moved call was a bug. This stage is *expected* to move
# calls, so the gate has to say in advance which ones, and how many, and check that
# they moved the way it said. Plan sec 4 Stage 2's wording is "calls may move, but
# only in a predicted direction"; below is the prediction.
#
# THE EXPECTATION, corrected 2026-08-02 by job 49502115 having actually run:
#
#   * ZERO calls move. 1493 in, 1493 out, per family 188/1223/82 in both arms.
#   * The reason is structural, not lucky. The barrier does not REPLACE the postF_hc
#     merge -- the merge still runs and still resolves family over whatever survives.
#     So the barrier changes output only where it removes a row that would have WON at
#     the merge, and on HG002 it never does: of the 16 conflicts it resolves, 7 are
#     also seen by the merge and all 7 get the identical winner/loser/rule, and 9 the
#     merge never sees at all.
#   * The merge's 22 conflicts become 15 in the ON arm (22 - 7). Those 15 are NOT
#     reachable from cns: the losing family has no cns row, because transduction and
#     sibling create its record after the barrier has already run. 12 of them are the
#     Alu-over-L1 tprt_override population.
#
# An earlier version of this header predicted 4 moved loci, from
# docs/handoff/runs/probe_early_assign.py. That prediction was wrong and the probe is
# fixed: it compared the two resolutions as ALTERNATIVES when both actually run, and it
# matched a late L1 locus to the nearest early locus of ANY family, reporting a "flip"
# at positions where L1 simply had no cns row yet. The barrier left L1 at 300 of 300.
#
# SO WHAT THIS GATE NOW CHECKS is that the relocation is output-neutral, and WHY:
# not "did the numbers stay the same" but "did the two decision points ever disagree".
# A disagreement is a moved call and must be justified against truth; agreement means
# the decision can be moved without changing answers, which is what licenses removing
# the second decision point later.
#
# TWO ARMS, for the same reason Stage 1 used two: BASELINE-0 was frozen at 9a43b80 and
# HEAD is later, so a lone ON arm conflates "the barrier did this" with "something else
# moved since". Arm OFF must reproduce BASELINE-0 exactly; only then does arm ON's delta
# mean what it claims.
#
# Frozen worktree, not the live checkout -- run_xtea spawns a fresh x_TEA_main per step,
# so a job sharing the checkout picks up mid-run edits.
# =============================================================================
set -euo pipefail

SCRATCH=/n/scratch/users/c/cos689/xtea2_testdata
WT="${XTEA2_WORKTREE:-$SCRATCH/stage2c_wt}"
OUT="${XTEA2_OUT:-$SCRATCH/stage2c_family_assign_${SLURM_JOB_ID:-manual}}"
BASE0="${XTEA2_BASELINE0:-$SCRATCH/baseline0_hg002_L1-Alu-SVA}"
TESTDATA=/n/data1/hms/dbmi/park/corinne/claude_sandbox/test_data/HG002
BAM="$TESTDATA/alignment/HG002.novaseq.pcr-free.30x.dedup.grch38.bam"
TRUTH="$TESTDATA/truth"
SAMPLE=HG002.wgs
TYPES="${1:-L1,Alu,SVA}"
REF=/n/data1/hms/dbmi/park-smaht_dac/ref/GRCh38_no_alt/hg38_no_alt.fa
CORES="${SLURM_CPUS_PER_TASK:-24}"

source /n/app/conda/miniforge3/24.11.3-0/etc/profile.d/conda.sh
conda activate xtea2
export PATH=/home/cos689/software/bin:$PATH
export XTEA_BASIC_INFO_SEED=42   # pin the insert-size draw, as BASELINE-0 did
test "$PYTHONNOUSERSITE" = 1     # the env's guard against ~/.local's onnxruntime

[ -d "$WT" ]                    || { echo "FATAL: worktree $WT missing"; exit 1; }
[ -x "$WT/src/csrc/xtea_scan" ] || { echo "FATAL: $WT/src/csrc/xtea_scan not built"; exit 1; }
[ -f "$BAM" ]                   || { echo "FATAL: $BAM missing"; exit 1; }
[ -f "$BAM.bai" ]               || { echo "FATAL: $BAM.bai missing"; exit 1; }
[ -d "$TRUTH" ]                 || { echo "FATAL: truth dir $TRUTH missing"; exit 1; }
[ -z "$(git -C "$WT" status --porcelain)" ] || { echo "FATAL: worktree is dirty"; exit 1; }
grep -q "family-assign-early" "$WT/xtea/run_xtea.py" || {
    echo "FATAL: worktree has no --family-assign-early; wrong commit"; exit 1; }
[ -f "$WT/benchmark/benchmark_truth_set.py" ] || {
    echo "FATAL: worktree has no benchmark/benchmark_truth_set.py; this gate scores"
    echo "       against truth and cannot fall back to a byte-diff"; exit 1; }
# The comparator has to exist BEFORE ~90 minutes of compute, not after.
[ -d "$BASE0/BASELINE-0" ] || { echo "FATAL: $BASE0/BASELINE-0 missing"; exit 1; }

mkdir -p "$OUT"
cd "$WT"

echo "=============================================================="
echo "STAGE 2c.3 GATE  --  family assignment at cns_filter, HG002 WGS, types=$TYPES"
echo "  job      : ${SLURM_JOB_ID:-<none>} on $(hostname)  cores=$CORES"
echo "  commit   : $(git rev-parse --short HEAD) $(git log -1 --format=%s | cut -c1-60)"
echo "  bam      : $BAM"
echo "  baseline : $BASE0  ($(grep -m1 '^commit:' "$BASE0/BASELINE-0/PROVENANCE.txt" || true))"
echo "  out      : $OUT"
echo "  started  : $(date '+%F %T')"
echo "=============================================================="

run_arm() {   # $1 = arm name, $2 = "on"|"off"
    local arm=$1 mode=$2 t0 wall armout
    armout="$OUT/$arm"
    rm -rf "$armout"; mkdir -p "$armout"
    echo
    echo "=== arm $arm: family-assign-early=$mode ==="
    t0=$(date +%s)
    # The OFF arm passes no flag because OFF is the default; that is also what makes the
    # "arms ran the paths they claim" check below non-negotiable -- if the flag ever
    # flips to default-on, this arm becomes a silent second ON arm.
    if [ "$mode" = on ]; then FA_FLAG=(--family-assign-early); else FA_FLAG=(); fi
    python xtea/run_xtea.py \
      -i "$BAM" -s "$SAMPLE" -o "$armout" \
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
      "${FA_FLAG[@]}" \
      --parallel-downstream \
      --keep-intermediates \
      -n "$CORES" \
      > "$armout/driver.log" 2>&1 || { echo "FAIL: arm $arm exited $?"; tail -40 "$armout/driver.log"; return 1; }
    wall=$(( $(date +%s) - t0 ))
    printf '  wall %02d:%02d:%02d\n' $((wall/3600)) $((wall%3600/60)) $((wall%60))
    echo "$wall" > "$armout/wall_s.txt"
}

# XTEA2_ANALYZE_ONLY=1 re-runs the checks below against arms that already exist in
# $OUT, without recomputing them. That is for iterating on THIS SCRIPT -- the checks
# had two path bugs that made job 49502115 misreport a run whose arms were fine -- not
# for skipping the compute on a real gate.
if [ "${XTEA2_ANALYZE_ONLY:-0}" = 1 ]; then
    echo; echo "*** ANALYZE-ONLY: reusing existing arms in $OUT (no pipeline run) ***"
    for arm in OFF ON; do
        [ -f "$OUT/$arm/$SAMPLE/$SAMPLE.vcf" ] || {
            echo "FATAL: XTEA2_ANALYZE_ONLY set but $OUT/$arm has no completed run"; exit 1; }
    done
else
    run_arm OFF off
    run_arm ON  on
fi

echo
echo "##############################################################"
echo "# STAGE 2c.3 GATE RESULT"
echo "##############################################################"

RC=0
# driver.log, NOT <sample>_xtea2.log: run_xtea's log() writes to stdout, which run_arm
# redirects here. The xtea2.log carries the per-STEP output. Grepping the wrong one made
# job 49502115 report "arm ON did not take the barrier" about a run whose barrier had in
# fact resolved 16 conflicts and written its sidecar.
LOG_OFF="$OUT/OFF/driver.log"
LOG_ON="$OUT/ON/driver.log"

# --- 0. the arms ran the paths they claim -------------------------------------------
# Nothing below can tell a real result from two OFF arms agreeing with each other.
echo
echo "--- 0. arms ran the paths they claim ---"
if grep -q "family assignment: EARLY" "$LOG_OFF"; then
    echo "  FATAL: arm OFF took the barrier"; RC=1
else
    echo "  arm OFF : no barrier (correct)"
fi
if grep -q "family assignment: EARLY" "$LOG_ON"; then
    echo "  arm ON  : barrier present (correct)"
    grep -A8 "family assignment barrier" "$LOG_ON" | sed 's/^/    /'
else
    echo "  FATAL: arm ON did not take the barrier"; RC=1
fi
# The barrier must also have had something to do. A run where it resolves 0 conflicts
# is indistinguishable from a run where the filtering silently no-op'd.
[ -s "$OUT/ON/$SAMPLE/cns_family_conflicts.tsv" ] || {
    echo "  FATAL: arm ON wrote no cns_family_conflicts.tsv"; RC=1; }

# --- 1. arm OFF vs frozen BASELINE-0 -- isolates commit drift ------------------------
echo
echo "--- 1. arm OFF vs frozen BASELINE-0 VCFs (must be IDENTICAL: isolates drift) ---"
n_same=0; n_diff=0
for f in "$BASE0"/BASELINE-0/*.vcf; do
    b=$(find "$OUT/OFF/$SAMPLE" -name "$(basename "$f")" -type f | head -1)
    if [ -z "$b" ]; then echo "  MISSING: $(basename "$f")"; n_diff=$((n_diff+1)); continue; fi
    if cmp -s "$f" "$b"; then n_same=$((n_same+1)); else
        echo "  DIFFERS: $(basename "$f")"; n_diff=$((n_diff+1)); fi
done
echo "  identical: $n_same    differing: $n_diff"
[ "$n_diff" -eq 0 ] || { echo "  FAIL: arm OFF is not BASELINE-0; arm ON's delta is uninterpretable"; RC=1; }

# --- 2. the movement, locus by locus -------------------------------------------------
# Call identity for this comparison is (CHROM, POS, family). Everything else in the
# record -- INFO, genotype -- can move for reasons this stage does not control.
# The per-type VCFs live in merged/ and are named after the BAM, not the sample, because
# move A.2 made gVCF fan out per record from one shared tail: merged/<bam>_LINE1.vcf etc.
# The old glob ("$d/$SAMPLE"/*/"$SAMPLE"_*.vcf) matched nothing and this leg silently
# compared two empty sets, reporting "0 moved" as a pass (job 49502115).
calls() {   # $1 = arm dir -> "chrom<TAB>pos<TAB>family"
    local d=$1 vcf fam n=0
    for fam in LINE1 ALU SVA; do
        for vcf in "$d/$SAMPLE"/merged/*_"$fam".vcf; do
            [ -f "$vcf" ] || continue
            n=$((n+1))
            awk -v f="$fam" -F'\t' '!/^#/{print $1"\t"$2"\t"f}' "$vcf"
        done
    done | sort -k1,1V -k2,2n
    [ "$n" -gt 0 ] || echo "  FATAL: no per-type VCFs found under $d/$SAMPLE/merged" >&2
}
calls "$OUT/OFF" > "$OUT/calls_off.tsv"
calls "$OUT/ON"  > "$OUT/calls_on.tsv"
echo
echo "--- 2. movement: arm ON vs arm OFF ---"
echo "  calls OFF: $(wc -l < "$OUT/calls_off.tsv")    calls ON: $(wc -l < "$OUT/calls_on.tsv")"
# An empty call set is not agreement. Without this, a broken path makes leg 2 pass.
if [ ! -s "$OUT/calls_off.tsv" ] || [ ! -s "$OUT/calls_on.tsv" ]; then
    echo "  FATAL: one or both arms produced no per-type calls -- leg 2 cannot compare"; RC=1
fi
comm -23 "$OUT/calls_off.tsv" "$OUT/calls_on.tsv" > "$OUT/only_off.tsv"
comm -13 "$OUT/calls_off.tsv" "$OUT/calls_on.tsv" > "$OUT/only_on.tsv"
N_OFF_ONLY=$(wc -l < "$OUT/only_off.tsv"); N_ON_ONLY=$(wc -l < "$OUT/only_on.tsv")
echo "  only in OFF: $N_OFF_ONLY      only in ON: $N_ON_ONLY"
echo "  --- every moved locus ---"
{ sed 's/^/    OFF-only  /' "$OUT/only_off.tsv"; sed 's/^/    ON-only   /' "$OUT/only_on.tsv"; } | head -80
# The bound: 22 is the whole contended population measured on BASELINE-0
# (merged/locus_family_conflicts.tsv). Movement beyond it is not family reassignment.
N_MOVED=$(( N_OFF_ONLY > N_ON_ONLY ? N_OFF_ONLY : N_ON_ONLY ))
if [ "$N_MOVED" -gt 22 ]; then
    echo "  FAIL: $N_MOVED loci moved, more than the 22 loci any family contends for."
    echo "        This is not family reassignment -- find out what else changed."
    RC=1
else
    echo "  OK: $N_MOVED moved, within the 22-locus contended population"
fi
# Every moved locus should appear in ONE of the two conflict sidecars. One that appears
# in neither was never contested, so the barrier is not what moved it.
echo "  --- moved loci NOT on either conflict list (each one is unexplained) ---"
cut -f1,2 "$OUT/only_off.tsv" "$OUT/only_on.tsv" | sort -u > "$OUT/moved_pos.tsv"
cat "$OUT/ON/$SAMPLE/cns_family_conflicts.tsv" \
    "$OUT/OFF/$SAMPLE/merged/locus_family_conflicts.tsv" 2>/dev/null \
  | awk -F'\t' '!/^#/{print $1"\t"$3; print $1"\t"$6}' | sort -u > "$OUT/conflict_pos.tsv"
# +-50bp, the merge window: a conflict is recorded at the winner's and loser's own
# positions, and the VCF position is derived from the surviving record.
join -t$'\t' -j1 <(awk -F'\t' '{print $1"@"int($2/50)"\t"$0}' "$OUT/moved_pos.tsv" | sort) \
                 <(awk -F'\t' '{print $1"@"int($2/50)"\t"$0}' "$OUT/conflict_pos.tsv" | sort) \
     > "$OUT/moved_explained.tsv" 2>/dev/null || true
N_UNEXPLAINED=$(( $(wc -l < "$OUT/moved_pos.tsv") - $(cut -f1 "$OUT/moved_explained.tsv" | sort -u | wc -l) ))
echo "    unexplained: $N_UNEXPLAINED (bin-matched at 50bp; inspect by hand, this is advisory)"

# --- 2b. do the two decision points ever disagree? -----------------------------------
# THE STRUCTURAL CHECK, and the one that licenses eventually deleting the merge's
# resolution. Leg 2 says the output did not change; this says WHY. A locus both stages
# resolve must get the same winner, loser and rule at each. Compared on the winner
# position, which both sidecars record identically.
echo
echo "--- 2b. barrier decisions vs merge decisions, on loci both stages resolve ---"
awk -F'\t' '!/^#/{print $1":"$3"\t"$2"/"$5"/"$8}' \
    "$OUT/ON/$SAMPLE/cns_family_conflicts.tsv" | sort > "$OUT/dec_barrier.tsv"
awk -F'\t' '!/^#/{print $1":"$3"\t"$2"/"$5"/"$8}' \
    "$OUT/OFF/$SAMPLE/merged/locus_family_conflicts.tsv" | sort > "$OUT/dec_merge.tsv"
join -t$'\t' -j1 "$OUT/dec_barrier.tsv" "$OUT/dec_merge.tsv" > "$OUT/dec_both.tsv" || true
N_BOTH=$(wc -l < "$OUT/dec_both.tsv")
N_DISAGREE=$(awk -F'\t' '$2!=$3' "$OUT/dec_both.tsv" | tee "$OUT/dec_disagree.tsv" | wc -l)
echo "  resolved at BOTH stages: $N_BOTH    disagreements: $N_DISAGREE"
awk -F'\t' '$2!=$3{print "    DISAGREE "$1"  barrier="$2"  merge="$3}' "$OUT/dec_both.tsv"
echo "  barrier-only (merge never sees these): $(( $(wc -l < "$OUT/dec_barrier.tsv") - N_BOTH ))"
echo "  merge-only   (no cns row for the loser; transduction/sibling create it): \
$(( $(wc -l < "$OUT/dec_merge.tsv") - N_BOTH ))"
if [ "$N_DISAGREE" -ne 0 ]; then
    echo "  Each disagreement is a call that MOVES. Justify every one against truth in leg 3"
    echo "  before treating this as a pass -- movement is allowed, unexplained movement is not."
fi
[ "$N_BOTH" -gt 0 ] || { echo "  FATAL: the two stages share no resolved locus -- one sidecar is empty or mis-parsed"; RC=1; }

# --- 3. the direction: score both arms against the HG002 truth set -------------------
# THE decisive leg. sec 2 says movement is bounded; this says it is movement toward
# truth. NOT comparable to the paper's numbers (hs38DH; plan sec 7.1) -- the report
# says so on every run, and the comparison that matters here is ON vs OFF, same
# reference, same everything else.
echo
echo "--- 3. truth-set score, both arms ---"
for arm in OFF ON; do
    echo
    echo "  === arm $arm ==="
    python benchmark/benchmark_truth_set.py score \
        --calls "$OUT/$arm/$SAMPLE/$SAMPLE.vcf" \
        --truth-dir "$TRUTH" \
        2>&1 | sed 's/^/    /' | tee "$OUT/score_$arm.txt"
done
echo
echo "  --- the comparison the gate turns on ---"
echo "  arm ON must not lose sensitivity or precision in ANY family. Read the two"
echo "  blocks above side by side; the predicted delta is +1 SVA TP and -3 L1 FP."

echo
echo "=============================================================="
[ "$RC" -eq 0 ] && echo "STAGE 2c.3 GATE: PASS (bounds held; read leg 3 for the direction)" \
                || echo "STAGE 2c.3 GATE: FAIL"
echo "  finished : $(date '+%F %T')"
echo "  wall OFF : $(cat "$OUT/OFF/wall_s.txt" 2>/dev/null || echo '?')s"
echo "  wall ON  : $(cat "$OUT/ON/wall_s.txt" 2>/dev/null || echo '?')s"
echo "=============================================================="
exit $RC
