#!/bin/bash
#SBATCH -J td_noregress
#SBATCH -c 8
#SBATCH -t 3:00:00
#SBATCH --mem=32G
#SBATCH -p short
#SBATCH --account=park
#SBATCH -o /n/scratch/users/c/cos689/xtea2_testdata/td_gate_no_regression_%j.out
#SBATCH --mail-type=NONE
# =============================================================================
# `HANDOFF_2026-08-04_td_gate_built.md` sec 8 step 3's PRECONDITION: measure no-regression
# BEFORE flipping `--td-decoy-structural` on by default. Not "flip and check" -- sec 10.1 is
# explicit that the order is the obligation, because a measurement taken after the flip cannot
# distinguish "the flip was safe" from "the flip moved the thing I am measuring with".
#
# NO xTea RUN HAPPENS HERE. Every arm's VCFs already exist. This job reads them, builds pools,
# and adjudicates from the BAM. That is why it is 8 cores on `short` and not 24 on `priopark`.
#
# THE TWO CHECKS ARE DIFFERENT QUESTIONS AND CONFLATING THEM IS THE TRAP THIS SCRIPT AVOIDS.
#
#   CHECK 1 -- REPRODUCTION, on the frozen 3-arm pool. Same pool TSV, same BAM, same operating
#   point, so the verdicts must come back IDENTICAL. This is not a statement about the caller at
#   all: it says the shape x sidedness restratification (sec 10.2) moved no verdict and that the
#   `--baseline` machinery agrees with a file it did not write. If it fails, nothing downstream
#   of it means anything.
#
#   CHECK 2 -- THE ONE-FLAG DELTA: two run trees that differ by exactly one flag, each as its
#   OWN single-arm pool. THIS is the no-regression measurement. The DEFAULT pair is the
#   `--td-decoy-structural` flip -- BASELINE-0 (OFF, then the default) against the decoy arm
#   (ON, the proposed default) -- which is what `gate_decoy_structural_wgs.sbatch` already
#   relies on when it calls BASELINE-0 "the frozen OFF arm". Any other pair goes in through
#   XTEA2_OFF_ARM / XTEA2_ON_ARM; see the USAGE block.
#
# WHY CHECK 2 IS NOT "RUN THE FROZEN 3-ARM BASELINE AGAIN AFTER THE FLIP". Because the frozen
# baseline's pool is the union over master + BASELINE-0 + decoy, 115 entries of which master
# contributes 21 uniquely. A post-flip single-arm pool would hold ~90 entries and would therefore
# "lose" verified entries in almost every stratum -- from ARM SCOPE, not from the flag. That
# comparison would report a large FAIL that is entirely an artifact of which arms are in the
# pool. **This is the same error the plan already records three times** (sec 10's two
# "read a set difference as a decision" entries, and the handoff's "241 measured at the wrong
# level"): a number taken from the wrong population. The frozen 3-arm baseline is the right
# reference for the gate's PRECISION RANGE and the wrong reference for a one-flag delta.
#
# WHAT A PASS DOES NOT MEAN (sec 10.4). It means no stratum lost a verified entry. It does NOT
# mean precision did not fall -- a stratum can hold its VERIFIED count while gaining REFUTED
# entries, and this gate says nothing about that by construction. And the margin is PROVISIONAL
# until an arm exists that is neither master nor xTea2, because two of the three arms here differ
# by one flag: three-arm agreement is a floor on independence, not a measure of it.
#
# CHECK 2's TWO ARMS ARE OVERRIDABLE, and the defaults are the flip's. `XTEA2_OFF_ARM` /
# `XTEA2_ON_ARM` name any two run trees that differ by ONE flag, so every future one-flag arm
# gets this instrument instead of a copy of it -- which is the point: the "single-arm pool vs
# 3-arm baseline" error above is easy to re-make by hand and impossible to make through here.
# `XTEA2_OFF_LEDGER` is the ledger the OFF arm must reproduce; set it to `none` for an arm that
# has no frozen ledger (only BASELINE-0 does), which weakens check 2's OFF leg to "the pool was
# built" and says so out loud rather than silently.
#
# USAGE  sbatch docs/handoff/runs/td_gate_no_regression.sbatch [<repo or worktree>]
#   two arbitrary arms, one flag apart:
#     XTEA2_OFF_ARM=<tree> XTEA2_ON_ARM=<tree> XTEA2_OFF_LEDGER=none XTEA2_ARMS_LABEL='...' \
#       sbatch docs/handoff/runs/td_gate_no_regression.sbatch <worktree>
# =============================================================================
set -euo pipefail

SCRATCH=/n/scratch/users/c/cos689/xtea2_testdata
WT="${1:-/n/data1/hms/dbmi/park/corinne/claude_sandbox/xTea}"
OUT="${XTEA2_OUT:-$SCRATCH/td_gate_no_regression_${SLURM_JOB_ID:-manual}}"
FROZEN="${XTEA2_TD_FROZEN:-$SCRATCH/td_gate_baseline_hg002}"
BASE0="${XTEA2_BASELINE0:-$SCRATCH/baseline0_hg002_L1-Alu-SVA}"
DECOY="${XTEA2_DECOY_ARM:-$SCRATCH/gate_decoy_structural_wgs_49724056/HG002.wgs}"
OFF_ARM="${XTEA2_OFF_ARM:-$BASE0/BASELINE-0}"
ON_ARM="${XTEA2_ON_ARM:-$DECOY}"
OFF_LEDGER="${XTEA2_OFF_LEDGER:-baseline0_hg002}"
ARMS_LABEL="${XTEA2_ARMS_LABEL:-the --td-decoy-structural flip: BASELINE-0 (OFF) vs the decoy arm (ON)}"
TESTDATA=/n/data1/hms/dbmi/park/corinne/claude_sandbox/test_data/HG002
BAM="$TESTDATA/alignment/HG002.novaseq.pcr-free.30x.dedup.grch38.bam"
REF=/n/data1/hms/dbmi/park-smaht_dac/ref/GRCh38_no_alt/hg38_no_alt.fa
REPLIB=/home/cos689/xTea/rep_lib_annotation
CORES="${SLURM_CPUS_PER_TASK:-8}"

source /n/app/conda/miniforge3/24.11.3-0/etc/profile.d/conda.sh
conda activate xtea2
export PATH=/home/cos689/software/bin:$PATH
test "$PYTHONNOUSERSITE" = 1     # the env's guard against ~/.local's onnxruntime

PROBE="$WT/docs/handoff/runs/probe_td_pool.py"
VERIFY="$WT/benchmark/verify_transduction.py"

[ -d "$WT" ]            || { echo "FATAL: $WT missing"; exit 1; }
[ -f "$BAM" ]           || { echo "FATAL: $BAM missing"; exit 1; }
[ -f "$BAM.bai" ]       || { echo "FATAL: $BAM.bai missing"; exit 1; }
[ -f "$REF" ]           || { echo "FATAL: $REF missing"; exit 1; }
[ -d "$OFF_ARM" ]       || { echo "FATAL: $OFF_ARM missing (the OFF arm)"; exit 1; }
[ -d "$ON_ARM" ]        || { echo "FATAL: $ON_ARM missing (the ON arm)"; exit 1; }
command -v bwa >/dev/null || { echo "FATAL: bwa not on PATH"; exit 1; }
# The code under test has to be present BEFORE the compute, not after.
grep -q 'shape_stratum' "$VERIFY" || { echo "FATAL: $VERIFY has no shape_stratum; wrong commit"; exit 1; }
grep -q '"--baseline"' "$VERIFY"  || { echo "FATAL: $VERIFY has no --baseline; wrong commit"; exit 1; }

mkdir -p "$OUT"
RC=0

echo "=============================================================="
echo "TD GATE: NO-REGRESSION, one flag apart"
echo "  job      : ${SLURM_JOB_ID:-<none>} on $(hostname)  cores=$CORES"
echo "  commit   : $(git -C "$WT" rev-parse --short HEAD) $(git -C "$WT" log -1 --format=%s | cut -c1-60)"
echo "  frozen   : $FROZEN"
echo "  arms     : $ARMS_LABEL"
echo "  OFF arm  : $OFF_ARM  (ledger: $OFF_LEDGER)"
echo "  ON  arm  : $ON_ARM"
echo "  out      : $OUT"
echo "=============================================================="

# --- check 0: the frozen baseline is what it was frozen as -------------------------------
# Before it is used as a reference, not after. A silently edited baseline turns every check
# below into a comparison against an unknown.
echo
echo "--- check 0: the frozen baseline's sha256s ---"
if ( cd "$FROZEN" && sha256sum -c SHA256SUMS ); then
    echo "  PASS: the frozen pool and verdicts are byte-identical to what was frozen"
else
    echo "  FATAL: the frozen baseline has CHANGED. Every check below would compare against an"
    echo "  unknown reference. Stopping."
    exit 1
fi

# --- check 1: reproduction on the frozen 3-arm pool --------------------------------------
echo
echo "--- check 1: REPRODUCTION -- the frozen pool, re-adjudicated ---"
echo "  Same pool, same BAM, same operating point => the verdicts must be IDENTICAL. This is a"
echo "  statement about the restratification and the --baseline reader, NOT about the caller."
R1="$OUT/repro_3arm"
mkdir -p "$R1"
set +e
python "$VERIFY" \
    --pool-tsv "$FROZEN/td_pool_hg002_3arm.tsv" \
    --bam "$BAM" --ref "$REF" --rep-lib-dir "$REPLIB" \
    --threads "$CORES" --out "$R1" \
    --baseline "$FROZEN/verdicts.tsv" > "$OUT/repro_3arm.report.txt" 2>"$OUT/repro_3arm.err.txt"
rc1=$?
set -e
sed 's/^/    /' "$OUT/repro_3arm.report.txt"
if [ "$rc1" -eq 3 ]; then
    echo "  FATAL: the operating point did not match. See the report above -- this is the"
    echo "  refusal, not a regression, and it means the sweep midpoints moved."
    RC=1
elif [ "$rc1" -ne 0 ]; then
    echo "  FAIL: reproduction reported a regression, which on an UNCHANGED pool and BAM can only"
    echo "  mean the verifier itself moved a verdict. Do not proceed to check 2."
    RC=1
else
    echo "  PASS: no stratum lost a verified entry on the frozen pool."
fi

# The stronger statement, and it is cheap: the per-entry file must be byte-identical to the
# frozen one once the two NEW columns are removed. `shape` and `sidedness` are inserted at
# columns 6-7, so 1-5,8- is exactly the frozen column set. Comment lines carry no tab and pass
# through `cut` unchanged, so the header comparison is included rather than skipped.
echo
echo "  --- byte-identity of the per-entry verdicts, minus the two new columns ---"
if [ -f "$R1/verdicts.tsv" ]; then
    cut -f1-5,8- "$R1/verdicts.tsv" > "$OUT/repro_stripped.tsv"
    if diff -q "$FROZEN/verdicts.tsv" "$OUT/repro_stripped.tsv" >/dev/null; then
        echo "  PASS: byte-identical. The restratification is REPORTING ONLY -- it moved no"
        echo "        verdict, no support count and no repeat-divergence number."
    else
        echo "  FAIL: the per-entry verdicts differ beyond the two new columns:"
        diff "$FROZEN/verdicts.tsv" "$OUT/repro_stripped.tsv" | head -30 | sed 's/^/      /'
        RC=1
    fi
else
    echo "  FAIL: $R1/verdicts.tsv was not written"; RC=1
fi

# --- check 2: the flip, one flag apart ---------------------------------------------------
echo
echo "--- check 2: THE ONE-FLAG DELTA -- $ARMS_LABEL ---"
echo "  Each arm is its OWN single-arm pool. Comparing a single-arm pool against the frozen"
echo "  3-arm baseline would 'lose' entries from ARM SCOPE and not from the flag."

echo
LEDGER_ARGS=()
if [ "$OFF_LEDGER" != none ]; then
    LEDGER_ARGS=(--ledger "$OFF_LEDGER")
    echo "  building the OFF arm's pool (with its ledger, so a mis-resolved tree says so)"
else
    # Only BASELINE-0 has a frozen ledger. An arm without one still builds a pool; what is lost
    # is the assertion that the tree resolved to the run it claims to be. Said out loud, because
    # a check that silently degrades is worse than one that is absent.
    echo "  building the OFF arm's pool (NO LEDGER for this arm -- the 'is this the run it says"
    echo "  it is' assertion is NOT made; only BASELINE-0 has a frozen ledger)"
fi
set +e
python "$PROBE" --run-dir "$OFF_ARM" "${LEDGER_ARGS[@]}" \
    --pool-tsv "$OUT/pool_off.tsv" > "$OUT/pool_off.report.txt" 2>&1
rc_off=$?
set -e
sed 's/^/      /' "$OUT/pool_off.report.txt"
[ "$rc_off" -eq 0 ] || { echo "  FAIL: the OFF arm's pool did not build (or does not reproduce"
                         echo "  its ledger), so the comparison is void."; RC=1; }

echo
echo "  building the ON arm's pool (no ledger exists for it -- reported, not asserted)"
python "$PROBE" --run-dir "$ON_ARM" --pool-tsv "$OUT/pool_on.tsv" \
    > "$OUT/pool_on.report.txt" 2>&1
sed 's/^/      /' "$OUT/pool_on.report.txt"

# --- check 2a: the flip is PURELY ADDITIVE at the pool level ----------------------------
# STATED AS A PREDICTION, and measured offline before any BAM work, because it makes check 2
# sharp. Measured 2026-08-05 on these two arms: the ON pool is a STRICT SUPERSET of the OFF
# pool -- all 90 OFF entries appear verbatim (same locus, same source, same TD_SRC), plus
# exactly the 4 SVA reference-copy entries the handoff sec 4 names.
#
# WHY THIS MATTERS FOR CHECK 2. If no hypothesis is removed, then a per-stratum VERIFIED loss
# cannot come from the flag at all -- it could only come from the same entry, in the same BAM,
# at the same operating point, getting a different verdict. So on a purely additive pool, check
# 2 is ALSO a determinism check on the verifier, and a FAIL should be read that way first.
#
# And the two failures are different in kind, which is why they are separate checks: an entry
# DISAPPEARING from the pool means the fix suppressed a call the OFF arm made, which is a
# behavioural regression the verdict comparison would never see (the entry is simply absent
# from both sides of a per-stratum count that only sums VERIFIED).
echo
echo "--- check 2a: is the flip purely additive at the pool level? ---"
for f in off on; do
    grep -v '^#' "$OUT/pool_$f.tsv" | tail -n +2 | cut -f1-8,11 | sort > "$OUT/key_$f.txt"
done
n_off=$(wc -l < "$OUT/key_off.txt"); n_on=$(wc -l < "$OUT/key_on.txt")
lost=$(comm -23 "$OUT/key_off.txt" "$OUT/key_on.txt")
added=$(comm -13 "$OUT/key_off.txt" "$OUT/key_on.txt")
echo "  OFF $n_off entries, ON $n_on entries"
if [ -n "$lost" ]; then
    echo "  FAIL: the ON arm DROPPED $(printf '%s\n' "$lost" | wc -l) hypothes(es) the OFF arm asserted."
    echo "  This is a behavioural regression the per-stratum VERIFIED comparison cannot see, because"
    echo "  a dropped entry is absent from both sides of it. The flip is not purely additive:"
    printf '%s\n' "$lost" | sed 's/^/      /'
    RC=1
else
    echo "  PASS: purely additive -- every OFF hypothesis survives, verbatim."
    echo "  added by the flag ($(printf '%s\n' "$added" | grep -c . ) entr(ies)):"
    printf '%s\n' "$added" | sed 's/^/      /'
fi

echo
echo "  adjudicating the OFF arm"
R_OFF="$OUT/verify_off"; mkdir -p "$R_OFF"
python "$VERIFY" --pool-tsv "$OUT/pool_off.tsv" \
    --bam "$BAM" --ref "$REF" --rep-lib-dir "$REPLIB" \
    --threads "$CORES" --out "$R_OFF" > "$OUT/verify_off.report.txt" 2>"$OUT/verify_off.err.txt"
sed -n '/PRIMARY STRATUM/,$p' "$OUT/verify_off.report.txt" | sed 's/^/      /'

echo
echo "  adjudicating the ON arm, against the OFF arm as baseline"
R_ON="$OUT/verify_on"; mkdir -p "$R_ON"
set +e
python "$VERIFY" --pool-tsv "$OUT/pool_on.tsv" \
    --bam "$BAM" --ref "$REF" --rep-lib-dir "$REPLIB" \
    --threads "$CORES" --out "$R_ON" \
    --baseline "$R_OFF/verdicts.tsv" > "$OUT/verify_on.report.txt" 2>"$OUT/verify_on.err.txt"
rc2=$?
set -e
sed -n '/PRIMARY STRATUM/,$p' "$OUT/verify_on.report.txt" | sed 's/^/      /'
if [ "$rc2" -eq 3 ]; then
    echo "  FATAL: the two arms were adjudicated at different operating points. Both runs use"
    echo "  this script's defaults, so this can only mean the verifier's sweep defaults changed"
    echo "  between them -- impossible in one job, so read it as a bug."
    RC=1
elif [ "$rc2" -ne 0 ]; then
    echo "  FAIL: flipping the flag LOSES verified entries in at least one stratum. Per sec 10.1"
    echo "  the flip's justification is structural correctness AND no precision regression --"
    echo "  so this is a blocker for the flip, not a number to weigh against the 0 -> 241."
    RC=1
else
    echo "  PASS: no stratum loses a verified entry when the flag goes on."
fi

# --- the precision ranges, side by side, REPORTED and not gated -------------------------
echo
echo "--- the precision RANGES, reported (sec 3.1: never quote a point) ---"
for n in repro_3arm verify_off verify_on; do
    printf '  %-12s %s\n' "$n" \
        "$(grep -m1 'precision over the FULL cross product' "$OUT/$n.report.txt" || echo 'not reported')"
done
echo "  These are NOT the gate. The gate is the per-stratum VERIFIED counts above. A range is"
echo "  quoted here because a reader will want it and because omitting it invites a point"
echo "  estimate to be quoted instead."

echo
echo "=============================================================="
if [ "$RC" -eq 0 ]; then
    echo "NO-REGRESSION: PASS -- the flip is clear to proceed (sec 8 step 3)"
    echo "  MARGIN: PROVISIONAL. Not validated until an arm exists that is neither master nor"
    echo "  xTea2. And PASS means 'no stratum lost a verified entry', which is weaker than"
    echo "  'precision did not fall'."
else
    echo "NO-REGRESSION: FAIL -- do NOT flip the default"
fi
echo "  artifacts: $OUT"
echo "=============================================================="
exit "$RC"
