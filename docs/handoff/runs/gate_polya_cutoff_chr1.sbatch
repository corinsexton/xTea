#!/bin/bash
#SBATCH -J polya_chr1
#SBATCH -c 12
#SBATCH -t 4:00:00
#SBATCH --mem=32G
#SBATCH -p short
#SBATCH --account=park
#SBATCH -o /n/scratch/users/c/cos689/xtea2_testdata/gate_polya_cutoff_chr1_%j.out
#SBATCH --mail-type=NONE
# =============================================================================
# Decoupling gate [A]'s polyA threshold from its clip threshold. FOUR arms, and the fourth is
# the one that makes the other three mean anything.
#
# THE CHANGE. `n_polya_cutoff` has always been DERIVED as `n_clip_cutoff / 2 - 1` inside
# TransductionParams, so the two thresholds at the gate that does 97% of the step's rejecting
# could not move independently -- a clip-support sweep is single-variable only up to 4, and at
# `--nclip 5` the polyA arm wakes and rejects 286 (L1) / 364 (SVA) loci the shipping
# configuration keeps. It is now `run_xtea --td-polya-cutoff` / `XTEA_TD_POLYA_CUTOFF`,
# **absent-means-DERIVED on both entry points**.
#
# THE PREDICTION, STATED BEFORE THE RUN, and it is unusually sharp for a call-moving flag:
#
#   1. **BEFORE == AFTER_DEFAULT, byte-identical over the whole output tree.** Containment: the
#      new parameter defaults to the frozen derivation, so a run with no flag is the old code.
#   2. **AFTER_DEFAULT == AFTER_EXPLICIT_HALF, byte-identical.** At HG002's operating point
#      (`n_clip_cutoff = 3`) the derived value is **0.5**, so passing `--td-polya-cutoff 0.5`
#      asks for exactly what the derivation produces. This is the decoupling's no-op arm.
#   3. **AFTER_RAISED MUST DIFFER, AND THE LEDGER SAYS WHERE.** `--td-polya-cutoff 2` must make
#      the `low_polya_clip` reason FIRE.
#
# **WHY ARM 4 EXISTS, AND WHY THIS GATE IS WORTHLESS WITHOUT IT.** Arms 1-3 predict
# byte-identity. Three byte-identical arms are equally consistent with "the flag works and is
# inert at this operating point" and with "**the flag never reached the gate at all**" -- a
# typo in the env var name, a resolver returning None, a params object built before the value
# was read. That is precisely the confusion the decoy default flip had to be rescued from, and
# the reason the concat step reports its append COUNT rather than only its policy banner. Arm 4
# is the positive control: it demands that the threshold BITE, at a named reason in the ledger.
#
# AND THE REASON IT BITES ON IS A DEAD PARAMETER COMING ALIVE. `low_polya_clip` is
# **unreachable** at the shipping operating point -- `m_polyA`'s record is created and one of its
# six counters incremented in the same block, so every entry has `sum(rcd[:6]) >= 1`, measured 0
# zero-sum entries of 4,604 (L1) and 5,261 (SVA). So its ledger count is **0 in arms 1-3 and must
# be nonzero in arm 4**, which is a cleaner positive control than a VCF delta: at chr1 the whole
# call set is L1 1 / Alu 14 / SVA 2, far too small to reliably register a gate [A] change.
#
# WHAT IS NOT CLAIMED. That raising the threshold is a good idea. This gate tests that the two
# thresholds are SEPARABLE and that separating them changes nothing until asked. Whether the
# polyA signature should be a veto at all is gate [A]'s other semantic decision and is not this
# change (`HANDOFF_2026-08-05_gate_a_measured.md` sec 2.2).
#
# Frozen worktrees, not the live checkout -- run_xtea spawns a fresh x_TEA_main per step, so a
# job sharing the checkout picks up mid-run edits.
#
# USAGE  sbatch gate_polya_cutoff_chr1.sbatch <before-wt> <after-wt>
# =============================================================================
set -euo pipefail

SCRATCH=/n/scratch/users/c/cos689/xtea2_testdata
WT_BEFORE="${1:?usage: sbatch gate_polya_cutoff_chr1.sbatch <before-wt> <after-wt>}"
WT_AFTER="${2:?usage: sbatch gate_polya_cutoff_chr1.sbatch <before-wt> <after-wt>}"
OUT="${XTEA2_OUT:-$SCRATCH/gate_polya_cutoff_chr1_${SLURM_JOB_ID:-manual}}"
TESTDATA=/n/data1/hms/dbmi/park/corinne/claude_sandbox/test_data/HG002
BAM="$TESTDATA/chr1/HG002.chr1.bam"
SAMPLE=HG002.chr1
TYPES="${XTEA2_TYPES:-L1,Alu,SVA}"
REF=/n/data1/hms/dbmi/park-smaht_dac/ref/GRCh38_no_alt/hg38_no_alt.fa
CORES="${SLURM_CPUS_PER_TASK:-12}"
SCAN_BIN="${XTEA2_SCAN_BIN:-/n/scratch/users/c/cos689/xtea2_wt/td_decoy_55055ba/src/csrc/xtea_scan}"
REPLIB=/home/cos689/xTea/rep_lib_annotation
# The value arm 4 raises the threshold to. Must exceed 1 to be reachable at all -- every polyA
# record carries a count of at least 1 by construction, so 0.5 and 1.0 cannot fire.
RAISED="${XTEA2_POLYA_RAISED:-2}"

source /n/app/conda/miniforge3/24.11.3-0/etc/profile.d/conda.sh
conda activate xtea2
export PATH=/home/cos689/software/bin:$PATH
export XTEA_BASIC_INFO_SEED=42
test "$PYTHONNOUSERSITE" = 1
# An INHERITED value would silently turn BEFORE and AFTER_DEFAULT into a raised arm and make
# this gate pass while proving nothing. Cleared explicitly, once; each arm passes its own.
unset XTEA_TD_POLYA_CUTOFF
unset XTEA_TD_DECOY_CONCAT
unset XTEA_TD_DECOY_STRUCTURAL

for wt in "$WT_BEFORE" "$WT_AFTER"; do
    [ -d "$wt" ] || { echo "FATAL: worktree $wt missing"; exit 1; }
    [ -z "$(git -C "$wt" status --porcelain)" ] || { echo "FATAL: $wt is dirty"; exit 1; }
done
[ -x "$SCAN_BIN" ] || { echo "FATAL: $SCAN_BIN not built"; exit 1; }
[ -f "$BAM" ]      || { echo "FATAL: $BAM missing"; exit 1; }
# The code under test must be present BEFORE the compute, not after it -- and its ABSENCE from
# the BEFORE arm asserted just as positively, or this gate compares the change against itself.
if ! grep -q -- 'td-polya-cutoff' "$WT_AFTER/xtea/run_xtea.py"; then
    echo "FATAL: AFTER has no --td-polya-cutoff; wrong commit"; exit 1
fi
if grep -q -- 'td-polya-cutoff' "$WT_BEFORE/xtea/run_xtea.py"; then
    echo "FATAL: BEFORE already has --td-polya-cutoff; wrong commit"; exit 1
fi
if ! grep -q 'ENV_TD_POLYA_CUTOFF' "$WT_AFTER/xtea/global_values.py"; then
    echo "FATAL: AFTER has no ENV_TD_POLYA_CUTOFF; wrong commit"; exit 1
fi

mkdir -p "$OUT"

echo "=============================================================="
echo "GATE: gate [A]'s polyA cutoff, decoupled  --  HG002 chr1, types=$TYPES"
echo "  job     : ${SLURM_JOB_ID:-<none>} on $(hostname)  cores=$CORES"
echo "  BEFORE  : $(git -C "$WT_BEFORE" rev-parse --short HEAD)  $WT_BEFORE"
echo "  AFTER   : $(git -C "$WT_AFTER"  rev-parse --short HEAD)  $WT_AFTER"
echo "  raised  : --td-polya-cutoff $RAISED (arm 4, the positive control)"
echo "  scan bin: $SCAN_BIN (shared by all four arms)"
echo "  out     : $OUT"
echo "  started : $(date '+%F %T')"
echo "=============================================================="

run_arm() {   # $1 = arm name, $2 = worktree, $3 = extra flags
    local arm=$1 wt=$2 extra=$3
    local t0 wall armout
    armout="$OUT/$arm"
    rm -rf "$armout"; mkdir -p "$armout"
    echo
    echo "=== arm $arm ($(git -C "$wt" rev-parse --short HEAD)) flags: ${extra:-<none>} ==="
    t0=$(date +%s)
    ( cd "$wt" && python xtea/run_xtea.py \
      -i "$BAM" -s "$SAMPLE" -o "$armout" \
      --reference "$REF" \
      --rep-lib-dir "$REPLIB" \
      --gene-gff3 /home/cos689/data1/corinne/ref/gencode.v45.basic.annotation.gff3 \
      --blacklist /home/cos689/xTea/rep_lib_annotation/blacklist/hg38/centromere.bed \
      --types "$TYPES" \
      --scan-bin "$SCAN_BIN" \
      --bwa /home/cos689/software/bin/bwa \
      --scan-min-clip 2 \
      --mate-in-rep-families alu,sva \
      --cns-emit-scan \
      --td-gather-scan \
      --td-decoy-structural \
      --parallel-downstream \
      --keep-intermediates \
      ${extra} \
      -n "$CORES" ) > "$armout/driver.log" 2>&1 || {
          echo "FAIL: arm $arm exited $?"; tail -60 "$armout/driver.log"; return 1; }
    wall=$(( $(date +%s) - t0 ))
    printf '  wall %02d:%02d:%02d\n' $((wall/3600)) $((wall%3600/60)) $((wall%60))
}

RC=0
# BEFORE passes nothing because its commit has no such flag to pass. AFTER_DEFAULT also passes
# nothing -- and that is the point, since "absent means derived" is the property under test;
# there is no `--no-td-polya-cutoff` to say it with, because the default is a VALUE and not a
# boolean. Check 0 distinguishes the two arms by their banners instead.
run_arm BEFORE              "$WT_BEFORE" ""                              || RC=1
run_arm AFTER_DEFAULT       "$WT_AFTER"  ""                              || RC=1
run_arm AFTER_EXPLICIT_HALF "$WT_AFTER"  "--td-polya-cutoff 0.5"         || RC=1
run_arm AFTER_RAISED        "$WT_AFTER"  "--td-polya-cutoff $RAISED"     || RC=1
[ "$RC" -eq 0 ] || { echo; echo "GATE polya (chr1): FAIL (an arm did not complete)"; exit 1; }

A="$OUT/BEFORE/$SAMPLE"
B="$OUT/AFTER_DEFAULT/$SAMPLE"
C="$OUT/AFTER_EXPLICIT_HALF/$SAMPLE"
D="$OUT/AFTER_RAISED/$SAMPLE"

echo
echo "--- check 0: the arms are the code and the configuration they claim ---"
DERIVED_BANNER='\[td-polya-cutoff\] gate \[A\].s polyA threshold is DERIVED'
SET_BANNER='\[td-polya-cutoff\] gate \[A\].s polyA threshold is DECOUPLED'
n=$(grep -ac "$DERIVED_BANNER" "$OUT/BEFORE/driver.log" || true)
echo "  DERIVED banner in BEFORE: $n (expect 0 -- its commit has no such banner)"
[ "$n" -eq 0 ] || { echo "    FATAL: BEFORE is not the parent commit"; RC=1; }
# Checked POSITIVELY: an arm whose only evidence is an absent line passes when the run dies.
n=$(grep -ac "$DERIVED_BANNER" "$OUT/AFTER_DEFAULT/driver.log" || true)
echo "  DERIVED banner in AFTER_DEFAULT: $n (expect 1)"
[ "$n" -eq 1 ] || { echo "    FATAL: AFTER_DEFAULT did not take the derived path"; RC=1; }
for arm in AFTER_EXPLICIT_HALF AFTER_RAISED; do
    n=$(grep -ac "$SET_BANNER" "$OUT/$arm/driver.log" || true)
    echo "  DECOUPLED banner in $arm: $n (expect 1)"
    [ "$n" -eq 1 ] || { echo "    FATAL: $arm's flag did not reach the driver"; RC=1; }
done
# And the STEP's own line -- the driver banner says the POLICY arrived, this says what the gate
# will actually compare against. Same "armed vs changed anything" split as the concat step.
echo "  the step's own threshold lines (what gate [A] compared against):"
grep -ah 'td-polya-cutoff' "$OUT"/AFTER_*/"$SAMPLE"_*.steps.log 2>/dev/null |
    sort -u | sed 's/^/    /' | head -12 || echo "    (none found)"

# ---------------------------------------------------------------------------------------
compare_tree() {   # $1 = label, $2 = ref sample dir, $3 = other sample dir, $4 = ref arm, $5 = other arm
    local label=$1 ra=$2 rb=$3 an=$4 bn=$5
    local n_same=0 n_diff=0 n_missing=0 n_tol=0 rel
    while IFS= read -r f; do
        rel=${f#"$ra/"}
        case "$rel" in
            tmp/*|*/tmp/*|tmp_basic/*|pub_clip/*) continue ;;
            *.log|step_timings.tsv) continue ;;
        esac
        if [ ! -f "$rb/$rel" ]; then echo "  MISSING in $bn: $rel"; n_missing=$((n_missing+1)); continue; fi
        if cmp -s "$f" "$rb/$rel"; then n_same=$((n_same+1)); continue; fi
        case "$rel" in
          cns_manifest.txt)
            if cmp -s <(sed "s|$OUT/$an/||g" "$f") <(sed "s|$OUT/$bn/||g" "$rb/$rel"); then
                n_tol=$((n_tol+1))
            else echo "  DIFFERS (beyond the path prefix): $rel"; n_diff=$((n_diff+1)); fi ;;
          scan_disc_pos.tsv)
            if cmp -s <(sort "$f") <(sort "$rb/$rel"); then n_tol=$((n_tol+1))
            else echo "  DIFFERS (not a permutation): $rel"; n_diff=$((n_diff+1)); fi ;;
          *)
            echo "  DIFFERS: $rel"; n_diff=$((n_diff+1))
            diff <(grep -v '^#' "$f" | cut -f1,2) <(grep -v '^#' "$rb/$rel" | cut -f1,2) | head -10 || true ;;
        esac
    done < <(find "$ra" -type f | sort)
    echo "  [$label] identical: $n_same   tolerated: $n_tol   DIFFERING: $n_diff   missing: $n_missing"
    [ "$n_same" -ge 95 ] || { echo "    FATAL: only $n_same files compared, expected >=95"; return 2; }
    [ "$n_diff" -eq 0 ] && [ "$n_missing" -eq 0 ]
}

echo
echo "--- check 1: CONTAINMENT -- AFTER_DEFAULT byte-identical to BEFORE ---"
compare_tree "BEFORE vs AFTER_DEFAULT" "$A" "$B" BEFORE AFTER_DEFAULT || RC=1

echo
echo "--- check 2: THE NO-OP ARM -- explicit 0.5 byte-identical to the derived default ---"
# 0.5 IS what n_clip_cutoff/2-1 produces at HG002's n_clip_cutoff=3, so this arm asks for
# exactly the derivation by hand. A difference here means the override path and the derivation
# path disagree about the same number, which would be the decoupling introducing a defect.
compare_tree "AFTER_DEFAULT vs AFTER_EXPLICIT_HALF" "$B" "$C" AFTER_DEFAULT AFTER_EXPLICIT_HALF || RC=1

echo
echo "--- check 3: THE POSITIVE CONTROL -- a dead reason must come alive at $RAISED ---"
# Without this, checks 1 and 2 pass just as happily when the flag never reaches the gate.
tot_raised=0
for fam in L1 SVA; do
    for arm in AFTER_DEFAULT AFTER_EXPLICIT_HALF AFTER_RAISED; do
        led="$OUT/$arm/$SAMPLE/$fam/candidate_disc_filtered_cns2.txt.all_non_sibling_td.txt.decision_ledger.tsv"
        if [ ! -f "$led" ]; then echo "  $fam/$arm: no ledger -- CANNOT CHECK"; RC=1; continue; fi
        n=$(awk -F'\t' 'NR>1 && $4=="low_polya_clip"' "$led" | wc -l)
        echo "  $fam $arm: low_polya_clip = $n"
        case "$arm" in
          AFTER_RAISED)
            tot_raised=$(( tot_raised + n )) ;;
          *)
            [ "$n" -eq 0 ] || {
                echo "    FAIL: low_polya_clip fired at the shipping operating point, where"
                echo "    every polyA record carries a count >= 1 and the cutoff is 0.5."; RC=1; } ;;
        esac
    done
done
echo "  total low_polya_clip across both families at --td-polya-cutoff $RAISED: $tot_raised"
if [ "$tot_raised" -le 0 ]; then
    echo "    FAIL: THE FLAG DID NOT REACH THE GATE. Checks 1 and 2 are therefore vacuous --"
    echo "    they would pass identically if the threshold were never read. Do not read this"
    echo "    run as evidence that the decoupling is contained."; RC=1
else
    echo "    PASS: the threshold is live, so checks 1 and 2 mean what they say"
fi

echo
echo "--- check 4: the raised arm's call sets, REPORTED not gated ---"
# chr1's whole call set is L1 1 / Alu 14 / SVA 2, far too small to gate a gate-[A] change on.
# Reported so the direction is on the record; the ledger above is the measurement.
for arm in BEFORE AFTER_DEFAULT AFTER_EXPLICIT_HALF AFTER_RAISED; do
    v="$OUT/$arm/$SAMPLE/$SAMPLE.vcf"
    if [ -f "$v" ]; then
        echo "  $arm: $(grep -vc '^#' "$v") calls"
    else
        echo "  $arm: no merged VCF"; RC=1
    fi
done
echo "  Alu must not move in ANY arm -- it never runs the transduction step:"
for arm in AFTER_DEFAULT AFTER_EXPLICIT_HALF AFTER_RAISED; do
    a_alu=$(find "$A" -name '*_ALU.vcf' -type f | head -1)
    b_alu=$(find "$OUT/$arm/$SAMPLE" -name '*_ALU.vcf' -type f | head -1)
    if [ -n "$a_alu" ] && [ -n "$b_alu" ] && cmp -s "$a_alu" "$b_alu"; then
        echo "    $arm: Alu IDENTICAL to BEFORE"
    else
        echo "    $arm: Alu DIFFERS -- a family that never runs this gate moved"; RC=1
    fi
done

echo
if [ "$RC" -eq 0 ]; then
    echo "GATE polya (chr1): PASS -- the thresholds are separable, the default is the frozen"
    echo "  derivation byte-for-byte, and the flag is demonstrably live."
else
    echo "GATE polya (chr1): FAIL"
fi
echo "  finished: $(date '+%F %T')"
exit $RC
