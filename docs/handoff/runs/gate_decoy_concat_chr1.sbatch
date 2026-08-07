#!/bin/bash
#SBATCH -J concat_chr1
#SBATCH -A park_contrib
#SBATCH -p priopark
#SBATCH -t 3:00:00
#SBATCH -c 12
#SBATCH --mem 64G
#SBATCH -o /n/scratch/users/c/cos689/xtea2_testdata/gate_decoy_concat_chr1_%j.out
# =============================================================================
# The decoy-consensus CONCATENATION (`xtea/td_decoy_library.py`, `run_xtea --td-decoy-concat`),
# gated on HG002 chr1. FOUR ARMS -- three from 2026-08-05, and a fourth added 2026-08-07 when
# the flag became the default:
#
#   BEFORE        the parent commit (457eb32).
#   AFTER_OFF     this commit with `--no-td-decoy-concat` -> must be BYTE-IDENTICAL to BEFORE.
#                 The containment claim. The seam is a `cat` argument list, and with the flag off
#                 it must build the same command it always did -- provable by inspection, pinned
#                 by a unit test, and still checked here at whole-output-tree scale because a
#                 predicate differential says nothing about the seam.
#   AFTER_ON      `--td-decoy-concat` -> THE MEASUREMENT.
#   AFTER_DEFAULT no flag at all -> THE SHIPPING PATH, and the only arm that exercises it.
#
# WHY THE FOURTH ARM IS NOT REDUNDANT, and it is the same argument `gate_polya_cutoff_chr1.sbatch`
# had to make for its own arm 4. Arms 1-3 all name a flag explicitly. If the 2026-08-07 default
# never reaches the step -- a resolver returning the wrong branch, an option parsed and dropped,
# a banner printed from one variable while the env is written from another -- **all three still
# pass, unchanged**. A gate whose every arm is explicit cannot see a broken default. So check 1b
# demands AFTER_DEFAULT == AFTER_ON byte-for-byte, and check 1c demands AFTER_DEFAULT != AFTER_OFF
# so that the agreement in 1b is not the vacuous kind that a wholly inert flag would also produce.
#
# THE PREDICTION THIS GATE EXISTS TO STATE BEFORE THE RUN, and it is sharper than the structural
# fix's was, because it is a claim about the DATA and not only about the code:
#
#   **L1 must be BYTE-IDENTICAL in AFTER_ON too, and not because L1 is insensitive to decoys.**
#
# `hg38_FL_L1_flanks_3k.fa` already carries all 45 decoy records as its final 45 -- the block is
# byte-identical to the one `td_decoy_library` rebuilds from `<rep_lib>/consensus/` -- so
# `write_missing_decoys` selects nothing for it and the bwa reference is the same file it was.
# hg38's SVA library carries **none**, so SVA gains all 45. Both halves are checked directly on
# the reference the step built (check 3a), which is a stronger statement than the VCF: it says
# the input to bwa was or was not the same bytes, before any of the step's own logic runs.
#
# WHAT chr1 CAN AND CANNOT SAY. chr1's gate B SVA survivors are 0, so AFTER_ON's SVA call delta
# here may well be zero -- the same limitation the structural gate's header records, for the same
# reason. What this leg proves is (a) the flag is inert when off over the whole tree, (b) it is
# WIRED -- it reaches the subprocess and changes the reference, (c) the L1 no-op is real at
# pipeline scale and not just in a unit test, and (d) nothing crashes on a bwa reference with 45
# extra records. The population is WGS's (`gate_decoy_concat_wgs.sbatch`). A zero SVA delta here
# is not evidence the concatenation does nothing.
#
# The two characterised exceptions from the parent gates are unchanged: `cns_manifest.txt` embeds
# each arm's output paths, `scan_disc_pos.tsv` is a permutation (the multithreaded scan writes in
# completion order).
#
# EVERY ARM RUNS WITH `--td-decoy-structural` EXPLICIT. It has been the default since 2026-08-05,
# so passing nothing would also arm it -- but then the banner cannot say which channel did, and
# these two decoy flags are easy to confuse in a log. Saying it costs nothing and makes each arm
# self-describing.
#
# USAGE  sbatch gate_decoy_concat_chr1.sbatch <before-worktree> <after-worktree>
#   e.g. /n/scratch/users/c/cos689/xtea2_wt/td_concat_457eb32 and .../td_concat_ba83d7e
#   Neither has a built xtea_scan; XTEA2_SCAN_BIN defaults to the td_decoy_55055ba worktree's,
#   whose src/ is byte-identical across this whole range (`git diff --name-only 11b0e61 HEAD --
#   src/` is empty), so one binary removes the C++ as a variable rather than introducing one.
# =============================================================================
set -euo pipefail

SCRATCH=/n/scratch/users/c/cos689/xtea2_testdata
WT_BEFORE="${1:?usage: sbatch gate_decoy_concat_chr1.sbatch <before-wt> <after-wt>}"
WT_AFTER="${2:?usage: sbatch gate_decoy_concat_chr1.sbatch <before-wt> <after-wt>}"
OUT="${XTEA2_OUT:-$SCRATCH/gate_decoy_concat_chr1_${SLURM_JOB_ID:-manual}}"
TESTDATA=/n/data1/hms/dbmi/park/corinne/claude_sandbox/test_data/HG002
BAM="$TESTDATA/chr1/HG002.chr1.bam"
SAMPLE=HG002.chr1
TYPES="${XTEA2_TYPES:-L1,Alu,SVA}"
REF=/n/data1/hms/dbmi/park-smaht_dac/ref/GRCh38_no_alt/hg38_no_alt.fa
CORES="${SLURM_CPUS_PER_TASK:-12}"
SCAN_BIN="${XTEA2_SCAN_BIN:-/n/scratch/users/c/cos689/xtea2_wt/td_decoy_55055ba/src/csrc/xtea_scan}"
REPLIB=/home/cos689/xTea/rep_lib_annotation

source /n/app/conda/miniforge3/24.11.3-0/etc/profile.d/conda.sh
conda activate xtea2
export PATH=/home/cos689/software/bin:$PATH
export XTEA_BASIC_INFO_SEED=42
test "$PYTHONNOUSERSITE" = 1
# An INHERITED opt-in would silently turn BEFORE and AFTER_OFF into the ON arm and make this gate
# pass while proving nothing. Both decoy vars cleared explicitly, once; each arm passes its own.
unset XTEA_TD_DECOY_CONCAT
unset XTEA_TD_DECOY_STRUCTURAL

for wt in "$WT_BEFORE" "$WT_AFTER"; do
    [ -d "$wt" ] || { echo "FATAL: worktree $wt missing"; exit 1; }
    [ -z "$(git -C "$wt" status --porcelain)" ] || { echo "FATAL: $wt is dirty"; exit 1; }
done
[ -x "$SCAN_BIN" ] || { echo "FATAL: $SCAN_BIN not built"; exit 1; }
[ -f "$BAM" ]      || { echo "FATAL: $BAM missing"; exit 1; }
# The code under test must be present BEFORE the compute, not after it.
if [ ! -f "$WT_AFTER/xtea/td_decoy_library.py" ]; then
    echo "FATAL: AFTER has no td_decoy_library; wrong commit"; exit 1
fi
if [ -f "$WT_BEFORE/xtea/td_decoy_library.py" ]; then
    echo "FATAL: BEFORE already has td_decoy_library; wrong commit"; exit 1
fi
for flag in 'td-decoy-concat' 'no-td-decoy-concat'; do
    if ! grep -q -- "$flag" "$WT_AFTER/xtea/run_xtea.py"; then
        echo "FATAL: AFTER has no --$flag; wrong commit"; exit 1
    fi
done
# The recipe must be reconstructible on THIS machine's rep_lib before an hour of compute, and
# the probe is the same code the step will run.
if ! python "$WT_AFTER/docs/handoff/runs/probe_decoy_library.py" --rep-lib-dir "$REPLIB" \
        > "${TMPDIR:-/tmp}/decoy_census.$$" 2>&1; then
    echo "FATAL: the decoy recipe does not reconstruct from $REPLIB:"; cat "${TMPDIR:-/tmp}/decoy_census.$$"
    exit 1
fi

mkdir -p "$OUT"
cp "${TMPDIR:-/tmp}/decoy_census.$$" "$OUT/decoy_census.txt"; rm -f "${TMPDIR:-/tmp}/decoy_census.$$"

echo "=============================================================="
echo "GATE: the decoy-consensus concatenation  --  HG002 chr1, types=$TYPES"
echo "  job     : ${SLURM_JOB_ID:-<none>} on $(hostname)  cores=$CORES"
echo "  BEFORE  : $(git -C "$WT_BEFORE" rev-parse --short HEAD)  $WT_BEFORE"
echo "  AFTER   : $(git -C "$WT_AFTER"  rev-parse --short HEAD)  $WT_AFTER"
echo "  scan bin: $SCAN_BIN (shared by all four arms)"
echo "  out     : $OUT"
echo "  started : $(date '+%F %T')"
echo "=============================================================="
echo
echo "--- the library census this gate's prediction rests on ---"
sed 's/^/  /' "$OUT/decoy_census.txt"

run_arm() {   # $1 = arm name, $2 = worktree, $3 = extra flags
    local arm=$1
    local wt=$2
    local extra=$3
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
# BEFORE passes nothing because its commit has no such flag to pass; AFTER_OFF says it
# explicitly, so the banner can tell "the default" from "the opt-out" and this gate cannot end
# up comparing the change against itself.
run_arm BEFORE        "$WT_BEFORE" ""                      || RC=1
run_arm AFTER_OFF     "$WT_AFTER"  "--no-td-decoy-concat"  || RC=1
run_arm AFTER_ON      "$WT_AFTER"  "--td-decoy-concat"     || RC=1
# ARM 4, added 2026-08-07 WITH THE DEFAULT FLIP, and it is the only arm that tests the thing
# that actually ships. Arms 1-3 all name a flag, so all three would pass unchanged if the new
# default never reached the step at all -- a typo in the resolver, an option consumed and
# dropped, a banner printed from one variable while the env is written from another. This is
# the same reason `gate_polya_cutoff_chr1.sbatch` needed its arm 4: a gate whose every arm is
# explicit cannot see a broken default. Check 1b demands DEFAULT == AFTER_ON byte-for-byte,
# and check 1c demands DEFAULT != AFTER_OFF so that agreement is not vacuous.
run_arm AFTER_DEFAULT "$WT_AFTER"  ""                      || RC=1
[ "$RC" -eq 0 ] || { echo; echo "GATE concat (chr1): FAIL (an arm did not complete)"; exit 1; }

A="$OUT/BEFORE/$SAMPLE"
B="$OUT/AFTER_OFF/$SAMPLE"
C="$OUT/AFTER_ON/$SAMPLE"
D="$OUT/AFTER_DEFAULT/$SAMPLE"

echo
echo "--- check 0: the arms are the code and the configuration they claim ---"
b_mod=$([ -f "$WT_BEFORE/xtea/td_decoy_library.py" ] && echo 1 || echo 0)
a_mod=$([ -f "$WT_AFTER/xtea/td_decoy_library.py" ] && echo 1 || echo 0)
echo "  td_decoy_library present: BEFORE=$b_mod AFTER=$a_mod (expect 0 and 1)"
{ [ "$b_mod" -eq 0 ] && [ "$a_mod" -eq 1 ]; } || { echo "    FATAL: wrong commits"; RC=1; }
ON_BANNER='\[decoy-concat\] the transduction bwa reference GAINS'
OFF_BANNER='\[decoy-concat\] the transduction bwa reference is the packaged flank library'
n=$(grep -ac "$ON_BANNER" "$OUT/BEFORE/driver.log" || true)
echo "  [decoy-concat] ON banner in BEFORE: $n (expect 0 -- its commit has no such banner)"
[ "$n" -eq 0 ] || { echo "    FATAL: BEFORE is not the parent commit"; RC=1; }
# Checked POSITIVELY as well as negatively: an off-arm gate whose only evidence is an absent
# line passes when the run dies before it gets there.
n=$(grep -ac "$OFF_BANNER" "$OUT/AFTER_OFF/driver.log" || true)
echo "  [decoy-concat] OFF banner in AFTER_OFF: $n (expect 1)"
[ "$n" -eq 1 ] || { echo "    FATAL: AFTER_OFF's opt-out did not reach the driver"; RC=1; }
n=$(grep -ac "$ON_BANNER" "$OUT/AFTER_ON/driver.log" || true)
echo "  [decoy-concat] ON banner in AFTER_ON: $n (expect 1)"
[ "$n" -eq 1 ] || { echo "    FATAL: AFTER_ON did not run with the concatenation ON"; RC=1; }
# And the STEP's own line, which is the one that separates "armed" from "changed anything".
echo "  the step's per-family append counts in AFTER_ON:"
grep -a 'td-decoy-concat' "$OUT/AFTER_ON"/*/*.steps.log "$OUT/AFTER_ON/driver.log" 2>/dev/null |
    sed 's/^/    /' | cut -c1-160 || echo "    (none found)"

# The whole-output-tree comparison, ONE definition. It was inline in check 1 until arm 4 needed
# the identical comparison for DEFAULT vs AFTER_ON; copying it would have been the fourth
# private copy of a shipping rule this project has paid for (plan sec 10), and the copy would
# have had to re-derive the two tolerated cases below to mean anything.
#   $1 arm name A   $2 arm name B      -> compares $OUT/$1/$SAMPLE against $OUT/$2/$SAMPLE
# Sets n_same / n_tol / n_diff / n_missing. Returns 1 if the trees differ.
compare_trees() {
    local arm_a=$1 arm_b=$2
    local dir_a="$OUT/$arm_a/$SAMPLE" dir_b="$OUT/$arm_b/$SAMPLE"
    local f rel
    n_same=0; n_diff=0; n_missing=0; n_tol=0
    while IFS= read -r f; do
        rel=${f#"$dir_a/"}
        case "$rel" in
            tmp/*|*/tmp/*|tmp_basic/*|pub_clip/*) continue ;;
            *.log|step_timings.tsv) continue ;;
        esac
        if [ ! -f "$dir_b/$rel" ]; then
            echo "  MISSING in $arm_b: $rel"; n_missing=$((n_missing+1)); continue; fi
        if cmp -s "$f" "$dir_b/$rel"; then n_same=$((n_same+1)); continue; fi
        case "$rel" in
          cns_manifest.txt)
            # Holds absolute paths, so it names its own arm. Identical modulo that prefix.
            if cmp -s <(sed "s|$OUT/$arm_a/||g" "$f") <(sed "s|$OUT/$arm_b/||g" "$dir_b/$rel"); then
                n_tol=$((n_tol+1)); echo "  $rel: identical after stripping the arm prefix"
            else
                echo "  DIFFERS (beyond the path prefix): $rel"; n_diff=$((n_diff+1))
            fi ;;
          scan_disc_pos.tsv)
            # Emitted in thread-completion order, so a permutation is not a difference.
            if cmp -s <(sort "$f") <(sort "$dir_b/$rel"); then
                n_tol=$((n_tol+1)); echo "  $rel: identical when sorted ($(wc -l < "$f") rows)"
            else
                echo "  DIFFERS (not a permutation): $rel"; n_diff=$((n_diff+1))
            fi ;;
          *)
            echo "  DIFFERS: $rel"; n_diff=$((n_diff+1))
            diff <(grep -v '^#' "$f" | cut -f1,2) <(grep -v '^#' "$dir_b/$rel" | cut -f1,2) \
                | head -20 || true ;;
        esac
    done < <(find "$dir_a" -type f | sort)
    echo "  identical: $n_same    tolerated: $n_tol    DIFFERING: $n_diff    missing: $n_missing"
    [ "$n_diff" -eq 0 ] && [ "$n_missing" -eq 0 ]
}

echo
echo "--- check 1: AFTER_OFF is byte-identical to BEFORE, over the whole output tree ---"
compare_trees BEFORE AFTER_OFF || RC=1
[ "$n_same" -ge 95 ] || { echo "    FATAL: only $n_same files compared, expected >=95"; RC=1; }

echo
echo "--- check 1b: AFTER_DEFAULT is byte-identical to AFTER_ON -- THE DEFAULT FLIP ITSELF ---"
echo "  Arms 1-3 each name a flag, so all three pass unchanged even if the 2026-08-07 default"
echo "  never reaches the step. This is the arm that fails when it does not."
compare_trees AFTER_ON AFTER_DEFAULT || {
    echo "    FATAL: passing NO flag does not reproduce --td-decoy-concat. The default flip"
    echo "    did not reach the step -- check resolve_decoy_concat and the env write."
    RC=1; }
[ "$n_same" -ge 95 ] || { echo "    FATAL: only $n_same files compared, expected >=95"; RC=1; }

echo
echo "--- check 1c: the DEFAULT arm's SVA bwa REFERENCE differs from AFTER_OFF's ---"
echo "  This is what stops 1b being vacuous, and it must be asked of the REFERENCE and not of"
echo "  the output tree. chr1 has 0 SVA gate-B survivors, so the concatenation is EXPECTED to"
echo "  move no output file here -- this header says so, check 3b measures it as L1=0 Alu=0"
echo "  SVA=0, and the 3-arm gate 49879284 passed with exactly that. An output-tree differential"
echo "  would therefore fail a working flag. The reference is the artifact the flag actually"
echo "  changes at this scale: 45 records appear in it whether or not any call moves."
r_off_sva="$B/SVA/tmp/transduction/all_with_polymerphic_flanks.fa"
r_def_sva="$D/SVA/tmp/transduction/all_with_polymerphic_flanks.fa"
if [ ! -f "$r_off_sva" ] || [ ! -f "$r_def_sva" ]; then
    echo "  FATAL: an SVA reference is missing (need --keep-intermediates) -- CANNOT CHECK"; RC=1
else
    n_off_sva=$(grep -c '^>' "$r_off_sva"); n_def_sva=$(grep -c '^>' "$r_def_sva")
    echo "  SVA reference records: --no-td-decoy-concat=$n_off_sva  DEFAULT=$n_def_sva  delta=$(( n_def_sva - n_off_sva ))"
    if [ "$(( n_def_sva - n_off_sva ))" -eq 45 ]; then
        echo "  PASS: the DEFAULT gains exactly the 45 decoy records the opt-out does not, so the"
        echo "  new default is live and check 1b has content."
    else
        echo "  FATAL: the DEFAULT arm's SVA reference is not 45 records longer than the opt-out's."
        echo "  Either the default did not arm the concatenation, or --no-td-decoy-concat no longer"
        echo "  disarms it. Check 1b's agreement with AFTER_ON would then prove nothing."
        RC=1
    fi
fi

echo
echo "--- check 2: no new or removed RESULT files in either after arm ---"
# The ON arm DOES write one new intermediate -- all_with_polymerphic_flanks.fa.decoy_only.fa --
# but it lands under tmp/, which this loop skips, so 0 remains the right expectation for both.
for arm in AFTER_OFF AFTER_ON; do
    n_new=0
    while IFS= read -r f; do
        rel=${f#"$OUT/$arm/$SAMPLE/"}
        case "$rel" in
            tmp/*|*/tmp/*|tmp_basic/*|pub_clip/*|*.log|step_timings.tsv) continue ;;
        esac
        [ -f "$A/$rel" ] && continue
        echo "  UNEXPECTED new file in $arm: $rel"; n_new=$((n_new+1))
    done < <(find "$OUT/$arm/$SAMPLE" -type f | sort)
    echo "  $arm: unexpected new result files: $n_new (expect 0)"
    [ "$n_new" -eq 0 ] || RC=1
done

echo
echo "--- check 3a: THE BWA REFERENCE ITSELF -- L1 unchanged, SVA longer by exactly 45 ---"
# The sharpest form of the prediction, and it is upstream of every gate: what did bwa actually
# index? Both L1 references must be the same BYTES with the flag on; SVA's must gain exactly the
# 45 records, in the block position the code claims (between the library and the polymorphic
# flanks). A VCF comparison cannot distinguish "the reference was identical" from "the reference
# differed and nothing survived to the call set".
for fam in L1 SVA; do
    r_off="$B/$fam/tmp/transduction/all_with_polymerphic_flanks.fa"
    r_on="$C/$fam/tmp/transduction/all_with_polymerphic_flanks.fa"
    if [ ! -f "$r_off" ] || [ ! -f "$r_on" ]; then
        echo "  $fam: reference not kept (need --keep-intermediates) -- CANNOT CHECK"; RC=1; continue
    fi
    n_off=$(grep -c '^>' "$r_off"); n_on=$(grep -c '^>' "$r_on")
    echo "  $fam: records OFF=$n_off ON=$n_on  delta=$(( n_on - n_off ))"
    if [ "$fam" = L1 ]; then
        if cmp -s "$r_off" "$r_on"; then
            echo "    PASS: byte-identical -- the packaged L1 library already carries all 45"
        else
            echo "    FAIL: L1's bwa reference MOVED. hg38_FL_L1_flanks_3k.fa carries the 45-record"
            echo "    block already, so write_missing_decoys must select nothing for it."
            diff <(grep '^>' "$r_off") <(grep '^>' "$r_on") | head -20 || true
            RC=1
        fi
    else
        [ "$(( n_on - n_off ))" -eq 45 ] || {
            echo "    FAIL: expected exactly 45 new records on hg38's SVA library"; RC=1; }
        # and they must be the recipe's, in order, in the block position
        python - "$r_on" "$WT_AFTER" "$REPLIB" <<'EOF' || RC=1
import os, sys
r_on, wt, replib = sys.argv[1], sys.argv[2], sys.argv[3]
sys.path.insert(0, os.path.join(wt, "xtea"))
import td_decoy_library as tdl
block = tdl.decoy_block_text(tdl.build_decoy_block(os.path.join(replib, "consensus")))
text = open(r_on).read()
if block not in text:
    print("    FAIL: the 45-record block is not present verbatim in the ON reference")
    raise SystemExit(1)
i = text.index(block)
tail = text[i + len(block):]
if not tail.lstrip().startswith(">"):
    print("    FAIL: the block is not followed by more records")
    raise SystemExit(1)
first_after = tail.split("\n", 1)[0]
print("    PASS: the block appears verbatim, followed by {0}".format(first_after[:60]))
print("    (expected the polymorphic flanks, i.e. a ~polymerphic~ record)")
EOF
    fi
done

echo
echo "--- check 3b: THE PREDICTION on the output tree -- AFTER_ON moves SVA and nothing else ---"
n_l1_diff=0; n_alu_diff=0; n_sva_diff=0; n_other_diff=0
while IFS= read -r f; do
    rel=${f#"$A/"}
    case "$rel" in
        tmp/*|*/tmp/*|tmp_basic/*|pub_clip/*|*.log|step_timings.tsv) continue ;;
        cns_manifest.txt|scan_disc_pos.tsv) continue ;;
    esac
    [ -f "$C/$rel" ] || continue
    cmp -s "$f" "$C/$rel" && continue
    case "$rel" in
        L1/*|*_LINE1.vcf|*LINE1*)  echo "  L1 MOVED: $rel";  n_l1_diff=$((n_l1_diff+1)) ;;
        Alu/*|*_ALU.vcf|*ALU*)     echo "  Alu MOVED: $rel"; n_alu_diff=$((n_alu_diff+1)) ;;
        SVA/*|*_SVA.vcf|*SVA*)     n_sva_diff=$((n_sva_diff+1)) ;;
        *)                         echo "  cross-family file moved: $rel"
                                   n_other_diff=$((n_other_diff+1)) ;;
    esac
done < <(find "$A" -type f | sort)
echo "  files moved with the concatenation ON: L1=$n_l1_diff Alu=$n_alu_diff SVA=$n_sva_diff other=$n_other_diff"
[ "$n_l1_diff" -eq 0 ] || { echo "    FAIL: L1 output moved, but its bwa reference did not (3a)."
                            echo "    Something other than the concatenation changed."; RC=1; }
[ "$n_alu_diff" -eq 0 ] || { echo "    FAIL: Alu moved, and it never runs the transduction step"; RC=1; }

echo
echo "--- check 4: the clip channel, per arm -- the 0-of-128,262 asymmetry, before and after ---"
# HANDOFF_2026-08-05_gate_a_measured.md §2.4: L1 counts 3,992 of 99,245 qualified clip records as
# decoy-consensus support and SVA counts 0 of 128,262, because its library has no decoy records.
# This is the direct before/after on that number, from the run's own surviving clip SAM.
for arm in BEFORE AFTER_ON; do
    echo "  --- $arm ---"
    python "$WT_AFTER/docs/handoff/runs/probe_clip_gate.py" \
        --run-dir "$OUT/$arm/$SAMPLE" > "$OUT/probe_clip_gate_$arm.txt" 2>&1 &&
        grep -aE 'decoy|family|survivors' "$OUT/probe_clip_gate_$arm.txt" | sed 's/^/    /' | head -20 ||
        { echo "    (the probe exited nonzero -- reporting, not failing)"
          tail -15 "$OUT/probe_clip_gate_$arm.txt" | sed 's/^/    /'; }
done

echo
echo "--- check 5: picked sources and ledger populations, per arm ---"
for arm in BEFORE AFTER_OFF AFTER_ON; do
    for fam in L1 SVA; do
        led="$OUT/$arm/$SAMPLE/$fam/candidate_disc_filtered_cns2.txt.all_non_sibling_td.txt.decision_ledger.tsv"
        src="$OUT/$arm/$SAMPLE/$fam/candidate_disc_filtered_cns2.txt.all_non_sibling_td.txt.unique_trsdct_disc_only"
        n_src=$([ -f "$src" ] && wc -l < "$src" || echo NA)
        if [ -f "$led" ]; then
            echo "  $arm/$fam: picked sources=$n_src  ledger rows=$(( $(wc -l < "$led") - 1 ))"
            echo "    by reason: $(awk -F'\t' 'NR>1 {c[$4]++} END {for (r in c) printf "%s=%d ", r, c[r]}' "$led")"
        else
            echo "  $arm/$fam: picked sources=$n_src  (no ledger)"
        fi
    done
done

echo
echo "--- check 6: call counts (BEFORE and AFTER_OFF must match; AFTER_ON is the delta) ---"
for arm in BEFORE AFTER_OFF AFTER_ON; do
    d="$OUT/$arm/$SAMPLE"
    l1=$(grep -vc '^#' "$d/merged/${SAMPLE}_LINE1.vcf" 2>/dev/null || echo 0)
    alu=$(grep -vc '^#' "$d/merged/${SAMPLE}_ALU.vcf" 2>/dev/null || echo 0)
    sva=$(grep -vc '^#' "$d/merged/${SAMPLE}_SVA.vcf" 2>/dev/null || echo 0)
    mrg=$(grep -vc '^#' "$d/${SAMPLE}.vcf" 2>/dev/null || echo 0)
    echo "  $arm: L1=$l1 Alu=$alu SVA=$sva merged=$mrg"
done

echo
echo "=============================================================="
if [ "$RC" -eq 0 ]; then
    echo "GATE concat (chr1): PASS"
else
    echo "GATE concat (chr1): FAIL"
fi
echo "  finished: $(date '+%F %T')"
echo "=============================================================="
exit "$RC"
