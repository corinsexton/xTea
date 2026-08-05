#!/bin/bash
#SBATCH -J decoy_chr1
#SBATCH -A park_contrib
#SBATCH -p priopark
#SBATCH -t 3:00:00
#SBATCH -c 12
#SBATCH --mem 64G
#SBATCH -o /n/scratch/users/c/cos689/xtea2_testdata/gate_decoy_structural_chr1_%j.out
# =============================================================================
# The SVA decoy fix (`transduction_evidence.decoy_tags`, `run_xtea --td-decoy-structural`),
# gated on HG002 chr1. THREE ARMS, because this change ships off and the two claims about it
# are different claims that need different arms:
#
#   BEFORE     the parent commit.
#   AFTER_OFF  the fix's commit with `--no-td-decoy-structural` -> must be BYTE-IDENTICAL to
#              BEFORE. This is the containment claim, and it is the one that keeps BASELINE-0 and
#              every handoff gate frozen against it valid. The unit differential already pins
#              `decoy_tags(name, False)` against a verbatim copy of the pre-fix expression, but
#              a predicate differential says nothing about the SEAM -- that is the lesson
#              written into tests/test_transduction_evidence.py, and this arm is the seam at
#              whole-output-tree scale.
#   AFTER_ON   the fix's commit, --td-decoy-structural -> the MEASUREMENT. Calls are expected
#              to move for SVA and NOT for L1 or Alu, and check 4 states that as a prediction
#              before the run rather than reading it off afterwards.
#
# **SINCE THE 2026-08-05 FLIP, AFTER_OFF PASSES `--no-td-decoy-structural` EXPLICITLY.** It used
# to pass nothing, because nothing WAS off. Now nothing is ON, so an arm that passes no flag is
# the AFTER_ON arm under another name -- and this gate would have compared the fix against itself
# and reported a byte-identity FAIL it had caused. The containment claim did not change; the
# channel that expresses it did. This is the concrete meaning of "the off arm becomes the one
# that needs a reason to exist".
#
# WHAT chr1 CAN AND CANNOT SAY -- read this before trusting a PASS.
# chr1's gate B survivors are 5 (L1) and **0** (SVA) on the parent commit, so AFTER_ON's SVA
# delta at chr1 may be zero or near it: the fix's whole point is that hits to the 5,454 hg38 SVA
# flank records become sources, and chr1 holds ~8% of the genome's disc reads. What this leg
# proves is (a) the flag is inert when off, over 99 files, and (b) the flag is WIRED -- it
# reaches the subprocess and changes the classification. The population is WGS's, which is
# `gate_decoy_structural_wgs.sbatch`. A zero delta here is not evidence the fix does nothing.
#
# The two characterised exceptions are the ones the parent gates established over a dry run of
# two arms of the SAME commit (job 49599837): `cns_manifest.txt` embeds each arm's own output
# paths, `scan_disc_pos.tsv` is a permutation because the multithreaded scan writes rows in
# completion order.
#
# USAGE  sbatch gate_decoy_structural_chr1.sbatch <before-worktree> <after-worktree>
#   BEFORE may be any worktree whose xtea/ is behaviourally the parent --
#   /n/scratch/users/c/cos689/xtea2_wt/td_b2_365f352 qualifies (365f352..c26f2f8 is docs plus
#   two comment-only line-number fixes), which is also why its already-built xtea_scan serves
#   all three arms: src/ is byte-identical across the range, so one binary removes the C++ as a
#   variable.
# =============================================================================
set -euo pipefail

SCRATCH=/n/scratch/users/c/cos689/xtea2_testdata
WT_BEFORE="${1:?usage: sbatch gate_decoy_structural_chr1.sbatch <before-wt> <after-wt>}"
WT_AFTER="${2:?usage: sbatch gate_decoy_structural_chr1.sbatch <before-wt> <after-wt>}"
OUT="${XTEA2_OUT:-$SCRATCH/gate_decoy_structural_chr1_${SLURM_JOB_ID:-manual}}"
TESTDATA=/n/data1/hms/dbmi/park/corinne/claude_sandbox/test_data/HG002
BAM="$TESTDATA/chr1/HG002.chr1.bam"
SAMPLE=HG002.chr1
TYPES="${XTEA2_TYPES:-L1,Alu,SVA}"
REF=/n/data1/hms/dbmi/park-smaht_dac/ref/GRCh38_no_alt/hg38_no_alt.fa
CORES="${SLURM_CPUS_PER_TASK:-12}"
SCAN_BIN="${XTEA2_SCAN_BIN:-$WT_BEFORE/src/csrc/xtea_scan}"
REPLIB=/home/cos689/xTea/rep_lib_annotation

source /n/app/conda/miniforge3/24.11.3-0/etc/profile.d/conda.sh
conda activate xtea2
export PATH=/home/cos689/software/bin:$PATH
export XTEA_BASIC_INFO_SEED=42
test "$PYTHONNOUSERSITE" = 1
# An INHERITED opt-in would silently turn the BEFORE and AFTER_OFF arms into the fixed
# behaviour and make this gate pass while proving nothing. Cleared explicitly, once.
unset XTEA_TD_DECOY_STRUCTURAL

for wt in "$WT_BEFORE" "$WT_AFTER"; do
    [ -d "$wt" ] || { echo "FATAL: worktree $wt missing"; exit 1; }
    [ -z "$(git -C "$wt" status --porcelain)" ] || { echo "FATAL: $wt is dirty"; exit 1; }
done
[ -x "$SCAN_BIN" ] || { echo "FATAL: $SCAN_BIN not built"; exit 1; }
[ -f "$BAM" ]      || { echo "FATAL: $BAM missing"; exit 1; }
# The code under test must be present BEFORE the compute, not after it. Written as `if`, not
# `grep && { ... }`: a failing AND-list is the whole statement's status, so under `set -e` the
# negative form aborts the script on the case it is meant to accept.
if ! grep -q 'def decoy_tags' "$WT_AFTER/xtea/transduction_evidence.py"; then
    echo "FATAL: AFTER has no decoy_tags; wrong commit"; exit 1
fi
if grep -q 'def decoy_tags' "$WT_BEFORE/xtea/transduction_evidence.py"; then
    echo "FATAL: BEFORE already has decoy_tags; wrong commit"; exit 1
fi
if ! grep -q 'td-decoy-structural' "$WT_AFTER/xtea/run_xtea.py"; then
    echo "FATAL: AFTER has no --td-decoy-structural flag; wrong commit"; exit 1
fi
# The OFF arm's channel has to exist before three arms of compute, not after. Without it,
# AFTER_OFF would run the fix ON and check 1 would fail for a reason that is this script's fault.
if ! grep -q 'no-td-decoy-structural' "$WT_AFTER/xtea/run_xtea.py"; then
    echo "FATAL: AFTER has no --no-td-decoy-structural opt-out; this script needs the post-flip"
    echo "commit. Before the flip, use the version of this script from that commit."; exit 1
fi

mkdir -p "$OUT"

echo "=============================================================="
echo "GATE: the SVA decoy fix  --  HG002 chr1, types=$TYPES"
echo "  job     : ${SLURM_JOB_ID:-<none>} on $(hostname)  cores=$CORES"
echo "  BEFORE  : $(git -C "$WT_BEFORE" rev-parse --short HEAD)  $WT_BEFORE"
echo "  AFTER   : $(git -C "$WT_AFTER"  rev-parse --short HEAD)  $WT_AFTER"
echo "  scan bin: $SCAN_BIN (shared by all three arms)"
echo "  out     : $OUT"
echo "  started : $(date '+%F %T')"
echo "=============================================================="

run_arm() {   # $1 = arm name, $2 = worktree, $3 = "" or "--td-decoy-structural"
    # Separate statements, not one `local a=$1 b=...` line: under `set -u` the whole line
    # expands before `local` assigns (jobs 49452897 / 49516776).
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
      --parallel-downstream \
      --keep-intermediates \
      ${extra} \
      -n "$CORES" ) > "$armout/driver.log" 2>&1 || {
          echo "FAIL: arm $arm exited $?"; tail -60 "$armout/driver.log"; return 1; }
    wall=$(( $(date +%s) - t0 ))
    printf '  wall %02d:%02d:%02d\n' $((wall/3600)) $((wall%3600/60)) $((wall%60))
}

RC=0
run_arm BEFORE    "$WT_BEFORE" ""                          || RC=1
# BEFORE passes nothing because its commit has no such flag to pass. AFTER_OFF must now say so
# explicitly -- see the header. Getting this wrong makes the gate compare the fix against itself.
run_arm AFTER_OFF "$WT_AFTER"  "--no-td-decoy-structural"  || RC=1
run_arm AFTER_ON  "$WT_AFTER"  "--td-decoy-structural"     || RC=1
[ "$RC" -eq 0 ] || { echo; echo "GATE decoy (chr1): FAIL (an arm did not complete)"; exit 1; }

A="$OUT/BEFORE/$SAMPLE"
B="$OUT/AFTER_OFF/$SAMPLE"
C="$OUT/AFTER_ON/$SAMPLE"

echo
echo "--- check 0: the arms are the code and the configuration they claim ---"
# An arm that silently ran the other configuration would pass every comparison below, and here
# the discriminator is not only the commit but the FLAG -- which is exactly what a shared
# worktree makes easy to get wrong.
b_fn=$(grep -c 'def decoy_tags' "$WT_BEFORE/xtea/transduction_evidence.py" || true)
a_fn=$(grep -c 'def decoy_tags' "$WT_AFTER/xtea/transduction_evidence.py" || true)
echo "  decoy_tags defined: BEFORE=$b_fn AFTER=$a_fn (expect 0 and 1)"
{ [ "$b_fn" -eq 0 ] && [ "$a_fn" -eq 1 ]; } || { echo "    FATAL: wrong commits"; RC=1; }
# Since the flip, run_xtea prints a banner on EVERY path and says which way it resolved. So the
# off arms are checked POSITIVELY as well as negatively: "no STRUCTURAL banner" alone cannot be
# told from "the driver died before it got there", and an off-arm gate whose only evidence is an
# absent line passes when the run fails early.
STRUCT='\[decoy\] the transduction decoy test is STRUCTURAL'
SUBSTR='\[decoy\] the transduction decoy test is the SUBSTRING match'
for arm in BEFORE AFTER_OFF; do
    n=$(grep -ac "$STRUCT" "$OUT/$arm/driver.log" || true)
    echo "  [decoy] STRUCTURAL banner in $arm: $n (expect 0)"
    [ "$n" -eq 0 ] || { echo "    FATAL: arm $arm ran with the fix ON"; RC=1; }
done
# BEFORE's commit predates the banner entirely, so it has neither line -- which is why only
# AFTER_OFF gets the positive check. For BEFORE, `decoy_tags` being absent (above) is the proof.
n=$(grep -ac "$SUBSTR" "$OUT/AFTER_OFF/driver.log" || true)
echo "  [decoy] SUBSTRING banner in AFTER_OFF: $n (expect 1)"
[ "$n" -eq 1 ] || { echo "    FATAL: AFTER_OFF's opt-out did not reach the driver -- this arm is"
                    echo "    NOT the off arm, and check 1 below would be the fix against itself"; RC=1; }
n=$(grep -ac "$STRUCT" "$OUT/AFTER_ON/driver.log" || true)
echo "  [decoy] STRUCTURAL banner in AFTER_ON: $n (expect 1)"
[ "$n" -eq 1 ] || { echo "    FATAL: AFTER_ON did not run with the fix ON"; RC=1; }

echo
echo "--- check 1: AFTER_OFF is byte-identical to BEFORE, over the whole output tree ---"
# THE CONTAINMENT CLAIM. Every file, including both decision ledgers -- unlike the gate B2
# ledger gates, nothing here is *supposed* to change, so the ledgers are compared outright.
n_same=0; n_diff=0; n_missing=0; n_tol=0
while IFS= read -r f; do
    rel=${f#"$A/"}
    case "$rel" in
        tmp/*|*/tmp/*|tmp_basic/*|pub_clip/*) continue ;;   # scratch, not results
        *.log|step_timings.tsv) continue ;;                 # wall clocks and stderr
    esac
    if [ ! -f "$B/$rel" ]; then echo "  MISSING in AFTER_OFF: $rel"; n_missing=$((n_missing+1)); continue; fi
    if cmp -s "$f" "$B/$rel"; then n_same=$((n_same+1)); continue; fi
    case "$rel" in
      cns_manifest.txt)
        if cmp -s <(sed "s|$OUT/BEFORE/||g" "$f") <(sed "s|$OUT/AFTER_OFF/||g" "$B/$rel"); then
            n_tol=$((n_tol+1)); echo "  $rel: identical after stripping the arm prefix"
        else
            echo "  DIFFERS (beyond the path prefix): $rel"; n_diff=$((n_diff+1))
        fi ;;
      scan_disc_pos.tsv)
        if cmp -s <(sort "$f") <(sort "$B/$rel"); then
            n_tol=$((n_tol+1)); echo "  $rel: identical when sorted ($(wc -l < "$f") rows)"
        else
            echo "  DIFFERS (not a permutation): $rel"; n_diff=$((n_diff+1))
        fi ;;
      *)
        echo "  DIFFERS: $rel"; n_diff=$((n_diff+1))
        diff <(grep -v '^#' "$f" | cut -f1,2) <(grep -v '^#' "$B/$rel" | cut -f1,2) | head -20 || true ;;
    esac
done < <(find "$A" -type f | sort)
echo "  identical: $n_same    tolerated: $n_tol    DIFFERING: $n_diff    missing: $n_missing"
[ "$n_diff" -eq 0 ] && [ "$n_missing" -eq 0 ] || RC=1
# A run that produced almost nothing would otherwise "pass" with n_diff 0. The parent gate
# measured 99 identical + 2 characterised at this sample and type set.
[ "$n_same" -ge 95 ] || { echo "    FATAL: only $n_same files compared, expected >=95"; RC=1; }

echo
echo "--- check 2: no new or removed files in EITHER after arm ---"
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
    echo "  $arm: unexpected new files: $n_new (expect 0 -- this fix adds no artifact)"
    [ "$n_new" -eq 0 ] || RC=1
done

echo
echo "--- check 3: THE PREDICTION -- AFTER_ON moves SVA and nothing else ---"
# Stated before the run. L1 flank ids collide with no decoy tag ("LINE1" is not "L1"), and Alu
# never runs the transduction step at all, so those two families' output must be byte-identical
# to BEFORE with the fix ON. If an L1 or Alu file moves, the classification is doing something
# other than what the probe measured and the SVA delta below cannot be attributed to the fix.
n_l1_diff=0; n_alu_diff=0; n_sva_diff=0; n_other_diff=0
while IFS= read -r f; do
    rel=${f#"$A/"}
    case "$rel" in
        tmp/*|*/tmp/*|tmp_basic/*|pub_clip/*|*.log|step_timings.tsv) continue ;;
        cns_manifest.txt|scan_disc_pos.tsv) continue ;;     # characterised in check 1
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
echo "  files moved with the fix ON: L1=$n_l1_diff Alu=$n_alu_diff SVA=$n_sva_diff other=$n_other_diff"
echo "  (the 'other' bucket is the merged VCF and the locus table, which contain SVA rows and"
echo "   so are EXPECTED to move if SVA moved -- listed above, not failed on)"
[ "$n_l1_diff" -eq 0 ] || { echo "    FAIL: the fix moved L1 output. L1 collides with no decoy"
                            echo "    tag, so this contradicts the classification's premise"; RC=1; }
[ "$n_alu_diff" -eq 0 ] || { echo "    FAIL: the fix moved Alu output, which never runs the"
                             echo "    transduction step at all"; RC=1; }

echo
echo "--- check 4: what the fix reached, per arm (the measurement) ---"
# `probe_td_source_classes.py`'s run half: of the sources the step PICKED, how many are
# reference copies (`chrom:start-end~rc` with start != end) and how many polymorphic. The whole
# defect is that SVA's reference-copy count is structurally 0, so this is the payoff column.
for arm in BEFORE AFTER_OFF AFTER_ON; do
    echo "  --- $arm ---"
    python "$WT_AFTER/docs/handoff/runs/probe_td_source_classes.py" \
        --rep-lib-dir "$REPLIB" --run-dir "$OUT/$arm/$SAMPLE" 2>&1 |
        sed -n '/a finished run/,$p' | sed 's/^/    /'
done

echo
echo "--- check 5: gate B's domain and the ledger populations, per arm ---"
# The fix adds hits that used to be discarded, so gate B sees more evidence per locus. Where
# that lands -- more picked sources, or more `no_dominant_source` rejections -- is the thing
# a WGS-scale arm has to answer, and this is the same measurement at chr1 scale.
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
    echo "GATE decoy (chr1): PASS"
else
    echo "GATE decoy (chr1): FAIL"
fi
echo "  finished: $(date '+%F %T')"
echo "=============================================================="
exit "$RC"
