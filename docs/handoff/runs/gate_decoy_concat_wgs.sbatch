#!/bin/bash
#SBATCH -J concat_wgs
#SBATCH -c 24
#SBATCH -t 6:00:00
#SBATCH --mem=64G
#SBATCH -p priopark
#SBATCH --account=park_contrib
#SBATCH -o /n/scratch/users/c/cos689/xtea2_testdata/gate_decoy_concat_wgs_%j.out
#SBATCH --mail-type=NONE
# =============================================================================
# The decoy-consensus concatenation at WGS scale. ONE arm, WITH `--td-decoy-concat` ON -- and
# this is the leg the change exists for, because chr1's gate B has 0 SVA survivors and therefore
# no population to move.
#
# WHY ONE ARM, AND WHICH RUN IS THE OTHER ONE. **Job 49792012** is the frozen OFF arm at no
# compute: same BAM, same seed, same flags, `--td-decoy-structural` on, `--td-decoy-concat`
# nonexistent -- which is the current default configuration exactly. Its tree, its VCFs, its
# ledgers, its probe outputs and its truth-set score are all on disk. **BASELINE-0 is NOT the
# comparator here** and using it would silently measure two changes at once: it predates the
# 2026-08-05 structural flip, so an SVA delta against it is the flip plus the concatenation.
# That is the same "a number taken from the wrong population" error plan §10 records twice.
#
# THE PREDICTION, STATED BEFORE THE RUN, and it is sharper than the structural fix's:
#
#   1. **L1 and Alu byte-identical to the OFF arm.** Not because L1 collides with no decoy tag
#      -- that is the OTHER fix's argument -- but because `hg38_FL_L1_flanks_3k.fa` ALREADY
#      CARRIES all 45 decoy records, so `write_missing_decoys` selects nothing and bwa indexes
#      the same bytes. Alu never runs the transduction step. Check 1a verifies this at the
#      reference file, upstream of every gate.
#   2. **SVA's bwa reference gains exactly 45 records** and the rebuilt block appears verbatim.
#   3. **SVA's clip-side decoy support goes 0 -> nonzero.** HANDOFF_2026-08-05 §2.4 measured L1
#      at 3,992 of 99,245 qualified clip records resting on decoy-consensus support and SVA at
#      **0 of 128,262**, for the packaging reason this change fixes. If it is still 0, the block
#      reached the file and not the aligner.
#   4. **SVA's disc side gains a live `decoy_cross_family` arm.** It fires 144x for L1 and 0 for
#      SVA today, and the SVA arm counting zero is not a measurement of zero competition -- it
#      is a gate with nothing to count. Whether it now REJECTS anything is the open question, and
#      it is the one plan §5's SVA-is-SINE-VNTR-Alu hazard has been waiting on.
#
# WHAT IS REPORTED AND NOT GATED, and read it that way: the truth-set score, the Stage 3
# transduction gate in `--baseline` no-regression mode, and the FN attribution. §7.1's obstacle
# applies unchanged -- the derived transduction denominator on HG002 is 0 (L1) and 2 (SVA), so
# the truth set scores POSITIONS, not transductions. With 90 SVA truth records and a 5-FP budget
# one SVA call is ~1.2 precision points (§10), so quote a small delta as the calls behind it and
# not as the number.
#
# **DIRECTION IS NOT PREDICTED.** The concatenation adds competition, so its expected effect on
# SVA is to REMOVE spurious sources (a read that came from SVA/Alu consensus body now maps to the
# consensus record instead of to some reference copy's flank), which is a precision-shaped change
# and could equally show up as lost calls. That is a legitimate outcome, not a failure: plan §4's
# FN attribution says 306 of BASELINE-0's 322 false negatives never reach this step at all, so
# sensitivity is not what is on the table. State what moved and why; do not root for a direction.
#
# Frozen worktree, not the live checkout -- run_xtea spawns a fresh x_TEA_main per step, so a job
# sharing the checkout picks up mid-run edits.
#
# USAGE  sbatch [--dependency=afterok:<chr1 job>] gate_decoy_concat_wgs.sbatch <worktree>
# =============================================================================
set -euo pipefail

SCRATCH=/n/scratch/users/c/cos689/xtea2_testdata
WT="${1:?usage: sbatch gate_decoy_concat_wgs.sbatch <worktree> [types]}"
TYPES="${2:-L1,Alu,SVA}"
OUT="${XTEA2_OUT:-$SCRATCH/gate_decoy_concat_wgs_${SLURM_JOB_ID:-manual}}"
# The frozen OFF arm: job 49792012, the decoy-structural WGS gate. Same BAM, same seed, same
# flags, and the current default configuration.
OFFARM="${XTEA2_OFF_ARM:-$SCRATCH/gate_decoy_structural_wgs_49792012/HG002.wgs}"
TESTDATA=/n/data1/hms/dbmi/park/corinne/claude_sandbox/test_data/HG002
BAM="$TESTDATA/alignment/HG002.novaseq.pcr-free.30x.dedup.grch38.bam"
SAMPLE=HG002.wgs
REF=/n/data1/hms/dbmi/park-smaht_dac/ref/GRCh38_no_alt/hg38_no_alt.fa
CORES="${SLURM_CPUS_PER_TASK:-24}"
MAX_DOWNSTREAM_S="${XTEA2_MAX_DOWNSTREAM_S:-900}"
TRUTH_DIR="${XTEA2_TRUTH_DIR:-$TESTDATA/truth}"
REPLIB=/home/cos689/xTea/rep_lib_annotation
VCF_PREFIX=HG002.novaseq.pcr-free.30x.dedup.grch38
SCAN_BIN="${XTEA2_SCAN_BIN:-/n/scratch/users/c/cos689/xtea2_wt/td_decoy_55055ba/src/csrc/xtea_scan}"

source /n/app/conda/miniforge3/24.11.3-0/etc/profile.d/conda.sh
conda activate xtea2
export PATH=/home/cos689/software/bin:$PATH
export XTEA_BASIC_INFO_SEED=42   # pin the insert-size draw, as the OFF arm did
test "$PYTHONNOUSERSITE" = 1
# Both flags are passed on the command line below; an inherited env var would make it impossible
# to tell which channel armed which, and the banners are check 0's evidence either way.
unset XTEA_TD_DECOY_CONCAT
unset XTEA_TD_DECOY_STRUCTURAL

[ -d "$WT" ]                 || { echo "FATAL: worktree $WT missing"; exit 1; }
[ -x "$SCAN_BIN" ]           || { echo "FATAL: $SCAN_BIN not built"; exit 1; }
[ -f "$BAM" ]                || { echo "FATAL: $BAM missing"; exit 1; }
[ -f "$BAM.bai" ]            || { echo "FATAL: $BAM.bai missing"; exit 1; }
[ -z "$(git -C "$WT" status --porcelain)" ] || { echo "FATAL: worktree is dirty"; exit 1; }
[ -f "$OFFARM/$SAMPLE.vcf" ] || { echo "FATAL: the frozen OFF arm $OFFARM is not there"; exit 1; }
# The OFF arm must be the configuration it is being used as. Checking it, rather than trusting
# the path name, is the whole lesson of plan §10's "a single-arm pool is not comparable" entry.
if ! grep -aq '\[decoy\] the transduction decoy test is STRUCTURAL' "$OFFARM/../driver.log"; then
    echo "FATAL: $OFFARM did not run with --td-decoy-structural, so an SVA delta against it"
    echo "would be two changes at once. Point XTEA2_OFF_ARM at a post-flip run."; exit 1
fi
if [ ! -f "$WT/xtea/td_decoy_library.py" ]; then
    echo "FATAL: worktree has no td_decoy_library; wrong commit"; exit 1
fi
if ! grep -q -- 'td-decoy-concat' "$WT/xtea/run_xtea.py"; then
    echo "FATAL: worktree has no --td-decoy-concat flag; wrong commit"; exit 1
fi

mkdir -p "$OUT"
# The recipe has to reconstruct on this machine's rep_lib BEFORE 40 minutes of compute.
python "$WT/docs/handoff/runs/probe_decoy_library.py" --rep-lib-dir "$REPLIB" \
    > "$OUT/decoy_census.txt" 2>&1 || {
        echo "FATAL: the decoy recipe does not reconstruct from $REPLIB:"
        cat "$OUT/decoy_census.txt"; exit 1; }

cd "$WT"

echo "=============================================================="
echo "GATE: the decoy-consensus concatenation, WGS, FLAG ON  --  types=$TYPES"
echo "  job      : ${SLURM_JOB_ID:-<none>} on $(hostname)  cores=$CORES"
echo "  commit   : $(git rev-parse --short HEAD) $(git log -1 --format=%s | cut -c1-60)"
echo "  worktree : $WT"
echo "  scan bin : $SCAN_BIN"
echo "  OFF arm  : $OFFARM  (job 49792012, frozen)"
echo "  out      : $OUT"
echo "  started  : $(date '+%F %T')"
echo "=============================================================="
echo
echo "--- the library census this gate's prediction rests on ---"
sed 's/^/  /' "$OUT/decoy_census.txt"

t0=$(date +%s)
python xtea/run_xtea.py \
  -i "$BAM" -s "$SAMPLE" -o "$OUT" \
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
  --td-decoy-concat \
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
echo "--- check 0: this arm really ran with BOTH flags on ---"
n=$(grep -ac '\[decoy-concat\] the transduction bwa reference GAINS' "$OUT/driver.log" || true)
echo "  [decoy-concat] ON banner: $n (expect 1)"
[ "$n" -eq 1 ] || { echo "  FATAL: the concat flag did not reach the driver -- this arm is a no-op"; RC=1; }
n=$(grep -ac '\[decoy\] the transduction decoy test is STRUCTURAL' "$OUT/driver.log" || true)
echo "  [decoy] STRUCTURAL banner: $n (expect 1 -- the OFF arm has it too)"
[ "$n" -eq 1 ] || { echo "  FATAL: this arm differs from the OFF arm by TWO flags, not one"; RC=1; }
echo "  the step's own per-family append counts (the line that separates armed from changed):"
grep -a 'td-decoy-concat' "$B"/*.steps.log "$OUT/driver.log" 2>/dev/null |
    sed 's/^/    /' | cut -c1-200 || echo "    (none found)"

echo
echo "--- check 1a: THE BWA REFERENCE -- L1 unchanged, SVA longer by exactly 45 ---"
# Upstream of every gate: what did bwa actually index? The OFF arm kept its references too.
for fam in L1 SVA; do
    r_off="$OFFARM/$fam/tmp/transduction/all_with_polymerphic_flanks.fa"
    r_on="$B/$fam/tmp/transduction/all_with_polymerphic_flanks.fa"
    if [ ! -f "$r_off" ] || [ ! -f "$r_on" ]; then
        echo "  $fam: a reference is missing -- CANNOT CHECK"; RC=1; continue
    fi
    n_off=$(grep -c '^>' "$r_off"); n_on=$(grep -c '^>' "$r_on")
    echo "  $fam: records OFF=$n_off ON=$n_on  delta=$(( n_on - n_off ))"
    if [ "$fam" = L1 ]; then
        cmp -s "$r_off" "$r_on" && echo "    PASS: byte-identical -- L1's library already carries all 45" || {
            echo "    FAIL: L1's bwa reference MOVED, but its packaged library carries the block"
            diff <(grep '^>' "$r_off") <(grep '^>' "$r_on") | head -20 || true; RC=1; }
    else
        [ "$(( n_on - n_off ))" -eq 45 ] || {
            echo "    FAIL: expected exactly 45 new records on hg38's SVA library"; RC=1; }
        python - "$r_on" "$WT" "$REPLIB" <<'EOF' || RC=1
import os, sys
r_on, wt, replib = sys.argv[1], sys.argv[2], sys.argv[3]
sys.path.insert(0, os.path.join(wt, "xtea"))
import td_decoy_library as tdl
block = tdl.decoy_block_text(tdl.build_decoy_block(os.path.join(replib, "consensus")))
text = open(r_on).read()
if block not in text:
    print("    FAIL: the 45-record block is not present verbatim in the ON reference")
    raise SystemExit(1)
tail = text[text.index(block) + len(block):]
print("    PASS: the block appears verbatim, followed by {0}".format(
    tail.split("\n", 1)[0][:60] or "<end of file>"))
EOF
    fi
done

echo
echo "--- check 1b: the call sets -- L1 and Alu identical to the OFF arm, SVA not ---"
n_same=0; n_diff=0
for f in "$OFFARM"/merged/*.vcf "$OFFARM/$SAMPLE.vcf"; do
    bn=$(basename "$f")
    b=$(find "$B" -name "$bn" -type f | head -1)
    if [ -z "$b" ]; then echo "  MISSING: $bn"; n_diff=$((n_diff+1)); continue; fi
    if cmp -s "$f" "$b"; then s=IDENTICAL; n_same=$((n_same+1)); else s=DIFFERS; n_diff=$((n_diff+1)); fi
    echo "  $bn: $s"
    case "$bn" in
      *_LINE1.vcf|*_ALU.vcf)
        [ "$s" = IDENTICAL ] || {
            echo "    FAIL: a family whose bwa reference did not change moved anyway."
            diff <(grep -v '^#' "$f" | cut -f1,2) <(grep -v '^#' "$b" | cut -f1,2) | head -20 || true
            RC=1; } ;;
      *_SVA.vcf)
        echo "    positions gained/lost vs the OFF arm:"
        diff <(grep -v '^#' "$f" | cut -f1,2) <(grep -v '^#' "$b" | cut -f1,2) |
            grep -E '^[<>]' | sed 's/^/      /' | head -40 || true
        # Deliberately NOT failed on in either direction. A concatenation that removes a
        # spurious source is doing its job; so is one that changes nothing at the call level
        # while changing what the disc parse sees. What would be a failure is check 1a.
        ;;
    esac
done
echo "  identical: $n_same    differing: $n_diff"

echo
echo "--- check 1c: TD_SRC conversions, which a position diff cannot see ---"
# plan §10: the decoy FLIP's call set grew by 3 while its verification pool grew by 4, because a
# position BOTH arms call had its TD_SRC rewritten. A field this change can rewrite in place must
# be diffed as a field.
python - "$OFFARM/$SAMPLE.vcf" "$B/$SAMPLE.vcf" <<'EOF' || RC=1
import re, sys
def td_src(path):
    m = {}
    for line in open(path):
        if line.startswith("#"):
            continue
        f = line.split("\t")
        hit = re.search(r"TD_SRC=([^;\t\n]*)", line)
        m[(f[0], f[1])] = hit.group(1) if hit else "<absent>"
    return m
a, b = td_src(sys.argv[1]), td_src(sys.argv[2])
both = sorted(set(a) & set(b))
changed = [k for k in both if a[k] != b[k]]
print("  positions in both arms: {0}   TD_SRC changed on: {1}".format(len(both), len(changed)))
for k in changed[:40]:
    print("    {0}:{1}  {2}  ->  {3}".format(k[0], k[1], a[k][:60], b[k][:60]))
print("  gained positions: {0}   lost positions: {1}".format(
    len(set(b) - set(a)), len(set(a) - set(b))))
EOF

echo
echo "--- check 2: call counts. L1 = 188 and Alu = 1223 are gated; SVA is reported ---"
l1=$(grep -vc '^#' "$B/merged/${VCF_PREFIX}_LINE1.vcf" 2>/dev/null || echo 0)
alu=$(grep -vc '^#' "$B/merged/${VCF_PREFIX}_ALU.vcf" 2>/dev/null || echo 0)
sva=$(grep -vc '^#' "$B/merged/${VCF_PREFIX}_SVA.vcf" 2>/dev/null || echo 0)
mrg=$(grep -vc '^#' "$B/$SAMPLE.vcf" 2>/dev/null || echo 0)
sva0=$(grep -vc '^#' "$OFFARM/merged/${VCF_PREFIX}_SVA.vcf" 2>/dev/null || echo 0)
echo "  L1=$l1 (188)  Alu=$alu (1223)  SVA=$sva (OFF arm $sva0)  merged=$mrg"
{ [ "$l1" -eq 188 ] && [ "$alu" -eq 1223 ]; } || { echo "  FAIL: L1 or Alu moved"; RC=1; }
echo "  SVA delta: $(( sva - sva0 )) calls"

echo
echo "--- check 3: THE CLIP CHANNEL -- SVA's 0 of 128,262 must stop being 0 ---"
# HANDOFF_2026-08-05 §2.4, measured on the run's own surviving clip SAM. This is the most direct
# statement that the block reached the ALIGNER and not merely the file.
for arm_name in OFF ON; do
    [ "$arm_name" = OFF ] && d="$OFFARM" || d="$B"
    echo "  --- $arm_name arm ($d) ---"
    if python "$WT/docs/handoff/runs/probe_clip_gate.py" --run-dir "$d" \
            > "$OUT/probe_clip_gate_$arm_name.txt" 2>&1; then
        grep -aE 'decoy-consensus|^  [A-Z]|survivors' "$OUT/probe_clip_gate_$arm_name.txt" |
            sed 's/^/    /' | head -30
    else
        echo "    (the probe exited nonzero -- reporting, not failing)"
        tail -20 "$OUT/probe_clip_gate_$arm_name.txt" | sed 's/^/    /'
    fi
done
n_sva_decoy=$(grep -aA4 'SVA' "$OUT/probe_clip_gate_ON.txt" 2>/dev/null |
              grep -aoE 'decoy-consensus clip records: [0-9]+' | head -1 |
              grep -oE '[0-9]+' || true)
echo "  SVA decoy-consensus clip records, ON arm: ${n_sva_decoy:-<unparsed>} (OFF arm: 0 of 128,262)"
if [ -z "${n_sva_decoy:-}" ]; then
    echo "    (could not parse -- read the probe output above; reporting, not failing)"
elif [ "$n_sva_decoy" -le 0 ]; then
    echo "    FAIL: still zero. The 45 records are in the reference (check 1a) but no SVA clip"
    echo "    read prefers one, which contradicts L1's 3,992 of 99,245 on the same channel."; RC=1
else
    echo "    PASS: SVA has a decoy channel on the clip side for the first time"
fi

echo
echo "--- check 4: the disc side -- decoy_cross_family and the picked sources ---"
for arm_name in OFF ON; do
    [ "$arm_name" = OFF ] && d="$OFFARM" || d="$B"
    echo "  --- $arm_name arm ---"
    for fam in L1 SVA; do
        led="$d/$fam/candidate_disc_filtered_cns2.txt.all_non_sibling_td.txt.decision_ledger.tsv"
        src="$d/$fam/candidate_disc_filtered_cns2.txt.all_non_sibling_td.txt.unique_trsdct_disc_only"
        n_src=$([ -f "$src" ] && wc -l < "$src" || echo NA)
        if [ -f "$led" ]; then
            echo "    $fam: picked sources=$n_src  ledger rows=$(( $(wc -l < "$led") - 1 ))"
            echo "      decoy_cross_family=$(awk -F'\t' 'NR>1 && $4=="decoy_cross_family"' "$led" | wc -l)"
            echo "      by reason: $(awk -F'\t' 'NR>1 {c[$4]++} END {for (r in c) printf "%s=%d ", r, c[r]}' "$led")"
        else
            echo "    $fam: picked sources=$n_src  (no ledger)"
        fi
    done
done

echo
echo "--- check 5: the source classes -- did reference-copy SVA sources survive the competition? ---"
# The structural flip took SVA reference-copy picked sources 0 -> 241 of 346. The concatenation
# adds competitors for exactly those hits, so this number MAY fall, and a fall is not a failure:
# a hit that prefers the SVA consensus over a specific reference copy was never evidence for that
# copy. What would be alarming is 241 -> 0, i.e. the block outcompeting the whole library.
for arm_name in OFF ON; do
    [ "$arm_name" = OFF ] && d="$OFFARM" || d="$B"
    echo "  --- $arm_name arm ---"
    python "$WT/docs/handoff/runs/probe_td_source_classes.py" \
        --rep-lib-dir "$REPLIB" --run-dir "$d" > "$OUT/probe_src_$arm_name.txt" 2>&1 &&
        sed -n '/a finished run/,$p' "$OUT/probe_src_$arm_name.txt" | sed 's/^/    /' ||
        { echo "    (the probe exited nonzero)"; tail -15 "$OUT/probe_src_$arm_name.txt" | sed 's/^/    /'; }
done

echo
echo "--- check 6: downstream_parallel, loosely bounded ---"
# Prior WGS runs at this configuration: 462.7s, 511.75s (a slower node), 463.3s, 485.8s. The
# concatenation makes SVA's bwa reference 45 records longer out of 5,454 -- immaterial to the
# index -- but hits that used to go to a flank record now go to a decoy and are skipped earlier,
# so if anything expect SVA's disc parse to get slightly cheaper. NOT a measurement: one arm on
# one node has no control.
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

echo
echo "--- check 7: the truth-set score, OFF arm vs this arm (REPORTED, not gated) ---"
# §7.1: the derived transduction denominator on HG002 is 0 (L1) and 2 (SVA), so this scores
# POSITIONS. SVA has 90 truth records and a 5-FP budget, so one call is ~1.2 precision points --
# quote the calls, not the percentage.
if [ -d "$TRUTH_DIR" ]; then
    for arm_name in OFF ON; do
        [ "$arm_name" = OFF ] && v="$OFFARM/$SAMPLE.vcf" || v="$B/$SAMPLE.vcf"
        echo "  === $arm_name arm ($v) ==="
        if python "$WT/benchmark/benchmark_truth_set.py" score \
                --calls "$v" --truth-dir "$TRUTH_DIR" --families L1,Alu,SVA \
                > "$OUT/score_$arm_name.txt" 2>&1; then
            sed -n '/--- L1 ---/,$p' "$OUT/score_$arm_name.txt" | sed 's/^/    /'
        else
            echo "    (the scorer exited nonzero -- reporting, not failing)"
            tail -20 "$OUT/score_$arm_name.txt" | sed 's/^/    /'
        fi
    done
else
    echo "  (no truth dir at $TRUTH_DIR -- skipped; set XTEA2_TRUTH_DIR)"
fi

echo
echo "--- check 8: Stage 3's transduction gate -- NOT run here, and that is deliberate ---"
# The instrument for "two arms, one flag apart" already exists and is generalised:
# `td_gate_no_regression.sbatch` with XTEA2_OFF_ARM / XTEA2_ON_ARM. Re-implementing its pool
# build and its `--baseline` adjudication inline would be a second copy of a rule whose whole
# value is that there is one -- and the specific error it exists to prevent (comparing a
# single-arm pool against the frozen 3-arm baseline) is easy to re-make by hand. It is also a
# different resource class: 8 cores on `short`, ~20 s, no xTea run.
cat <<CMD
  run it as a dependent job:

    XTEA2_OFF_ARM=$OFFARM \\
    XTEA2_ON_ARM=$B \\
    XTEA2_OFF_LEDGER=none \\
    XTEA2_ARMS_LABEL='the --td-decoy-concat arm: 49792012 (OFF) vs ${SLURM_JOB_ID:-this run} (ON)' \\
      sbatch $WT/docs/handoff/runs/td_gate_no_regression.sbatch $WT

  XTEA2_OFF_LEDGER=none because only BASELINE-0 has a frozen ledger; job 49792012 does not, so
  the "is this the run it says it is" assertion is not available for the OFF arm here. Its
  identity is instead pinned by check 0's banner test above, which this gate already ran.
CMD

echo
echo "--- check 9: the funnel probe still reconciles ---"
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
    echo "GATE concat (WGS): PASS"
else
    echo "GATE concat (WGS): FAIL"
fi
echo "  finished: $(date '+%F %T')"
echo "=============================================================="
exit "$RC"
