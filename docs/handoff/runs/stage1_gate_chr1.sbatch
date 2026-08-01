#!/bin/bash
#SBATCH -c 12
#SBATCH -t 6:00:00
#SBATCH --mem=48G
#SBATCH -p priopark
#SBATCH --account=park_contrib
#SBATCH -o /n/scratch/users/c/cos689/xtea2_testdata/stage1_gate_chr1_%j.out
#SBATCH --mail-type=NONE

# =============================================================================
# STAGE 1 GATE -- DOWNSTREAM_OVERHAUL_PLAN.md Stage 1 ("mechanics only, agreement-gated").
#
#     sbatch docs/handoff/runs/stage1_gate_chr1.sbatch
#
# What the plan asks for, and how each part is produced here:
#
#   "a chr1 diff vs BASELINE-0 -- byte-identical calls required"
#       BASELINE-0 is a WGS artifact (three arms, HG002 whole genome), so it cannot be
#       diffed against a chr1 run directly. The equivalent, and the stronger test, is an
#       A/B on the chr1 slice at ONE commit: arm A runs the legacy Python gather, arm B
#       the native one, everything else identical. Byte-identical per-type call lists and
#       VCFs is then a statement about the port alone, with no confound from the sample
#       or the commit. (BASELINE-0 stays the Stage 2 comparator, where calls may move.)
#
#   "plus a two-run determinism check"
#       arm C repeats arm B. B == C byte for byte, or the native gather is not
#       deterministic and no A/B result means anything.
#
#   "tests/test_downstream_from_cns.py (VCF cols 1-8)"
#       run LAST, with the native gather forced on, against the NA12878 chr1 golden --
#       a second sample, and the only arm here that compares against a frozen external
#       baseline rather than a same-job control.
#
# The three pipeline arms run SEQUENTIALLY on purpose: they are being compared, and
# co-scheduling them on one node would let a noisy neighbour arm change the timing
# numbers reported at the end. Correctness would survive it; the speedup figure would not.
#
# Frozen worktree, not the live checkout -- run_xtea spawns a fresh x_TEA_main per step,
# so a job sharing the checkout picks up mid-run edits.
# =============================================================================
set -euo pipefail

SCRATCH=/n/scratch/users/c/cos689/xtea2_testdata
WT="${XTEA2_WORKTREE:-$SCRATCH/stage1_wt}"
OUT="${XTEA2_OUT:-$SCRATCH/stage1_gate_chr1_${SLURM_JOB_ID:-manual}}"
TESTDATA=/n/data1/hms/dbmi/park/corinne/claude_sandbox/test_data/HG002
BAM="$TESTDATA/chr1/HG002.chr1.bam"
SAMPLE=HG002.chr1
TYPES="${1:-L1,Alu,SVA}"
REF=/n/data1/hms/dbmi/park-smaht_dac/ref/GRCh38_no_alt/hg38_no_alt.fa
CORES="${SLURM_CPUS_PER_TASK:-12}"

source /n/app/conda/miniforge3/24.11.3-0/etc/profile.d/conda.sh
conda activate xtea2
export PATH=/home/cos689/software/bin:$PATH
export XTEA_BASIC_INFO_SEED=42   # pin the insert-size draw, as every comparator run did
test "$PYTHONNOUSERSITE" = 1     # the env's guard against ~/.local's onnxruntime

[ -d "$WT" ]                    || { echo "FATAL: worktree $WT missing"; exit 1; }
[ -x "$WT/src/csrc/xtea_scan" ] || { echo "FATAL: $WT/src/csrc/xtea_scan not built"; exit 1; }
[ -f "$BAM" ]                   || { echo "FATAL: $BAM missing"; exit 1; }
[ -f "$BAM.bai" ]               || { echo "FATAL: $BAM.bai missing"; exit 1; }
[ -z "$(git -C "$WT" status --porcelain)" ] || { echo "FATAL: worktree is dirty"; exit 1; }
# A gate that silently ran the legacy path in both arms would PASS and mean nothing.
grep -q "td-gather-scan" "$WT/xtea/run_xtea.py" || {
    echo "FATAL: worktree has no --td-gather-scan; wrong commit"; exit 1; }
"$WT/src/csrc/xtea_scan" --help 2>&1 | grep -q -- "--td-sites" || {
    echo "FATAL: $WT/src/csrc/xtea_scan predates --td-sites; rebuild it"; exit 1; }

mkdir -p "$OUT"
cd "$WT"

echo "=============================================================="
echo "STAGE 1 GATE  --  HG002 chr1, types=$TYPES"
echo "  job    : ${SLURM_JOB_ID:-<none>} on $(hostname)  cores=$CORES"
echo "  commit : $(git rev-parse --short HEAD) $(git log -1 --format=%s | cut -c1-66)"
echo "  bam    : $BAM"
echo "  out    : $OUT"
echo "  started: $(date '+%F %T')"
echo "=============================================================="

run_arm() {   # $1 = arm name, $2 = "native"|"legacy"
    # Separate statements, not `local arm=$1 armout="$OUT/$arm"`: the shell expands the WHOLE
    # command line before `local` assigns anything, so $arm is still unset there -- which under
    # `set -u` kills the job before the first arm runs.
    local arm=$1 mode=$2 t0 wall armout
    armout="$OUT/$arm"
    rm -rf "$armout"; mkdir -p "$armout"
    echo
    echo "=== arm $arm: gather=$mode ==="
    t0=$(date +%s)
    # The native gather is selected by an env var because the transduction step is its own
    # subprocess; --td-gather-scan is what sets it. Unset it explicitly for the legacy arm
    # so a stray value in the job environment cannot silently make both arms native.
    if [ "$mode" = native ]; then TD_FLAG=(--td-gather-scan); else TD_FLAG=(); unset XTEA_SCAN_TD_GATHER || true; fi
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
      "${TD_FLAG[@]}" \
      --parallel-downstream \
      --keep-intermediates \
      -n "$CORES" \
      > "$armout/driver.log" 2>&1 || { echo "FAIL: arm $arm exited $?"; tail -40 "$armout/driver.log"; return 1; }
    wall=$(( $(date +%s) - t0 ))
    printf '  wall %02d:%02d:%02d\n' $((wall/3600)) $((wall%3600/60)) $((wall%60))
    echo "$wall" > "$armout/wall_s.txt"
    # [TD-PHASE] lands in the STEP's stderr, never the driver's stdout -- and WHICH log it lands
    # in depends on --parallel-downstream: the serial path collects every step into
    # <sample>_xtea2.log, while the parallel path gives each type its own
    # <sample>_<TYPE>.steps.log and the timers go there. Glob both, or the timing section is
    # silently empty (which is how this was found).
    grep -ah '^\[TD-PHASE\]' "$armout/$SAMPLE/${SAMPLE}_xtea2.log" \
        "$armout/$SAMPLE/${SAMPLE}"_*.steps.log 2>/dev/null \
        > "$armout/td_phases.txt" || true
    awk -F'\t' '{printf "  rtype=%-3s %-40s %8.1fs\n",$3,$4,$5}' "$armout/td_phases.txt" || true
}

run_arm A legacy
run_arm B native
run_arm C native

# ---- the comparison ----
# Per-type call lists AND per-type VCFs AND the merged VCF. The call list is the transduction
# step's own output (candidate_disc_filtered_cns2.txt) and is what moves first if the gather
# moved; the VCFs are what the plan actually gates on.
echo
echo "##############################################################"
echo "# STAGE 1 GATE RESULT"
echo "##############################################################"

RC=0
compare_arms() {   # $1 = label, $2 = arm x, $3 = arm y
    # `local` with no cross-references, unlike run_arm's -- see the note there.
    local label=$1 x=$2 y=$3 n_same=0 n_diff=0 rel f a b
    echo
    echo "--- $label ($x vs $y) ---"
    # Enumerate from arm x, so a file arm y failed to produce shows up as a difference.
    while IFS= read -r a; do
        rel=${a#"$OUT/$x/"}
        b="$OUT/$y/$rel"
        if [ ! -f "$b" ]; then
            echo "  MISSING in $y: $rel"; n_diff=$((n_diff+1)); continue
        fi
        if cmp -s "$a" "$b"; then n_same=$((n_same+1)); else
            echo "  DIFFERS: $rel"; n_diff=$((n_diff+1))
        fi
    done < <(find "$OUT/$x" \( -name '*.vcf' -o -name 'candidate_disc_filtered_cns2.txt' \
                               -o -name 'candidate_sibling_transduction2.txt' \) -type f | sort)
    echo "  identical: $n_same    differing: $n_diff"
    [ "$n_diff" -eq 0 ] || RC=1
}

# FIRST: prove the arms actually differ in the way the gate claims. Nothing else here can tell
# a real agreement result from three legacy runs agreeing with each other -- the flag reaches the
# transduction step through an INHERITED ENV VAR, so a wiring break is silent and passes.
echo
echo "--- arms ran the paths they claim (checked before believing any comparison) ---"
for arm in A B C; do
    n_legacy=$(awk -F'\t' '$4=="clip_gather"' "$OUT/$arm/td_phases.txt" 2>/dev/null | wc -l)
    n_native=$(awk -F'\t' '$4=="native_gather"' "$OUT/$arm/td_phases.txt" 2>/dev/null | wc -l)
    echo "  arm $arm: native_gather phases=$n_native  clip_gather phases=$n_legacy"
    case "$arm" in
      A) [ "$n_legacy" -gt 0 ] && [ "$n_native" -eq 0 ] || { echo "    FATAL: arm A is not the legacy path"; RC=1; } ;;
      *) [ "$n_native" -gt 0 ] && [ "$n_legacy" -eq 0 ] || { echo "    FATAL: arm $arm is not the native path"; RC=1; } ;;
    esac
done

compare_arms "AGREEMENT  legacy vs native" A B
compare_arms "DETERMINISM native vs native" B C

echo
echo "--- call counts per arm (sanity: the arms must not all be empty) ---"
n_records() { [ -f "$1" ] || { echo 0; return 0; }; grep -vc '^#' "$1" || true; }
first_match() { local m; for m in $1; do [ -e "$m" ] && { printf '%s' "$m"; return 0; }; done; return 0; }
for arm in A B C; do
    S="$OUT/$arm/$SAMPLE"
    line="  arm $arm:"
    for pair in "L1:LINE1" "Alu:ALU" "SVA:SVA"; do
        f=$(first_match "$S/merged/*_${pair##*:}.vcf")
        [ -n "$f" ] && line="$line ${pair%%:*}=$(n_records "$f")"
    done
    mf=$(first_match "$S/*.vcf")
    [ -n "$mf" ] && line="$line merged=$(n_records "$mf")"
    echo "$line"
done
# All-empty arms would compare "identical" and prove nothing.
TOTAL_B=$(first_match "$OUT/B/$SAMPLE/*.vcf")
if [ -z "$TOTAL_B" ] || [ "$(n_records "$TOTAL_B")" -eq 0 ]; then
    echo "FATAL: arm B produced no calls -- the agreement result above is vacuous"; RC=1
fi

echo
echo "--- gather cost, legacy vs native ---"
# clip_gather + disc_gather (arm A) against native_gather (arms B/C), per rep type.
for arm in A B C; do
    [ -f "$OUT/$arm/td_phases.txt" ] || continue
    awk -v arm="$arm" -F'\t' '
        $4=="clip_gather"||$4=="disc_gather"||$4=="native_gather" {g[$3]+=$5}
        $4~/^TOTAL/ {t[$3]=$5}
        END {for (r in t) printf "  arm %s rtype=%-3s gather %7.1fs of %7.1fs (%4.1f%%)\n",
              arm, r, g[r], t[r], (t[r]>0 ? g[r]/t[r]*100 : 0)}' "$OUT/$arm/td_phases.txt"
done
echo
for arm in A B C; do
    [ -f "$OUT/$arm/wall_s.txt" ] && echo "  arm $arm whole-run wall: $(cat "$OUT/$arm/wall_s.txt")s"
done

# ---- the downstream-from-cns gate, native gather forced on ----
echo
echo "--- tests/test_downstream_from_cns.py (NA12878 chr1 golden, native gather) ---"
set +e
XTEA_SCAN_TD_GATHER="$WT/src/csrc/xtea_scan" \
XTEA2_SCAN_BIN="$WT/src/csrc/xtea_scan" \
    python -m pytest tests/test_downstream_from_cns.py -q > "$OUT/downstream_from_cns.log" 2>&1
PYRC=$?
set -e
tail -25 "$OUT/downstream_from_cns.log"
[ "$PYRC" -eq 0 ] || { echo "  FAIL (exit $PYRC) -- full log: $OUT/downstream_from_cns.log"; RC=1; }

echo
if [ "$RC" -eq 0 ]; then echo "STAGE 1 GATE: PASS"; else echo "STAGE 1 GATE: FAIL"; fi
echo "finished $(date '+%F %T')"
exit "$RC"
