#!/bin/bash
#SBATCH -c 12
#SBATCH -t 6:00:00
#SBATCH --mem=48G
#SBATCH -p priopark
#SBATCH --account=park_contrib
#SBATCH -o /n/scratch/users/c/cos689/xtea2_testdata/refactor_gate_chr1_%j.out
#SBATCH --mail-type=NONE

# =============================================================================
# REFACTOR AGREEMENT GATE -- chr1, both gather paths, against a PRIOR gate's arms.
#
#     sbatch docs/handoff/runs/refactor_gate_chr1.sbatch <worktree> [reference_run_dir]
#
# For a change that is supposed to move nothing: run HG002 chr1 twice at the new
# commit (legacy gather, native gather) and cmp every call list, sibling list and
# VCF against the corresponding arm of an earlier PASSING gate run.
#
# Why compare against a stored run instead of running a before-arm here. The
# reference arms already exist -- stage1_gate_chr1_49453032/{A,B}, the Stage 1 gate
# -- and re-running them would cost two more chr1 arms to reproduce bytes that are
# already on disk under a recorded commit. The cost is that the comparison is only
# valid while the Python pipeline between the two commits is otherwise unchanged,
# which is not something a script can check for you. IT IS CHECKED BY HAND, and the
# check belongs in the commit message that cites this gate:
#
#     git diff <ref_commit>..HEAD -- xtea/ src/
#
# must contain nothing but the change under test. If it doesn't, run the full
# stage1_gate_chr1.sbatch instead -- it builds its own before-arm.
#
# BOTH gather paths, not just one. A refactor inside XClipDiscFilter/XTransduction
# is reached by both, and the native path skips code the legacy path runs.
# =============================================================================
set -euo pipefail

SCRATCH=/n/scratch/users/c/cos689/xtea2_testdata
WT="${1:?usage: sbatch refactor_gate_chr1.sbatch <worktree> [reference_run_dir]}"
REF_RUN="${2:-$SCRATCH/stage1_gate_chr1_49453032}"
OUT="${XTEA2_OUT:-$SCRATCH/refactor_gate_chr1_${SLURM_JOB_ID:-manual}}"
TESTDATA=/n/data1/hms/dbmi/park/corinne/claude_sandbox/test_data/HG002
BAM="$TESTDATA/chr1/HG002.chr1.bam"
SAMPLE=HG002.chr1
TYPES="${XTEA2_TYPES:-L1,Alu,SVA}"
REF=/n/data1/hms/dbmi/park-smaht_dac/ref/GRCh38_no_alt/hg38_no_alt.fa
CORES="${SLURM_CPUS_PER_TASK:-12}"

source /n/app/conda/miniforge3/24.11.3-0/etc/profile.d/conda.sh
conda activate xtea2
export PATH=/home/cos689/software/bin:$PATH
export XTEA_BASIC_INFO_SEED=42
test "$PYTHONNOUSERSITE" = 1

[ -d "$WT" ]                    || { echo "FATAL: worktree $WT missing"; exit 1; }
[ -x "$WT/src/csrc/xtea_scan" ] || { echo "FATAL: $WT/src/csrc/xtea_scan not built"; exit 1; }
[ -f "$BAM" ]                   || { echo "FATAL: $BAM missing"; exit 1; }
[ -z "$(git -C "$WT" status --porcelain)" ] || { echo "FATAL: worktree is dirty"; exit 1; }
[ -d "$REF_RUN/A" ] && [ -d "$REF_RUN/B" ] || {
    echo "FATAL: $REF_RUN needs arms A (legacy) and B (native)"; exit 1; }

mkdir -p "$OUT"
cd "$WT"

echo "=============================================================="
echo "REFACTOR GATE  --  HG002 chr1, types=$TYPES"
echo "  job    : ${SLURM_JOB_ID:-<none>} on $(hostname)  cores=$CORES"
echo "  commit : $(git rev-parse --short HEAD) $(git log -1 --format=%s | cut -c1-60)"
echo "  ref run: $REF_RUN"
echo "  out    : $OUT"
echo "  started: $(date '+%F %T')"
echo "=============================================================="

run_arm() {   # $1 = arm name, $2 = "native"|"legacy"
    # Separate statements, not `local arm=$1 armout=...`: under `set -u` the whole line
    # expands before `local` assigns, and $arm is still unset (job 49452897).
    local arm=$1 mode=$2 t0 wall armout
    armout="$OUT/$arm"
    rm -rf "$armout"; mkdir -p "$armout"
    echo
    echo "=== arm $arm: gather=$mode ==="
    t0=$(date +%s)
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
    grep -ah '^\[TD-PHASE\]' "$armout/$SAMPLE/${SAMPLE}_xtea2.log" \
        "$armout/$SAMPLE/${SAMPLE}"_*.steps.log 2>/dev/null \
        > "$armout/td_phases.txt" || true
}

run_arm A legacy
run_arm B native

echo
echo "##############################################################"
echo "# REFACTOR GATE RESULT"
echo "##############################################################"
RC=0

# The arms must have run different paths, or "A matches the reference A" could be a
# native run matching a native reference and the legacy leg would be untested.
echo
echo "--- arms ran the paths they claim ---"
for arm in A B; do
    n_legacy=$(awk -F'\t' '$4=="clip_gather"' "$OUT/$arm/td_phases.txt" 2>/dev/null | wc -l)
    n_native=$(awk -F'\t' '$4=="native_gather"' "$OUT/$arm/td_phases.txt" 2>/dev/null | wc -l)
    echo "  arm $arm: native_gather=$n_native  clip_gather=$n_legacy"
    case "$arm" in
      A) [ "$n_legacy" -gt 0 ] && [ "$n_native" -eq 0 ] || { echo "    FATAL: arm A is not legacy"; RC=1; } ;;
      B) [ "$n_native" -gt 0 ] && [ "$n_legacy" -eq 0 ] || { echo "    FATAL: arm B is not native"; RC=1; } ;;
    esac
done

compare_to_ref() {   # $1 = arm
    local arm=$1 n_same=0 n_diff=0 n_new=0 rel a b
    echo
    echo "--- arm $arm vs $REF_RUN/$arm ---"
    while IFS= read -r a; do
        rel=${a#"$REF_RUN/$arm/"}
        b="$OUT/$arm/$rel"
        if [ ! -f "$b" ]; then
            echo "  MISSING now: $rel"; n_diff=$((n_diff+1)); continue
        fi
        if cmp -s "$a" "$b"; then n_same=$((n_same+1)); else
            echo "  DIFFERS: $rel"; n_diff=$((n_diff+1))
            diff <(grep -v '^#' "$a" | cut -f1,2) <(grep -v '^#' "$b" | cut -f1,2) | head -20 || true
        fi
    done < <(find "$REF_RUN/$arm" \( -name '*.vcf' -o -name 'candidate_disc_filtered_cns2.txt' \
                                     -o -name 'candidate_sibling_transduction2.txt' \) -type f | sort)
    # A file the new run produces and the reference does not is also a change.
    while IFS= read -r b; do
        rel=${b#"$OUT/$arm/"}
        [ -f "$REF_RUN/$arm/$rel" ] || { echo "  NEW file: $rel"; n_new=$((n_new+1)); }
    done < <(find "$OUT/$arm" \( -name '*.vcf' -o -name 'candidate_disc_filtered_cns2.txt' \
                                 -o -name 'candidate_sibling_transduction2.txt' \) -type f | sort)
    echo "  identical: $n_same    differing: $n_diff    new: $n_new"
    [ "$n_diff" -eq 0 ] && [ "$n_new" -eq 0 ] || RC=1
    # An empty comparison passes vacuously; refuse it.
    [ "$n_same" -gt 0 ] || { echo "  FATAL: nothing was compared"; RC=1; }
}

compare_to_ref A
compare_to_ref B

echo
echo "--- call counts ---"
n_records() { [ -f "$1" ] || { echo 0; return 0; }; grep -vc '^#' "$1" || true; }
first_match() { local m; for m in $1; do [ -e "$m" ] && { printf '%s' "$m"; return 0; }; done; return 0; }
for d in "$REF_RUN/A" "$REF_RUN/B" "$OUT/A" "$OUT/B"; do
    S="$d/$SAMPLE"; line="  $d:"
    for pair in "L1:LINE1" "Alu:ALU" "SVA:SVA"; do
        f=$(first_match "$S/merged/*_${pair##*:}.vcf")
        [ -n "$f" ] && line="$line ${pair%%:*}=$(n_records "$f")"
    done
    mf=$(first_match "$S/*.vcf"); [ -n "$mf" ] && line="$line merged=$(n_records "$mf")"
    echo "$line"
done

echo
if [ "$RC" -eq 0 ]; then echo "REFACTOR GATE: PASS"; else echo "REFACTOR GATE: FAIL"; fi
echo "finished $(date '+%F %T')"
exit "$RC"
