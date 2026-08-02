#!/bin/bash
#SBATCH -c 24
#SBATCH -t 10:00:00
#SBATCH --mem=64G
#SBATCH -p priopark
#SBATCH --account=park_contrib
#SBATCH -o /n/scratch/users/c/cos689/xtea2_testdata/stage1_flip_wgs_%j.out
#SBATCH --mail-type=NONE

# =============================================================================
# STAGE 1 DEFAULT-FLIP GATE -- the one thing HANDOFF_2026-08-01_stage1_gather.md
# leaves open: "what remains before flipping [--td-gather-scan to default] is one
# WGS run whose calls match a BASELINE-0 arm -- the chr1 gate does not substitute,
# being 17 calls."
#
#     sbatch docs/handoff/runs/stage1_flip_wgs.sbatch
#
# TWO arms at HEAD, not one, and the reason is a confound the single-arm version
# cannot remove. BASELINE-0 was frozen at 9a43b80; HEAD is five commits later. A
# lone native-at-HEAD arm compared against BASELINE-0 conflates "the gather is
# faithful" with "nothing else moved in those five commits". So:
#
#   arm LEGACY  at HEAD, legacy gather  -- vs BASELINE-0: isolates commit drift
#   arm NATIVE  at HEAD, native gather  -- vs LEGACY:     isolates the gather
#                                       -- vs BASELINE-0: the flip criterion itself
#
# The five commits touch only src/csrc/td_gather.hpp (new), the xtea_scan CLI, and
# additive gated hooks in run_xtea.py / x_clip_disc_filter.py / x_transduction.py,
# so LEGACY == BASELINE-0 is EXPECTED. Expected is not verified, and this is the
# last stage where a byte-diff is a valid gate -- spend it.
#
# No determinism arm here: determinism was established at WGS scale on the gather's
# own artifacts (job 49453857, two native runs byte-identical on all three files in
# both families) and on the chr1 pipeline (job 49453032, arms B vs C). A third
# ~50-minute WGS arm would re-measure what is already measured.
#
# Frozen worktree, not the live checkout -- run_xtea spawns a fresh x_TEA_main per
# step, so a job sharing the checkout picks up mid-run edits.
# =============================================================================
set -euo pipefail

SCRATCH=/n/scratch/users/c/cos689/xtea2_testdata
WT="${XTEA2_WORKTREE:-$SCRATCH/stage1_flip_wt}"
OUT="${XTEA2_OUT:-$SCRATCH/stage1_flip_wgs_${SLURM_JOB_ID:-manual}}"
BASE0="${XTEA2_BASELINE0:-$SCRATCH/baseline0_hg002_L1-Alu-SVA}"
TESTDATA=/n/data1/hms/dbmi/park/corinne/claude_sandbox/test_data/HG002
BAM="$TESTDATA/alignment/HG002.novaseq.pcr-free.30x.dedup.grch38.bam"
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
[ -z "$(git -C "$WT" status --porcelain)" ] || { echo "FATAL: worktree is dirty"; exit 1; }
grep -q "td-gather-scan" "$WT/xtea/run_xtea.py" || {
    echo "FATAL: worktree has no --td-gather-scan; wrong commit"; exit 1; }
"$WT/src/csrc/xtea_scan" --help 2>&1 | grep -q -- "--td-sites" || {
    echo "FATAL: $WT/src/csrc/xtea_scan predates --td-sites; rebuild it"; exit 1; }
# The comparator has to exist BEFORE ~2 hours of compute, not after.
[ -d "$BASE0/BASELINE-0" ] || { echo "FATAL: $BASE0/BASELINE-0 missing"; exit 1; }
[ -d "$BASE0/$SAMPLE/L1" ] || {
    echo "FATAL: $BASE0/$SAMPLE/L1 missing -- BASELINE-0's per-type intermediates were"
    echo "       cleaned up. Only the VCF comparison would be possible; re-run"
    echo "       baseline0_hg002.sbatch with --keep-intermediates first."; exit 1; }

mkdir -p "$OUT"
cd "$WT"

echo "=============================================================="
echo "STAGE 1 DEFAULT-FLIP GATE  --  HG002 WGS, types=$TYPES"
echo "  job      : ${SLURM_JOB_ID:-<none>} on $(hostname)  cores=$CORES"
echo "  commit   : $(git rev-parse --short HEAD) $(git log -1 --format=%s | cut -c1-60)"
echo "  bam      : $BAM"
echo "  baseline : $BASE0  ($(grep -m1 '^commit:' "$BASE0/BASELINE-0/PROVENANCE.txt" || true))"
echo "  out      : $OUT"
echo "  started  : $(date '+%F %T')"
echo "=============================================================="

run_arm() {   # $1 = arm name, $2 = "native"|"legacy"
    # Separate statements, not `local arm=$1 armout="$OUT/$arm"`: the shell expands the
    # WHOLE command line before `local` assigns anything, so $arm is still unset there --
    # which under `set -u` kills the job before the first arm runs (job 49452897).
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
    # [TD-PHASE] lands in the STEP's stderr. Under --parallel-downstream that is
    # <sample>_<TYPE>.steps.log, NOT <sample>_xtea2.log -- glob both.
    grep -ah '^\[TD-PHASE\]' "$armout/$SAMPLE/${SAMPLE}_xtea2.log" \
        "$armout/$SAMPLE/${SAMPLE}"_*.steps.log 2>/dev/null \
        > "$armout/td_phases.txt" || true
    awk -F'\t' '{printf "  rtype=%-3s %-40s %8.1fs\n",$3,$4,$5}' "$armout/td_phases.txt" || true
}

run_arm LEGACY legacy
run_arm NATIVE native

echo
echo "##############################################################"
echo "# STAGE 1 DEFAULT-FLIP GATE RESULT"
echo "##############################################################"

RC=0

# FIRST: prove the arms actually ran different paths. Nothing below can tell a real
# agreement result from two legacy runs agreeing with each other -- the flag reaches the
# transduction step through an INHERITED ENV VAR, so a wiring break is silent and passes.
echo
echo "--- arms ran the paths they claim (checked before believing any comparison) ---"
for arm in LEGACY NATIVE; do
    n_legacy=$(awk -F'\t' '$4=="clip_gather"' "$OUT/$arm/td_phases.txt" 2>/dev/null | wc -l)
    n_native=$(awk -F'\t' '$4=="native_gather"' "$OUT/$arm/td_phases.txt" 2>/dev/null | wc -l)
    echo "  arm $arm: native_gather phases=$n_native  clip_gather phases=$n_legacy"
    case "$arm" in
      LEGACY) [ "$n_legacy" -gt 0 ] && [ "$n_native" -eq 0 ] || { echo "    FATAL: arm LEGACY is not the legacy path"; RC=1; } ;;
      *)      [ "$n_native" -gt 0 ] && [ "$n_legacy" -eq 0 ] || { echo "    FATAL: arm NATIVE is not the native path"; RC=1; } ;;
    esac
done

# --- comparison 1 & 2: the two arms' full artifact sets, against each other ---
compare_arms() {   # $1 = label, $2 = arm x dir, $3 = arm y dir
    local label=$1 x=$2 y=$3 n_same=0 n_diff=0 rel a b
    echo
    echo "--- $label ---"
    while IFS= read -r a; do
        rel=${a#"$x/"}
        b="$y/$rel"
        if [ ! -f "$b" ]; then
            echo "  MISSING in $(basename "$y"): $rel"; n_diff=$((n_diff+1)); continue
        fi
        if cmp -s "$a" "$b"; then n_same=$((n_same+1)); else
            echo "  DIFFERS: $rel"; n_diff=$((n_diff+1))
        fi
    done < <(find "$x" \( -name '*.vcf' -o -name 'candidate_disc_filtered_cns2.txt' \
                          -o -name 'candidate_sibling_transduction2.txt' \
                          -o -name '*.high_confident.post_filtering.txt' \) -type f | sort)
    echo "  identical: $n_same    differing: $n_diff"
    [ "$n_diff" -eq 0 ] || RC=1
}

compare_arms "AGREEMENT  arm LEGACY vs arm NATIVE (isolates the gather)" \
             "$OUT/LEGACY" "$OUT/NATIVE"

# --- comparison 3: each arm's VCFs against the frozen BASELINE-0 VCFs ---
# BASELINE-0/ holds only VCFs (that is what was frozen), so this leg is VCF-only by
# construction. The per-type call lists are covered against BASELINE-0's surviving
# intermediates separately below.
compare_vs_baseline0() {   # $1 = arm
    local arm=$1 n_same=0 n_diff=0 f b
    echo
    echo "--- arm $arm vs frozen BASELINE-0 VCFs ---"
    for f in "$BASE0"/BASELINE-0/*.vcf; do
        b=$(find "$OUT/$arm/$SAMPLE" -name "$(basename "$f")" -type f | head -1)
        if [ -z "$b" ]; then
            echo "  MISSING in arm $arm: $(basename "$f")"; n_diff=$((n_diff+1)); continue
        fi
        if cmp -s "$f" "$b"; then n_same=$((n_same+1)); else
            echo "  DIFFERS: $(basename "$f")"; n_diff=$((n_diff+1))
            # A call-level summary, because "the VCF differs" on 1493 calls is not a finding.
            # Compare on CHROM/POS/ALT only: everything else (INFO, genotype) can move for
            # reasons this gate is not about.
            comm -3 <(grep -v '^#' "$f" | cut -f1,2,5 | sort) \
                    <(grep -v '^#' "$b" | cut -f1,2,5 | sort) \
                | head -40
            echo "    baseline-only: $(comm -23 <(grep -v '^#' "$f" | cut -f1,2,5 | sort) <(grep -v '^#' "$b" | cut -f1,2,5 | sort) | wc -l)"
            echo "    arm-only:      $(comm -13 <(grep -v '^#' "$f" | cut -f1,2,5 | sort) <(grep -v '^#' "$b" | cut -f1,2,5 | sort) | wc -l)"
        fi
    done
    echo "  identical: $n_same    differing: $n_diff"
    [ "$n_diff" -eq 0 ] || RC=1
}

compare_vs_baseline0 LEGACY     # isolates drift across the five commits since 9a43b80
compare_vs_baseline0 NATIVE     # THE FLIP CRITERION

# --- comparison 4: per-type call lists vs BASELINE-0's surviving intermediates ---
# The transduction step's own output moves first if the gather moved, well before a VCF
# does; a VCF-only comparison can hide a difference that post-filtering happened to absorb.
echo
echo "--- per-type call lists vs BASELINE-0 intermediates ---"
for arm in LEGACY NATIVE; do
    n_same=0; n_diff=0
    for t in L1 Alu SVA; do
        for fn in candidate_disc_filtered_cns2.txt \
                  candidate_disc_filtered_cns.txt.high_confident.post_filtering.txt; do
            a="$BASE0/$SAMPLE/$t/$fn"; b="$OUT/$arm/$SAMPLE/$t/$fn"
            [ -f "$a" ] || continue
            if [ ! -f "$b" ]; then echo "  MISSING in $arm: $t/$fn"; n_diff=$((n_diff+1)); continue; fi
            if cmp -s "$a" "$b"; then n_same=$((n_same+1)); else
                echo "  DIFFERS ($arm): $t/$fn  [$(wc -l < "$a") vs $(wc -l < "$b") lines]"
                n_diff=$((n_diff+1))
            fi
        done
    done
    echo "  arm $arm: identical $n_same, differing $n_diff"
    [ "$n_diff" -eq 0 ] || RC=1
done

echo
echo "--- call counts ---"
n_records() { [ -f "$1" ] || { echo 0; return 0; }; grep -vc '^#' "$1" || true; }
first_match() { local m; for m in $1; do [ -e "$m" ] && { printf '%s' "$m"; return 0; }; done; return 0; }
for arm in LEGACY NATIVE; do
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
line="  BASELINE-0:"
for pair in "L1:LINE1" "Alu:ALU" "SVA:SVA"; do
    f=$(first_match "$BASE0/BASELINE-0/*_${pair##*:}.vcf")
    [ -n "$f" ] && line="$line ${pair%%:*}=$(n_records "$f")"
done
mf=$(first_match "$BASE0/BASELINE-0/$SAMPLE.vcf")
[ -n "$mf" ] && line="$line merged=$(n_records "$mf")"
echo "$line"
# Two empty arms compare "identical" and prove nothing.
TOTAL_N=$(first_match "$OUT/NATIVE/$SAMPLE/*.vcf")
if [ -z "$TOTAL_N" ] || [ "$(n_records "$TOTAL_N")" -eq 0 ]; then
    echo "FATAL: arm NATIVE produced no calls -- every agreement result above is vacuous"; RC=1
fi

echo
echo "--- the speedup this flip is for (WGS, both arms at n=$CORES on one node) ---"
for arm in LEGACY NATIVE; do
    [ -f "$OUT/$arm/td_phases.txt" ] || continue
    awk -v arm="$arm" -F'\t' '
        $4=="clip_gather"||$4=="disc_gather"||$4=="native_gather" {g[$3]+=$5}
        $4~/^TOTAL/ {t[$3]=$5}
        END {for (r in t) printf "  arm %-6s rtype=%-3s gather %8.1fs of %8.1fs (%4.1f%%)\n",
              arm, r, g[r], t[r], (t[r]>0 ? g[r]/t[r]*100 : 0)}' "$OUT/$arm/td_phases.txt"
done
echo
for arm in LEGACY NATIVE; do
    awk -F'\t' '$1=="phase"{printf "  '"$arm"' %-22s %8.1fs\n",$2,$3}' \
        "$OUT/$arm/$SAMPLE/step_timings.tsv" 2>/dev/null || true
done
echo
for arm in LEGACY NATIVE; do
    [ -f "$OUT/$arm/wall_s.txt" ] && echo "  arm $arm whole-run wall: $(cat "$OUT/$arm/wall_s.txt")s"
done

echo
if [ "$RC" -eq 0 ]; then
    echo "STAGE 1 FLIP GATE: PASS -- --td-gather-scan may become the default."
    echo "  Note what flipping ALSO does (handoff section 4): it makes the multi-BAM"
    echo "  divergence reachable, so the mosaic path would run the fixed behaviour with"
    echo "  no baseline of its own. That is a Stage 4 concern, not a blocker here."
else
    echo "STAGE 1 FLIP GATE: FAIL -- do not flip the default."
fi
echo "finished $(date '+%F %T')"
exit "$RC"
