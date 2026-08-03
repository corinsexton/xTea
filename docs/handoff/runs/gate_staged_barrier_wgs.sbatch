#!/bin/bash
#SBATCH -c 24
#SBATCH -t 10:00:00
#SBATCH --mem=64G
#SBATCH -p priopark
#SBATCH --account=park_contrib
#SBATCH -o /n/scratch/users/c/cos689/xtea2_testdata/gate_staged_barrier_wgs_%j.out
#SBATCH --mail-type=NONE
# =============================================================================
# THE STAGED cns_filter -> transduction BARRIER at WGS scale. THREE arms, ONE commit.
#
#     sbatch docs/handoff/runs/gate_staged_barrier_wgs.sbatch <worktree> [types]
#
#   OFF     no barrier at all -- the production default, and the reference for both
#           correctness (vs frozen BASELINE-0) and timing.
#   GLOBAL  --td-gather-driver --global-barrier: the join, i.e. exactly what jobs 49548949
#           and 49550925 gated.
#   STAGED  --td-gather-driver: the same work with per-family waits (run_xtea.StagedBarrier).
#
# WHY WGS AND NOT ONLY chr1. Two independent reasons, and they are different reasons:
#
#   a. **The last bug in this wave split was invisible on chr1.** Job 49498881 died at
#      `SVA:transduction` with a bare KeyError because a step opened an input it was never
#      handed, and chr1's smaller candidate sets never reached it. This change rewrites how
#      the split is *scheduled*, so it earns that leg the same way the flags themselves did.
#   b. **chr1 cannot show the saving.** The whole cost being recovered is that every family
#      waited for Alu's cns_filter, which at WGS is 124.1s against L1's 57.3 and SVA's 62.5.
#      At chr1 those are seconds and the delta is inside the noise. The chr1 gate therefore
#      asserts the *structure* (the flankless family stopped waiting) and this one adds the
#      *magnitude*.
#
# WHAT MUST BE TRUE:
#   0. the arms ran what they claim. Nothing below distinguishes a real result from three
#      arms that all quietly did the same thing.
#   1. arm OFF == frozen BASELINE-0. Without it a clean STAGED-vs-GLOBAL diff could be two
#      equally-wrong arms -- the reasoning that earned this leg in job 49550925.
#   2. GLOBAL == OFF and STAGED == OFF on every VCF / cns2 / sibling list.
#   3. STAGED == GLOBAL on the universe artifacts themselves. Check 2 compares CALLS, and
#      the universe pass is diagnostic -- it cannot move one, so a staged enumeration that
#      produced a different or empty table would pass check 2 in silence.
#   4. in BOTH barrier arms, the driver's replica == the step's own `.slct`, and the step
#      skipped its own gather. Alu has no step `.slct` (no transduction flank), so the
#      expected count is the FLANKED count -- printed, because a gate whose only check runs
#      at comparison time cannot tell "the code is wrong" from "nothing was compared"
#      (job 49516776).
#   5. the flankless family stopped waiting for a gather it never reads -- structural, see
#      the chr1 gate's check 5 for why it is measured this way and not as "did anything
#      overlap".
#   6. **the magnitude**: STAGED's downstream_parallel must beat GLOBAL's by more than
#      noise. Predicted from job 49550925's decomposition: OFF 462.7s, GLOBAL 505.3s,
#      STAGED ~449s (OFF minus the hoist's own -13.7s), i.e. a ~55s recovery. Gated at >20s
#      so a genuine regression fails and a few seconds of jitter does not.
#
# Frozen worktree, not the live checkout -- run_xtea spawns a fresh x_TEA_main per step, so
# a job sharing the checkout picks up mid-run edits.
# =============================================================================
set -euo pipefail

SCRATCH=/n/scratch/users/c/cos689/xtea2_testdata
WT="${1:?usage: sbatch gate_staged_barrier_wgs.sbatch <worktree> [types]}"
TYPES="${2:-L1,Alu,SVA}"
OUT="${XTEA2_OUT:-$SCRATCH/gate_staged_barrier_wgs_${SLURM_JOB_ID:-manual}}"
BASE0="${XTEA2_BASELINE0:-$SCRATCH/baseline0_hg002_L1-Alu-SVA}"
TESTDATA=/n/data1/hms/dbmi/park/corinne/claude_sandbox/test_data/HG002
BAM="$TESTDATA/alignment/HG002.novaseq.pcr-free.30x.dedup.grch38.bam"
SAMPLE=HG002.wgs
REF=/n/data1/hms/dbmi/park-smaht_dac/ref/GRCh38_no_alt/hg38_no_alt.fa
CORES="${SLURM_CPUS_PER_TASK:-24}"

source /n/app/conda/miniforge3/24.11.3-0/etc/profile.d/conda.sh
conda activate xtea2
export PATH=/home/cos689/software/bin:$PATH
export XTEA_BASIC_INFO_SEED=42   # pin the insert-size draw, as BASELINE-0 did. NOT optional
                                 # here: the driver's gather window and the step's evaluation
                                 # are computed from separate basic-info files that agree only
                                 # because the draw is seeded (downstream_unified.
                                 # transduction_extnd documents this).
test "$PYTHONNOUSERSITE" = 1     # the env's guard against ~/.local's onnxruntime

[ -d "$WT" ]                    || { echo "FATAL: worktree $WT missing"; exit 1; }
[ -x "$WT/src/csrc/xtea_scan" ] || { echo "FATAL: $WT/src/csrc/xtea_scan not built"; exit 1; }
[ -f "$BAM" ]                   || { echo "FATAL: $BAM missing"; exit 1; }
[ -f "$BAM.bai" ]               || { echo "FATAL: $BAM.bai missing"; exit 1; }
[ -z "$(git -C "$WT" status --porcelain)" ] || { echo "FATAL: worktree is dirty"; exit 1; }
# The comparator and the flags have to exist BEFORE ~2 hours of compute, not after.
grep -q "global-barrier" "$WT/xtea/run_xtea.py" || {
    echo "FATAL: worktree has no --global-barrier; wrong commit"; exit 1; }
grep -q "class StagedBarrier" "$WT/xtea/run_xtea.py" || {
    echo "FATAL: worktree has no StagedBarrier; wrong commit"; exit 1; }
[ -d "$BASE0/BASELINE-0" ]      || { echo "FATAL: $BASE0/BASELINE-0 missing"; exit 1; }

mkdir -p "$OUT"
cd "$WT"

echo "=============================================================="
echo "STAGED BARRIER GATE  --  HG002 WGS, types=$TYPES"
echo "  job      : ${SLURM_JOB_ID:-<none>} on $(hostname)  cores=$CORES"
echo "  commit   : $(git rev-parse --short HEAD) $(git log -1 --format=%s | cut -c1-60)"
echo "  baseline : $BASE0"
echo "  out      : $OUT"
echo "  started  : $(date '+%F %T')"
echo "=============================================================="

run_arm() {   # $1 = arm name, $2... = extra driver flags
    # Separate statements, not `local arm=$1 armout=...`: bash expands the whole `local`
    # line before assigning any of it, so $arm would be a leaked global (jobs 49452897,
    # 49516776 -- the second made a passing change report FAIL).
    local arm=$1; shift
    local t0 wall armout
    armout="$OUT/$arm"
    rm -rf "$armout"; mkdir -p "$armout"
    echo
    echo "=== arm $arm: extra flags: $* ==="
    t0=$(date +%s)
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
      --td-gather-scan \
      --parallel-downstream \
      --keep-intermediates \
      -n "$CORES" \
      "$@" \
      > "$armout/driver.log" 2>&1 || {
          echo "FAIL: arm $arm exited $?"; tail -60 "$armout/driver.log"; return 1; }
    wall=$(( $(date +%s) - t0 ))
    printf '  wall %02d:%02d:%02d\n' $((wall/3600)) $((wall%3600/60)) $((wall%60))
    echo "$wall" > "$armout/wall_s.txt"
}

# XTEA2_ANALYZE_ONLY=1 re-runs the checks against arms already in $OUT. For iterating on
# THIS SCRIPT, not for skipping the compute on a real gate.
if [ "${XTEA2_ANALYZE_ONLY:-0}" = 1 ]; then
    echo; echo "*** ANALYZE-ONLY: reusing existing arms in $OUT ***"
    for arm in OFF GLOBAL STAGED; do
        [ -f "$OUT/$arm/$SAMPLE/$SAMPLE.vcf" ] || {
            echo "FATAL: $OUT/$arm has no completed run"; exit 1; }
    done
else
    run_arm OFF                                       || exit 1
    run_arm GLOBAL --td-gather-driver --global-barrier || exit 1
    run_arm STAGED --td-gather-driver                  || exit 1
fi

echo
echo "##############################################################"
echo "# STAGED BARRIER GATE RESULT"
echo "##############################################################"
RC=0
FLANKED=""
for t in ${TYPES//,/ }; do
    case "$t" in L1|SVA) FLANKED="$FLANKED $t" ;; esac
done
n_exp=$(echo $FLANKED | wc -w)

# --- 0. the arms ran what they claim -------------------------------------------------
echo
echo "--- 0. arms ran what they claim ---"
for arm in OFF GLOBAL STAGED; do
    L="$OUT/$arm/driver.log"
    n_uni=$(grep -ac 'locus universe (cns_filter' "$L" || true)
    n_gth=$(grep -ac 'td gather (driver):' "$L" || true)
    n_stg=$(grep -ac 'barrier: STAGED' "$L" || true)
    n_glb=$(grep -ac 'barrier: GLOBAL' "$L" || true)
    echo "  $arm: universe=$n_uni gather=$n_gth staged-line=$n_stg global-line=$n_glb"
    case "$arm" in
      OFF)    [ "$n_uni" -eq 0 ] && [ "$n_gth" -eq 0 ] && [ "$n_stg" -eq 0 ] && [ "$n_glb" -eq 0 ] || {
                  echo "    FATAL: OFF opened a barrier"; RC=1; } ;;
      GLOBAL) [ "$n_uni" -eq 1 ] && [ "$n_gth" -eq 1 ] && [ "$n_stg" -eq 0 ] && [ "$n_glb" -eq 1 ] || {
                  echo "    FATAL: GLOBAL arm did not run the global barrier"; RC=1; } ;;
      STAGED) [ "$n_uni" -eq 1 ] && [ "$n_gth" -eq 1 ] && [ "$n_stg" -eq 1 ] && [ "$n_glb" -eq 0 ] || {
                  echo "    FATAL: STAGED arm did not run the staged barrier"; RC=1; } ;;
    esac
done

# --- 1. arm OFF vs frozen BASELINE-0 -- isolates commit drift -------------------------
echo
echo "--- 1. arm OFF vs frozen BASELINE-0 VCFs (must be IDENTICAL) ---"
n_same=0; n_diff=0
for f in "$BASE0"/BASELINE-0/*.vcf; do
    b=$(find "$OUT/OFF/$SAMPLE" -name "$(basename "$f")" -type f | head -1)
    if [ -z "$b" ]; then echo "  MISSING: $(basename "$f")"; n_diff=$((n_diff+1)); continue; fi
    if cmp -s "$f" "$b"; then n_same=$((n_same+1)); else
        echo "  DIFFERS: $(basename "$f")"; n_diff=$((n_diff+1)); fi
done
echo "  identical: $n_same    differing: $n_diff"
[ "$n_diff" -eq 0 ] && [ "$n_same" -gt 0 ] || {
    echo "  FAIL: arm OFF is not BASELINE-0; every delta below is uninterpretable"; RC=1; }

# --- 2. both barrier arms vs OFF: no call may move -----------------------------------
echo
echo "--- 2. GLOBAL and STAGED vs OFF: every VCF / cns2 / sibling list ---"
compare_to_off() {
    local arm=$1
    local n_same=0 n_diff=0 n_new=0 rel a b
    echo "  -- $arm vs OFF --"
    while IFS= read -r a; do
        rel=${a#"$OUT/OFF/"}
        b="$OUT/$arm/$rel"
        if [ ! -f "$b" ]; then echo "    MISSING in $arm: $rel"; n_diff=$((n_diff+1)); continue; fi
        if cmp -s "$a" "$b"; then n_same=$((n_same+1)); else
            echo "    DIFFERS: $rel"; n_diff=$((n_diff+1))
            diff <(grep -v '^#' "$a" | cut -f1,2) <(grep -v '^#' "$b" | cut -f1,2) | head -20 || true
        fi
    done < <(find "$OUT/OFF" \( -name '*.vcf' -o -name 'candidate_disc_filtered_cns2.txt' \
                                -o -name 'candidate_sibling_transduction2.txt' \) -type f | sort)
    # A file an arm produces that OFF does not is expected ONLY under merged/locus_universe/.
    while IFS= read -r b; do
        rel=${b#"$OUT/$arm/"}
        case "$rel" in */merged/locus_universe/*) continue ;; esac
        [ -f "$OUT/OFF/$rel" ] || { echo "    NEW file: $rel"; n_new=$((n_new+1)); }
    done < <(find "$OUT/$arm" \( -name '*.vcf' -o -name 'candidate_disc_filtered_cns2.txt' \
                                 -o -name 'candidate_sibling_transduction2.txt' \) -type f | sort)
    echo "    identical: $n_same    differing: $n_diff    unexpected new: $n_new"
    [ "$n_diff" -eq 0 ] && [ "$n_new" -eq 0 ] || RC=1
    [ "$n_same" -gt 0 ] || { echo "    FATAL: nothing was compared"; RC=1; }
}
compare_to_off GLOBAL
compare_to_off STAGED

# --- 3. STAGED vs GLOBAL on the universe artifacts themselves ------------------------
echo
echo "--- 3. STAGED vs GLOBAL: the universe table and the per-family replicas ---"
n_u=0; n_u_ok=0
for f in "locus_universe.tsv" $(for t in $FLANKED; do echo "$t.slct"; done); do
    g="$OUT/GLOBAL/$SAMPLE/merged/locus_universe/$f"
    s="$OUT/STAGED/$SAMPLE/merged/locus_universe/$f"
    n_u=$((n_u+1))
    if [ ! -f "$g" ] || [ ! -f "$s" ]; then
        echo "  $f: FATAL: missing (GLOBAL=$([ -f "$g" ] && echo y || echo n) STAGED=$([ -f "$s" ] && echo y || echo n))"
        continue
    fi
    if cmp -s "$g" "$s"; then
        n_u_ok=$((n_u_ok+1))
        echo "  $f: IDENTICAL ($(wc -c < "$s") bytes, $(wc -l < "$s") rows)"
    else
        echo "  $f: DIFFERS"; diff "$g" "$s" | head -10 || true
    fi
done
echo "  universe artifacts identical: $n_u_ok / $n_u"
[ "$n_u" -gt 0 ] && [ "$n_u_ok" -eq "$n_u" ] || RC=1
# The offline gate's own numbers (docs/handoff/runs/gate_locus_universe.py): 41,604 loci,
# 38,815 created-only. Reported, and the row count gated -- an empty table is the failure
# that every byte comparison above would call a pass.
echo "  the universe, as the STAGED arm reported it:"
grep -aA 8 'locus universe: ' "$OUT/STAGED/driver.log" | head -6 | sed 's/^/    /' || true
n_rows=$(grep -vc '^#' "$OUT/STAGED/$SAMPLE/merged/locus_universe/locus_universe.tsv" || echo 0)
echo "  table rows: $n_rows (offline gate: 41604)"
[ "$n_rows" -gt 0 ] || { echo "  FATAL: empty universe table"; RC=1; }

# --- 4. the replica == the step's own list, and the step skipped its own gather -------
echo
echo "--- 4. in BOTH barrier arms: replica == step .slct, and the bypass fired ---"
n_cmp=0; n_ident=0
for arm in GLOBAL STAGED; do
    S="$OUT/$arm/$SAMPLE"
    for t in $FLANKED; do
        step="$S/$t/candidate_list_from_disc.txt.clip_sites_raw_disc.txt.slct"
        drv="$S/merged/locus_universe/$t.slct"
        if [ ! -f "$step" ] || [ ! -f "$drv" ]; then
            echo "  $arm/$t: FATAL: missing (step=$([ -f "$step" ] && echo y || echo n) drv=$([ -f "$drv" ] && echo y || echo n))"
            continue
        fi
        n_cmp=$((n_cmp+1))
        if cmp -s "$step" "$drv"; then
            n_ident=$((n_ident+1))
            echo "  $arm/$t: IDENTICAL ($(wc -c < "$step") bytes, $(wc -l < "$step") rows)"
        else
            echo "  $arm/$t: DIFFERS"; diff "$step" "$drv" | head -20 || true
        fi
    done
    if [ -f "$S/Alu/candidate_list_from_disc.txt.clip_sites_raw_disc.txt.slct" ]; then
        echo "  FATAL: Alu grew a step .slct in $arm; the flankless assumption is wrong"; RC=1
    fi
done
echo "  replicas compared: $n_cmp (expected $((n_exp * 2)))   identical: $n_ident"
[ "$n_cmp" -eq $((n_exp * 2)) ] && [ "$n_ident" -eq "$n_cmp" ] && [ "$n_exp" -gt 0 ] || RC=1

n_ok=0; n_bad=0
for t in $FLANKED; do
    for arm in OFF GLOBAL STAGED; do
        slog="$OUT/$arm/$SAMPLE/${SAMPLE}_${t}.steps.log"
        if [ ! -f "$slog" ]; then echo "  $arm/$t: FATAL: no steps log"; n_bad=$((n_bad+1)); continue; fi
        n_pre=$(grep -ac 'precomputed by the driver' "$slog" || true)
        n_own=$(grep -ac -- '--td-sites' "$slog" || true)
        echo "  $arm/$t: precomputed=$n_pre  own-gather=$n_own"
        case "$arm" in
          OFF) [ "$n_pre" -eq 0 ] && [ "$n_own" -eq 1 ] && n_ok=$((n_ok+1)) || {
                   echo "    FATAL: OFF arm must gather in the step"; n_bad=$((n_bad+1)); } ;;
          *)   [ "$n_pre" -eq 1 ] && [ "$n_own" -eq 0 ] && n_ok=$((n_ok+1)) || {
                   echo "    FATAL: $arm arm must bypass the step's gather"; n_bad=$((n_bad+1)); } ;;
        esac
    done
done
echo "  bypass checks passed: $n_ok (expected $((n_exp * 3)))   failed: $n_bad"
[ "$n_bad" -eq 0 ] && [ "$n_ok" -eq $((n_exp * 3)) ] || RC=1

# --- 5. the flankless family stopped waiting ------------------------------------------
echo
echo "--- 5. the flankless family stopped waiting for a gather it never reads ---"
# Structural, and deliberately not "did anything overlap" -- see the chr1 gate's check 5
# for why that version passes or fails on which family happened to be slowest.
stamp_of() {
    grep -aoE "^\[xtea2 \+ *[0-9.]+s\].*$2" "$1" 2>/dev/null | head -1 \
        | sed -E 's/^\[xtea2 \+ *([0-9.]+)s\].*/\1/' || true
}
case ",$TYPES," in
  *,Alu,*)
    gap_of() {
        local L t_enum t_td
        L="$OUT/$1/driver.log"
        t_enum=$(stamp_of "$L" 'Alu +[0-9]+ cns candidates')
        t_td=$(stamp_of "$L" '=== type Alu .*from transduction')
        [ -n "$t_enum" ] && [ -n "$t_td" ] || { printf '?'; return 0; }
        awk -v a="$t_enum" -v b="$t_td" 'BEGIN{printf "%.1f", b - a}'
    }
    g_global=$(gap_of GLOBAL); g_staged=$(gap_of STAGED)
    d_gather=$(awk -F'\t' '$2=="td_gather_driver"{print $3}' \
        "$OUT/GLOBAL/$SAMPLE/step_timings.tsv" 2>/dev/null | head -1)
    echo "  Alu enumeration -> Alu transduction:  GLOBAL=${g_global}s  STAGED=${g_staged}s"
    echo "  driver gather duration (GLOBAL arm):  ${d_gather:-?}s"
    if [ "$g_global" = "?" ] || [ "$g_staged" = "?" ] || [ -z "$d_gather" ]; then
        echo "    FATAL: could not read the stamps (log format changed?)"; RC=1
    else
        awk -v s="$g_staged" -v g="$g_global" -v d="$d_gather" 'BEGIN{exit !(s < g && s < d)}' || {
            echo "    FATAL: the flankless family is still waiting for the driver gather"; RC=1; }
        awk -v g="$g_global" -v d="$d_gather" 'BEGIN{exit !(g >= d * 0.5)}' || {
            echo "    FATAL: the GLOBAL arm's gap does not contain the gather -- the control"
            echo "           is not measuring what this check assumes"; RC=1; }
    fi
    ;;
  *) echo "  SKIPPED: no flankless family in types=$TYPES" ;;
esac

# --- 6. the magnitude ------------------------------------------------------------------
echo
echo "--- 6. STAGED recovers the join's cost (the reason the change exists) ---"
# downstream_parallel, not the total run wall: job 49550925 moved 2263 -> 2228s only
# because scan_pass varied by 66s, which nothing here can touch. The phase row is the
# quantity the critical path decomposes to (462.6 predicted vs 462.7 measured, OFF).
phase_s() {  # $1 = arm -> the downstream_parallel phase seconds
    awk -F'\t' '$1=="phase" && $2=="downstream_parallel"{print $3}' \
        "$OUT/$1/$SAMPLE/step_timings.tsv" 2>/dev/null | head -1
}
p_off=$(phase_s OFF); p_glb=$(phase_s GLOBAL); p_stg=$(phase_s STAGED)
echo "  downstream_parallel:  OFF=${p_off:-?}s  GLOBAL=${p_glb:-?}s  STAGED=${p_stg:-?}s"
echo "  (job 49550925 predicted OFF 462.7 / GLOBAL 505.3; STAGED ~449 = OFF - the hoist's 13.7)"
if [ -z "$p_glb" ] || [ -z "$p_stg" ]; then
    echo "    FATAL: could not read downstream_parallel from step_timings.tsv"; RC=1
else
    awk -v g="$p_glb" -v s="$p_stg" 'BEGIN{printf "  recovered: %.1fs\n", g - s}'
    # >20s, not >0: the predicted recovery is ~55s, and a threshold at zero would call
    # jitter a pass. A real regression -- staging that reintroduced a wait -- fails here.
    awk -v g="$p_glb" -v s="$p_stg" 'BEGIN{exit !(g - s > 20)}' || {
        echo "    FATAL: STAGED did not beat GLOBAL by more than noise"; RC=1; }
    [ -n "$p_off" ] && awk -v o="$p_off" -v s="$p_stg" \
        'BEGIN{printf "  vs OFF (no barrier at all): %+.1fs\n", s - o}'
fi
echo "  per-family cns_filter (what the join waited on):"
for arm in OFF GLOBAL STAGED; do
    for t in ${TYPES//,/ }; do
        v=$(awk -F'\t' -v k="$t:cns_filter" '$2==k{print $3}' \
            "$OUT/$arm/$SAMPLE/step_timings.tsv" 2>/dev/null | head -1)
        printf '    %-7s %-4s %ss\n' "$arm" "$t" "${v:-?}"
    done
done
echo "  total run wall (NOT the place to read this -- scan_pass varies by ~66s):"
for arm in OFF GLOBAL STAGED; do
    printf '    %-7s %ss\n' "$arm" "$(cat "$OUT/$arm/wall_s.txt" 2>/dev/null || echo '?')"
done

# --- call counts ----------------------------------------------------------------------
echo
echo "--- call counts, all three arms ---"
n_records() { [ -f "$1" ] || { echo 0; return 0; }; grep -vc '^#' "$1" || true; }
first_match() { local m; for m in $1; do [ -e "$m" ] && { printf '%s' "$m"; return 0; }; done; return 0; }
counts=""
for arm in OFF GLOBAL STAGED; do
    d="$OUT/$arm/$SAMPLE"; line=""
    for pair in "L1:LINE1" "Alu:ALU" "SVA:SVA"; do
        f=$(first_match "$d/merged/*_${pair##*:}.vcf")
        [ -n "$f" ] && line="$line ${pair%%:*}=$(n_records "$f")"
    done
    mf=$(first_match "$d/*.vcf"); [ -n "$mf" ] && line="$line merged=$(n_records "$mf")"
    echo "  $arm:$line"
    counts="$counts|$line"
done
n_uniq=$(printf '%s' "$counts" | tr '|' '\n' | grep -v '^$' | sort -u | wc -l)
[ "$n_uniq" -eq 1 ] || { echo "  FATAL: the arms produced different call counts"; RC=1; }

echo
if [ "$RC" -eq 0 ]; then echo "STAGED BARRIER GATE (WGS): PASS"; else echo "STAGED BARRIER GATE (WGS): FAIL"; fi
echo "finished $(date '+%F %T')"
exit "$RC"
