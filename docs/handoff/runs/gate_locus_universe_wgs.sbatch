#!/bin/bash
#SBATCH -c 24
#SBATCH -t 10:00:00
#SBATCH --mem=64G
#SBATCH -p priopark
#SBATCH --account=park_contrib
#SBATCH -o /n/scratch/users/c/cos689/xtea2_testdata/gate_locus_universe_wgs_%j.out
#SBATCH --mail-type=NONE
# =============================================================================
# LOCUS UNIVERSE + DRIVER GATHER at WGS scale.
#
#     sbatch docs/handoff/runs/gate_locus_universe_wgs.sbatch <worktree> [types]
#
# WHY WGS AND NOT ONLY chr1, which already passed (job 49547574). Both flags open the
# `cns_filter -> transduction` wave split, and **the last bug in that split was invisible
# on chr1**: job 49498881 died at `SVA:transduction` with a bare `KeyError` because a step
# opened an input it was never handed (`x_TEA_main.py:551` derives `.high_confident` from
# its own `-i`), and chr1's smaller candidate sets never reached it. Anything that changes
# the wave topology earns a WGS leg for that reason alone.
#
# ONE ARM COVERS BOTH FLAGS: --td-gather-driver IMPLIES --locus-universe (it consumes the
# site lists that pass writes), so the GATHER arm exercises the strictly larger change.
# A separate ON arm would add ~90 minutes to re-test a subset.
#
# WHAT MUST BE TRUE. Both flags are answer-preserving by construction, so this is a
# byte-diff gate, not a prediction gate:
#
#   0. the arms ran what they claim -- OFF has no barrier and no driver gather, GATHER has
#      exactly one of each. Nothing below can tell a real result from two OFF arms.
#   1. arm OFF == frozen BASELINE-0 (4/4 VCFs). Isolates commit drift: without it, a clean
#      GATHER-vs-OFF diff could be two equally-wrong arms.
#   2. GATHER == OFF on every VCF / cns2 / sibling list.
#   3. per flanked family, the driver's replica == the step's own `.slct`. Alu has no step
#      `.slct` (no transduction flank), so the expected count is the FLANKED count -- and it
#      is printed, because a gate whose only check runs at comparison time cannot tell "the
#      code is wrong" from "the harness compared nothing" (job 49516776).
#   4. the GATHER arm's steps SKIPPED their own gather. Both gathering is silently correct
#      and buys nothing; neither gathering is silently empty. Read off the step logs.
#   5. the universe numbers reproduce the offline gate: 41,604 loci, 38,815 created-only.
#
# Frozen worktree, not the live checkout -- run_xtea spawns a fresh x_TEA_main per step, so
# a job sharing the checkout picks up mid-run edits.
# =============================================================================
set -euo pipefail

SCRATCH=/n/scratch/users/c/cos689/xtea2_testdata
WT="${1:?usage: sbatch gate_locus_universe_wgs.sbatch <worktree> [types]}"
TYPES="${2:-L1,Alu,SVA}"
OUT="${XTEA2_OUT:-$SCRATCH/gate_locus_universe_wgs_${SLURM_JOB_ID:-manual}}"
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
# The comparator and the flags have to exist BEFORE ~3 hours of compute, not after.
grep -q "td-gather-driver" "$WT/xtea/run_xtea.py" || {
    echo "FATAL: worktree has no --td-gather-driver; wrong commit"; exit 1; }
[ -d "$BASE0/BASELINE-0" ]      || { echo "FATAL: $BASE0/BASELINE-0 missing"; exit 1; }

mkdir -p "$OUT"
cd "$WT"

echo "=============================================================="
echo "LOCUS UNIVERSE + DRIVER GATHER GATE  --  HG002 WGS, types=$TYPES"
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
    for arm in OFF GATHER; do
        [ -f "$OUT/$arm/$SAMPLE/$SAMPLE.vcf" ] || {
            echo "FATAL: $OUT/$arm has no completed run"; exit 1; }
    done
else
    run_arm OFF                       || exit 1
    run_arm GATHER --td-gather-driver || exit 1
fi

echo
echo "##############################################################"
echo "# LOCUS UNIVERSE + DRIVER GATHER GATE RESULT"
echo "##############################################################"
RC=0
LOG_OFF="$OUT/OFF/driver.log"
LOG_GTH="$OUT/GATHER/driver.log"

# --- 0. the arms ran what they claim -------------------------------------------------
echo
echo "--- 0. arms ran what they claim ---"
n_u_off=$(grep -c 'locus universe (cns_filter -> transduction' "$LOG_OFF" || true)
n_u_gth=$(grep -c 'locus universe (cns_filter -> transduction' "$LOG_GTH" || true)
n_g_off=$(grep -c 'td gather (driver)' "$LOG_OFF" || true)
n_g_gth=$(grep -c 'td gather (driver)' "$LOG_GTH" || true)
echo "  universe barrier: OFF=$n_u_off  GATHER=$n_u_gth   (want 0 / 1)"
echo "  driver gather   : OFF=$n_g_off  GATHER=$n_g_gth   (want 0 / 1)"
[ "$n_u_off" -eq 0 ] && [ "$n_g_off" -eq 0 ] || { echo "  FATAL: arm OFF is not OFF"; RC=1; }
# The universe must be present in GATHER: --td-gather-driver implies --locus-universe, and
# a GATHER arm without it would mean the implication silently stopped working.
[ "$n_u_gth" -eq 1 ] && [ "$n_g_gth" -eq 1 ] || {
    echo "  FATAL: arm GATHER did not take both the universe pass and the driver gather"; RC=1; }

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
    echo "  FAIL: arm OFF is not BASELINE-0; the GATHER delta is uninterpretable"; RC=1; }

# --- 2. GATHER vs OFF: no call may move ----------------------------------------------
echo
echo "--- 2. GATHER vs OFF: every VCF / cns2 / sibling list ---"
n_same=0; n_diff=0; n_new=0
while IFS= read -r a; do
    rel=${a#"$OUT/OFF/"}
    b="$OUT/GATHER/$rel"
    if [ ! -f "$b" ]; then echo "  MISSING in GATHER: $rel"; n_diff=$((n_diff+1)); continue; fi
    if cmp -s "$a" "$b"; then n_same=$((n_same+1)); else
        echo "  DIFFERS: $rel"; n_diff=$((n_diff+1))
        diff <(grep -v '^#' "$a" | cut -f1,2) <(grep -v '^#' "$b" | cut -f1,2) | head -20 || true
    fi
done < <(find "$OUT/OFF" \( -name '*.vcf' -o -name 'candidate_disc_filtered_cns2.txt' \
                            -o -name 'candidate_sibling_transduction2.txt' \) -type f | sort)
# A file GATHER produces that OFF does not is expected ONLY under merged/locus_universe/.
while IFS= read -r b; do
    rel=${b#"$OUT/GATHER/"}
    case "$rel" in */merged/locus_universe/*) continue ;; esac
    [ -f "$OUT/OFF/$rel" ] || { echo "  NEW file: $rel"; n_new=$((n_new+1)); }
done < <(find "$OUT/GATHER" \( -name '*.vcf' -o -name 'candidate_disc_filtered_cns2.txt' \
                               -o -name 'candidate_sibling_transduction2.txt' \) -type f | sort)
echo "  identical: $n_same    differing: $n_diff    unexpected new: $n_new"
[ "$n_diff" -eq 0 ] && [ "$n_new" -eq 0 ] && [ "$n_same" -gt 0 ] || RC=1

# --- 3. the driver's replica == the step's own site list -------------------------------
echo
echo "--- 3. driver replica vs step .slct (flanked families only) ---"
S="$OUT/GATHER/$SAMPLE"
FLANKED=""
for t in ${TYPES//,/ }; do
    case "$t" in L1|SVA) FLANKED="$FLANKED $t" ;; esac
done
n_exp=$(echo $FLANKED | wc -w)
n_cmp=0; n_ident=0
for t in $FLANKED; do
    step="$S/$t/candidate_list_from_disc.txt.clip_sites_raw_disc.txt.slct"
    drv="$S/merged/locus_universe/$t.slct"
    [ -f "$step" ] || { echo "  $t: FATAL: step .slct missing ($step)"; continue; }
    [ -f "$drv" ]  || { echo "  $t: FATAL: driver replica missing ($drv)"; continue; }
    n_cmp=$((n_cmp+1))
    if cmp -s "$step" "$drv"; then
        n_ident=$((n_ident+1))
        echo "  $t: IDENTICAL ($(wc -c < "$step") bytes, $(wc -l < "$step") rows)"
    else
        echo "  $t: DIFFERS"; diff "$step" "$drv" | head -20 || true
    fi
done
echo "  products compared: $n_cmp (expected $n_exp)   identical: $n_ident"
[ "$n_cmp" -eq "$n_exp" ] && [ "$n_ident" -eq "$n_exp" ] && [ "$n_exp" -gt 0 ] || RC=1
[ -f "$S/Alu/candidate_list_from_disc.txt.clip_sites_raw_disc.txt.slct" ] && {
    echo "  FATAL: Alu grew a step .slct; the flankless assumption is wrong"; RC=1; } || true

# --- 4. the GATHER arm's steps skipped their own gather --------------------------------
echo
echo "--- 4. the bypass actually fired ---"
n_ok=0; n_bad=0
for t in $FLANKED; do
    for arm in OFF GATHER; do
        slog="$OUT/$arm/$SAMPLE/${SAMPLE}_${t}.steps.log"
        [ -f "$slog" ] || { echo "  $arm/$t: FATAL: no steps log"; n_bad=$((n_bad+1)); continue; }
        n_pre=$(grep -c 'precomputed by the driver' "$slog" || true)
        n_own=$(grep -c -- '--td-sites' "$slog" || true)
        echo "  $arm/$t: precomputed=$n_pre  own-gather=$n_own"
        case "$arm" in
          OFF)    [ "$n_pre" -eq 0 ] && [ "$n_own" -eq 1 ] && n_ok=$((n_ok+1)) || {
                      echo "    FATAL: OFF arm must gather in the step"; n_bad=$((n_bad+1)); } ;;
          GATHER) [ "$n_pre" -eq 1 ] && [ "$n_own" -eq 0 ] && n_ok=$((n_ok+1)) || {
                      echo "    FATAL: GATHER arm must bypass the step's gather"; n_bad=$((n_bad+1)); } ;;
        esac
    done
done
echo "  bypass checks passed: $n_ok (expected $((n_exp * 2)))   failed: $n_bad"
[ "$n_bad" -eq 0 ] && [ "$n_ok" -eq $((n_exp * 2)) ] || RC=1

# --- 5. the universe reproduces the offline gate ---------------------------------------
echo
echo "--- 5. the universe, as the run reported it (offline gate: 41604 loci, 38815 created-only) ---"
grep -A 8 'locus universe (cns_filter' "$LOG_GTH" | sed 's/^/  /' || true
n_rows=$(grep -vc '^#' "$S/merged/locus_universe/locus_universe.tsv" || echo 0)
echo "  table rows: $n_rows"
[ "$n_rows" -gt 0 ] || { echo "  FATAL: empty universe table"; RC=1; }
echo
echo "--- the driver gather, as the GATHER arm reported it ---"
grep -A 6 'td gather (driver)' "$LOG_GTH" | sed 's/^/  /' || true

# --- timings, reported not gated ------------------------------------------------------
echo
echo "--- wall clock (the saving is CPU, not wall -- see run_td_gather_driver) ---"
for arm in OFF GATHER; do
    printf '  %-7s total %ss   ' "$arm" "$(cat "$OUT/$arm/wall_s.txt" 2>/dev/null || echo '?')"
    grep -aoE 'phase (downstream_parallel|locus_universe|td_gather_driver): [0-9.]+s' \
        "$OUT/$arm/driver.log" 2>/dev/null | tr '\n' ' ' || true
    echo
done
echo "  transduction gather CPU, per family:"
for arm in OFF GATHER; do
    for t in $FLANKED; do
        v=$(grep -ah '^\[TD-PHASE\]' "$OUT/$arm/$SAMPLE/${SAMPLE}_${t}.steps.log" 2>/dev/null \
            | awk -F'\t' '$4=="native_gather"{print $5}' | head -1)
        printf '    %-7s %-4s native_gather=%s\n' "$arm" "$t" "${v:-none}"
    done
done

echo
echo "--- call counts ---"
n_records() { [ -f "$1" ] || { echo 0; return 0; }; grep -vc '^#' "$1" || true; }
first_match() { local m; for m in $1; do [ -e "$m" ] && { printf '%s' "$m"; return 0; }; done; return 0; }
for arm in OFF GATHER; do
    d="$OUT/$arm/$SAMPLE"; line="  $arm:"
    for pair in "L1:LINE1" "Alu:ALU" "SVA:SVA"; do
        f=$(first_match "$d/merged/*_${pair##*:}.vcf")
        [ -n "$f" ] && line="$line ${pair%%:*}=$(n_records "$f")"
    done
    mf=$(first_match "$d/*.vcf"); [ -n "$mf" ] && line="$line merged=$(n_records "$mf")"
    echo "$line"
done

echo
if [ "$RC" -eq 0 ]; then echo "LOCUS UNIVERSE WGS GATE: PASS"; else echo "LOCUS UNIVERSE WGS GATE: FAIL"; fi
echo "finished $(date '+%F %T')"
exit "$RC"
