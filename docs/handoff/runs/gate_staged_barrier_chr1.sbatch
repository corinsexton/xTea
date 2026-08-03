#!/bin/bash
#SBATCH -J sb_chr1
#SBATCH -A park_contrib
#SBATCH -p priopark
#SBATCH -t 2:00:00
#SBATCH -c 12
#SBATCH --mem 64G
#SBATCH -o /n/scratch/users/c/cos689/xtea2_testdata/gate_staged_barrier_chr1_%j.out
# =============================================================================
# The STAGED cns_filter -> transduction barrier, end to end on HG002 chr1. THREE arms,
# ONE commit: OFF (no barrier at all), GLOBAL (--td-gather-driver --global-barrier, the
# form jobs 49548949 and 49550925 gated), STAGED (--td-gather-driver, the new default).
#
# WHAT CHANGED AND WHY IT NEEDS A GATE. Until now, anything that ran at this point made
# EVERY family wait for the SLOWEST family's cns_filter. At WGS that is Alu's, ~2x L1's
# and SVA's, and job 49550925 measured the join at +52.9s -- none of it work. The unsplit
# chain overlaps L1/SVA transduction with the tail of Alu's cns_filter, and a join throws
# that overlap away. `StagedBarrier` gives each family only the waits it actually has:
# enumeration is per-family (own cns list, own raw_disc, nobody else's anything), the
# driver gather needs the FLANKED families and is needed by only them, and the union table
# and --family-assign-early are the only things that need all three.
#
# So this is a CONCURRENCY change to code that must not move a call, which is the
# combination that produces silent damage. Three things can go wrong and only one of them
# shows up as a wrong number:
#
#   1. A family runs before something it reads. Fatal, and fatal LATE -- job 49498881 died
#      inside SVA:transduction at WGS, invisible on chr1. Here the risk is a flanked
#      family entering transduction before the driver gather wrote the FASTQ its bypass
#      tells it to skip to. Checks 1 and 4.
#   2. The staged enumeration produces different files from the all-at-once one -- it runs
#      the same function per family, but concurrently, against a working folder each step
#      also uses. Check 3 compares the two arms' replicas to each other AND to the steps'.
#   3. It changes an answer. Checks 1 and 6.
#
# AND THE POINT OF THE CHANGE, which needs its own evidence or the whole thing is a
# refactor with no payoff: check 5 reads the driver log's elapsed stamps and asserts that
# in the STAGED arm work genuinely OVERLAPS the barrier -- some family is past it while
# another is still short of it -- and that in the GLOBAL arm nothing does. That is the
# ordering claim stated as a fact about the run rather than as a wall-clock delta, which
# at chr1 scale would be inside the noise.
#
# USAGE  sbatch gate_staged_barrier_chr1.sbatch <worktree>
# =============================================================================
set -euo pipefail

SCRATCH=/n/scratch/users/c/cos689/xtea2_testdata
WT="${1:?usage: sbatch gate_staged_barrier_chr1.sbatch <worktree>}"
OUT="${XTEA2_OUT:-$SCRATCH/gate_staged_barrier_chr1_${SLURM_JOB_ID:-manual}}"
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

mkdir -p "$OUT"
cd "$WT"

echo "=============================================================="
echo "STAGED BARRIER GATE  --  HG002 chr1, types=$TYPES"
echo "  job    : ${SLURM_JOB_ID:-<none>} on $(hostname)  cores=$CORES"
echo "  commit : $(git rev-parse --short HEAD) $(git log -1 --format=%s | cut -c1-60)"
echo "  out    : $OUT"
echo "  started: $(date '+%F %T')"
echo "=============================================================="

run_arm() {   # $1 = arm name, $2... = extra driver flags
    # Separate statements, not `local arm=$1 armout=...`: under `set -u` the whole line
    # expands before `local` assigns and $arm is still unset (job 49452897 / 49516776).
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
}

RC=0
run_arm OFF || RC=1
run_arm GLOBAL --td-gather-driver --global-barrier || RC=1
run_arm STAGED --td-gather-driver || RC=1
[ "$RC" -eq 0 ] || { echo; echo "STAGED BARRIER GATE: FAIL (an arm did not complete)"; exit 1; }

FLANKED=""
for t in ${TYPES//,/ }; do
    case "$t" in L1|SVA) FLANKED="$FLANKED $t" ;; esac
done
n_exp=$(echo $FLANKED | wc -w)

echo
echo "--- check 0: the arms ran the topologies they claim ---"
# Both barrier arms must do the SAME WORK -- one universe pass, one driver gather -- and
# differ only in how they wait. An arm that silently lost the gather would sail through
# every byte comparison below, because not gathering is answer-preserving too.
for arm in OFF GLOBAL STAGED; do
    n_uni=$(grep -ac 'locus universe (cns_filter' "$OUT/$arm/driver.log" || true)
    n_gth=$(grep -ac 'td gather (driver):' "$OUT/$arm/driver.log" || true)
    n_stg=$(grep -ac 'barrier: STAGED' "$OUT/$arm/driver.log" || true)
    n_glb=$(grep -ac 'barrier: GLOBAL' "$OUT/$arm/driver.log" || true)
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

echo
echo "--- check 1: every arm vs OFF, no call may move ---"
compare_to_off() {
    # Separate `local` statements: bash expands the whole line before assigning (job 49516776).
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

echo
echo "--- check 2: STAGED vs GLOBAL, the universe artifacts themselves ---"
# Check 1 compares CALLS, and the universe pass is diagnostic -- it cannot move one. So a
# staged enumeration that produced a different (or empty) table would pass check 1 in
# silence. These are the files the change actually rewrites the production of.
n_u=0; n_u_ok=0
for f in "locus_universe.tsv" $(for t in $FLANKED; do echo "$t.slct"; done) ; do
    g="$OUT/GLOBAL/$SAMPLE/merged/locus_universe/$f"
    s="$OUT/STAGED/$SAMPLE/merged/locus_universe/$f"
    if [ ! -f "$g" ] || [ ! -f "$s" ]; then
        echo "  $f: FATAL: missing (GLOBAL=$([ -f "$g" ] && echo y || echo n) STAGED=$([ -f "$s" ] && echo y || echo n))"
        n_u=$((n_u+1)); continue
    fi
    n_u=$((n_u+1))
    if cmp -s "$g" "$s"; then
        n_u_ok=$((n_u_ok+1))
        echo "  $f: IDENTICAL ($(wc -c < "$s") bytes, $(wc -l < "$s") rows)"
    else
        echo "  $f: DIFFERS"; diff "$g" "$s" | head -10 || true
    fi
done
echo "  universe artifacts identical: $n_u_ok / $n_u"
[ "$n_u" -gt 0 ] && [ "$n_u_ok" -eq "$n_u" ] || RC=1

echo
echo "--- check 3: in BOTH arms, the driver's replica == the step's own site list ---"
# Run per arm, not just once: the staged arm enumerates concurrently with other families'
# steps, and XTransduction builds its annotation index inside a working folder the step
# also uses. If that collided, this is where it shows.
n_cmp=0; n_ident=0
for arm in GLOBAL STAGED; do
    S="$OUT/$arm/$SAMPLE"
    for t in $FLANKED; do
        step="$S/$t/candidate_list_from_disc.txt.clip_sites_raw_disc.txt.slct"
        drv="$S/merged/locus_universe/$t.slct"
        if [ ! -f "$step" ]; then echo "  $arm/$t: FATAL: step .slct missing"; continue; fi
        if [ ! -f "$drv" ]; then echo "  $arm/$t: FATAL: driver replica missing"; continue; fi
        n_cmp=$((n_cmp+1))
        if cmp -s "$step" "$drv"; then
            n_ident=$((n_ident+1))
            echo "  $arm/$t: IDENTICAL ($(wc -c < "$step") bytes, $(wc -l < "$step") rows)"
        else
            echo "  $arm/$t: DIFFERS"; diff "$step" "$drv" | head -20 || true
        fi
    done
    # Alu must NOT have a step .slct -- if it grows one, the flankless assumption is wrong.
    if [ -f "$S/Alu/candidate_list_from_disc.txt.clip_sites_raw_disc.txt.slct" ]; then
        echo "  FATAL: Alu grew a step .slct in $arm; the flankless assumption is wrong"; RC=1
    fi
done
echo "  products compared: $n_cmp (expected $((n_exp * 2)))   identical: $n_ident"
[ "$n_cmp" -eq $((n_exp * 2)) ] && [ "$n_ident" -eq "$n_cmp" ] && [ "$n_exp" -gt 0 ] || RC=1

echo
echo "--- check 4: in BOTH barrier arms, the steps skipped their own gather ---"
# The bypass is where failure mode 1 lands. Both gathering is silently correct and buys
# nothing; neither gathering is silently empty. Per flanked family: one 'precomputed by
# the driver' (the bypass, on stderr) and zero '--td-sites' (the step's own gather). The
# OFF arm is the control and must be the exact reverse.
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

echo
echo "--- check 5: the flankless family stopped waiting for a gather it never reads ---"
# The reason the change exists, asserted rather than inferred from a wall-clock delta that
# chr1 is too small to show. Every driver.log line now carries seconds since driver start.
#
# WHAT IS MEASURED, and why it is this and not "did anything overlap". The tempting check
# -- did some family cross the barrier before the slowest cns_filter landed -- is real but
# ORDER-DEPENDENT: whether it holds depends on which family happens to be slowest and on
# whether the gather outruns the gap behind it, neither of which is a property of the
# code. A gate that passes for a reason the change does not control is worse than no gate.
#
# This one is structural. Alu has no transduction flank, so it never reads the gather, so
# under the staged barrier it waits for NOTHING: the time from its own enumeration line to
# its wave-2 announcement is just the union-table write. Under the global barrier the same
# gap must contain the entire driver gather, whichever family was slowest. So:
#
#     gap(STAGED)  <  gap(GLOBAL),  and  gap(STAGED)  <  the gather's own duration
#
# holds on any ordering, and fails the moment the flankless family is put back in the
# queue behind work it does not consume.
stamp_of() {  # $1 = log, $2 = grep -E pattern -> the elapsed seconds on the FIRST match
    grep -aoE "^\[xtea2 \+ *[0-9.]+s\].*$2" "$1" 2>/dev/null | head -1 \
        | sed -E 's/^\[xtea2 \+ *([0-9.]+)s\].*/\1/' || true
}
case ",$TYPES," in
  *,Alu,*)
    gap_of() {  # $1 = arm -> seconds from Alu's enumeration line to Alu's transduction
        local L t_enum t_td
        L="$OUT/$1/driver.log"
        t_enum=$(stamp_of "$L" 'Alu +[0-9]+ cns candidates')
        t_td=$(stamp_of "$L" '=== type Alu .*from transduction')
        [ -n "$t_enum" ] && [ -n "$t_td" ] || { printf '?'; return 0; }
        awk -v a="$t_enum" -v b="$t_td" 'BEGIN{printf "%.1f", b - a}'
    }
    g_global=$(gap_of GLOBAL); g_staged=$(gap_of STAGED)
    # The gather's own duration, from its timing row -- "phase" in GLOBAL, "step" in STAGED
    # (it overlaps other families' steps there, so recording it as a phase would break the
    # timing file's rule that phase rows are disjoint).
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
            echo "    FATAL: the GLOBAL arm's gap does not contain the gather -- the"
            echo "           control is not measuring what this check assumes"; RC=1; }
    fi
    ;;
  *)
    echo "  SKIPPED: no flankless family in types=$TYPES; nothing waits on a gather it"
    echo "           does not read, so there is no gap to measure."
    ;;
esac

echo
echo "--- reported, not gated: did anything cross the barrier early ---"
# Order-dependent (see check 5), so observed rather than asserted -- but it is the shape
# of the saving, so it is worth having in the record for the WGS arm.
last_stamp_of() {
    grep -aoE "^\[xtea2 \+ *[0-9.]+s\].*$2" "$1" 2>/dev/null | tail -1 \
        | sed -E 's/^\[xtea2 \+ *([0-9.]+)s\].*/\1/' || true
}
for arm in GLOBAL STAGED; do
    L="$OUT/$arm/driver.log"
    echo "  $arm: last enumeration=+$(last_stamp_of "$L" 'cns candidates \+')s" \
         " first transduction=+$(stamp_of "$L" 'from transduction')s" \
         " gather=+$(stamp_of "$L" 'td gather \(driver\):')s"
done

echo
echo "--- check 6: call counts, all three arms ---"
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
# Three identical strings, or the arms disagree on the call set.
n_uniq=$(printf '%s' "$counts" | tr '|' '\n' | grep -v '^$' | sort -u | wc -l)
[ "$n_uniq" -eq 1 ] || { echo "  FATAL: the arms produced different call counts"; RC=1; }

echo
echo "--- reported, not gated: where the time went ---"
for arm in OFF GLOBAL STAGED; do
    printf '  %-7s' "$arm"
    grep -aoE '(phase|step) (downstream_parallel|locus_universe|td_gather_driver|[A-Za-z0-9]+:locus_universe): [0-9.]+s' \
        "$OUT/$arm/driver.log" 2>/dev/null | tr '\n' ' ' || true
    echo
done
echo "  per-family cns_filter (the thing the join waited on):"
for arm in OFF GLOBAL STAGED; do
    for t in ${TYPES//,/ }; do
        v=$(awk -F'\t' -v k="$t:cns_filter" '$2==k{print $3}' \
            "$OUT/$arm/$SAMPLE/step_timings.tsv" 2>/dev/null | head -1)
        printf '    %-7s %-4s %ss\n' "$arm" "$t" "${v:-?}"
    done
done

echo
if [ "$RC" -eq 0 ]; then echo "STAGED BARRIER GATE: PASS"; else echo "STAGED BARRIER GATE: FAIL"; fi
echo "finished $(date '+%F %T')"
exit "$RC"
