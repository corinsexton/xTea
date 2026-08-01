#!/bin/bash
#SBATCH -c 24
#SBATCH -t 3:00:00
#SBATCH --mem=64G
#SBATCH -p priopark
#SBATCH --account=park_contrib
#SBATCH -o /n/scratch/users/c/cos689/xtea2_testdata/transduction_phase_timing_%j.out
#SBATCH --mail-type=NONE

# DOWNSTREAM_OVERHAUL_PLAN.md Stage 0, measure item 2: split transduction's ~1481s (L1) /
# ~1772s (SVA) between read gather, the three bwa passes, and SAM parsing.  This is what
# unblocks the Stage 1 scope decision -- 2x off 26% of the step is a different call from
# 2x off 80% -- so the plan forbids pre-committing the C++ gather boundary until these
# numbers exist.
#
# Differences from its parent, docs/handoff/runs/transduction_attribution.sbatch:
#   * runs the LIVE tree, not the frozen $SCRATCH/port_wt worktree.  The phase timers
#     ([TD-PHASE] lines, x_transduction.py:td_phase) only exist in the live tree, and the
#     live tree is the code the plan rewrites.
#   * drops the master arm entirely.  Master carries no timers, and the plan retires
#     master agreement as a gate (§2).  Halving the work buys the second repeat instead.
#   * two repeats, because the gate is "reproducible within ~5%", which one run cannot show.
#
# Shape matches production: L1 and SVA concurrent at -n 8 each.  Same staged C++-run
# cns-filter and raw-disc inputs both repeats, read-only; every output tree is job-local
# and deliberately retained as diagnostic.
set -euo pipefail

SCRATCH=/n/scratch/users/c/cos689/xtea2_testdata
INPUT_ROOT=$SCRATCH/port_wgs_out/NA12878.wgs
REPO=/n/data1/hms/dbmi/park/corinne/claude_sandbox/xTea
# Run against a SNAPSHOT of the tree, not the live checkout.  Repeat 2 starts ~30 min after
# repeat 1 and re-imports the package, so a live path would silently let concurrent editing
# of the working tree change what the second repeat measures -- and the reproducibility gate
# would then be comparing two different programs.  Pass TD_PHASE_XTEA=<snapshot>/xtea.
LIVE_XTEA=${TD_PHASE_XTEA:-$REPO/xtea}
OUT=${TD_PHASE_OUT:-$SCRATCH/transduction_phase_timing_${SLURM_JOB_ID}}
REF=/n/data1/hms/dbmi/park-smaht_dac/ref/GRCh38_no_alt/hg38_no_alt.fa
REP=/home/cos689/xTea/rep_lib_annotation

# xtea2 is the repository's supported pipeline/test environment.  Activation supplies its
# registered PYTHONNOUSERSITE=1 guard against ~/.local's incompatible onnxruntime; do not
# replace this with a bare interpreter-path invocation.
source /n/app/conda/miniforge3/24.11.3-0/etc/profile.d/conda.sh
conda activate xtea2
PY=$(command -v python)
test "$PYTHONNOUSERSITE" = 1
export PATH=/home/cos689/software/bin:$PATH

for required in "$INPUT_ROOT/bam_list.txt" "$LIVE_XTEA/x_TEA_main.py" "$REF"; do
    [ -e "$required" ] || { echo "FATAL: missing $required"; exit 1; }
done
# The whole point is measuring instrumented code; fail loud rather than produce a
# timing run with no timings in it.
grep -q "TD-PHASE" "$LIVE_XTEA/x_transduction.py" || {
    echo "FATAL: $LIVE_XTEA/x_transduction.py carries no phase timers"; exit 1; }

mkdir -p "$OUT"
printf 'repeat\ttype\tthreads\telapsed_s\tuser_s\tsystem_s\tmax_rss_kb\texit\n' > "$OUT/results.tsv"
printf 'repeat\trep_type\tpid\tphase\twall_s\n' > "$OUT/phases.tsv"

family_args() {
    case "$1" in
        L1)
            printf '%s\n' 1 \
                "$REP/LINE/hg38/hg38_FL_L1_flanks_3k.fa" \
                "$REP/consensus/LINE1.fa" \
                "$REP/LINE/hg38/hg38_L1.fa.out"
            ;;
        SVA)
            printf '%s\n' 4 \
                "$REP/SVA/hg38/hg38_FL_SVA_flanks_3k.fa" \
                "$REP/consensus/SVA.fa" \
                "$REP/SVA/hg38/hg38_SVA.out"
            ;;
        *) echo "unknown family: $1" >&2; return 2 ;;
    esac
}

run_one() {
    local rep=$1 family=$2 threads=$3
    local rtype flank cns annot work result raw_disc timefile rc
    mapfile -t args < <(family_args "$family")
    rtype=${args[0]}; flank=${args[1]}; cns=${args[2]}; annot=${args[3]}
    work="$OUT/rep${rep}/$family/n${threads}"
    result="$work/candidate_disc_filtered_cns2.txt"
    raw_disc="$INPUT_ROOT/$family/candidate_list_from_disc.txt.clip_sites_raw_disc.txt"
    timefile="$work/time.txt"
    mkdir -p "$work/tmp"

    echo "=== repeat $rep $family: -n $threads ==="
    set +e
    /usr/bin/time -v -o "$timefile" "$PY" "$LIVE_XTEA/x_TEA_main.py" \
        --transduction --cr 4 --nd 5 -b "$INPUT_ROOT/bam_list.txt" -p "$work/tmp/" \
        --fflank "$flank" --flklen 3000 -n "$threads" \
        -i "$INPUT_ROOT/$family/candidate_disc_filtered_cns.txt" -r "$cns" \
        --ref "$REF" --input2 "$raw_disc" --rtype "$rtype" -a "$annot" -o "$result" \
        > "$work/stdout.txt" 2> "$work/stderr.txt"
    rc=$?
    set -e
    # NOT -F: here. /usr/bin/time -v prints the wall clock as "h:mm:ss or m:ss", so the
    # VALUE contains colons and a colon-split field grab returns "mm" -- which is what the
    # parent script (transduction_attribution.sbatch) does, and why its elapsed column is
    # useless. Take everything after the last ": " on the matched line instead.
    awk -v rep="$rep" -v family="$family" -v threads="$threads" -v rc="$rc" '
        function val(line) { sub(/^.*: /, "", line); gsub(/^[ \t]+|[ \t]+$/, "", line); return line }
        /Elapsed \(wall clock\) time/      {elapsed=val($0)}
        /User time \(seconds\)/            {user=val($0)}
        /System time \(seconds\)/          {sys_time=val($0)}
        /Maximum resident set size/        {rss=val($0)}
        END {printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n", rep, family, threads, elapsed, user, sys_time, rss, rc}
    ' "$timefile" >> "$OUT/results.tsv"
    # The timers write to stderr, mixed into a very noisy step log; pull them out here.
    awk -v rep="$rep" -F'\t' '$1=="[TD-PHASE]" {print rep"\t"$3"\t"$2"\t"$4"\t"$5}' \
        "$work/stderr.txt" >> "$OUT/phases.tsv"
    [ "$rc" -eq 0 ] || { echo "FAIL: repeat $rep $family (exit $rc)" >&2; return "$rc"; }
}

{
    echo "code run: $LIVE_XTEA"
    echo "repo commit at job start: $(git -C "$REPO" rev-parse --short HEAD)"
    echo "repo dirty at job start:"
    git -C "$REPO" status --short
    echo "input source: $INPUT_ROOT"
    [ -f "$LIVE_XTEA/../SNAPSHOT_PROVENANCE.txt" ] && cat "$LIVE_XTEA/../SNAPSHOT_PROVENANCE.txt"
} | tee "$OUT/provenance.txt"

for rep in 1 2; do
    run_one "$rep" L1 8 & p1=$!
    run_one "$rep" SVA 8 & p2=$!
    wait "$p1"
    wait "$p2"
done

echo
echo "--- step wall time ---"
column -t -s $'\t' "$OUT/results.tsv"

echo
echo "--- phase breakdown, and repeat-to-repeat spread (gate: within ~5%) ---"
awk -F'\t' 'NR>1 {t[$2"\t"$4"\t"$1]=$5}
    END {
        printf "%-12s %-38s %10s %10s %8s\n", "rep_type", "phase", "rep1_s", "rep2_s", "delta%"
        for (k in t) {
            split(k, a, "\t")
            if (a[3]==1) {
                k2=a[1]"\t"a[2]"\t2"
                d = (k2 in t && t[k]>0) ? (t[k2]-t[k])/t[k]*100 : "NA"
                printf "%-12s %-38s %10.2f %10s %8s\n", a[1], a[2], t[k],
                    (k2 in t ? sprintf("%.2f", t[k2]) : "NA"),
                    (d=="NA" ? "NA" : sprintf("%+.1f", d))
            }
        }
    }' "$OUT/phases.tsv" | sort -k1,1 -k2,2

echo
echo "--- share of TOTAL, repeat 1 (this is the Stage 1 scope decision) ---"
awk -F'\t' 'NR>1 && $1==1 {
        if ($4=="TOTAL_call_candidate_transduction_v3") tot[$2]=$5; else t[$2"\t"$4]=$5
    }
    END {
        for (k in t) { split(k, a, "\t")
            printf "%-6s %-38s %9.1fs  %5.1f%%\n", a[1], a[2], t[k],
                (tot[a[1]]>0 ? t[k]/tot[a[1]]*100 : 0) }
        for (r in tot) printf "%-6s %-38s %9.1fs  %5.1f%%\n", r, "TOTAL", tot[r], 100
    }' "$OUT/phases.tsv" | sort -k1,1 -k2,2

echo
echo "--- output equivalence across repeats (determinism check) ---"
for family in L1 SVA; do
    a="$OUT/rep1/$family/n8/candidate_disc_filtered_cns2.txt"
    b="$OUT/rep2/$family/n8/candidate_disc_filtered_cns2.txt"
    printf '%s\t' "$family"
    if diff -q "$a" "$b" >/dev/null 2>&1; then echo "IDENTICAL"; else echo "DIFFER"; fi
done
