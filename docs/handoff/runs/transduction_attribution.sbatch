#!/bin/bash
#SBATCH -c 24
#SBATCH -t 6:00:00
#SBATCH --mem=64G
#SBATCH -p priopark
#SBATCH --account=park_contrib
#SBATCH -o /n/scratch/users/c/cos689/xtea2_testdata/transduction_attribution_%j.out
#SBATCH --mail-type=NONE

# Controlled answer to the WGS tail question.  The master and xTea2 invocations
# below receive the *same* C++-run cns-filter input and raw-disc input, use the
# same annotations/BAM/reference.  It measures both the production tail shape
# (-n 8 per family, L1/SVA concurrent) and -n 24 serially.  L1 and SVA are the
# only families measured because they dominate the observed downstream critical
# path.
#
# Each implementation gets its own output/temp directory; source inputs are
# read-only.  Do not use the normal run_xtea cleanup policy here: every artifact
# is diagnostic and deliberately retained under this job-specific output tree.
set -euo pipefail

SCRATCH=/n/scratch/users/c/cos689/xtea2_testdata
INPUT_ROOT=$SCRATCH/port_wgs_out/NA12878.wgs
MASTER_XTEA=$SCRATCH/master_cmp_wt/xtea
PORT_XTEA=$SCRATCH/port_wt/xtea
OUT=${TD_ATTR_OUT:-$SCRATCH/transduction_attribution_${SLURM_JOB_ID}}
PY_MASTER=/n/data1/hms/dbmi/park/corinne/.conda/envs/xTea/bin/python
REF=/n/data1/hms/dbmi/park-smaht_dac/ref/GRCh38_no_alt/hg38_no_alt.fa
REP=/home/cos689/xTea/rep_lib_annotation

# xtea2 is the repository's supported pipeline/test environment.  Activation
# supplies its registered PYTHONNOUSERSITE=1 guard against ~/.local's incompatible
# onnxruntime; do not replace this with a bare interpreter-path invocation.
source /n/app/conda/miniforge3/24.11.3-0/etc/profile.d/conda.sh
conda activate xtea2
PY_PORT=$(command -v python)
test "$PYTHONNOUSERSITE" = 1
export PATH=/home/cos689/software/bin:$PATH

for required in "$INPUT_ROOT/bam_list.txt" "$MASTER_XTEA/x_TEA_main.py" \
                "$PORT_XTEA/x_TEA_main.py" "$REF"; do
    [ -e "$required" ] || { echo "FATAL: missing $required"; exit 1; }
done

mkdir -p "$OUT"
printf 'implementation\ttype\tthreads\telapsed_s\tuser_s\tsystem_s\tmax_rss_kb\texit\n' > "$OUT/results.tsv"

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
    local impl=$1 family=$2 py=$3 code=$4 threads=$5
    local rtype flank cns annot work result raw_disc timefile rc
    mapfile -t args < <(family_args "$family")
    rtype=${args[0]}; flank=${args[1]}; cns=${args[2]}; annot=${args[3]}
    work="$OUT/$impl/$family/n${threads}"
    result="$work/candidate_disc_filtered_cns2.txt"
    raw_disc="$INPUT_ROOT/$family/candidate_list_from_disc.txt.clip_sites_raw_disc.txt"
    timefile="$work/time.txt"
    mkdir -p "$work/tmp"

    echo "=== $impl $family: same staged inputs, -n $threads ==="
    set +e
    /usr/bin/time -v -o "$timefile" "$py" "$code/x_TEA_main.py" \
        --transduction --cr 4 --nd 5 -b "$INPUT_ROOT/bam_list.txt" -p "$work/tmp/" \
        --fflank "$flank" --flklen 3000 -n "$threads" \
        -i "$INPUT_ROOT/$family/candidate_disc_filtered_cns.txt" -r "$cns" \
        --ref "$REF" --input2 "$raw_disc" --rtype "$rtype" -a "$annot" -o "$result" \
        > "$work/stdout.txt" 2> "$work/stderr.txt"
    rc=$?
    set -e
    awk -F: -v impl="$impl" -v family="$family" -v threads="$threads" -v rc="$rc" '
        /Elapsed \(wall clock\) time/ {gsub(/^[ \t]+/, "", $2); elapsed=$2}
        /User time \(seconds\)/ {gsub(/^[ \t]+/, "", $2); user=$2}
        /System time \(seconds\)/ {gsub(/^[ \t]+/, "", $2); sys_time=$2}
        /Maximum resident set size/ {gsub(/^[ \t]+/, "", $2); rss=$2}
        END {printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n", impl, family, threads, elapsed, user, sys_time, rss, rc}
    ' "$timefile" >> "$OUT/results.tsv"
    [ "$rc" -eq 0 ] || { echo "FAIL: $impl $family (exit $rc)" >&2; return "$rc"; }
}

echo "master commit: $(git -C "$SCRATCH/master_cmp_wt" rev-parse --short HEAD)" | tee "$OUT/provenance.txt"
echo "port commit: $(git -C "$SCRATCH/port_wt" rev-parse --short HEAD)" | tee -a "$OUT/provenance.txt"
echo "input source: $INPUT_ROOT" | tee -a "$OUT/provenance.txt"

# Run the two 8-core family jobs concurrently per implementation.  This matches
# the meaningful part of xTea2's parallel tail without oversubscribing the node.
for impl_spec in "master:$PY_MASTER:$MASTER_XTEA" "port:$PY_PORT:$PORT_XTEA"; do
    IFS=: read -r impl py code <<< "$impl_spec"
    run_one "$impl" L1 "$py" "$code" 8 & p1=$!
    run_one "$impl" SVA "$py" "$code" 8 & p2=$!
    wait "$p1"
    wait "$p2"
done

# Now give one family the whole 24-core allocation at a time.  This isolates
# each implementation's within-step scaling without concurrent-family effects.
for impl_spec in "master:$PY_MASTER:$MASTER_XTEA" "port:$PY_PORT:$PORT_XTEA"; do
    IFS=: read -r impl py code <<< "$impl_spec"
    run_one "$impl" L1 "$py" "$code" 24
    run_one "$impl" SVA "$py" "$code" 24
done

echo
column -t -s $'\t' "$OUT/results.tsv"
echo
echo "--- output equivalence (master vs port, same staged inputs) ---"
for threads in 8 24; do
    for family in L1 SVA; do
        m="$OUT/master/$family/n${threads}/candidate_disc_filtered_cns2.txt"
        p="$OUT/port/$family/n${threads}/candidate_disc_filtered_cns2.txt"
        printf '%s n=%s\t' "$family" "$threads"
        sha256sum "$m" "$p"
        diff -q "$m" "$p" || true
    done
done
