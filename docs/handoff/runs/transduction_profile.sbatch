#!/bin/bash
#SBATCH -c 16
#SBATCH -t 3:00:00
#SBATCH --mem=48G
#SBATCH -p priopark
#SBATCH --account=park_contrib
#SBATCH -o /n/scratch/users/c/cos689/xtea2_testdata/transduction_profile_%j.out
#SBATCH --mail-type=NONE

# Per-candidate phase profiling for parse_td_sibling_by_site (2026-07-29 handoff item 2,
# "transduction refinement"). Runs master and port against the SAME staged cns-filter inputs
# used by job 49075238, but against separate copies (master_wt_profile / port_wt_profile) with
# XTEA_TD_PROFILE=1 instrumentation, so as not to disturb 49075238's own master_cmp_wt/port_wt
# trees while it's still mid-flight. Each candidate's phase split (open/fetch/cluster) is
# written to stderr as a TD_PROFILE line; this script summarizes those after each run.
set -euo pipefail

SCRATCH=/n/scratch/users/c/cos689/xtea2_testdata
INPUT_ROOT=$SCRATCH/port_wgs_out/NA12878.wgs
MASTER_XTEA=$SCRATCH/master_wt_profile/xtea
PORT_XTEA=$SCRATCH/port_wt_profile/xtea
OUT=${TD_PROFILE_OUT:-$SCRATCH/transduction_profile_${SLURM_JOB_ID}}
PY_MASTER=/n/data1/hms/dbmi/park/corinne/.conda/envs/xTea/bin/python
REF=/n/data1/hms/dbmi/park-smaht_dac/ref/GRCh38_no_alt/hg38_no_alt.fa
REP=/home/cos689/xTea/rep_lib_annotation

source /n/app/conda/miniforge3/24.11.3-0/etc/profile.d/conda.sh
conda activate xtea2
PY_PORT=$(command -v python)
test "$PYTHONNOUSERSITE" = 1
export PATH=/home/cos689/software/bin:$PATH
export XTEA_TD_PROFILE=1

for required in "$INPUT_ROOT/bam_list.txt" "$MASTER_XTEA/x_TEA_main.py" "$PORT_XTEA/x_TEA_main.py" "$REF"; do
    [ -e "$required" ] || { echo "FATAL: missing $required"; exit 1; }
done

mkdir -p "$OUT"

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
    work="$OUT/$impl/$family"
    result="$work/candidate_disc_filtered_cns2.txt"
    raw_disc="$INPUT_ROOT/$family/candidate_list_from_disc.txt.clip_sites_raw_disc.txt"
    timefile="$work/time.txt"
    mkdir -p "$work/tmp"

    echo "=== $impl $family: profiling run, -n $threads ==="
    set +e
    /usr/bin/time -v -o "$timefile" "$py" "$code/x_TEA_main.py" \
        --transduction --cr 4 --nd 5 -b "$INPUT_ROOT/bam_list.txt" -p "$work/tmp/" \
        --fflank "$flank" --flklen 3000 -n "$threads" \
        -i "$INPUT_ROOT/$family/candidate_disc_filtered_cns.txt" -r "$cns" \
        --ref "$REF" --input2 "$raw_disc" --rtype "$rtype" -a "$annot" -o "$result" \
        > "$work/stdout.txt" 2> "$work/stderr.txt"
    rc=$?
    set -e
    grep -c '^TD_PROFILE' "$work/stderr.txt" > "$work/n_profiled.txt" || echo 0 > "$work/n_profiled.txt"
    [ "$rc" -eq 0 ] || { echo "FAIL: $impl $family (exit $rc)" >&2; return "$rc"; }
}

echo "master commit: $(git -C "$SCRATCH/master_wt_profile" rev-parse --short HEAD)" | tee "$OUT/provenance.txt"
echo "port commit: $(git -C "$SCRATCH/port_wt_profile" rev-parse --short HEAD)" | tee -a "$OUT/provenance.txt"
echo "input source: $INPUT_ROOT" | tee -a "$OUT/provenance.txt"

for impl_spec in "master:$PY_MASTER:$MASTER_XTEA" "port:$PY_PORT:$PORT_XTEA"; do
    IFS=: read -r impl py code <<< "$impl_spec"
    run_one "$impl" L1 "$py" "$code" 8 & p1=$!
    run_one "$impl" SVA "$py" "$code" 8 & p2=$!
    wait "$p1"
    wait "$p2"
done

echo
echo "--- per-phase summary (mean seconds/candidate, n candidates) ---"
printf 'impl\tfamily\tn\tmean_open_s\tmean_fetch_s\tmean_cluster_s\tmean_total_s\tmean_n_reads\n'
for impl in master port; do
    for family in L1 SVA; do
        f="$OUT/$impl/$family/stderr.txt"
        [ -f "$f" ] || continue
        awk -F'\t' -v impl="$impl" -v family="$family" '
            $1=="TD_PROFILE" {n++; nreads+=$4; open+=$5; fetch+=$6; cluster+=$7; total+=$8}
            END {
                if (n>0) printf "%s\t%s\t%d\t%.4f\t%.4f\t%.4f\t%.4f\t%.1f\n", \
                    impl, family, n, open/n, fetch/n, cluster/n, total/n, nreads/n
            }' "$f"
    done
done
