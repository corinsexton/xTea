#!/bin/bash
#SBATCH -c 8
#SBATCH -t 2:00:00
#SBATCH --mem=32G
#SBATCH -p priopark
#SBATCH --account=park_contrib
#SBATCH -o /n/scratch/users/c/cos689/xtea2_testdata/derive_truth_td_%j.out
#SBATCH --mail-type=NONE

# =============================================================================
# The truth side's TRANSDUCTION annotation, derived -- because the HG002 truth VCFs
# do not carry one (verified: no transduct*/TD=/SOURCE=/ORPHAN/SIBLING token in any
# of the five files). DOWNSTREAM_OVERHAUL_PLAN.md section 7.1.
#
#     sbatch docs/handoff/runs/derive_truth_td.sbatch [outdir]
#
# Depends only on the truth set + the reference, never on any call set, so its output
# is a CACHED ARTIFACT: regenerate it when the reference or the truth set changes, not
# once per scored run. It is committed to the plan's record by sha, not re-derived.
#
# Runs bwa twice per family: every truth ALT against the family consensus (to find the
# 3' soft-clip residual), then the long residuals against the whole genome (to name a
# source locus). Both are a few hundred sequences; the wall time is dominated by loading
# the hg38 bwa index, which is also why this is 32G and not 8G.
#
# It is submitted rather than run interactively for the usual reason: the session has
# 2 cores, and nothing that loads a 5.5 GB index belongs there.
#
# ONE CORRECTION TO THE PLAN, MEASURED HERE (job 49481569). Plan section 7.1 reports the
# L1 ALT-vs-consensus residual as "median 18 bp, >=100 bp: 60, >=200 bp: 51, max 4620".
# Those are TOTAL soft clip, LEADING PLUS TRAILING -- reproduced exactly (median 17,
# 60, 51, max 4619). The plan's own criterion asks for the residual at the 3' END only,
# which is a different and much smaller number: trailing-only is median 7, >=100 bp: 22,
# >=200 bp: 15. The leading clip is the larger of the two (40 records >=100 bp, max
# 1745) and is 5' inversion or genomic flank in the ALT -- the very confounders the plan
# names as its reason for requiring 3'. So the transduction denominator here is ~20 per
# family, not ~60, and nothing sized off the 60 is right.
# =============================================================================
set -euo pipefail

SCRATCH=/n/scratch/users/c/cos689/xtea2_testdata
REPO=/n/data1/hms/dbmi/park/corinne/claude_sandbox/xTea
OUT="${1:-$SCRATCH/truth_td_$(git -C "$REPO" rev-parse --short HEAD)}"
TRUTH=/n/data1/hms/dbmi/park/corinne/claude_sandbox/test_data/HG002/truth
REF=/n/data1/hms/dbmi/park-smaht_dac/ref/GRCh38_no_alt/hg38_no_alt.fa
REPLIB=/home/cos689/xTea/rep_lib_annotation
BWA=/home/cos689/software/bin/bwa
CORES="${SLURM_CPUS_PER_TASK:-8}"

source /n/app/conda/miniforge3/24.11.3-0/etc/profile.d/conda.sh
conda activate xtea2
export PATH=/home/cos689/software/bin:$PATH

[ -f "$REF.bwt" ] || { echo "FATAL: $REF is not bwa-indexed"; exit 1; }
[ -x "$BWA" ]     || { echo "FATAL: $BWA missing"; exit 1; }
[ -d "$TRUTH" ]   || { echo "FATAL: $TRUTH missing"; exit 1; }

mkdir -p "$OUT"
echo "=============================================================="
echo "DERIVE TRUTH TRANSDUCTIONS"
echo "  job    : ${SLURM_JOB_ID:-<none>} on $(hostname)  cores=$CORES"
echo "  commit : $(git -C "$REPO" rev-parse --short HEAD)"
echo "  truth  : $TRUTH"
echo "  ref    : $REF"
echo "  out    : $OUT"
echo "  started: $(date '+%F %T')"
echo "=============================================================="

# Alu is deliberately absent: a non-autonomous SINE has no internal Pol II promoter and
# cannot 3'-read-through, so there is nothing to derive. The tool refuses --family Alu
# rather than emitting an all-negative table that would read as a measurement.
for FAM in L1 SVA; do
    echo
    echo "=== $FAM ==="
    python "$REPO/benchmark/benchmark_truth_set.py" derive-td \
        --truth-dir "$TRUTH" --family "$FAM" \
        --reference "$REF" --rep-lib-dir "$REPLIB" \
        --bwa "$BWA" --threads "$CORES" \
        --workdir "$OUT/work_$FAM" \
        --out "$OUT/truth_td_$FAM.tsv"
done

# One table for the scorer: it keys on (chrom, pos), which is unique across families here.
cat "$OUT/truth_td_L1.tsv" > "$OUT/truth_td.tsv"
grep -hv '^#' "$OUT/truth_td_SVA.tsv" >> "$OUT/truth_td.tsv"

echo
echo "--- combined table ---"
awk -F'\t' '!/^#/{n[$3"\t"$5]++} END{for(k in n) printf "  %-30s %d\n", k, n[k]}' \
    "$OUT/truth_td.tsv" | sort
echo
echo "--- why records were left INDETERMINATE ---"
awk -F'\t' '!/^#/ && $5=="TRANSDUCTION_INDETERMINATE"{
        r=$10; sub(/_[0-9.]+$/,"_N",r); n[$3"\t"r]++}
     END{for(k in n) printf "  %-52s %d\n", k, n[k]}' "$OUT/truth_td.tsv" | sort
echo
sha256sum "$OUT"/truth_td*.tsv

# Score BASELINE-0 with it immediately: a table nobody has scored against is a table
# whose column meanings have not been checked.
echo
echo "=============================================================="
echo "BASELINE-0 scored against the derived table (a smoke test of the table itself)"
echo "=============================================================="
B0="$SCRATCH/baseline0_hg002_L1-Alu-SVA/BASELINE-0/HG002.wgs.vcf"
if [ -f "$B0" ]; then
    python "$REPO/benchmark/benchmark_truth_set.py" score \
        --calls "$B0" --truth-dir "$TRUTH" --td-table "$OUT/truth_td.tsv" \
        | tee "$OUT/baseline0_score.txt"
else
    echo "  (BASELINE-0 VCF not found at $B0 -- skipped)"
fi

echo
echo "finished $(date '+%F %T')"
