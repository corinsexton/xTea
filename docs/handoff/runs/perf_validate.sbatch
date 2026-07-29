#!/bin/bash
#SBATCH -c 24
#SBATCH -t 8:00:00
#SBATCH --mem=64G
#SBATCH -p priopark
#SBATCH --account=park_contrib
#SBATCH -o /n/scratch/users/c/cos689/xtea2_testdata/perf_validate_%j.out
#SBATCH --mail-type=NONE

# =============================================================================
# End-to-end validation of the transduction fixes at fb1d766: timing, memory, disk.
#
# What is being validated (both landed in fb1d766, both projected but NOT yet measured
# end-to-end -- that is the whole point of this run):
#
#   1. Per-process BAM handle/meta cache (x_alignments.get_cached_bam). The transduction
#      step was paying 2 BAM opens + a 3,366-contig get_tid loop PER candidate site, over
#      the ~54K raw disc sites. Projected saving: ~5,112s of worker CPU per family, which
#      is ~10.7 min of wall at -n 8. Derived from a microbenchmark (0.0948s/site) that
#      agrees with the in-situ profiling job 49078859 (0.0928s) to within 4%.
#
#   2. Transduction-only core split. Alu has no transduction flank, so transduction now
#      splits over the flanked types (L1, SVA -> 12 cores each at -c 24) instead of all
#      three (8 each). Projected on the measured -n 8 vs -n 24 scaling: ~1.1x on that step.
#
# The number that matters: `transduction` in step_timings.tsv for L1 and SVA. Baselines to
# beat, from the completed cpp-cns WGS run 49070030 (same node class, -n 8 each):
#
#      type   cns_filter  transduction  sibling  postF_all  postF_hc   type total
#      Alu       116s          1.5s       38s        7s        6s         169s
#      L1         55s         1481s       55s        5s        6s        1602s
#      SVA        57s         1772s       52s        1s        1s        1884s   (= whole phase)
#
# and the whole-run three-way table (49070030 is the arm this run replaces):
#      master (pre-B10)  50:20   40.6 GB   6.0 G
#      python cns        60:37   21.99 GB  7.1 G
#      cpp cns           54:39   20.56 GB  7.0 G
#
# Run chr1 first as a cheap correctness+direction check, then WGS for the real numbers:
#   sbatch --export=ALL,XTEA2_SCALE=chr1 -t 3:00:00 docs/handoff/runs/perf_validate.sbatch
#   sbatch --export=ALL,XTEA2_SCALE=wgs        docs/handoff/runs/perf_validate.sbatch
# The type list is a POSITIONAL arg, never --export (see TYPES below):
#   sbatch --export=ALL,XTEA2_SCALE=chr1,XTEA2_OUT=/path/out ... perf_validate.sbatch L1,SVA
#
# --keep-intermediates is passed deliberately: disk is one of the numbers being collected,
# and the 6.0-7.1 G figures above are all pre-cleanup, so cleaning would make this run
# incomparable. The post-cleanup size is reported separately at the end.
#
# Frozen worktree, not the live checkout: run_xtea spawns a fresh x_TEA_main per step, so a
# job sharing the checkout would pick up mid-run edits.
# =============================================================================
set -euo pipefail

SCRATCH=/n/scratch/users/c/cos689/xtea2_testdata
WT="${XTEA2_WORKTREE:-$SCRATCH/perf_wt}"
SCALE="${XTEA2_SCALE:-chr1}"

case "$SCALE" in
  chr1) BAM="$SCRATCH/NA12878.chr1.bam"; SAMPLE=NA12878.chr1 ;;
  wgs)  BAM="$SCRATCH/NA12878.wgs.bam";  SAMPLE=NA12878.wgs  ;;
  *)    echo "FATAL: XTEA2_SCALE must be chr1 or wgs, got '$SCALE'"; exit 1 ;;
esac
OUT="${XTEA2_OUT:-$SCRATCH/perf_${SCALE}_out}"

# The type list comes in as $1, NOT as an env var. `sbatch --export=ALL,FOO=a,b,c` splits on
# every comma, so --export=...,XTEA2_TYPES=L1,Alu,SVA silently sets XTEA2_TYPES=L1 and then tries
# to export "Alu" and "SVA" as separate names. That cost a round of three jobs -- all three ran
# L1-only while reporting success-shaped output. Positional args are not comma-parsed.
#
# Why a type-list knob exists at all: it drives the A/B that isolates the BAM cache from the core
# split. The cache must be output-neutral; the core split changes thread counts, and this pipeline
# is KNOWN to move calls with thread count (job 49075238: master's own -n 8 output differed from
# master's own -n 24 output on identical inputs). Comparing fb1d766 against its parent with all
# three types therefore confounds the two. Running "L1,SVA" makes transduction's -n identical in
# both arms -- parent: 24//2 types = 12; fb1d766: 24//2 flanked types = 12 -- so any remaining
# diff is attributable to the cache alone.
TYPES="${1:-L1,Alu,SVA}"

source /n/app/conda/miniforge3/24.11.3-0/etc/profile.d/conda.sh
conda activate xtea2
export PATH=/home/cos689/software/bin:$PATH
export XTEA_BASIC_INFO_SEED=42   # pin the insert-size draw, as every comparator run did

[ -d "$WT" ] || { echo "FATAL: worktree $WT missing"; exit 1; }
[ -x "$WT/src/csrc/xtea_scan" ] || { echo "FATAL: $WT/src/csrc/xtea_scan missing"; exit 1; }
[ -f "$BAM" ] || { echo "FATAL: bam $BAM missing"; exit 1; }

rm -rf "$OUT"; mkdir -p "$OUT"
cd "$WT"

# --keep-intermediates landed in fb1d766. A/B arms run against the PARENT worktree, which does
# not have the flag and exits 2 on it ("unrecognized arguments"), so probe rather than assume.
if python xtea/run_xtea.py --help 2>&1 | grep -q -- --keep-intermediates; then
  KEEP_FLAG=--keep-intermediates
else
  KEEP_FLAG=
  echo "NOTE: this worktree predates --keep-intermediates; intermediates are retained anyway"
  echo "      because cleanup only runs when the flag exists to be omitted."
fi

echo "=============================================================="
echo "xTea2 $SCALE -- transduction fixes validation"
echo "  job    : ${SLURM_JOB_ID:-<none>} on $(hostname)  cores=${SLURM_CPUS_PER_TASK:-24}"
echo "  commit : $(git rev-parse --short HEAD) $(git log -1 --format=%s | cut -c1-66)"
echo "  dirty  : $(git status --porcelain | wc -l) file(s) (expect 0)"
echo "  bam    : $BAM"
echo "  types  : $TYPES"
echo "  pysam  : $(python -c 'import pysam; print(pysam.__version__)')"
echo "  started: $(date '+%F %T')"
echo "=============================================================="

T0=$(date +%s)
/usr/bin/time -v python xtea/run_xtea.py \
  -i "$BAM" \
  -s "$SAMPLE" \
  -o "$OUT" \
  --reference /home/cos689/hg38_no_alt.fa \
  --rep-lib-dir /home/cos689/xTea/rep_lib_annotation \
  --gene-gff3 /home/cos689/data1/corinne/ref/gencode.v45.basic.annotation.gff3 \
  --blacklist /home/cos689/xTea/rep_lib_annotation/blacklist/hg38/centromere.bed \
  --types "$TYPES" \
  --scan-bin "$WT/src/csrc/xtea_scan" \
  --bwa /home/cos689/software/bin/bwa \
  --scan-min-clip 2 \
  --col5-families alu,sva \
  --cns-emit-scan \
  --parallel-downstream \
  $KEEP_FLAG \
  -n "${SLURM_CPUS_PER_TASK:-24}" \
  2>&1 | tee "$OUT/time_v.txt"
WALL=$(( $(date +%s) - T0 ))
S="$OUT/$SAMPLE"

echo
echo "##############################################################"
echo "# TIMING / MEMORY / DISK  ($SCALE)"
echo "##############################################################"
printf 'wall %02d:%02d:%02d (%ds)\n' $((WALL/3600)) $((WALL%3600/60)) $((WALL%60)) "$WALL"
grep -E "Percent of CPU|Elapsed \(wall|Maximum resident" "$OUT/time_v.txt" 2>/dev/null || true

echo
echo "--- phases ---"
awk -F'\t' '$1=="phase"{printf "  %-22s %8.1fs  childRSS %5.2f GB\n",$2,$3,$4}' \
    "$S/step_timings.tsv" 2>/dev/null || echo "  MISSING step_timings.tsv"

echo
echo "--- THE NUMBER THIS RUN IS ABOUT: per-type steps ---"
echo "  (WGS baselines from 49070030: L1 transduction 1481s, SVA 1772s, Alu 1.5s)"
awk -F'\t' '$1 ~ /^step/ {printf "  %-34s %8.1fs\n",$2,$3}' "$S/step_timings.tsv" 2>/dev/null || true

echo
echo "--- transduction only, side by side with the 49070030 WGS baseline ---"
awk -F'\t' -v scale="$SCALE" '
  $1 ~ /^step/ && $2 ~ /transduction/ {
     split($2, a, ":"); t=a[1]; if (t=="") t=$2;
     base["L1"]=1481; base["SVA"]=1772; base["Alu"]=1.5;
     b = (t in base) ? base[t] : 0;
     if (scale=="wgs" && b>0)
        printf "  %-12s %8.1fs   baseline %7.1fs   %+6.1f%%  (%.2fx)\n", t,$3,b,($3-b)/b*100,b/$3;
     else
        printf "  %-12s %8.1fs   (chr1: no WGS baseline to compare)\n", t,$3;
  }' "$S/step_timings.tsv" 2>/dev/null || true

echo
echo "--- the core split actually took effect? (expect L1/SVA n=12, Alu n=1) ---"
grep -aE "transduction n=|per-type downstreams CONCURRENT" "$OUT/time_v.txt" 2>/dev/null | sed 's/^/  /' || echo "  (not found)"

echo
echo "--- disk ---"
du -sh "$OUT" 2>/dev/null | sed 's/^/  total (with intermediates)  /'
mf=$(ls "$S"/*.vcf 2>/dev/null | head -1) || true
[ -n "$mf" ] && du -sh "$mf" | sed 's/^/  merged VCF (what survives cleanup)  /' || true
echo "  largest intermediates:"
du -ah "$S" 2>/dev/null | sort -rh | sed -n '2,9p' | sed 's/^/    /'

echo
echo "--- call counts ---"
for pair in "L1:LINE1" "Alu:ALU" "SVA:SVA"; do
  t=${pair%%:*}; alt=${pair##*:}
  f=$(ls "$S"/merged/*_${alt}.vcf 2>/dev/null | head -1) || true   # non-matching glob + pipefail
  [ -n "$f" ] && printf '  %-4s %s\n' "$t" "$(grep -vc '^#' "$f")" || true
done
[ -n "$mf" ] && echo "  merged $(grep -vc '^#' "$mf")" || true

echo
echo "--- raw disc site counts (the real transduction denominator) ---"
for t in L1 SVA; do
  f="$S/$t/candidate_list_from_disc.txt.clip_sites_raw_disc.txt"
  [ -f "$f" ] && printf '  %-4s %s sites\n' "$t" "$(wc -l < "$f")" || true
done

echo
echo "finished $(date '+%F %T')"
