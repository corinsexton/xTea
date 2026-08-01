#!/bin/bash
#SBATCH -c 24
#SBATCH -t 8:00:00
#SBATCH --mem=64G
#SBATCH -p priopark
#SBATCH --account=park_contrib
#SBATCH -o /n/scratch/users/c/cos689/xtea2_testdata/port_wgs_%j.out
#SBATCH --mail-type=NONE

# =============================================================================
# WGS run of the new implementation: xTea2 with the C++ cns collection (--cns-emit-scan).
#
# This is the only run needed for a three-way whole-genome comparison, because the other
# two arms already exist on disk and are intact (VCFs, temp SAMs, step timings):
#
#   master_serial_out    job 48825977  master 46520ab, serial -n 24   50:20   40.6 GB
#   b10_out              job 49065573  xTea2, Python cns, B10-fixed   59:59   18.59 GB
#   <this>                             xTea2, C++ cns collection
#
# Re-running master would cost 50 min of node time to reproduce a result already on file
# from the same node class, so it is not re-run. Note what master carries that the other
# two do not: the B10 disc double-count, since 46520ab predates that fix. That is the
# point of keeping it as an arm -- it is the reference implementation as published.
#
# Frozen worktree, not the live checkout: run_xtea spawns a fresh x_TEA_main per step, so
# a job sharing the checkout picks up mid-run edits.
# =============================================================================
set -euo pipefail

SCRATCH=/n/scratch/users/c/cos689/xtea2_testdata
WT="${XTEA2_WORKTREE:-$SCRATCH/port_wt}"
OUT="${XTEA2_OUT:-$SCRATCH/port_wgs_out}"

source /n/app/conda/miniforge3/24.11.3-0/etc/profile.d/conda.sh
conda activate xtea2
export PATH=/home/cos689/software/bin:$PATH
export XTEA_BASIC_INFO_SEED=42   # pin the insert-size draw, as every comparator run did

[ -d "$WT" ] || { echo "FATAL: worktree $WT missing"; exit 1; }
[ -x "$WT/src/csrc/xtea_scan" ] || { echo "FATAL: $WT/src/csrc/xtea_scan missing"; exit 1; }

rm -rf "$OUT"; mkdir -p "$OUT"
cd "$WT"

echo "=============================================================="
echo "xTea2 WGS -- C++ cns collection (--cns-emit-scan)"
echo "  job    : ${SLURM_JOB_ID:-<none>} on $(hostname)  cores=${SLURM_CPUS_PER_TASK:-24}"
echo "  commit : $(git rev-parse --short HEAD) $(git log -1 --format=%s | cut -c1-66)"
echo "  dirty  : $(git status --porcelain | wc -l) file(s) (expect 0)"
echo "  started: $(date '+%F %T')"
echo "=============================================================="

T0=$(date +%s)
/usr/bin/time -v python xtea/run_xtea.py \
  -i "$SCRATCH/NA12878.wgs.bam" \
  -s NA12878.wgs \
  -o "$OUT" \
  --reference /home/cos689/hg38_no_alt.fa \
  --rep-lib-dir /home/cos689/xTea/rep_lib_annotation \
  --gene-gff3 /home/cos689/data1/corinne/ref/gencode.v45.basic.annotation.gff3 \
  --blacklist /home/cos689/xTea/rep_lib_annotation/blacklist/hg38/centromere.bed \
  --types L1,Alu,SVA \
  --scan-bin "$WT/src/csrc/xtea_scan" \
  --bwa /home/cos689/software/bin/bwa \
  --scan-min-clip 2 \
  --col5-families alu,sva \
  --cns-emit-scan \
  --parallel-downstream \
  -n "${SLURM_CPUS_PER_TASK:-24}" \
  2>&1 | tee "$OUT/time_v.txt"
WALL=$(( $(date +%s) - T0 ))
S="$OUT/NA12878.wgs"

echo
echo "##############################################################"
echo "# TIMING / MEMORY / DISK"
echo "##############################################################"
printf 'wall %02d:%02d:%02d (%ds)\n' $((WALL/3600)) $((WALL%3600/60)) $((WALL%60)) "$WALL"
grep -E "Percent of CPU|Elapsed \(wall|Maximum resident" "$OUT/time_v.txt" 2>/dev/null || true
echo
echo "--- phases ---"
awk -F'\t' '$1=="phase"{printf "  %-22s %8.1fs  childRSS %5.2f GB\n",$2,$3,$4}' "$S/step_timings.tsv" 2>/dev/null || echo "  MISSING"
echo
echo "--- the cns step should now be on the bypass (parse+filter only) ---"
awk -F'\t' '$1 ~ /^step/ && $2 ~ /cns/ {printf "  %-22s %8.1fs\n",$2,$3}' "$S/step_timings.tsv" 2>/dev/null || true
echo
echo "--- C++ collection internals ---"
grep -a '^\[cns_fan\]' "$OUT/time_v.txt" 2>/dev/null | sed 's/^/  /' || echo "  (none)"
grep -ac '^\[bwa_sam\]' "$OUT/time_v.txt" 2>/dev/null | sed 's/^/  bwa_sam calls: /' || true
echo
echo "--- output size / no spill expected ---"
du -sh "$OUT" 2>/dev/null | sed 's/^/  total  /'
ls -la "$S"/cns_spill_*.txt 2>/dev/null || echo "  no spill files (expected -- the spill path is deleted)"
echo
echo "--- call counts (comparison done separately, against master + b10) ---"
# `f=$(ls MISSING 2>/dev/null | head -1)` USED TO BE HERE and is a landmine under
# `set -euo pipefail`: 2>/dev/null hides ls's message but not its exit code 2, pipefail carries
# that through the pipe into the assignment, and set -e kills the script. It never fired here
# only because merged/ always existed -- but the identical line in baseline0_hg002.sbatch killed
# job 49440169 AFTER a fully successful 22-minute pipeline run. Same fix as there: a glob that
# tolerates no match, and grep -c's exit 1 on a zero count swallowed rather than doubled.
first_match() {  # first existing path matching $1, or empty
    local m
    for m in $1; do [ -e "$m" ] && { printf '%s' "$m"; return 0; }; done
    return 0
}
n_records() { [ -f "$1" ] || { echo 0; return 0; }; grep -vc '^#' "$1" || true; }
for pair in "L1:LINE1" "Alu:ALU" "SVA:SVA"; do
  t=${pair%%:*}; alt=${pair##*:}
  f=$(first_match "$S/merged/*_${alt}.vcf")
  [ -n "$f" ] && printf '  %-4s %s\n' "$t" "$(n_records "$f")"
done
mf=$(first_match "$S/*.vcf")
[ -n "$mf" ] && echo "  merged $(n_records "$mf")"
echo
echo "finished $(date '+%F %T')"
