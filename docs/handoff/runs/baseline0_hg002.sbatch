#!/bin/bash
#SBATCH -c 24
#SBATCH -t 10:00:00
#SBATCH --mem=64G
#SBATCH -p priopark
#SBATCH --account=park_contrib
#SBATCH -o /n/scratch/users/c/cos689/xtea2_testdata/baseline0_hg002_%j.out
#SBATCH --mail-type=NONE

# =============================================================================
# BASELINE-0 -- DOWNSTREAM_OVERHAUL_PLAN.md Stage 0, Freeze item 1.
#
# The frozen reference call set the whole downstream rewrite is measured against.
# Stage 1's gate is "HG002 before vs after, byte-identical calls", and that only means
# something if the before is pinned to a commit and a command line.
#
# Run TWICE, once per type set -- pass the type list as $1:
#     sbatch docs/handoff/runs/baseline0_hg002.sbatch L1,Alu,SVA
#     sbatch docs/handoff/runs/baseline0_hg002.sbatch L1
# The single-type arm is not optional. Section 6.2 makes single-type behaviour an explicit
# design requirement, and a single-type regression is invisible to a three-type baseline --
# which is exactly how the mate_in_rep type-set-dependence bug survived until 2026-07-29.
#
# Type list goes as a POSITIONAL argument, never through --export: sbatch --export splits
# on every comma, so "L1,Alu,SVA" would silently arrive as three broken variables and the
# run would look successful.
#
# Taken at 9a43b80 or later, i.e. AFTER both Stage 0 standalone fixes (de96727
# type-set-dependence, 39dd582 the rename). That ordering is deliberate: every later
# baseline is then against correct behaviour rather than against a known bug.
#
# Frozen worktree, not the live checkout -- run_xtea spawns a fresh x_TEA_main per step,
# so a job sharing the checkout picks up mid-run edits.
# =============================================================================
set -euo pipefail

TYPES="${1:?usage: sbatch baseline0_hg002.sbatch <TYPE_LIST>   e.g. L1,Alu,SVA or L1}"
TAG=$(echo "$TYPES" | tr ',' '-')          # L1,Alu,SVA -> L1-Alu-SVA, for the output path

SCRATCH=/n/scratch/users/c/cos689/xtea2_testdata
WT="${XTEA2_WORKTREE:-$SCRATCH/baseline0_wt}"
OUT="${XTEA2_OUT:-$SCRATCH/baseline0_hg002_${TAG}}"
TESTDATA=/n/data1/hms/dbmi/park/corinne/claude_sandbox/test_data/HG002
BAM="$TESTDATA/alignment/HG002.novaseq.pcr-free.30x.dedup.grch38.bam"
SAMPLE=HG002.wgs

source /n/app/conda/miniforge3/24.11.3-0/etc/profile.d/conda.sh
conda activate xtea2
export PATH=/home/cos689/software/bin:$PATH
export XTEA_BASIC_INFO_SEED=42   # pin the insert-size draw, as every comparator run did

[ -d "$WT" ]                  || { echo "FATAL: worktree $WT missing"; exit 1; }
[ -x "$WT/src/csrc/xtea_scan" ] || { echo "FATAL: $WT/src/csrc/xtea_scan not built"; exit 1; }
[ -f "$BAM" ]                 || { echo "FATAL: $BAM missing"; exit 1; }
[ -f "$BAM.bai" ]             || { echo "FATAL: $BAM.bai missing"; exit 1; }
# A baseline taken from a dirty tree cannot be reproduced; refuse rather than mislabel it.
[ -z "$(git -C "$WT" status --porcelain)" ] || { echo "FATAL: worktree is dirty"; exit 1; }

# The mate-in-rep family flag was renamed from --col5-families in 39dd582. This same script
# has to run against PRE-rename worktrees too -- that is the whole point of the Stage 0
# byte-identity gate ("three-type calls byte-identical across the rename+fix"), which needs
# a before and an after produced by an otherwise identical command line. Detect rather than
# hardcode, and fail loud if neither flag exists rather than silently dropping the option
# and changing the signal being measured.
# Absolute path, because this runs before the cd into $WT below.
if python "$WT/xtea/run_xtea.py" --help 2>&1 | grep -q -- "--mate-in-rep-families"; then
    MATE_IN_REP_FLAG=--mate-in-rep-families
elif python "$WT/xtea/run_xtea.py" --help 2>&1 | grep -q -- "--col5-families"; then
    MATE_IN_REP_FLAG=--col5-families
else
    echo "FATAL: worktree exposes neither --mate-in-rep-families nor --col5-families"; exit 1
fi
echo "  mate-in-rep flag for this worktree: $MATE_IN_REP_FLAG"

rm -rf "$OUT"; mkdir -p "$OUT"
cd "$WT"

echo "=============================================================="
echo "BASELINE-0  --  HG002, types=$TYPES"
echo "  job    : ${SLURM_JOB_ID:-<none>} on $(hostname)  cores=${SLURM_CPUS_PER_TASK:-24}"
echo "  commit : $(git rev-parse --short HEAD) $(git log -1 --format=%s | cut -c1-66)"
echo "  dirty  : $(git status --porcelain | wc -l) file(s) (must be 0)"
echo "  bam    : $BAM"
echo "  md5    : 7378e27ba9a623ccb5eb1a4353aaebd9 (per test_data/HG002/PROVENANCE.md)"
echo "  NOTE   : this BAM is hs38DH-aligned, NOT no_alt -- see plan section 7.1. Absolute"
echo "           sensitivity/FDR may sit below the paper's for alignment reasons alone."
echo "  started: $(date '+%F %T')"
echo "=============================================================="

T0=$(date +%s)
/usr/bin/time -v python xtea/run_xtea.py \
  -i "$BAM" \
  -s "$SAMPLE" \
  -o "$OUT" \
  --reference /n/data1/hms/dbmi/park-smaht_dac/ref/GRCh38_no_alt/hg38_no_alt.fa \
  --rep-lib-dir /home/cos689/xTea/rep_lib_annotation \
  --gene-gff3 /home/cos689/data1/corinne/ref/gencode.v45.basic.annotation.gff3 \
  --blacklist /home/cos689/xTea/rep_lib_annotation/blacklist/hg38/centromere.bed \
  --types "$TYPES" \
  --scan-bin "$WT/src/csrc/xtea_scan" \
  --bwa /home/cos689/software/bin/bwa \
  --scan-min-clip 2 \
  "$MATE_IN_REP_FLAG" alu,sva \
  --cns-emit-scan \
  --parallel-downstream \
  --keep-intermediates \
  -n "${SLURM_CPUS_PER_TASK:-24}" \
  2>&1 | tee "$OUT/time_v.txt"
WALL=$(( $(date +%s) - T0 ))
S="$OUT/$SAMPLE"

echo
echo "##############################################################"
echo "# BASELINE-0 SUMMARY -- types=$TYPES"
echo "##############################################################"
printf 'wall %02d:%02d:%02d (%ds)\n' $((WALL/3600)) $((WALL%3600/60)) $((WALL%60)) "$WALL"
grep -E "Percent of CPU|Elapsed \(wall|Maximum resident" "$OUT/time_v.txt" 2>/dev/null || true

echo
echo "--- phases ---"
awk -F'\t' '$1=="phase"{printf "  %-22s %8.1fs  childRSS %5.2f GB\n",$2,$3,$4}' \
    "$S/step_timings.tsv" 2>/dev/null || echo "  MISSING"

echo
echo "--- transduction phase breakdown ([TD-PHASE], unconditional since 6ce36c7) ---"
# NOT time_v.txt. run_xtea spawns each step as a subprocess writing to the driver's own
# log, so the timers -- which go to the STEP's stderr -- never reach the parent's stdout.
# They land in <sample>_xtea2.log.
{ grep -ah '^\[TD-PHASE\]' "$S/${SAMPLE}_xtea2.log" 2>/dev/null \
  | awk -F'\t' '{printf "  rtype=%-3s %-38s %8.1fs\n",$3,$4,$5}'; } || true
{ grep -ah '^\[TD-PHASE\]' "$S/${SAMPLE}_xtea2.log" 2>/dev/null | awk -F'\t' '
    $4~/gather/{g[$3]+=$5} $4~/bwa/{b[$3]+=$5} $4~/sam_parse/{p[$3]+=$5} $4~/TOTAL/{t[$3]=$5}
    END{for (r in t) printf "  rtype=%-3s gather %.1fs (%.1f%%)  bwa %.1fs (%.1f%%)  parse %.1fs (%.1f%%)  TOTAL %.1fs\n",
          r, g[r], g[r]/t[r]*100, b[r], b[r]/t[r]*100, p[r], p[r]/t[r]*100, t[r]}'; } || true

echo
echo "--- call counts ---"
# `f=$(ls MISSING 2>/dev/null | head -1)` is a trap under `set -euo pipefail`, and it is
# the one that killed job 49440169: 2>/dev/null hides ls's MESSAGE but not its exit code 2,
# pipefail propagates that through the pipe into the assignment, and set -e kills the run
# -- AFTER the pipeline had succeeded and written its VCF. Inherited from port_wgs.sbatch,
# where it never fired only because merged/ always existed. Use a glob that tolerates no
# match instead of shelling out to ls.
first_match() {  # first existing path matching $1, or empty
    local m
    for m in $1; do [ -e "$m" ] && { printf '%s' "$m"; return 0; }; done
    return 0
}
# grep -c exits 1 when the count is 0 while still PRINTING 0, so `grep -vc ... || echo 0`
# would emit "0" twice. Swallow the status, keep the number.
n_records() { [ -f "$1" ] || { echo 0; return 0; }; grep -vc '^#' "$1" || true; }

for pair in "L1:LINE1" "Alu:ALU" "SVA:SVA"; do
  t=${pair%%:*}; alt=${pair##*:}
  f=$(first_match "$S/merged/*_${alt}.vcf")
  if [ -n "$f" ]; then printf '  %-4s %s\n' "$t" "$(n_records "$f")"; fi
done
mf=$(first_match "$S/*.vcf")
if [ -n "$mf" ]; then echo "  merged $(n_records "$mf")"; fi

# The baseline artifact itself: the VCFs plus enough provenance to reproduce the command.
echo
echo "--- freezing the baseline ---"
FROZEN="$OUT/BASELINE-0"
mkdir -p "$FROZEN"
find "$S" -name '*.vcf' -print0 | while IFS= read -r -d '' v; do
    cp "$v" "$FROZEN/$(basename "$v")"
done
{
  echo "BASELINE-0"
  echo "types:        $TYPES"
  echo "commit:       $(git rev-parse HEAD)"
  echo "job:          ${SLURM_JOB_ID:-<none>}"
  echo "bam:          $BAM"
  echo "bam md5:      7378e27ba9a623ccb5eb1a4353aaebd9"
  echo "reference:    /n/data1/hms/dbmi/park-smaht_dac/ref/GRCh38_no_alt/hg38_no_alt.fa"
  echo "seed:         XTEA_BASIC_INFO_SEED=$XTEA_BASIC_INFO_SEED"
  echo "wall_s:       $WALL"
  echo "generated:    $(date '+%F %T')"
  echo
  echo "sha256 of every frozen VCF:"
  (cd "$FROZEN" && sha256sum ./*.vcf 2>/dev/null) || echo "  (none)"
} > "$FROZEN/PROVENANCE.txt"
cat "$FROZEN/PROVENANCE.txt"

echo
echo "finished $(date '+%F %T')"
