#!/bin/bash
#SBATCH -c 24
#SBATCH -t 12:00:00
#SBATCH --mem=192G
#SBATCH -p priopark
#SBATCH --account=park_contrib
#SBATCH -o /n/scratch/users/c/cos689/xtea2_testdata/mosaic_hapmap_chr1_%j.out
#SBATCH --mail-type=NONE

# =============================================================================
# The FIRST mosaic run: --mosaic end to end on the HapMap mixture, BCM arm, chr1.
#
#     sbatch docs/handoff/runs/mosaic_hapmap_chr1.sbatch
#
# Input: 4 BCM BAMs sliced to chr1 (120 GB, ~466x pooled), one SM tag so NO
# --ignore-rg. Truth: 240 tier1 somatic loci on chr1 (L1 32, Alu 196, SVA 12) at
# known VAF, 0.25%-8.25%.
#
# -c 24 REQUIRES priopark: short/medium cap at 20 CPUs.
#
# --mem=192G rather than the germline 64G, and the reason is not caution. xtea-scan's
# peak is O(distinct candidate sites), not O(reads), and this run multiplies that term
# twice over: ~466x depth against germline WGS's ~30x, AND clip/disc cutoffs dropped to
# 2/2, which is what admits the low-AF sites this mode exists for. Germline WGS peaked
# at 18 GB across all three families; chr1 is ~8% of the genome but the site count per
# base is what scales here, so the germline number is not a useful anchor. If this
# lands far under 192G, cut it -- but a first run of an untested memory profile is the
# wrong place to discover an OOM after 6 hours.
#
# WHY chr1 AND NOT THE FULL 1443 GB: MOSAIC_PLAN calls a chr1 slice mandatory for
# iteration. The full-depth timing run is M.3 and comes after this one produces a
# scoreable call set.
# =============================================================================
set -euo pipefail

REPO=/n/data1/hms/dbmi/park/corinne/claude_sandbox/xTea
SLICE=/n/data1/hms/dbmi/park/corinne/claude_sandbox/test_data/HapMap_Mixture/bcm_chr1
REF=/n/data1/hms/dbmi/park-smaht_dac/ref/GRCh38_no_alt/hg38_no_alt.fa
REPLIB=/home/cos689/xTea/rep_lib_annotation
GFF3=/home/cos689/data1/corinne/ref/gencode.v45.basic.annotation.gff3
# Default OUT carries the job id, so a rerun cannot clobber the previous run's outputs.
# It did: relaunching chr1 as 50440585 overwrote 50434768's step_timings.tsv mid-analysis,
# destroying the WITH half of the --cns-emit-scan A/B while it was being read. Pass OUT
# explicitly to reuse a tree on purpose.
OUT=${OUT:-/n/scratch/users/c/cos689/xtea2_testdata/mosaic_hapmap_bcm_chr1_${SLURM_JOB_ID:-manual}}
# CNS_EMIT=1 (default) passes --cns-emit-scan; CNS_EMIT=0 omits it. An A/B on chr1 is the
# cheapest way to settle which cns collection path is faster AT MOSAIC CANDIDATE COUNTS --
# the only direct measurement we have (HANDOFF_2026-07-26) says the C++ path is 31 min
# SLOWER, but that was at germline's ~1,500 candidates, where there is nothing to merge and
# its whole-genome mate-gather is pure overhead. chr1 mosaic has ~82,000. The crossover is
# candidate count, and nobody has measured it.
CNS_EMIT=${CNS_EMIT:-1}
if [ "$CNS_EMIT" = "1" ]; then CNS_EMIT_ARG="--cns-emit-scan"; SAMPLE=hapmap_bcm_chr1
else CNS_EMIT_ARG=""; SAMPLE=hapmap_bcm_chr1_nocnsemit; fi

source /n/app/conda/miniforge3/24.11.3-0/etc/profile.d/conda.sh
conda activate xtea2
export PATH=/home/cos689/software/bin:$PATH

echo "job ${SLURM_JOB_ID:-<none>} on $(hostname), ${SLURM_CPUS_PER_TASK:-24} cpus, $(date '+%F %T')"
echo "commit: $(git -C $REPO rev-parse --short HEAD)"

# The scan binary must be newer than its sources or the run is testing a stale engine --
# the exact trap tests/test_scan_bin_freshness.py exists for. Assert it here too, because
# a batch job does not run pytest.
# PKG_CONFIG_PATH must be EXPORTED, not passed as a make variable: the Makefile calls
# pkg-config through $(shell ...), which does not inherit make-level variable
# assignments. Passing it the wrong way printed "htslib 1.22 via " with an empty
# libdir and carried on -- a build that looks fine and would link nothing.
export PKG_CONFIG_PATH=/home/cos689/software/lib/pkgconfig
make -C "$REPO/src/csrc"
echo "scan bin: $(ls -l --time-style=+%F_%T $REPO/src/csrc/xtea_scan | awk '{print $6, $7}')"

for f in "$SLICE"/*.chr1.bam; do [ -f "$f" ] || { echo "FATAL: no slices in $SLICE"; exit 1; }; done
[ -f "$SLICE/bam_list.txt" ] || { echo "FATAL: no bam_list.txt"; exit 1; }

mkdir -p "$OUT"
echo "=== inputs ==="; cat "$SLICE/bam_list.txt"

ARGS=(
    --mosaic
    -i "$(tr '\n' ' ' < "$SLICE/bam_list.txt")"
    -s "$SAMPLE" -o "$OUT"
    --types L1,Alu,SVA
    ${CNS_EMIT_ARG}
    # CNS_EMIT=1 (default) adds --cns-emit-scan; CNS_EMIT=0 omits it. Toggleable because
    # which path is FASTER at mosaic candidate counts is genuinely unknown -- see below. baseline0_hg002.sbatch and the
    --reference "$REF" --rep-lib-dir "$REPLIB" --gene-gff3 "$GFF3"
    --cores "${SLURM_CPUS_PER_TASK:-24}"
)

# PREFLIGHT. --dry-run resolves every path and builds every step's argv without touching a
# BAM, so an argv mistake fails here in seconds instead of after the job has queued. Added
# because the first attempt (50427531) died on `--types L1 Alu SVA` -- --types takes ONE
# comma-separated value, and a single-type dry-run had not exercised it.
echo "=== preflight (--dry-run) ==="
python3 "$REPO/xtea/run_xtea.py" "${ARGS[@]}" --dry-run > "$OUT/preflight.log" 2>&1 || {
    echo "FATAL: preflight failed; argv is wrong. Last lines:"; tail -20 "$OUT/preflight.log"; exit 1; }
echo "  preflight OK: $(grep -cE '^# ' "$OUT/preflight.log") steps resolved"

# No --ignore-rg: all four BCM files share SM:SMACUIX7SY21 (asserted by chr1_slices.sbatch).
# Omitting it deliberately -- if BCM ever ships a mismatched file the run should FAIL, not
# silently pool two samples.
/usr/bin/time -v python3 "$REPO/xtea/run_xtea.py" "${ARGS[@]}" \
    2>&1

echo
echo "=== mosaic call sets ==="
for t in L1 Alu SVA; do
    f="$OUT/$SAMPLE/$t/candidate_disc_filtered_cns_post_filtering_mosaic.txt"
    printf '  %-4s %s\n' "$t" "$([ -f "$f" ] && wc -l < "$f" || echo MISSING)"
done
echo "=== VCF ==="; ls -la "$OUT/$SAMPLE/"*.vcf 2>/dev/null || echo "  none"
echo "finished $(date '+%F %T')"
