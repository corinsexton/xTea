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
# CNS_EMIT=1 passes --cns-emit-scan; CNS_EMIT=0 passes --no-cns-emit-scan. BOTH arms name
# the flag explicitly, and the OFF arm MUST: as of 2026-08-13 --cns-emit-scan defaults ON
# for every mode, so the old "omit the flag" spelling of OFF would now silently run the ON
# path and the A/B would compare a run against itself.
#
# THE A/B THIS SCRIPT EXISTS FOR IS ANSWERED (50434768 ON vs 50434907 OFF, chr1 mosaic,
# 83,775 candidates): L1:cns_filter 19.0 s ON vs 1833.5 s OFF, ~96x, and ON still wins by
# ~14 min net after its own 954.8 s collection phase. The prior contrary measurement
# (HANDOFF_2026-07-26, "31 min slower") was taken against the candidate-agnostic
# implementation that has since been deleted, so it does not transfer. Kept runnable
# because a default this load-bearing should stay falsifiable, not because it is open.
CNS_EMIT=${CNS_EMIT:-1}
if [ "$CNS_EMIT" = "1" ]; then CNS_EMIT_ARG="--cns-emit-scan"; SAMPLE=hapmap_bcm_chr1
else CNS_EMIT_ARG="--no-cns-emit-scan"; SAMPLE=hapmap_bcm_chr1_nocnsemit; fi

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
    # --blacklist: baseline0_hg002.sbatch:121 and repro_master_vs_xtea2.sbatch:166 both pass
    # centromere.bed, so every germline result was produced with it and omitting it here was
    # the same class of error as omitting --cns-emit-scan. It also matters more for mosaic:
    # x_mosaic_calling's own header lists "add blacklist, like gnomAD-SV" as design intent.
    #
    # NOT the truth set's own hg005_blacklist.bed, deliberately. That file is HG005's
    # insertions -- the 83.5% background donor -- and filtering CALLS with it before scoring
    # against the same truth set would be circular. The scorer already excludes HG005 at
    # SCORING time (benchmark_hapmap_mosaic.py), which is the right place for it.
    --blacklist "$REPLIB/blacklist/hg38/centromere.bed"
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
