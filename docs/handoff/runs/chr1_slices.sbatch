#!/bin/bash
#SBATCH -c 12
#SBATCH -t 12:00:00
#SBATCH --mem=32G
#SBATCH -p priopark
#SBATCH --account=park_contrib
#SBATCH -o /n/scratch/users/c/cos689/xtea2_testdata/chr1_slice_%j.out
#SBATCH --mail-type=NONE

# =============================================================================
# chr1 slices -- DOWNSTREAM_OVERHAUL_PLAN.md Stage 0, Freeze item 2.
#
# "Without it every gate is a ~50 min job and iteration dies" (section 7.1). The rewrite
# spans five stages of edit-measure-edit; a ~5 min loop is the difference between the
# plan being executable and not.
#
#     sbatch docs/handoff/runs/chr1_slices.sbatch hg002
#     sbatch docs/handoff/runs/chr1_slices.sbatch smht005
#
# Two separate jobs, not one: the HG002 slice is small and unblocks the germline loop
# immediately, while the SMHT005 slice reads chr1 out of ~561 GB of CRAM across five
# files and is much longer. Chaining them would hold the cheap one hostage to the
# expensive one.
#
# WHY chr1: it is what every existing comparator used (NA12878.chr1.bam, the chr1
# goldens), so the slice is comparable to the frozen chr1 history, and it is the largest
# chromosome -- the most MEI content per slice.
# =============================================================================
set -euo pipefail

WHICH="${1:?usage: sbatch chr1_slices.sbatch <hg002|smht005|hapmap>}"

SCRATCH=/n/scratch/users/c/cos689/xtea2_testdata
TESTDATA=/n/data1/hms/dbmi/park/corinne/claude_sandbox/test_data
REF=/n/data1/hms/dbmi/park-smaht_dac/ref/GRCh38_no_alt/hg38_no_alt.fa
THREADS="${SLURM_CPUS_PER_TASK:-12}"

source /n/app/conda/miniforge3/24.11.3-0/etc/profile.d/conda.sh
conda activate xtea2
export PATH=/home/cos689/software/bin:$PATH

echo "job ${SLURM_JOB_ID:-<none>} on $(hostname), $THREADS threads, started $(date '+%F %T')"
echo "samtools: $(samtools --version | head -1)"

case "$WHICH" in

hg002)
    SRC="$TESTDATA/HG002/alignment/HG002.novaseq.pcr-free.30x.dedup.grch38.bam"
    DSTDIR="$TESTDATA/HG002/chr1"
    DST="$DSTDIR/HG002.chr1.bam"
    [ -f "$SRC" ] || { echo "FATAL: $SRC missing"; exit 1; }
    mkdir -p "$DSTDIR"

    echo "=== HG002 chr1 slice ==="
    echo "  src: $SRC"
    echo "  dst: $DST"
    # -b BAM out (source is BAM; no reference needed for decode).
    # NOT filtering flags: the pipeline's own MAPQ/clip logic must see the same read
    # population it would see genome-wide, or the slice stops predicting the full run.
    samtools view -@ "$THREADS" -b -o "$DST" "$SRC" chr1
    samtools index -@ "$THREADS" "$DST"
    samtools quickcheck -v "$DST" && echo "  quickcheck OK"

    echo "  size:  $(du -h "$DST" | cut -f1)"
    echo "  reads: $(samtools view -@ "$THREADS" -c "$DST")"
    # The slice keeps the FULL @SQ dictionary (3366 contigs) -- samtools view does not
    # prune the header. That is wanted: the pipeline is proven on contig-rich references
    # and pruning would make the slice unrepresentative of the full run.
    echo "  @SQ:   $(samtools view -H "$DST" | grep -c '^@SQ') (expect 3366, unpruned)"
    ;;

smht005)
    BROAD=/n/data1/hms/dbmi/park-smaht_dac/DATA/GCC_Broad/SMHT005/SMHT005-3G/illuminaNovaseq_bulkWgs/seq_data
    UW=/n/data1/hms/dbmi/park-smaht_dac/DATA/GCC_UW/SMHT005/SMHT005-3G/illuminaNovaseq_bulkWgs/seq_data
    DSTDIR="$TESTDATA/SMHT005-3G/chr1"
    mkdir -p "$DSTDIR"

    # All five, pooled -- section 7.2 says use all of them via the multi-BAM bam_list.
    CRAMS=(
        "$BROAD/SMAFI5HKGSER.cram"
        "$BROAD/SMAFI6WNA1XM.cram"
        "$UW/SMAFI4852JH2.cram"
        "$UW/SMAFI4D3DFD4.cram"
        "$UW/SMAFIH9MAGML.cram"
    )
    for c in "${CRAMS[@]}"; do
        [ -f "$c" ]        || { echo "FATAL: $c missing"; exit 1; }
        [ -f "$c.crai" ]   || { echo "FATAL: $c.crai missing"; exit 1; }
    done
    [ -f "$REF" ] || { echo "FATAL: reference $REF missing"; exit 1; }

    echo "=== SMHT005-3G chr1 slices (5 CRAMs, ~561 GB total source) ==="
    LIST="$DSTDIR/bam_list.txt"
    : > "$LIST"
    for c in "${CRAMS[@]}"; do
        base=$(basename "$c" .cram)
        dst="$DSTDIR/${base}.chr1.cram"
        echo "--- $base"
        # Stay in CRAM: the mosaic input IS CRAM, and the per-file ftag machinery that
        # keeps same-named reads from different files apart (section 7.2 watch-item) has
        # never run at five-file scale. Converting to BAM here would quietly change the
        # thing under test.
        samtools view -@ "$THREADS" -C -T "$REF" -o "$dst" "$c" chr1
        samtools index -@ "$THREADS" "$dst"
        samtools quickcheck -v "$dst" || { echo "FATAL: quickcheck failed on $dst"; exit 1; }
        printf '%s\n' "$dst" >> "$LIST"
        echo "    -> $dst  $(du -h "$dst" | cut -f1)"
    done

    echo
    echo "--- bam_list for the pipeline ---"
    cat "$LIST"
    echo
    echo "REMINDER: pass --ignore-rg. The five CRAMs carry different @RG SM values"
    echo "(SMATSBL826FN vs SMATSKLUWAXC); @SQ is byte-identical once UR: is stripped, so"
    echo "validate_headers passes and --ignore-rg is the only flag needed (plan 7.2)."
    echo
    echo "--- @RG SM per slice (should still differ; slicing does not touch @RG) ---"
    while read -r f; do
        printf '  %-22s %s\n' "$(basename "$f")" \
            "$(samtools view -H "$f" | awk '/^@RG/{for(i=1;i<=NF;i++) if($i ~ /^SM:/) print $i}' | sort -u | tr '\n' ' ')"
    done < "$LIST"
    ;;

hapmap)
    # HapMap mixture -- the mosaic benchmark as of 2026-08-12, REPLACING SMHT005-3G
    # (plan sec 7.2: SMHT005-3G's 15 sMEI loci carry no short-read signal above background).
    #
    # Arm: BCM, 4 BAMs, ~466x pooled, 1443 GB source. Corinne's call 2026-08-12.
    # All four share SM:SMACUIX7SY21, so unlike smht005 this arm needs NO --ignore-rg.
    # Source is BAM, so slice to BAM -- no reference needed for decode.
    BCM=/n/data1/hms/dbmi/park-smaht_dac/DATA/GCC_BCM/HapMap_Mixture/illuminaNovaseq_bulkWgs/seq_data
    DSTDIR="$TESTDATA/HapMap_Mixture/bcm_chr1"
    mkdir -p "$DSTDIR"

    BAMS=(
        "$BCM/SMAFI3EQANWU.bam"
        "$BCM/SMAFIB1SYCDQ.bam"
        "$BCM/SMAFIVF1NGTJ.bam"
        "$BCM/SMAFIXLDO5QN.bam"
    )
    for b in "${BAMS[@]}"; do
        [ -f "$b" ]     || { echo "FATAL: $b missing"; exit 1; }
        [ -f "$b.bai" ] || { echo "FATAL: $b.bai missing"; exit 1; }
    done

    echo "=== HapMap mixture, BCM chr1 slices (4 BAMs, ~1443 GB source, ~466x pooled) ==="
    LIST="$DSTDIR/bam_list.txt"
    : > "$LIST"
    for b in "${BAMS[@]}"; do
        base=$(basename "$b" .bam)
        dst="$DSTDIR/${base}.chr1.bam"
        echo "--- $base"
        samtools view -@ "$THREADS" -b -o "$dst" "$b" chr1
        samtools index -@ "$THREADS" "$dst"
        samtools quickcheck -v "$dst" || { echo "FATAL: quickcheck failed on $dst"; exit 1; }
        printf '%s\n' "$dst" >> "$LIST"
        echo "    -> $dst  $(du -h "$dst" | cut -f1)"
    done

    echo
    echo "--- bam_list for the pipeline ---"
    cat "$LIST"
    echo
    # Assert the no-ignore-rg claim rather than restating it: if BCM ever ships a file
    # with a different SM this must fail loudly, not silently need a flag nobody passed.
    NSM=$(while read -r f; do
              samtools view -H "$f" | awk '/^@RG/{for(i=1;i<=NF;i++) if($i ~ /^SM:/) print $i}'
          done < "$LIST" | sort -u | wc -l)
    echo "--- distinct @RG SM across the 4 slices: $NSM ---"
    if [ "$NSM" -eq 1 ]; then
        echo "    one SM -> --ignore-rg NOT required for this arm (as expected)"
    else
        echo "FATAL: expected 1 distinct SM, found $NSM -- --ignore-rg assumption broken"; exit 1
    fi
    echo
    echo "Truth set: /n/data1/hms/dbmi/park/corinne/hapmap_mei_calls/hapmap_mix_calls/truth_set_10_18_2024"
    echo "  chr1 tier1 somatic loci: L1 32, Alu 196, SVA 12"
    ;;

*)
    echo "FATAL: unknown target '$WHICH' (expected hg002, smht005 or hapmap)"; exit 2
    ;;
esac

echo
echo "finished $(date '+%F %T')"
