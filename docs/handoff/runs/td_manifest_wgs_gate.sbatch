#!/bin/bash
#SBATCH -c 12
#SBATCH -t 2:00:00
#SBATCH --mem=64G
#SBATCH -p priopark
#SBATCH --account=park_contrib
#SBATCH -o /n/scratch/users/c/cos689/xtea2_testdata/td_manifest_wgs_%j.out
#SBATCH --mail-type=NONE

# =============================================================================
# --td-manifest AGREEMENT GATE -- WGS scale, real site lists.
#
#     sbatch docs/handoff/runs/td_manifest_wgs_gate.sbatch <worktree>
#
# WHY THIS AND NOT JUST THE UNIT TESTS. tests/test_td_gather.py already runs
# N-entries-together against each-alone and cmps every product, but on a chr1
# SLICE with ~40 sites. Every risk the shared decode introduces is
# scale-dependent and invisible there:
#
#   * the union fetch plan merges 34,197 L1 + 39,452 SVA sites at extnd=954 --
#     overlapping windows, thousands of merged regions, an ownership rule that
#     has to keep partitioning reads exactly at boundaries the solo runs never
#     produce;
#   * the mate gather's per-entry part files carry ~1 GB per family, striped
#     across 12 threads and concatenated -- an ordering bug here is a few
#     records out of millions, not a visibly wrong file;
#   * peak memory now holds the UNION's collected events plus one `needs` map
#     per entry.
#
# So: run L1 alone, SVA alone, then both from one manifest, and require `cmp`
# on all six products. Bytes, not records -- bwa's tie-break among equal-scoring
# placements depends on a read's index in its input, so a reordered FASTQ moves
# borderline calls (see the header of src/csrc/td_gather.hpp).
#
# The site lists are the REAL ones a WGS run produced (sf_raw_slct from the
# Stage 2c gate's OFF arm) at that run's own extnd, so this is the production
# shape, not a constructed one.
# =============================================================================
set -euo pipefail

SCRATCH=/n/scratch/users/c/cos689/xtea2_testdata
WT="${1:?usage: sbatch td_manifest_wgs_gate.sbatch <worktree>}"
OUT="${XTEA2_OUT:-$SCRATCH/td_manifest_wgs_${SLURM_JOB_ID:-manual}}"
SRC=$SCRATCH/stage2c_family_assign_49502115/OFF/HG002.wgs
BAM=/n/data1/hms/dbmi/park/corinne/claude_sandbox/test_data/HG002/alignment/HG002.novaseq.pcr-free.30x.dedup.grch38.bam
REF=/n/data1/hms/dbmi/park-smaht_dac/ref/GRCh38_no_alt/hg38_no_alt.fa
EXTND=954          # the value the real run used (grepped from its [td-gather] command line)
SLCT=candidate_list_from_disc.txt.clip_sites_raw_disc.txt.slct
CORES="${SLURM_CPUS_PER_TASK:-12}"
SCAN="$WT/src/csrc/xtea_scan"

[ -x "$SCAN" ]  || { echo "FATAL: $SCAN not built"; exit 1; }
[ -f "$BAM" ]   || { echo "FATAL: $BAM missing"; exit 1; }
for fam in L1 SVA; do
    [ -f "$SRC/$fam/$SLCT" ] || { echo "FATAL: $SRC/$fam/$SLCT missing"; exit 1; }
done
mkdir -p "$OUT"

echo "=============================================================="
echo "--td-manifest WGS GATE"
echo "  job    : ${SLURM_JOB_ID:-<none>} on $(hostname)  cores=$CORES"
echo "  commit : $(git -C "$WT" rev-parse --short HEAD) $(git -C "$WT" log -1 --format=%s | cut -c1-60)"
echo "  sites  : L1 $(wc -l < "$SRC/L1/$SLCT")  SVA $(wc -l < "$SRC/SVA/$SLCT")  extnd=$EXTND"
echo "  out    : $OUT"
echo "  started: $(date '+%F %T')"
echo "=============================================================="

solo() {   # $1 = family
    # SEPARATE STATEMENTS, and this is not style. `local t=$1 d="$OUT/solo_$t"` expands the whole
    # line BEFORE `local` assigns anything, so $t there is the CALLER's t, not $1 — the same trap
    # refactor_gate_chr1.sbatch documents from job 49452897. Job 49516776 walked into it: the
    # precondition loop above had left a global t=SVA, so `solo L1` wrote into solo_SVA and
    # `solo SVA` then rm -rf'd L1's 242 MB of output before writing its own. The gate reported
    # "solo L1 is empty" and read as a FAIL of the code under test, which was fine.
    local t=$1
    local d="$OUT/solo_$t"
    local t0
    rm -rf "$d"; mkdir -p "$d"
    echo; echo "=== solo $t ==="
    t0=$(date +%s)
    "$SCAN" --bam "$BAM" --reference "$REF" \
        --td-sites "$SRC/$t/$SLCT" --td-extnd "$EXTND" \
        --td-clip-fastq "$d/clip.fq" --td-disc-fa "$d/disc.fa" \
        --td-disc-pos "$d/all_disc_pos.disc_pos" --threads "$CORES"
    echo "  wall $(( $(date +%s) - t0 ))s"
    # Checked HERE, not at comparison time: an arm that wrote nothing — or wrote somewhere else —
    # must fail at the point of the mistake, where the path is still in hand.
    local f
    for f in clip.fq disc.fa all_disc_pos.disc_pos; do
        [ -s "$d/$f" ] || { echo "FATAL: solo $t produced no $d/$f"; exit 1; }
    done
}

solo L1
solo SVA

echo; echo "=== shared: both families, one manifest, one decode ==="
D="$OUT/shared"; rm -rf "$D"; mkdir -p "$D"
{
    echo "# tag sites extnd clip_fastq disc_fa disc_pos"
    for fam in L1 SVA; do
        printf '%s\t%s\t%s\t%s\t%s\t%s\n' "$fam" "$SRC/$fam/$SLCT" "$EXTND" \
            "$D/$fam.clip.fq" "$D/$fam.disc.fa" "$D/$fam.all_disc_pos.disc_pos"
    done
} > "$D/manifest.txt"
cat "$D/manifest.txt"
T0=$(date +%s)
"$SCAN" --bam "$BAM" --reference "$REF" --td-manifest "$D/manifest.txt" --threads "$CORES"
echo "  wall $(( $(date +%s) - T0 ))s"

echo
echo "##############################################################"
echo "# --td-manifest WGS GATE RESULT"
echo "##############################################################"
RC=0
N_CMP=0
for fam in L1 SVA; do
    for pair in "clip.fq:$fam.clip.fq" "disc.fa:$fam.disc.fa" \
                "all_disc_pos.disc_pos:$fam.all_disc_pos.disc_pos"; do
        a="$OUT/solo_$fam/${pair%%:*}"
        b="$D/${pair##*:}"
        sz=$(stat -c %s "$a" 2>/dev/null || echo 0)
        # An empty product would make cmp pass vacuously; that is a FAIL, not a pass.
        if [ "$sz" -eq 0 ]; then
            echo "  FATAL: solo $fam ${pair%%:*} is empty — nothing was compared"; RC=1; continue
        fi
        if cmp -s "$a" "$b"; then
            echo "  identical  $fam ${pair%%:*}  ($sz bytes)"
            N_CMP=$((N_CMP + 1))
        else
            echo "  DIFFERS    $fam ${pair%%:*}  (solo $sz vs shared $(stat -c %s "$b" 2>/dev/null || echo 0))"
            cmp "$a" "$b" | head -3 || true
            RC=1
        fi
    done
done

# Leftover part files mean a stripe died without being noticed.
if ls "$D"/*.part* >/dev/null 2>&1; then
    echo "  FATAL: part files were left behind:"; ls "$D"/*.part* | head; RC=1
fi

echo
echo "  products compared: $N_CMP (expected 6)"
[ "$N_CMP" -eq 6 ] || { echo "  FATAL: not every product was compared"; RC=1; }

echo
if [ "$RC" -eq 0 ]; then
    echo "TD-MANIFEST WGS GATE: PASS -- one decode serves both families byte for byte."
else
    echo "TD-MANIFEST WGS GATE: FAIL"
fi
echo "finished: $(date '+%F %T')"
exit $RC
