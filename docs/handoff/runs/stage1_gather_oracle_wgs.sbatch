#!/bin/bash
#SBATCH -c 12
#SBATCH -t 3:00:00
#SBATCH --mem=64G
#SBATCH -p priopark
#SBATCH --account=park_contrib
#SBATCH -o /n/scratch/users/c/cos689/xtea2_testdata/stage1_oracle_wgs_%j.out
#SBATCH --mail-type=NONE

# =============================================================================
# STAGE 1 -- the gather artifacts, byte-compared at WHOLE-GENOME scale.
#
#     sbatch docs/handoff/runs/stage1_gather_oracle_wgs.sbatch
#
# WHY, given that stage1_gate_chr1.sbatch already passes. The chr1 gate compares CALLS,
# which is the plan's gate, but chr1 cannot exercise the two things that only exist at
# genome scale -- and one of them is an ordering rule:
#
#   * MULTI-CONTIG emission order. The clip FASTQ follows the SITE FILE's contig order;
#     the disc FASTA follows the BAM HEADER's. Those are different orders, and on a
#     single-contig slice they are indistinguishable -- a port that used the wrong one for
#     either file passes every chr1 comparison. The real site list here spans 24 contigs.
#   * The ~50K-site profile: 49,042 sites (L1) / 53,569 (SVA) rather than ~3.7K, which is
#     also the only check on the collector's memory behaviour.
#
# It compares the three FILES the transduction step consumes rather than the calls, so a
# failure names the artifact and the byte offset instead of a moved VCF record.
#
# Inputs are the staged NA12878 WGS run's own outputs, read-only -- including its real
# `.slct` site list, so the comparison runs on exactly the list the step would build
# (which matters: re_slct_with_clip_raw_disc_sites is what keeps a duplicate position out
# of the list, and a duplicate CRASHES the legacy gather).
#
# extnd = 753 = int(mean_is + 3*std) = int(444.615 + 3*103.053) for this BAM, from the
# same run's stdout. It is x_TEA_main's iextnd, and it must match or the windows are not
# the ones production uses.
# =============================================================================
set -euo pipefail

SCRATCH=/n/scratch/users/c/cos689/xtea2_testdata
WT="${XTEA2_WORKTREE:-$SCRATCH/stage1_wt}"
OUT="${XTEA2_OUT:-$SCRATCH/stage1_oracle_wgs_${SLURM_JOB_ID:-manual}}"
INPUT_ROOT=$SCRATCH/port_wgs_out/NA12878.wgs
BAM=$SCRATCH/NA12878.wgs.bam
REF=/n/data1/hms/dbmi/park-smaht_dac/ref/GRCh38_no_alt/hg38_no_alt.fa
EXTND=753
THREADS="${SLURM_CPUS_PER_TASK:-12}"

source /n/app/conda/miniforge3/24.11.3-0/etc/profile.d/conda.sh
conda activate xtea2
export PATH=/home/cos689/software/bin:$PATH
test "$PYTHONNOUSERSITE" = 1

[ -d "$WT" ]                    || { echo "FATAL: worktree $WT missing"; exit 1; }
[ -x "$WT/src/csrc/xtea_scan" ] || { echo "FATAL: $WT/src/csrc/xtea_scan not built"; exit 1; }
[ -f "$BAM" ]                   || { echo "FATAL: $BAM missing"; exit 1; }
"$WT/src/csrc/xtea_scan" --help 2>&1 | grep -q -- "--td-sites" || {
    echo "FATAL: $WT/src/csrc/xtea_scan predates --td-sites; rebuild it"; exit 1; }

mkdir -p "$OUT"
cd "$WT"

echo "=============================================================="
echo "STAGE 1 GATHER ORACLE (WGS)  --  NA12878, extnd=$EXTND, threads=$THREADS"
echo "  job    : ${SLURM_JOB_ID:-<none>} on $(hostname)"
echo "  commit : $(git rev-parse --short HEAD)"
echo "  started: $(date '+%F %T')"
echo "=============================================================="

RC=0
for family in L1 SVA; do
    sites="$INPUT_ROOT/$family/candidate_list_from_disc.txt.clip_sites_raw_disc.txt.slct"
    [ -f "$sites" ] || { echo "FATAL: $sites missing"; exit 1; }
    echo
    echo "############### $family ###############"
    set +e
    python docs/handoff/runs/td_gather_oracle.py \
        --bam "$BAM" --ref "$REF" --scan-bin "$WT/src/csrc/xtea_scan" \
        --out "$OUT/$family" --extnd "$EXTND" --threads "$THREADS" --sites "$sites"
    rc=$?
    set -e
    [ "$rc" -eq 0 ] || { echo "  $family: MISMATCH (exit $rc)"; RC=1; }

    # Determinism: a second native pass must reproduce its own output byte for byte. The
    # collector stripes regions across threads, so this is the check that the stripe
    # boundaries do not reach the file.
    echo "--- $family: native determinism repeat ---"
    python docs/handoff/runs/td_gather_oracle.py \
        --bam "$BAM" --ref "$REF" --scan-bin "$WT/src/csrc/xtea_scan" \
        --out "$OUT/${family}_rep2" --extnd "$EXTND" --threads "$THREADS" --sites "$sites" \
        --skip-legacy
    for f in clip.fq disc.fa all_disc_pos.disc_pos; do
        if cmp -s "$OUT/$family/cpp/$f" "$OUT/${family}_rep2/cpp/$f"; then
            echo "  IDENTICAL  $f (two native runs)"
        else
            echo "  DIFFERS    $f (two native runs) -- the gather is NOT deterministic"; RC=1
        fi
    done
    # These are ~1 GB each; keep the comparison, not the bytes.
    rm -rf "$OUT/${family}_rep2" "$OUT/$family/py" "$OUT/$family/cpp"
done

echo
if [ "$RC" -eq 0 ]; then echo "WGS GATHER ORACLE: PASS"; else echo "WGS GATHER ORACLE: FAIL"; fi
echo "finished $(date '+%F %T')"
exit "$RC"
