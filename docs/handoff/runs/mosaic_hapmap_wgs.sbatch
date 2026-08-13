#!/bin/bash
#SBATCH -c 36
#SBATCH -t 96:00:00
#SBATCH --mem=256G
#SBATCH -p park
#SBATCH --account=park
#SBATCH -o /n/scratch/users/c/cos689/xtea2_testdata/mosaic_wgs_%x_%j.out
#SBATCH --mail-type=NONE

# =============================================================================
# Full-WGS mosaic on ONE HapMap-mixture sequencing arm.
#
#     sbatch -J bcm docs/handoff/runs/mosaic_hapmap_wgs.sbatch bcm
#     ... one job per arm: bcm broad nygc uwsc washu
#
# THE EXPERIMENT: the five GCC arms are not redundant replicates -- they are a DEPTH
# LADDER over the same contrived mixture:
#
#     broad ~176x   nygc ~313x   uwsc ~350x   bcm ~466x   washu ~513x
#
# Crossed with the truth set's VAF ladder (0.25%-8.25%, 2,996 tier1 somatic loci), that
# gives sensitivity as a function of BOTH allele fraction and depth -- which is the
# question a mosaic caller actually gets asked ("how deep must I sequence to find a 2%
# insertion"), and is a far stronger result than one sensitivity number.
#
# The all-17-pooled arm (~1819x, 5982 GB) is deliberately NOT here. It is the only
# configuration that could reach the 0.25% stratum, and also the single most expensive
# run; it is worth doing only once the per-arm runs show the mode works.
#
# CAVEAT ON THE LADDER'S TOP RUNG: washu is tissue_code H001; the other four arms are
# A001. That difference is unexplained and comes from the sample manifest, not from
# anything here. If H001 is a different aliquot then washu is not a clean point on a
# depth curve -- treat its number as provisional until that is settled.
#
# Resources are the user's explicit call (2026-08-12): -p park, 36 cores, 192G, 48h,
# submitted WITHOUT waiting for the chr1 run to validate the mode. park nodes are
# 96-core/500G with unlimited CPUs-per-node, so all three fit.
#
# MEMORY RAISED 192G -> 256G (user's call, 2026-08-13). 192G was not a near miss: two of
# the five arms were OOM-killed at essentially the same place -- uwsc (50435715) at 199.4G
# after 4:51, nygc (50435714) at 198.4G after 7:08, both during/just after the scan pass,
# and nygc is the SINGLE-BAM arm, so this is not a per-BAM effect. Both landed ~4% over the
# limit rather than wildly over, which is why 256G (33% headroom) is the bump rather than
# something larger. The observed high-water marks on the arms that survived to downstream
# were 71-101G, so the scan pass is the peak, not the type loop.
#
# WALL TIME RAISED 48h -> 96h (user's call, 2026-08-13), after the first submission
# (jobs 50434763-67) showed bcm scanning at ~3% of the genome per half hour -- a ~15h
# scan pass before any downstream step, so 48h left no useful margin on the deepest
# arms. park's MaxTime is 30 days, so 96h is well inside the partition limit. A running
# job's TimeLimit cannot be raised by the job owner (scontrol update -> "Access/permission
# denied"), so this change only takes effect via resubmission.
# =============================================================================
set -euo pipefail

ARM="${1:?usage: sbatch -J <arm> mosaic_hapmap_wgs.sbatch <bcm|broad|nygc|uwsc|washu>}"

REPO=/n/data1/hms/dbmi/park/corinne/claude_sandbox/xTea
MANIFEST="$REPO/hapmapmix_mosaictruthtest.tsv"
REF=/n/data1/hms/dbmi/park-smaht_dac/ref/GRCh38_no_alt/hg38_no_alt.fa
REPLIB=/home/cos689/xTea/rep_lib_annotation
GFF3=/home/cos689/data1/corinne/ref/gencode.v45.basic.annotation.gff3
OUT=${OUT:-/n/scratch/users/c/cos689/xtea2_testdata/mosaic_hapmap_wgs_$ARM}
SAMPLE="hapmap_${ARM}_wgs"

source /n/app/conda/miniforge3/24.11.3-0/etc/profile.d/conda.sh
conda activate xtea2
export PATH=/home/cos689/software/bin:$PATH
# EXPORTED, not passed as a make variable: the Makefile reads it through $(shell), which
# does not inherit make-level assignments (this silently half-worked in run 50427531).
export PKG_CONFIG_PATH=/home/cos689/software/lib/pkgconfig

echo "job ${SLURM_JOB_ID:-<none>} arm=$ARM on $(hostname), ${SLURM_CPUS_PER_TASK:-36} cpus, $(date '+%F %T')"
echo "commit: $(git -C $REPO rev-parse --short HEAD)"

make -C "$REPO/src/csrc"

# Paths come from the manifest, never retyped -- so an arm cannot silently run on the
# wrong files, and adding an arm to the TSV is all it takes to make it runnable.
mapfile -t BAMS < <(awk -F'\t' -v a="$ARM" 'NR>1 && $2==a {print $4}' "$MANIFEST")
[ ${#BAMS[@]} -gt 0 ] || { echo "FATAL: no BAMs for arm '$ARM' in $MANIFEST"; exit 1; }
for b in "${BAMS[@]}"; do
    [ -r "$b" ]      || { echo "FATAL: unreadable $b"; exit 1; }
    [ -r "$b.bai" ]  || { echo "FATAL: missing index $b.bai"; exit 1; }
done
GB=$(awk -F'\t' -v a="$ARM" 'NR>1 && $2==a {s+=$6} END {printf "%.0f", s}' "$MANIFEST")
echo "=== arm $ARM: ${#BAMS[@]} BAM(s), ${GB} GB ==="; printf '  %s\n' "${BAMS[@]}"

# Every arm is a single sequencing center, and each center registered ONE sample id, so
# no --ignore-rg is needed. Asserted rather than assumed: if an arm ever spans two SM
# tags the run must fail here, not silently pool two samples for 48 hours.
NSM=$(for b in "${BAMS[@]}"; do
          samtools view -H "$b" | awk '/^@RG/{for(i=1;i<=NF;i++) if($i ~ /^SM:/) print $i}'
      done | sort -u | wc -l)
[ "$NSM" -eq 1 ] || { echo "FATAL: arm $ARM spans $NSM SM tags; --ignore-rg would be required"; exit 1; }
echo "  single SM across the arm -> no --ignore-rg"

mkdir -p "$OUT"

ARGS=(
    --mosaic
    -i "${BAMS[*]}"
    -s "$SAMPLE" -o "$OUT"
    --types L1,Alu,SVA
    # --cns-emit-scan: the VALIDATED collection path. baseline0_hg002.sbatch and the
    # bit-reproducibility experiment (repro_master_vs_xtea2.sbatch) both pass it, so every
    # recent germline result was produced with it -- omitting it, as the first mosaic runs
    # did, silently selects a slower path nothing else uses. It matters far more here than
    # for germline: without it cns_filter does an indexed fetch and a temp file PER
    # CANDIDATE, and mosaic's 2/2 cutoffs admit 21,731 L1 and 60,052 Alu candidates on chr1
    # alone. Run 50427774 spent 93 minutes in L1:cns_filter and had not finished.
    --cns-emit-scan
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
    # --keep-intermediates (user's call, 2026-08-13): without it cleanup_successful_run
    # keeps ONLY the merged VCF and deletes the whole sample tree. That cost us on chr1 --
    # the per-family candidate_disc_filtered_cns_post_filtering_mosaic.txt files and
    # step_timings.tsv were deleted the moment the run succeeded, and they are exactly what
    # the truth-set matching and the sibling-timing comparison read. This round has to be
    # measurable: it is the first with sibling sharded by position range, so step_timings
    # is the evidence the fix worked.
    #
    # WATCH THE INODE COUNT, not the byte count. chr1 alone peaked at 42,304 files for ~19
    # GB -- this workload is small-file-heavy, and the scratch quota is 2.5M inodes against
    # 25 TiB. Five whole-genome arms retaining everything is the first configuration that
    # could plausibly reach that limit; bytes are not the binding constraint, files are.
    --keep-intermediates
    # --resume: reuse a completed scan_pass / cns_emit from a previous attempt on this same
    # OUT dir. These are the two phases that cost hours (broad 1.8h + 0.5h; bcm 3.0h + 2.0h;
    # uwsc was still in scan_pass at 6h), and every failure this fleet has actually hit
    # landed AFTER them: the 192G OOMs during/just after the scan, the `| head` pipefail bug
    # at the very end with the VCF already written, and code fixes that force a rerun.
    #
    # Safe by construction, not by trust: the marker under $OUT/<sample>/.resume is written
    # only after a phase RETURNS, so an OOM mid-scan leaves no marker and is redone; and the
    # recorded outputs are re-checked for existence and size before anything is reused, so a
    # marker that outlived its files is rejected. The fingerprint includes the xtea_scan
    # binary's size+mtime, so `make -C src/csrc` above invalidates the scan automatically.
    #
    # Pairs with --keep-intermediates rather than fighting it: cleanup_successful_run removes
    # .resume along with everything else on a SUCCESSFUL run, which is what you want (nothing
    # to resume). On a failed run cleanup never executes, which is the case this is for.
    --resume
    --cores "${SLURM_CPUS_PER_TASK:-36}"
)

# Preflight: --dry-run resolves every path and builds every step's argv without opening a
# BAM, so an argv mistake costs seconds rather than a 48-hour allocation.
echo "=== preflight (--dry-run) ==="
python3 "$REPO/xtea/run_xtea.py" "${ARGS[@]}" --dry-run > "$OUT/preflight.log" 2>&1 || {
    echo "FATAL: preflight failed. Last lines:"; tail -20 "$OUT/preflight.log"; exit 1; }
echo "  preflight OK: $(grep -cE '^# ' "$OUT/preflight.log") steps resolved"

# --- disk sampler ------------------------------------------------------------------
# How much scratch does one mosaic run need? Nobody knows: mosaic drops the clip/disc
# cutoffs to 2/2, which admits far more candidate sites than any germline run, and the
# intermediates scale with sites rather than with input size. Five of these run
# concurrently, so the answer decides whether they fit on scratch at all.
#
# Sampled over time rather than measured at the end, because the two differ: the PEAK is
# what has to fit, and intermediates are removed as steps complete, so a final-size
# measurement can be a small fraction of the high-water mark. Per-step timings in
# step_timings.tsv line up with the elapsed column, so a spike can be attributed to the
# step that caused it.
DISK_LOG="$OUT/disk_usage.tsv"
printf 'elapsed_s\ttotal_bytes\ttotal_h\tn_files\n' > "$DISK_LOG"
(
  T0=$SECONDS
  while true; do
      B=$(du -sb "$OUT" 2>/dev/null | cut -f1); B=${B:-0}
      N=$(find "$OUT" -type f 2>/dev/null | wc -l)
      printf '%s\t%s\t%s\t%s\n' "$((SECONDS-T0))" "$B" "$(numfmt --to=iec "$B" 2>/dev/null)" "$N" >> "$DISK_LOG"
      sleep 60
  done
) &
SAMPLER=$!
trap 'kill $SAMPLER 2>/dev/null || true' EXIT

/usr/bin/time -v python3 "$REPO/xtea/run_xtea.py" "${ARGS[@]}" 2>&1

kill $SAMPLER 2>/dev/null || true

echo
echo "=== DISK FOOTPRINT ($ARM, input ${GB} GB) ==="
awk -F'\t' 'NR>1 && $2+0>peak {peak=$2; at=$1; nf=$4}
     END {if (peak) printf "  peak      %.1f GB at t+%ds (%d files)\n", peak/1e9, at, nf}' "$DISK_LOG"
echo "  final     $(du -sh "$OUT" 2>/dev/null | cut -f1)"
echo "  samples   $(( $(wc -l < "$DISK_LOG") - 1 )) (60 s interval) -> $DISK_LOG"
echo "  --- largest surviving files ---"
# awk truncates, NOT `| head -12`. head exits as soon as it has its 12 lines, the upstream
# find/sort then take SIGPIPE (141), and `set -o pipefail` turns that into a script abort --
# AFTER the pipeline has already printed exactly what it was supposed to. That is what killed
# broad (job 50447194, FAILED at 05:54:32) with the run itself complete, "run complete" logged
# and the 1230-call VCF written: the science was fine, only this reporting tail died, and it
# took the call-set counts and the benchmark scoring below down with it. awk reads its input
# to EOF, so nothing upstream is ever signalled.
find "$OUT" -type f -printf '%s\t%p\n' 2>/dev/null \
  | sort -rn | awk -F'\t' 'NR<=12 {printf "  %8.2f GB  %s\n", $1/1e9, $2}'
echo "  --- by top-level subdir ---"
# Same truncation fix, plus `|| true`: if the run died before creating any per-type subdir the
# glob matches nothing and du exits non-zero, which pipefail would again treat as fatal. A
# disk-usage summary must never be able to fail the job whose output it is summarising.
{ du -sh "$OUT"/*/ 2>/dev/null | sort -rh | awk 'NR<=8 {print "  " $0}'; } || true

echo
echo "=== mosaic call sets ($ARM) ==="
for t in L1 Alu SVA; do
    f="$OUT/$SAMPLE/$t/candidate_disc_filtered_cns_post_filtering_mosaic.txt"
    printf '  %-4s %s\n' "$t" "$([ -f "$f" ] && wc -l < "$f" || echo MISSING)"
done
VCF="$OUT/$SAMPLE/$SAMPLE.vcf"
echo "=== score (VAF-stratified) ==="
if [ -f "$VCF" ]; then
    python3 "$REPO/benchmark/benchmark_hapmap_mosaic.py" --calls "$VCF" || true
else
    echo "  no VCF at $VCF"
fi
echo "finished $(date '+%F %T')"
