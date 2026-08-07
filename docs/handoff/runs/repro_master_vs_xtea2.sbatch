#!/bin/bash
#SBATCH -J repro_cmp
#SBATCH --array=1-6
#SBATCH -c 24
#SBATCH -t 12:00:00
#SBATCH --mem=64G
#SBATCH -p priopark
#SBATCH --account=park_contrib
#SBATCH -o /n/scratch/users/c/cos689/xtea2_testdata/repro_cmp_%A_%a.out
#SBATCH --mail-type=NONE
# =============================================================================
# MASTER vs xTea2, THREE REPLICATES EACH -- speed, memory, precision/recall.
# Corinne's direction, 2026-08-07.
#
# WHY THIS EXISTS. Every master-vs-xTea2 figure in this project is n=1 ON THE MASTER SIDE.
# `HANDOFF_2026-08-07_objective_met.md` §1 quotes "1.66-1.82x faster on 51-61% of master's peak
# memory" and says outright that the range's precision is "bounded by a variance nobody has
# measured": three xTea2 WGS runs exist (49792012 / 49879677 / 50038678) and exactly one master
# run does (49763841). Two of those xTea2 runs are the SAME configuration and came out 199 s and
# 3.41 GB apart, so the xTea2 noise floor is known and master's is not.
#
# And the accuracy side was worse than n=1: master's HG002 VCFs have sat on disk since 2026-08-04
# and were never scored against the truth set at all. Scored offline on 2026-08-07 they came to
# 80.6% sensitivity / 87.8% precision against xTea2's 80.4% / 88.3% -- a 4-true-positive,
# 9-false-positive difference over 1642 truth records. THAT DIFFERENCE IS UNINTERPRETABLE AT n=1,
# and not for a generic statistical reason: job 49075238 measured master's own `-n 8` output
# differing from master's own `-n 24` output ON IDENTICAL INPUTS. Master is known to be
# thread-scheduling sensitive, so 4 calls may be entirely its own run-to-run spread. This array
# is what decides that.
#
# THE DESIGN, and the two choices that are not defaults.
#
#   1. THE SEED IS NOT PINNED ON THE xTea2 ARM (Corinne, 2026-08-07). Every prior xTea2 gate
#      exported XTEA_BASIC_INFO_SEED=42 to pin the insert-size draw. Master has no such knob.
#      Pinning one arm and not the other would hold a source of variation fixed on the xTea2 side
#      only, and xTea2 would then look more reproducible for a reason that is not the pipeline.
#      So this harness UNSETS it and asserts it is unset. The cost is deliberate: these three runs
#      are therefore NOT drop-in comparable with 49792012/49879677/50038678, which were pinned.
#
#   2. ALL SIX RUN CONCURRENTLY, as one unthrottled array (Corinne, 2026-08-07). Both arms then
#      see identical cluster conditions and identical contention on the one 44 GB BAM they all
#      stream, so contention lands INSIDE the measured spread rather than beside it -- which is
#      what a noise floor is supposed to include. The cost is that these wall times are
#      under-load, not solo, so they are not directly comparable to the solo n=1 numbers either.
#      Node and CPU model are recorded per task so a node effect stays visible instead of being
#      averaged into the arm difference.
#
# WHAT IS HELD FIXED ACROSS ALL SIX: the BAM, `-c 24`, `--mem=64G`, priopark, hg38 no_alt
# reference, rep_lib, L1+Alu+SVA, and the 12 h limit. Nothing here passes a cutoff: master sets
# `b_automatic` and derives its own from coverage (see master_wgs_hg002.sbatch's header), and
# xTea2 does the same.
#
# WHAT EACH ARM RUNS.
#   master  -- `master_wgs_hg002.sbatch` invoked as a plain script with MASTER_OUT repointed. It
#              is reused rather than reimplemented ON PURPOSE: it drives the golden's own
#              per-type `run_xTEA_pipeline.sh` command lines, which are the output of master's own
#              `gnrt_pipeline_local_v38.py`, so the only thing differing from a run master
#              generated itself is the BAM. Retyping ~30 invocations is exactly the silent
#              divergence that file's header warns about.
#   xTea2   -- `run_xtea.py` at the frozen worktree, WITH NEITHER DECOY FLAG, so both flips come
#              from their 2026-08-05 / 2026-08-07 defaults. Check 0 below asserts both banners
#              name the DEFAULT as their cause, which is the same "armed vs changed anything"
#              split arm 4 of the polyA gate established.
#
# THIS HARNESS DOES NOT SCORE ANYTHING. Scoring six trees needs all six present, so it lives in
# `probe_repro_comparison.py` -- offline, seconds, re-runnable without re-running any of this.
#
# USAGE  sbatch docs/handoff/runs/repro_master_vs_xtea2.sbatch [<xtea2-worktree>]
#        then: python3 docs/handoff/runs/probe_repro_comparison.py --array-job <jobid>
# =============================================================================
set -euo pipefail

SCRATCH=/n/scratch/users/c/cos689/xtea2_testdata
WT="${1:-/n/scratch/users/c/cos689/xtea2_wt/repro_ea08f45}"
TESTDATA=/n/data1/hms/dbmi/park/corinne/claude_sandbox/test_data/HG002
BAM="$TESTDATA/alignment/HG002.novaseq.pcr-free.30x.dedup.grch38.bam"
REF=/n/data1/hms/dbmi/park-smaht_dac/ref/GRCh38_no_alt/hg38_no_alt.fa
REPLIB=/home/cos689/xTea/rep_lib_annotation
CORES="${SLURM_CPUS_PER_TASK:-24}"
IDX="${SLURM_ARRAY_TASK_ID:?this script must run as a job array (--array=1-6)}"

# 1-3 master, 4-6 xTea2. Kept as arithmetic rather than a case block so the replicate count is
# one number to change.
NREPS=3
if [ "$IDX" -le "$NREPS" ]; then ARM=master; REP="$IDX"; else ARM=xtea2; REP=$((IDX - NREPS)); fi
OUT="$SCRATCH/repro_${ARM}_rep${REP}"

echo "=============================================================="
echo "REPRO COMPARISON  --  arm=$ARM  replicate=$REP of $NREPS"
echo "  array      : ${SLURM_ARRAY_JOB_ID:-?}_${IDX}   job ${SLURM_JOB_ID:-?}"
echo "  node       : $(hostname)"
echo "  cpu        : $(grep -m1 'model name' /proc/cpuinfo | cut -d: -f2- | sed 's/^ //')"
echo "  cores      : $CORES     mem: ${SLURM_MEM_PER_NODE:-?} MB"
echo "  out        : $OUT"
echo "  started    : $(date '+%F %T')"
echo "=============================================================="

# The BAM is the one input both arms share; a silent substitution here would invalidate the whole
# comparison, so check it rather than trust the path.
[ -f "$BAM" ]     || { echo "FATAL: $BAM missing"; exit 1; }
[ -f "$BAM.bai" ] || { echo "FATAL: $BAM.bai missing"; exit 1; }
# Refuse to write over any of the evidence this project's frozen gates cite.
case "$OUT" in
  *master_wgs_hg002*|*master_serial_out*|*master_golden_wgs*|*baseline0*|*gate_*)
    echo "FATAL: OUT=$OUT names an existing frozen tree; refusing"; exit 1 ;;
esac

# THE SEED IS UNPINNED ON BOTH ARMS -- design note 1. Asserted, not assumed: an inherited
# XTEA_BASIC_INFO_SEED from a shell or a sourced profile would silently pin the xTea2 arm and
# there would be no trace of it in the output.
unset XTEA_BASIC_INFO_SEED || true
[ -z "${XTEA_BASIC_INFO_SEED:-}" ] || { echo "FATAL: XTEA_BASIC_INFO_SEED is set ($XTEA_BASIC_INFO_SEED)"; exit 1; }
echo "[seed] XTEA_BASIC_INFO_SEED is UNSET on both arms (design note 1)"

t0=$(date +%s)

if [ "$ARM" = master ]; then
    # ---------------------------------------------------------------------------------------
    # MASTER. Reuse the existing harness wholesale; repoint only its output tree. XTEA2_REPO is
    # pointed at the FROZEN worktree rather than the live checkout so its two report probes read
    # code that cannot change mid-run.
    # ---------------------------------------------------------------------------------------
    MASTER_SH="$WT/docs/handoff/runs/master_wgs_hg002.sbatch"
    [ -f "$MASTER_SH" ] || { echo "FATAL: $MASTER_SH missing"; exit 1; }
    MWT="$SCRATCH/master_cmp_wt"
    [ -d "$MWT/xtea" ] || { echo "FATAL: master worktree $MWT missing"; exit 1; }
    echo "  master code: $MWT ($(git -C "$MWT" rev-parse --short HEAD))"
    echo
    MASTER_OUT="$OUT" MASTER_WORKTREE="$MWT" XTEA2_REPO="$WT" NTHREADS="$CORES" \
        bash "$MASTER_SH"
else
    # ---------------------------------------------------------------------------------------
    # xTea2, SHIPPING DEFAULTS. Neither --td-decoy-structural nor --td-decoy-concat is passed.
    # ---------------------------------------------------------------------------------------
    SAMPLE=HG002.wgs
    SCAN_BIN="$WT/src/csrc/xtea_scan"
    [ -d "$WT" ]       || { echo "FATAL: worktree $WT missing"; exit 1; }
    [ -x "$SCAN_BIN" ] || { echo "FATAL: $SCAN_BIN not built -- make -C $WT/src/csrc"; exit 1; }
    [ -z "$(git -C "$WT" status --porcelain)" ] || { echo "FATAL: worktree is dirty"; exit 1; }

    source /n/app/conda/miniforge3/24.11.3-0/etc/profile.d/conda.sh
    conda activate xtea2
    export PATH=/home/cos689/software/bin:$PATH
    # ~/.local off sys.path. Named explicitly rather than `test "$PYTHONNOUSERSITE" = 1`, which
    # under `set -u` dies with "unbound variable" and tells the reader nothing. It cost a SIGSEGV
    # once (job 48777874).
    [ "${PYTHONNOUSERSITE:-}" = 1 ] || { echo "FATAL: PYTHONNOUSERSITE is not 1; the xtea2 env"; \
        echo "did not activate, and ~/.local on sys.path segfaulted job 48777874"; exit 1; }
    # Both decoy policies must come from their DEFAULTS. An inherited env var would arm one
    # without the command line saying so, and the banners could then not be read as evidence.
    unset XTEA_TD_DECOY_CONCAT
    unset XTEA_TD_DECOY_STRUCTURAL

    rm -rf "$OUT"; mkdir -p "$OUT"
    echo "  xtea2 code : $WT ($(git -C "$WT" rev-parse --short HEAD))"
    echo "  scan bin   : $SCAN_BIN"
    echo "  python     : $(python --version 2>&1)  [$(which python)]"
    echo

    cd "$WT"
    python xtea/run_xtea.py \
      -i "$BAM" -s "$SAMPLE" -o "$OUT" \
      --reference "$REF" \
      --rep-lib-dir "$REPLIB" \
      --gene-gff3 /home/cos689/data1/corinne/ref/gencode.v45.basic.annotation.gff3 \
      --blacklist "$REPLIB/blacklist/hg38/centromere.bed" \
      --types L1,Alu,SVA \
      --scan-bin "$SCAN_BIN" \
      --bwa /home/cos689/software/bin/bwa \
      --scan-min-clip 2 \
      --mate-in-rep-families alu,sva \
      --cns-emit-scan \
      --td-gather-scan \
      --parallel-downstream \
      --keep-intermediates \
      -n "$CORES" \
      > "$OUT/driver.log" 2>&1 || {
          rc=$?; echo "FAIL: the run exited $rc"; tail -80 "$OUT/driver.log"; exit $rc; }

    # --- check 0: BOTH decoy policies came from their DEFAULTS -------------------------------
    # Not decoration. This arm is the shipping configuration only if the defaults are what armed
    # it; a resolver returning the wrong branch, or an option parsed and dropped, would leave a
    # run that looks like the default path and is not.
    echo "--- check 0: both decoy policies are ON, and the DEFAULT is what turned them on ---"
    RC=0
    n=$(grep -ac '\[decoy\] the transduction decoy test is STRUCTURAL' "$OUT/driver.log" || true)
    d=$(grep -ac 'STRUCTURAL (DEFAULT since 2026-08-05)' "$OUT/driver.log" || true)
    echo "  structural ON: $n (expect 1)   attributed to the DEFAULT: $d (expect 1)"
    [ "$n" -eq 1 ] && [ "$d" -eq 1 ] || { echo "  FATAL: structural is not on via the default"; RC=1; }
    n=$(grep -ac '\[decoy-concat\] the transduction bwa reference GAINS' "$OUT/driver.log" || true)
    d=$(grep -ac 'lacks it (DEFAULT)' "$OUT/driver.log" || true)
    echo "  concat ON    : $n (expect 1)   attributed to the DEFAULT: $d (expect 1)"
    [ "$n" -eq 1 ] && [ "$d" -eq 1 ] || { echo "  FATAL: concat is not on via the default"; RC=1; }
    echo "  the step's own per-family append counts:"
    grep -a '\[td-decoy-concat\]' "$OUT/driver.log" | sed 's/^/    /' || true
    [ "$RC" -eq 0 ] || exit 1

    echo
    echo "--- call counts ---"
    for pair in L1:LINE1 Alu:ALU SVA:SVA; do
        t="${pair%%:*}"; alt="${pair#*:}"
        f=$(ls "$OUT/$SAMPLE/merged"/*_${alt}.vcf 2>/dev/null | head -1)
        printf '  %-4s %s\n' "$t" "$([ -n "$f" ] && grep -vc '^#' "$f" || echo NA)"
    done
    printf '  merged %s\n' "$(grep -vc '^#' "$OUT/$SAMPLE/$SAMPLE.vcf" 2>/dev/null || echo NA)"
fi

wall=$(( $(date +%s) - t0 ))
echo
echo "=============================================================="
printf 'arm=%s rep=%s WALL %02d:%02d:%02d   node=%s\n' \
    "$ARM" "$REP" $((wall/3600)) $((wall%3600/60)) $((wall%60)) "$(hostname)"
echo "finished $(date '+%F %T')"
# Written as a machine-readable sidecar so the aggregator needs neither this log nor sacct to
# know which tree came from which task -- sacct's retention is not forever and %A_%a logs get
# tidied, but the run tree is the artifact that persists.
printf 'arm\trep\tarray_job\ttask\tjob\tnode\tcores\twall_s\n' > "$OUT/repro_run.tsv"
printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%d\n' \
    "$ARM" "$REP" "${SLURM_ARRAY_JOB_ID:-}" "$IDX" "${SLURM_JOB_ID:-}" \
    "$(hostname)" "$CORES" "$wall" >> "$OUT/repro_run.tsv"
echo "  sidecar: $OUT/repro_run.tsv"
echo "  memory : sacct -j ${SLURM_JOB_ID:-<jobid>} --format=JobID%20,Elapsed,MaxRSS%12,AllocCPUS,State"
echo "=============================================================="
