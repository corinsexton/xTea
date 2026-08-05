#!/bin/bash
#SBATCH -c 24
#SBATCH -t 12:00:00
#SBATCH --mem=64G
#SBATCH -p priopark
#SBATCH --account=park_contrib
#SBATCH -o /n/scratch/users/c/cos689/xtea2_testdata/master_wgs_hg002_%j.out
#SBATCH --mail-type=NONE

# =============================================================================
# MASTER xTea, whole genome, on HG002 -- the Stage 3 transduction gate's FIRST POOL ARM.
#
# WHY THIS RUN EXISTS, and it is not the timing comparator.
# `master_serial_wgs.sbatch` runs the same code on NA12878 to measure WALL CLOCK.
# This run exists to produce a CALL SET: the transduction hypotheses that
# `HANDOFF_2026-08-04_stage3_td_gate_spec.md` §2 puts in the gate's pool, on the
# sample the truth set is for. Corinne's direction, 2026-08-04.
#
# It is a separate file rather than a flag on that one for two reasons that both
# bite: that script opens with `rm -rf "$OUT"`, and its OUT default is the
# NA12878 tree the gate spec's §2.1 table and the decoy probe both cite -- so a
# repointed invocation that forgot XTEA2_OUT would delete the evidence. And its
# report is entirely about timing (three-way wall-clock comparison, per-step
# totals, golden call-count deltas for NA12878) which says nothing here.
#
# WHY THE GOLDEN'S CHAIN SCRIPTS AND NOT A REGENERATED PIPELINE. Spec §7 step 1
# requires "master's OWN pipeline end to end" -- plan §10 records that master's
# transduction figures come in two non-interchangeable flavours, master's own
# pipeline on its own inputs versus master's code on the port's inputs, and only
# the former samples master's call set. The golden's per-type
# `run_xTEA_pipeline.sh` files ARE the output of master's own
# `gnrt_pipeline_local_v38.py`, so transforming them in place by sed keeps
# exactly one thing different from a run master generated itself: the BAM.
# Retyping ~30 invocations would be a source of silent divergence.
#
# WHAT DIFFERS FROM THE NA12878 SERIAL RUN, EXHAUSTIVELY:
#   1. BAM: HG002.novaseq.pcr-free.30x.dedup.grch38.bam instead of NA12878.wgs.bam
#   2. OUT tree and PREFIX/pub_clip paths follow from (1)
#   3. the report is the transduction pool, not wall clock
# Everything else -- cutoffs as written, annotations, reference, the py3.8 xTea
# env with deep-forest, step order, --sva placement -- is the golden's, verbatim.
#
# THE CUTOFFS ARE NOT A PROBLEM, AND THIS IS THE ONE THING WORTH CHECKING BEFORE
# DOUBTING THIS RUN. The golden's command lines carry `--lc 4 --rc 4 --cr 2`
# and `--nfclip 4 --nfdisc 5`, derived from NA12878's coverage. They do not
# transfer to HG002 -- and they do not have to: `x_TEA_main` sets
# `b_automatic = True` unless `--user` is passed (nothing passes it), which
# OVERRIDES the command-line cutoffs with `Parameters.get_par_by_cov(coverage)`
# computed from the BAM in hand. That mechanism is measured, not assumed: it is
# what `probe_decoy_fix_delta.py` had to model to reproduce 94 of a run's 105
# picked sources, and HG002 WGS at 31.816x resolves to `m_cov_par[30] = (3,4,1)`
# -- n_clip_cutoff 3, n_disc_cutoff 4 (commit e38d0f5). So master self-adjusts
# to HG002 and the golden's numbers on the command line are inert.
#
# REFERENCE SHAPE. HG002's BAM is hs38DH (3,366 @SQ) against a 195-contig
# no_alt reference, which is the SAME configuration NA12878.wgs.bam has already
# run in -- known-good, not a new risk (test_data/HG002/PROVENANCE.md §2). It is
# also the configuration every xTea2 HG002 gate has used, so the two pool arms
# see the same reads placed the same way.
#
# WHAT THIS RUN IS NOT. Not truth -- master enumerates hypotheses, the reads
# return the verdict (spec §6.4). Not a timing figure -- it shares a node class
# with nothing in particular and its call set is the point. Not comparable to
# the paper (hs38DH vs no_alt, plan §7.1).
#
# USAGE  sbatch docs/handoff/runs/master_wgs_hg002.sbatch
#        MASTER_OUT=... NTHREADS=8 sbatch ...      (both overridable)
# =============================================================================
set -euo pipefail

SCRATCH=/n/scratch/users/c/cos689/xtea2_testdata
WT="${MASTER_WORKTREE:-$SCRATCH/master_cmp_wt}"
OUT="${MASTER_OUT:-$SCRATCH/master_wgs_hg002}"
NTHREADS="${NTHREADS:-24}"
# The templates. NA12878's golden is the source of the command lines ONLY; no
# NA12878 output is read and none is written.
GOLDEN="$SCRATCH/master_golden_wgs/NA12878"
TESTDATA=/n/data1/hms/dbmi/park/corinne/claude_sandbox/test_data/HG002
BAM="$TESTDATA/alignment/HG002.novaseq.pcr-free.30x.dedup.grch38.bam"
SAMPLE=HG002
# master's gVCF step names its output from the BAM's basename, so this is
# derived, not chosen -- same string every HG002 gate script already uses.
VCF_PREFIX=HG002.novaseq.pcr-free.30x.dedup.grch38
REPLIB=/home/cos689/xTea/rep_lib_annotation
# This repo, for the probe -- NOT the master worktree. probe_td_source_classes
# imports xTea2's `transduction_evidence`, which is the point of it: the
# classification it applies is the shipping one.
REPO="${XTEA2_REPO:-/n/data1/hms/dbmi/park/corinne/claude_sandbox/xTea}"

# The golden's env, verbatim: the py3.8 "magical" xTea env (numpy 1.19.5 +
# deepforest). NOT the xtea2 env -- master's gntp_classify needs DF21.
export PATH=/n/data1/hms/dbmi/park/corinne/.conda/envs/xTea/bin:/home/cos689/software/bin:$PATH
export PYTHONNOUSERSITE=1     # keep ~/.local off sys.path (cost a SIGSEGV once: job 48777874)

[ -d "$WT/xtea" ] || { echo "FATAL: master worktree $WT/xtea missing"; exit 1; }
[ -d "$WT/xtea/genotyping/DF21_model_1_2" ] || { echo "FATAL: DF21 model missing in $WT"; exit 1; }
[ -f "$BAM" ] || { echo "FATAL: $BAM missing"; exit 1; }
[ -f "$BAM.bai" ] || [ -f "${BAM%.bam}.bai" ] || { echo "FATAL: no index beside $BAM"; exit 1; }
[ -d "$GOLDEN" ] || { echo "FATAL: template dir $GOLDEN missing"; exit 1; }
# Refuse to write anywhere that could be another run's evidence.
case "$OUT" in
  *master_serial_out*|*master_golden_wgs*|*baseline0*)
    echo "FATAL: OUT=$OUT names an existing frozen tree; refusing"; exit 1 ;;
esac

rm -rf "$OUT"; mkdir -p "$OUT/$SAMPLE/pub_clip"

echo "=============================================================="
echo "MASTER xTea -- whole genome, HG002 -- Stage 3 td-gate pool arm 1"
echo "  job        : ${SLURM_JOB_ID:-<none>}  on $(hostname)"
echo "  master code: $WT ($(cd "$WT" && git rev-parse --short HEAD) $(cd "$WT" && git log -1 --format=%s | cut -c1-50))"
echo "  python     : $(python --version 2>&1)  [$(which python)]"
echo "  bam        : $BAM"
echo "  templates  : $GOLDEN/{L1,Alu,SVA}/run_xTEA_pipeline.sh  (command lines only)"
echo "  out        : $OUT/$SAMPLE"
echo "  -n per step: $NTHREADS"
echo "  cutoffs    : as written in the templates, then OVERRIDDEN by"
echo "               get_par_by_cov(HG002 coverage) -- b_automatic, see header"
echo "  started    : $(date '+%F %T')"
echo "=============================================================="

TIMINGS="$OUT/step_timings.tsv"
printf 'type\tstep_idx\tstep_label\twall_s\n' > "$TIMINGS"

label_of() {   # derive a step name from the invocation's mode flag
  case "$1" in
    *" -C "*)              echo clip ;;
    *" -D "*)              echo disc ;;
    *" -N "*)              echo cns_filter ;;
    *--transduction*)      echo transduction ;;
    *--sibling*)           echo sibling ;;
    *--postF*high_confident*) echo postF_hc ;;
    *--postF*)             echo postF_all ;;
    *--gene*)              echo gene ;;
    *--gntp_classify*)     echo genotype ;;
    *--gVCF*)              echo gvcf ;;
    *)                     echo unknown ;;
  esac
}

T_RUN_START=$(date +%s)
declare -A TYPE_WALL

for t in L1 Alu SVA; do
  SRC="$GOLDEN/$t/run_xTEA_pipeline.sh"
  [ -f "$SRC" ] || { echo "FATAL: template script $SRC missing"; exit 1; }
  P="$OUT/$SAMPLE/$t"
  mkdir -p "$P/tmp/clip" "$P/tmp/cns" "$P/tmp/transduction"
  printf '%s\tillumina\n' "$BAM" > "$P/bam_list.txt"   # the ONLY place the BAM enters

  GEN="$OUT/$t.chain.sh"
  sed -e "s|^PREFIX=.*|PREFIX=$P/|" \
      -e "s|^XTEA_PATH=.*|XTEA_PATH=$WT/xtea/|" \
      -e "s|-n 8 |-n $NTHREADS |g" \
      -e "s|-n 8$|-n $NTHREADS|g" \
      -e "s|/home/cos689/scratch/xtea2_testdata/master_golden_wgs/NA12878/pub_clip/|$OUT/$SAMPLE/pub_clip/|g" \
      -e "/^#SBATCH/d" -e "/^export PATH=/d" -e "/^#!/d" \
      "$SRC" > "$GEN"

  # Fail loudly rather than silently reusing NA12878 paths: after the sed, no
  # invocation may still mention the template sample.
  if grep -q 'NA12878' "$GEN"; then
    echo "FATAL: $GEN still references NA12878 after the transform:"
    grep -n 'NA12878' "$GEN"
    exit 1
  fi

  # source the variable block (paths/cutoffs), then run the python lines one at a
  # time so each step is timed. eval is required: the invocations reference ${VARS}
  # defined in that same block.
  # shellcheck disable=SC1090
  source <(grep -v '^python ' "$GEN" | grep -v '^#')

  mapfile -t PYLINES < <(grep '^python ' "$GEN")
  echo
  echo "=============== $t : ${#PYLINES[@]} steps, -n $NTHREADS ==============="
  T_TYPE_START=$(date +%s)
  idx=0
  for line in "${PYLINES[@]}"; do
    idx=$((idx + 1))
    lbl=$(label_of "$line")
    echo "--- [$t $idx/${#PYLINES[@]}] $lbl  $(date '+%T')"
    t0=$(date +%s)
    eval "$line"
    t1=$(date +%s)
    printf '%s\t%d\t%s\t%d\n' "$t" "$idx" "$lbl" "$((t1 - t0))" >> "$TIMINGS"
    printf '    %s %s done in %d s (%.1f min)\n' "$t" "$lbl" "$((t1 - t0))" \
      "$(echo "scale=2; ($t1-$t0)/60" | bc)"
  done
  T_TYPE_END=$(date +%s)
  TYPE_WALL[$t]=$((T_TYPE_END - T_TYPE_START))
  printf '=== %s TOTAL %d s (%02d:%02d:%02d) ===\n' "$t" "${TYPE_WALL[$t]}" \
    $((TYPE_WALL[$t]/3600)) $((TYPE_WALL[$t]%3600/60)) $((TYPE_WALL[$t]%60))
done

T_RUN_END=$(date +%s)
WALL=$((T_RUN_END - T_RUN_START))

# -----------------------------------------------------------------------------
# report. Reporting failures must not mask the run: || true throughout, and the
# artifacts are all on disk for offline re-analysis either way.
# -----------------------------------------------------------------------------
echo
echo "##############################################################"
echo "# 1. THE CALL SET -- this is what the run is for"
echo "##############################################################"
printf 'wall %02d:%02d:%02d (informational; NA12878 serial was 50:18 at job 48825977)\n' \
  $((WALL/3600)) $((WALL%3600/60)) $((WALL%60))
echo
echo "NA12878 is shown only to say whether this run is in the same regime. HG002 is a"
echo "different sample at different coverage with different cutoffs, so a delta here is"
echo "NOT an error signal -- see the header on b_automatic."
printf '%-6s %10s %14s\n' type HG002 'NA12878(ref)'
for pair in "L1:LINE1:203" "Alu:ALU:1262" "SVA:SVA:73"; do
  t="${pair%%:*}"; rest="${pair#*:}"; alt="${rest%%:*}"; ref="${rest##*:}"
  f=$(ls "$OUT/$SAMPLE/$t"/*_${alt}.vcf 2>/dev/null | head -1)
  n=$([ -n "$f" ] && grep -vc '^#' "$f" || echo NA)
  printf '%-6s %10s %14s\n' "$t" "$n" "$ref"
done

echo
echo "##############################################################"
echo "# 2. THE TRANSDUCTION POOL -- spec §2.1 / §7 step 2"
echo "##############################################################"
echo "The number the gate spec calls unknown: HG002's pool size. NA12878's is 115"
echo "transduction-asserting calls (86 L1 + 29 SVA, 0 Alu by construction)."
echo
if [ -f "$REPO/docs/handoff/runs/probe_td_pool.py" ]; then
  python3 "$REPO/docs/handoff/runs/probe_td_pool.py" \
      --run-dir "$OUT/$SAMPLE" --vcf-prefix "$VCF_PREFIX" \
      --pool-tsv "$OUT/td_pool_master_hg002.tsv" || true
else
  echo "  probe_td_pool.py not present in $REPO -- raw TD_SRC counts only, no strata."
  for pair in "L1:LINE1" "Alu:ALU" "SVA:SVA"; do
    t="${pair%%:*}"; alt="${pair#*:}"
    f=$(ls "$OUT/$SAMPLE/$t"/*_${alt}.vcf 2>/dev/null | head -1)
    [ -n "$f" ] || { echo "  $t: no VCF"; continue; }
    awk -F'\t' -v t="$t" '
      /^#/ {next}
      {calls++
       if (match($8, /TD_SRC=[^;]*/)) {
         v = substr($8, RSTART + 7, RLENGTH - 7)
         if (v != "" && v != "not_transduction") td++
       }}
      END {printf "  %-4s calls %5d   assert a transduction %5d\n", t, calls, td}' "$f" || true
  done
fi

echo
echo "##############################################################"
echo "# 3. THE DECOY SIGNATURE -- is HG002's master run SVA-blind too"
echo "##############################################################"
echo "master shares the decoy defect (its NA12878 run: SVA 92 picked sources, 0"
echo "reference copies, against L1 135/84). If HG002's master run shows the same"
echo "signature then spec §5's rule -- score the decoy fix as GROWTH, never against"
echo "master's calls -- is established on the gate's OWN sample, not by transfer."
echo
python3 "$REPO/docs/handoff/runs/probe_td_source_classes.py" \
    --rep-lib-dir "$REPLIB" --genome hg38 --run-dir "$OUT/$SAMPLE" 2>&1 || true

echo
echo "##############################################################"
echo "# 4. TIMINGS (informational) and WHERE THINGS ARE"
echo "##############################################################"
sort -t$'\t' -k4 -rn "$TIMINGS" | grep -v '^type' | head -8 \
  | awk -F'\t' '{printf "  %-4s %-14s %8.2f min\n", $1, $3, $4/60}' || true
echo
du -sh "$OUT" 2>/dev/null || true
echo "  call set   : $OUT/$SAMPLE/{L1,Alu,SVA}/${VCF_PREFIX}_{LINE1,ALU,SVA}.vcf"
echo "  pool tsv   : $OUT/td_pool_master_hg002.tsv"
echo "  sources    : $OUT/$SAMPLE/{L1,SVA}/candidate_disc_filtered_cns2.txt.all_non_sibling_td.txt.unique_trsdct_disc_only"
echo "  memory     : sacct -j ${SLURM_JOB_ID:-<jobid>} --format=JobID%16,Elapsed,MaxRSS%12,AllocCPUS,State"

echo
echo "=============================================================="
echo "finished $(date '+%F %T')"
printf 'total wall %02d:%02d:%02d\n' $((WALL/3600)) $((WALL%3600/60)) $((WALL%60))
echo "=============================================================="
