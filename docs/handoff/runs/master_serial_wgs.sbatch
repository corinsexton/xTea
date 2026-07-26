#!/bin/bash
#SBATCH -c 24
#SBATCH -t 24:00:00
#SBATCH --mem=192G
#SBATCH -p priopark
#SBATCH --account=park_contrib
#SBATCH -o /n/scratch/users/c/cos689/xtea2_testdata/master_serial_wgs_%j.out
#SBATCH --mail-type=NONE

# =============================================================================
# TRUE COMPARATOR: master xTea, whole genome, NON-CONCURRENT, one job, 24 cores.
#
# Why this run exists. The frozen golden WGS run is NOT a fair comparator for
# xTea2's single-job runtime: it was three per-type jobs (48013142/43/44) that the
# scheduler ran CONCURRENTLY on 8 cores each, so its 38:31 "wall clock" is really
# max(37:15, 29:06, 38:31) across a 24-core total allocation. xTea2 is one process
# with one 24-core allocation. Comparing 38:31 to xTea2's 1:18:23 therefore charges
# xTea2 for concurrency master got free from SLURM, not for being slower code.
#
# This run removes that asymmetry: same master code, same input, same total
# allocation, but the three types run one after another in a single job -- the same
# shape xTea2 has. That makes wall clock directly comparable.
#
# The command lines are NOT retyped. They are the golden's own
# run_xTEA_pipeline.sh files, transformed in place by sed (paths, thread count,
# XTEA_PATH). Retyping 30 invocations would have been a source of silent
# divergence, and the whole value of a comparator is that only ONE thing differs.
#
# What differs from the golden, exhaustively:
#   1. the three chains run sequentially in one job instead of as three jobs
#   2. -n 8 -> -n 24, so each step gets the whole allocation (see NTHREADS below)
#   3. XTEA_PATH points at a git worktree of THIS repo's `master` branch, because
#      the golden's XTEA_PATH (/home/cos689/src/xtea/xtea/) no longer exists
#   4. output goes to a fresh dir; --cp pub_clip is per-run, not the golden's
# Everything else -- cutoffs, annotations, reference, the py3.8 xTea env with
# deep-forest, step order, --sva placement -- is byte-identical to the golden.
#
# On NTHREADS. -n 24 measures "master given the whole allocation per step", which is
# the number that is NOT already derivable. The faithful -n 8 serial figure needs no
# run: it is 37:15 + 29:06 + 38:31 = 1:44:52, by summing the golden's own three jobs.
# Set NTHREADS=8 to reproduce that arithmetic empirically instead.
#
# NOTE: this is a TIMING comparator. It writes its own output tree and does not
# touch the goldens. Its call set should reproduce the golden (same code, same
# cutoffs) modulo the documented -n-dependent spill nondeterminism; the report
# checks that as a sanity signal, not as a gate.
# =============================================================================
set -euo pipefail

WT="${MASTER_WORKTREE:-/n/scratch/users/c/cos689/xtea2_testdata/master_cmp_wt}"
OUT="${MASTER_OUT:-/n/scratch/users/c/cos689/xtea2_testdata/master_serial_out}"
NTHREADS="${NTHREADS:-24}"
SCRATCH=/n/scratch/users/c/cos689/xtea2_testdata
GOLDEN="$SCRATCH/master_golden_wgs/NA12878"
BAM="$SCRATCH/NA12878.wgs.bam"

# The golden's env, verbatim: the py3.8 "magical" xTea env (numpy 1.19.5 +
# deepforest). NOT the xtea2 env -- master's gntp_classify needs DF21, which cannot
# even be installed on py3.10 (deepforest publishes no cp310 wheel and no sdist).
export PATH=/n/data1/hms/dbmi/park/corinne/.conda/envs/xTea/bin:/home/cos689/software/bin:$PATH
export PYTHONNOUSERSITE=1     # keep ~/.local off sys.path (cost a SIGSEGV once: job 48777874)

[ -d "$WT/xtea" ] || { echo "FATAL: master worktree $WT/xtea missing"; exit 1; }
[ -d "$WT/xtea/genotyping/DF21_model_1_2" ] || { echo "FATAL: DF21 model missing in $WT"; exit 1; }
[ -f "$BAM" ] || { echo "FATAL: $BAM missing"; exit 1; }

rm -rf "$OUT"; mkdir -p "$OUT/NA12878/pub_clip"

echo "=============================================================="
echo "MASTER xTea -- whole genome, NON-CONCURRENT (true comparator)"
echo "  job        : ${SLURM_JOB_ID:-<none>}  on $(hostname)"
echo "  master code: $WT ($(cd "$WT" && git rev-parse --short HEAD) $(cd "$WT" && git log -1 --format=%s | cut -c1-50))"
echo "  python     : $(python --version 2>&1)  [$(which python)]"
echo "  -n per step: $NTHREADS   (golden used 8, on 8 cores, x3 concurrent jobs)"
echo "  cores      : ${SLURM_CPUS_PER_TASK:-24}"
echo "  types      : L1 -> Alu -> SVA, SEQUENTIALLY"
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
  [ -f "$SRC" ] || { echo "FATAL: golden script $SRC missing"; exit 1; }
  P="$OUT/NA12878/$t"
  mkdir -p "$P/tmp/clip" "$P/tmp/cns" "$P/tmp/transduction"
  printf '%s\tillumina\n' "$BAM" > "$P/bam_list.txt"   # same format as the golden's

  GEN="$OUT/$t.chain.sh"
  sed -e "s|^PREFIX=.*|PREFIX=$P/|" \
      -e "s|^XTEA_PATH=.*|XTEA_PATH=$WT/xtea/|" \
      -e "s|-n 8 |-n $NTHREADS |g" \
      -e "s|-n 8$|-n $NTHREADS|g" \
      -e "s|/home/cos689/scratch/xtea2_testdata/master_golden_wgs/NA12878/pub_clip/|$OUT/NA12878/pub_clip/|g" \
      -e "/^#SBATCH/d" -e "/^export PATH=/d" -e "/^#!/d" \
      "$SRC" > "$GEN"

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
# report (reporting failures must not mask the run: || true throughout)
# -----------------------------------------------------------------------------
echo
echo "##############################################################"
echo "# 1. TIMING -- master, serial, one 24-core job"
echo "##############################################################"
printf 'TOTAL WALL: %02d:%02d:%02d  (%d s)\n' $((WALL/3600)) $((WALL%3600/60)) $((WALL%60)) "$WALL"
for t in L1 Alu SVA; do
  w=${TYPE_WALL[$t]:-0}
  printf '  %-4s %02d:%02d:%02d   golden same type: ' "$t" $((w/3600)) $((w%3600/60)) $((w%60))
  case $t in L1) echo "37:15 (8 cores)";; Alu) echo "29:06 (8 cores)";; SVA) echo "38:31 (8 cores)";; esac
done
echo
echo "--- slowest steps (all types) ---"
sort -t$'\t' -k4 -rn "$TIMINGS" | grep -v '^type' | head -12 \
  | awk -F'\t' '{printf "  %-4s %-14s %8.2f min\n", $1, $3, $4/60}' || true
echo
echo "--- per-step-label totals across the three types ---"
grep -v '^type' "$TIMINGS" | awk -F'\t' '{s[$3]+=$4} END {for (k in s) printf "  %-14s %8.2f min\n", k, s[k]/60}' \
  | sort -k2 -rn || true

cat <<REF

--- THE COMPARISON (all same NA12878.wgs.bam, same node class) -----------------
  master, 3 CONCURRENT jobs, 8 cores each   38:31 wall   (the frozen golden;
                                                          1:44:52 summed job time)
  master, SERIAL, one 24-core job           $(printf '%02d:%02d:%02d' $((WALL/3600)) $((WALL%3600/60)) $((WALL%60))) wall   <-- THIS RUN
  xtea2,  SERIAL, one 24-core job           1:18:23 wall  (best known, job 48110538)
  xtea2,  after trim steps 1-4              see job 48824433

Read it this way: THIS RUN is the honest denominator for xtea2. The golden's 38:31
was three jobs deep; xtea2 has never had that luxury. If this run lands near or
above xtea2's 1:18:23, then xtea2 is at parity or ahead on equal footing and the
remaining gap to 38:31 is a CONCURRENCY deficit (xtea2 peaked at 925% of an
available 2400%), not slower code. If it lands well below, the trim has a real
regression and the per-step-label totals above say where.
REF

echo
echo "##############################################################"
echo "# 2. MEMORY"
echo "##############################################################"
echo "  This job runs one step at a time, so its peak RSS is the peak of the single"
echo "  heaviest step -- NOT comparable to a concurrent run's peak. Golden per-job"
echo "  peaks were 9.17 / 11.30 / 13.68 GB (L1/Alu/SVA) at -n 8; -n $NTHREADS should"
echo "  raise it (more workers, more per-worker annotation copies)."
echo "  xtea2 has sat at ~18.7 GB on every WGS run."
echo
echo "  read it from the scheduler once this finishes:"
echo "    sacct -j ${SLURM_JOB_ID:-<jobid>} --format=JobID%16,Elapsed,MaxRSS%12,MaxVMSize%12,AllocCPUS,State"

echo
echo "##############################################################"
echo "# 3. INTERMEDIATE FILE SIZES (master's own footprint)"
echo "##############################################################"
du -sh "$OUT" 2>/dev/null || true
echo
for t in L1 Alu SVA; do
  d="$OUT/NA12878/$t"
  [ -d "$d" ] || continue
  echo "[$t] $(du -sh "$d" 2>/dev/null | cut -f1) in $(find "$d" -type f | wc -l) files"
  find "$d" -type f -printf '%s\t%p\n' 2>/dev/null | sort -rn | head -10 \
    | awk -v r="$d/" '{s=$1;p=$2;sub(r,"",p); printf "    %8.1f MB  %s\n", s/1048576, p}' || true
done
echo
echo "--- shared pub_clip (master's cross-type clip sharing) ---"
du -sh "$OUT/NA12878/pub_clip" 2>/dev/null || true
echo
echo "  compare against xtea2's tree in the step45 run's section 3. xtea2 adds"
echo "  scan-pass artefacts (scan_clip_realign.tsv, scan_disc_pos.tsv, cns_spill_*)"
echo "  that master has no equivalent of -- that is a real disk cost of the C++ path."

echo
echo "##############################################################"
echo "# 4. CALL SET -- sanity signal, not a gate"
echo "##############################################################"
echo "Same code and cutoffs as the golden, so this SHOULD reproduce it; -n 8 -> $NTHREADS"
echo "can move things via the documented spill nondeterminism. A large gap means the"
echo "comparator is not measuring what the golden measured -- check before trusting"
echo "the timing."
printf '%-6s %10s %10s %8s\n' type this golden delta
for pair in "L1:LINE1:203" "Alu:ALU:1262" "SVA:SVA:73"; do
  t="${pair%%:*}"; rest="${pair#*:}"; alt="${rest%%:*}"; exp="${rest##*:}"
  f=$(ls "$OUT/NA12878/$t"/*_${alt}.vcf 2>/dev/null | head -1)
  n=$([ -n "$f" ] && grep -vc '^#' "$f" || echo NA)
  d="NA"; [ "$n" != "NA" ] && d=$((n - exp))
  printf '%-6s %10s %10s %8s\n' "$t" "$n" "$exp" "$d"
done

echo
echo "=============================================================="
echo "finished $(date '+%F %T')"
printf 'total wall %02d:%02d:%02d\n' $((WALL/3600)) $((WALL%3600/60)) $((WALL%60))
echo "=============================================================="
