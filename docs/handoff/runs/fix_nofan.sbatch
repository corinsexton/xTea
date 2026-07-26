#!/bin/bash
#SBATCH -c 24
#SBATCH -t 8:00:00
#SBATCH --mem=64G
#SBATCH -p priopark
#SBATCH --account=park_contrib
#SBATCH -o /n/scratch/users/c/cos689/xtea2_testdata/fix_nofan_%j.out
#SBATCH --mail-type=NONE

# =============================================================================
# xTea2 full-genome NA12878 run after trim steps 4 + 5, INSTRUMENTED.
#
# Deliverables this run exists to produce (user request, 2026-07-25 late night):
#   1. wall/CPU timing, compared against the naive original golden run
#   2. peak memory
#   3. intermediate file sizes
#
# NOT a baseline-conformance run. The B8 fix (trim step 2) deliberately moves calls
# in the more-conservative direction and the goldens predate it; the user explicitly
# accepted a shifted baseline for this first pass. So call-set differences are
# REPORTED and never asserted -- this script's exit status reflects whether the
# pipeline ran, not whether it matched.
#
# Runs from a git WORKTREE, not the live checkout. run_xtea.py spawns a fresh
# x_TEA_main.py subprocess per step, so a job sharing the checkout picks up
# half-applied edits mid-run -- which is why job 48765024 had to be cancelled. The
# worktree is created at submit time by step45_submit.sh and is frozen for the run's
# duration; the live repo can be edited freely while this runs.
# =============================================================================
set -euo pipefail

WT="${XTEA2_WORKTREE:-/n/scratch/users/c/cos689/xtea2_testdata/step45_wt}"
OUT="${XTEA2_OUT:-/n/scratch/users/c/cos689/xtea2_testdata/fix_nofan_out}"
SCRATCH=/n/scratch/users/c/cos689/xtea2_testdata

# One lever: activating xtea2 is the whole environment setup. The env itself carries
# PYTHONNOUSERSITE=1 (registered with `conda env config vars set`), so ~/.local stays
# off sys.path and the env's own onnxruntime 1.22.0 wins. run_xtea.py hands its own env
# down to every step subprocess, so this covers all 18 steps, genotype included.
source /n/app/conda/miniforge3/24.11.3-0/etc/profile.d/conda.sh
conda activate xtea2
# samtools 1.22 / bwa 0.7.17 are not conda-managed; same binaries as the golden runs.
export PATH=/home/cos689/software/bin:$PATH
export XTEA_BASIC_INFO_SEED=42   # pin the insert-size draw (matches the seeded golden)

[ -d "$WT" ] || { echo "FATAL: worktree $WT missing -- run step45_submit.sh, not sbatch directly"; exit 1; }
[ -x "$WT/src/csrc/xtea_scan" ] || { echo "FATAL: $WT/src/csrc/xtea_scan missing or not executable"; exit 1; }

rm -rf "$OUT"; mkdir -p "$OUT"
cd "$WT"

echo "=============================================================="
echo "xTea2 WGS FIX-A: --cns-emit-scan REMOVED (cns_filter does indexed per-site fetch)"
echo "  job          : ${SLURM_JOB_ID:-<none>}  on $(hostname)"
echo "  worktree     : $WT"
echo "  commit       : $(git rev-parse --short HEAD) $(git log -1 --format=%s | cut -c1-70)"
echo "  worktree diff: $(git status --porcelain | wc -l) file(s) (expect 0)"
echo "  xtea/ size   : $(find xtea -name '*.py' | wc -l) modules, $(find xtea -name '*.py' -exec cat {} + | wc -l) lines"
echo "  cores        : ${SLURM_CPUS_PER_TASK:-24}"
echo "  started      : $(date '+%F %T')"
echo "=============================================================="

T_START=$(date +%s)

/usr/bin/time -v python xtea/run_xtea.py \
  -i "$SCRATCH/NA12878.wgs.bam" \
  -s NA12878.wgs \
  -o "$OUT" \
  --reference /home/cos689/hg38_no_alt.fa \
  --rep-lib-dir /home/cos689/xTea/rep_lib_annotation \
  --gene-gff3 /home/cos689/data1/corinne/ref/gencode.v45.basic.annotation.gff3 \
  --blacklist /home/cos689/xTea/rep_lib_annotation/blacklist/hg38/centromere.bed \
  --types L1,Alu,SVA \
  --scan-bin "$WT/src/csrc/xtea_scan" \
  --bwa /home/cos689/software/bin/bwa \
  --scan-min-clip 2 \
  --col5-families alu,sva \
  --parallel-downstream \
  -n 24 \
  2>&1 | tee "$OUT/time_v.txt"

T_END=$(date +%s)
WALL=$((T_END - T_START))

# -----------------------------------------------------------------------------
# Everything below is reporting only. `|| true` throughout: a reporting failure
# must not mask a successful pipeline run (or hide a failed one behind a
# confusing error).
# -----------------------------------------------------------------------------
S="$OUT/NA12878.wgs"

echo
echo "##############################################################"
echo "# 1. TIMING"
echo "##############################################################"
printf 'this run: %02d:%02d:%02d wall (%ds)\n' $((WALL/3600)) $((WALL%3600/60)) $((WALL%60)) "$WALL"
grep -E "Percent of CPU|Elapsed \(wall clock\)|Maximum resident" "$OUT/time_v.txt" 2>/dev/null || true

cat <<'REF'

--- reference points (measured, same NA12878.wgs.bam, same node class) ----------
NAIVE ORIGINAL GOLDEN (master xTea, its own SLURM pipeline, py3.8 + deep-forest):
  three per-type jobs submitted together, -c 8 --mem=192G each = 24 cores total
    48013142  L1   37:15   MaxRSS  9.17 GB
    48013143  Alu  29:06   MaxRSS 11.30 GB
    48013144  SVA  38:31   MaxRSS 13.68 GB
  => WALL CLOCK 38:31 (they ran concurrently); 1:44:52 summed job time
  NOTE the shape difference: the golden got 3 processes x 8 cores; xtea2 is ONE
  process given 24 cores. Same total allocation, so wall clock is comparable, but
  the golden's parallelism was free (3 independent jobs) where xtea2 has to
  manufacture it (--parallel-downstream). Read the gap as "xtea2 has not yet
  recovered the concurrency master got from the job scheduler", not as pure
  slowdown.

xTEA2 HISTORY (single job, -c 24, /usr/bin/time -v):
  48026097  a1  clip swap only            2:08:42   650% CPU   18.41 GB
  48034244  a2  clip+disc swap            2:00:15   626% CPU   18.51 GB
  48051732  a2  clip+disc swap (repeat)   2:00:21   640% CPU   18.75 GB
  48075090  a2  + cns-emit                2:27:12   491% CPU   18.68 GB
  48096921  a2  speed2                    3:31:50   348% CPU   18.67 GB
  48110538  a2  speed3 (fastest known)    1:18:23   925% CPU   18.68 GB
  48228181  moveA2 unified downstream     1:22:06   895% CPU   18.67 GB
  <this>    trim steps 1-5                 see above

Interpreting CPU%: 24 cores == 2400%. Every xtea2 run so far has left most of the
allocation idle; speed3/moveA2 at ~900% are the best achieved. If this run's CPU%
is well under 900% the regression is in parallelism, not in the trim.
REF

echo
echo "--- per-phase and per-step wall time (from run_xtea's own instrumentation) ---"
if [ -f "$S/step_timings.tsv" ]; then
  python - "$S/step_timings.tsv" <<'PYEOF' || true
import sys, collections
rows = [l.rstrip("\n").split("\t") for l in open(sys.argv[1]).read().splitlines()[1:] if l.strip()]
if not rows:
    print("  (empty)"); raise SystemExit
phases = [(n, float(s), float(g)) for k, n, s, g in rows if k == "phase"]
steps  = [(n, float(s), float(g)) for k, n, s, g in rows if k.startswith("step")]
tot = sum(s for _, s, _ in phases)
print(f"  PHASES (these do not overlap; they sum to the run: {tot/60:.1f} min)")
for n, s, g in phases:
    pct = 100.0 * s / tot if tot else 0
    print(f"    {n:<24} {s/60:8.2f} min  {pct:5.1f}%   children peak RSS {g:5.2f} GB")
print(f"  STEPS ({len(steps)}; under --parallel-downstream the three types OVERLAP --")
print( "         do not sum these, and compare only within one type)")
by_type = collections.defaultdict(list)
for n, s, g in steps:
    t = n.split(":")[0] if ":" in n else "shared"
    by_type[t].append((n, s, g))
for t in sorted(by_type):
    sub = by_type[t]
    print(f"    [{t}] {sum(s for _, s, _ in sub)/60:.1f} min over {len(sub)} steps")
    for n, s, g in sorted(sub, key=lambda r: -r[1]):
        print(f"      {n:<28} {s/60:8.2f} min   peak-so-far {g:5.2f} GB")
slow = sorted(steps, key=lambda r: -r[1])[:5]
print("  SLOWEST 5 STEPS OVERALL")
for n, s, g in slow:
    print(f"    {n:<30} {s/60:8.2f} min")
PYEOF
else
  echo "  MISSING: $S/step_timings.tsv (run_xtea writes it unless --dry-run)"
fi

echo
echo "##############################################################"
echo "# 2. MEMORY"
echo "##############################################################"
grep -E "Maximum resident set size" "$OUT/time_v.txt" 2>/dev/null \
  | awk '{printf "  /usr/bin/time -v peak RSS: %.2f GB\n", $NF/1048576}' || true
echo "  golden per-job peak RSS was 9.17 / 11.30 / 13.68 GB (L1/Alu/SVA)."
echo "  xtea2 has sat at ~18.7 GB on every prior WGS run -- one process holding what"
echo "  the golden split across three, so a like-for-like read is 18.7 vs 13.68 (the"
echo "  largest single golden job), not vs their sum."
echo "  Which phase set the high-water mark is in the children-peak column above."
echo
echo "  after the job finishes, for the scheduler's own accounting:"
echo "    sacct -j ${SLURM_JOB_ID:-<jobid>} --format=JobID%16,Elapsed,MaxRSS%12,MaxVMSize%12,AllocCPUS,State"

echo
echo "##############################################################"
echo "# 3. INTERMEDIATE FILE SIZES"
echo "##############################################################"
echo "--- total output tree ---"
du -sh "$OUT" 2>/dev/null || true
echo
echo "--- per top-level component ---"
du -sh "$S"/* 2>/dev/null | sort -h || true
echo
echo "--- per type: total, file count, and the 12 largest files ---"
for t in L1 Alu SVA; do
  [ -d "$S/$t" ] || continue
  echo
  echo "[$t] $(du -sh "$S/$t" 2>/dev/null | cut -f1) in $(find "$S/$t" -type f | wc -l) files"
  find "$S/$t" -type f -printf '%s\t%p\n' 2>/dev/null | sort -rn | head -12 \
    | awk -v r="$S/$t/" '{s=$1; p=$2; sub(r,"",p); printf "    %8.1f MB  %s\n", s/1048576, p}' || true
done
echo
echo "--- scan-pass artefacts (shared across types; the new cost xtea2 adds) ---"
find "$S" -maxdepth 1 -type f -printf '%s\t%p\n' 2>/dev/null | sort -rn \
  | awk -v r="$S/" '{s=$1; p=$2; sub(r,"",p); if (s>1048576) printf "    %8.1f MB  %s\n", s/1048576, p}' || true
echo
echo "--- what is safe to delete afterwards (tmp scratch, not results) ---"
du -sh "$S"/*/tmp 2>/dev/null | sort -h || true

echo
echo "##############################################################"
echo "# 4. CALL SET -- reported, NOT gated (baseline shift is expected)"
echo "##############################################################"
echo "The goldens predate the trim step 2 B8 fix, which deliberately drops slightly"
echo "more calls (min-divergence rule, strictly more conservative). Differences here"
echo "are the expected consequence, not a regression. Read the numbers; the exit"
echo "status of this job says nothing about them."
echo
echo "--- per-type call counts vs the seeded golden ---"
printf '%-6s %10s %10s %8s\n' type this golden delta
for pair in "L1:LINE1" "Alu:ALU" "SVA:SVA"; do
  t="${pair%%:*}"; alt="${pair##*:}"
  f=$(ls "$S"/merged/*_${alt}.vcf 2>/dev/null | head -1)
  g=$(ls "$SCRATCH/master_golden_wgs_seeded/NA12878/$t"/*.vcf 2>/dev/null | head -1)
  n=$([ -n "$f" ] && grep -vc '^#' "$f" || echo NA)
  m=$([ -n "$g" ] && grep -vc '^#' "$g" || echo NA)
  d="NA"; [ "$n" != "NA" ] && [ "$m" != "NA" ] && d=$((n - m))
  printf '%-6s %10s %10s %8s\n' "$t" "$n" "$m" "$d"
done
echo
echo "  (frozen UNSEEDED golden counts for reference: L1 203, Alu 1262, SVA 73)"
echo
echo "--- merged single VCF ---"
mf=$(ls "$S"/*.vcf 2>/dev/null | head -1)
[ -n "$mf" ] && echo "  $(grep -vc '^#' "$mf") calls -> $mf" || echo "  MISSING merged VCF"
echo
echo "--- positions gained / lost per type vs the seeded golden ---"
for pair in "L1:LINE1" "Alu:ALU" "SVA:SVA"; do
  t="${pair%%:*}"; alt="${pair##*:}"
  f=$(ls "$S"/merged/*_${alt}.vcf 2>/dev/null | head -1)
  g=$(ls "$SCRATCH/master_golden_wgs_seeded/NA12878/$t"/*.vcf 2>/dev/null | head -1)
  if [ -n "$f" ] && [ -n "$g" ]; then
    lost=$(comm -23 <(grep -v '^#' "$g" | cut -f1,2 | sort -u) <(grep -v '^#' "$f" | cut -f1,2 | sort -u) | wc -l)
    gain=$(comm -13 <(grep -v '^#' "$g" | cut -f1,2 | sort -u) <(grep -v '^#' "$f" | cut -f1,2 | sort -u) | wc -l)
    echo "  $t: -$lost lost / +$gain gained (exact position match; B8 predicts lost >> gained)"
  else
    echo "  $t: could not locate both VCFs"
  fi
done

echo
echo "--- concordance gate, for information only (asserts >= 1.0; will 'fail') ---"
XTEA2_OUTPUT_DIR="$S" \
XTEA2_GOLDEN_DIR="$SCRATCH/master_golden_wgs_seeded/NA12878" \
XTEA2_SCAN_BIN="$WT/src/csrc/xtea_scan" \
python -m pytest tests/test_e2e_concordance.py -q 2>&1 | tail -20 || true

echo
echo "=============================================================="
echo "finished $(date '+%F %T')"
printf 'total wall %02d:%02d:%02d\n' $((WALL/3600)) $((WALL%3600/60)) $((WALL%60))
echo "=============================================================="
