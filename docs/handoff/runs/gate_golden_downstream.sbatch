#!/bin/bash
#SBATCH -J golden_down
#SBATCH -A park_contrib
#SBATCH -p priopark
#SBATCH -t 1:00:00
#SBATCH -c 8
#SBATCH --mem 32G
#SBATCH -o /n/scratch/users/c/cos689/xtea2_testdata/gate_golden_downstream_%j.out
#SBATCH --mail-type=NONE
# =============================================================================
# `tests/test_downstream_from_cns.py` -- the NA12878 chr1 golden downstream gate, as a job.
#
# WHY THIS SCRIPT EXISTS. This is **the cheapest end-to-end signal in the repo** (~10 minutes,
# all three families) and it is the gate that caught the return-arity bug the whole fast suite
# was blind to (job 49633492: `filter_transduction` returning three values where its callers
# unpack two, found ten minutes into a cluster job). Until now it was submitted with throwaway
# scripts -- jobs 49634291 and others -- so the invocation lived only in a scrollback. It is a
# committed script now so it can be the FIRST thing run on any downstream change.
#
# WHAT MAKES IT INDEPENDENT of the HG002 A/B gates, which is why it is worth its own job:
#   * a DIFFERENT sample (NA12878, not HG002), and
#   * a frozen EXTERNAL baseline (master xTea's own outputs), not this repo's own prior run.
# It stages the golden's own post-cns outputs and re-runs only the tail (transduction -> sibling
# -> postF x2 -> gene -> genotype -> gVCF), byte-comparing every intermediate and **gating on
# VCF columns 1-8**.
#
# TWO CLASSES OF EXPECTED NOTE, and neither is a failure -- do not "fix" them here:
#   * L1's four "intermediate moved (not gated)" lines are the pre-existing B8 carry-forward the
#     test's own docstring names: everything downstream of `…clip_sites_raw_disc.txt.slct` being
#     3,272 rows where the stale golden has 3,302.
#   * Alu's three genotype flips are the documented DF21-golden-vs-ONNX class.
#
# USAGE  sbatch gate_golden_downstream.sbatch <worktree> [scan-bin]
#   The conftest defaults resolve on this machine for the golden dir, reference, BAM,
#   rep_lib_annotation, GFF3 and blacklist, so only the worktree is required. Override any of
#   them with the XTEA2_* env vars conftest.py documents.
# =============================================================================
set -euo pipefail

SCRATCH=/n/scratch/users/c/cos689/xtea2_testdata
WT="${1:?usage: sbatch gate_golden_downstream.sbatch <worktree> [scan-bin]}"
SCAN_BIN="${2:-$WT/src/csrc/xtea_scan}"
OUT="${XTEA2_OUT:-$SCRATCH/gate_golden_downstream_${SLURM_JOB_ID:-manual}}"

source /n/app/conda/miniforge3/24.11.3-0/etc/profile.d/conda.sh
conda activate xtea2
export PATH=/home/cos689/software/bin:$PATH
test "$PYTHONNOUSERSITE" = 1     # the env's guard against ~/.local's onnxruntime

[ -d "$WT" ] || { echo "FATAL: worktree $WT missing"; exit 1; }
[ -z "$(git -C "$WT" status --porcelain)" ] || { echo "FATAL: $WT is dirty"; exit 1; }
[ -x "$SCAN_BIN" ] || { echo "FATAL: $SCAN_BIN not built"; exit 1; }

mkdir -p "$OUT"
cd "$WT"

echo "=============================================================="
echo "GOLDEN DOWNSTREAM GATE  --  NA12878 chr1, native gather forced on"
echo "  job     : ${SLURM_JOB_ID:-<none>} on $(hostname)"
echo "  commit  : $(git rev-parse --short HEAD) $(git log -1 --format=%s | cut -c1-60)"
echo "  worktree: $WT"
echo "  scan bin: $SCAN_BIN"
echo "  out     : $OUT"
echo "  started : $(date '+%F %T')"
echo "=============================================================="

# `set +e` around pytest: its nonzero exit IS the result, and under `set -e` the script would
# die before printing the log that says which family failed and how.
set +e
t0=$(date +%s)
XTEA_SCAN_TD_GATHER="$SCAN_BIN" \
XTEA2_SCAN_BIN="$SCAN_BIN" \
    python -m pytest tests/test_downstream_from_cns.py -q -s \
    > "$OUT/downstream_from_cns.log" 2>&1
PYRC=$?
wall=$(( $(date +%s) - t0 ))
set -e
printf '  wall %02d:%02d:%02d\n' $((wall/3600)) $((wall%3600/60)) $((wall%60))

echo
cat "$OUT/downstream_from_cns.log"

echo
echo "--- the notes, collected (expected: 4 L1 carry-forward, 3 Alu genotype flips) ---"
grep -c 'intermediate moved (not gated)' "$OUT/downstream_from_cns.log" |
    sed 's/^/  "intermediate moved" lines: /' || true
grep -c 'genotype flips' "$OUT/downstream_from_cns.log" |
    sed 's/^/  "genotype flips" lines:     /' || true

echo
echo "=============================================================="
if [ "$PYRC" -eq 0 ]; then
    echo "GOLDEN DOWNSTREAM GATE: PASS"
else
    echo "GOLDEN DOWNSTREAM GATE: FAIL (pytest exit $PYRC)"
    echo "  full log: $OUT/downstream_from_cns.log"
fi
echo "  finished: $(date '+%F %T')"
echo "=============================================================="
exit "$PYRC"
