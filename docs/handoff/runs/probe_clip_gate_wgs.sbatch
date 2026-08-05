#!/bin/bash
#SBATCH -J clipgate
#SBATCH -c 2
#SBATCH -t 6:00:00
#SBATCH --mem=96G
#SBATCH -p short
#SBATCH --account=park
#SBATCH -o /n/scratch/users/c/cos689/xtea2_testdata/probe_clip_gate_wgs_%j.out
#SBATCH --mail-type=NONE
# =============================================================================
# Gate [A] at WGS scale, swept offline. Plan §5's semantic half, on the gate that does 97% of
# the rejecting -- and it needs NO pipeline run, because a finished run's
# `tmp/transduction/temp_transduction_clip.sam` IS gate [A]'s entire input.
#
# WHY THIS IS A JOB AND NOT AN INLINE COMMAND, which is the only reason it is here: the L1 clip
# SAM is **2.46 GB** and the SVA one 1.6 GB, and the parse holds one string per qualified clip
# record. chr1's 147 MB equivalent runs inline in a couple of minutes; WGS does not belong on a
# login node or in a 2-core session. Everything else about it is identical to the chr1 invocation
# that already PASSED both of check 0's halves.
#
# WHAT IT MEASURES -- see `probe_clip_gate.py`'s docstring for the reasoning behind each:
#   1. the `n_clip_cutoff` sweep, and the fact that ONE swept value moves BOTH of gate [A]'s
#      thresholds because `n_polya_cutoff = n_clip_cutoff/2 - 1` is derived from it
#   2. whether `low_polya_clip` is reachable at all at the shipping cutoff (chr1: NO -- it first
#      fires at 5, so the polyA decision is a MEMBERSHIP test and the named threshold decides
#      nothing at this operating point)
#   3. the `n_clip` distribution of the `no_polya_clip` population, which the LEDGER CANNOT GIVE
#      -- `_drop(R_NO_POLYA_CLIP)` writes no detail
#   4. how much of gate [A]'s evidence is hits to the DECOY CONSENSUS, which the clip path does
#      not skip (the skip is commented out at `x_transduction.py:872-875`) while the disc path
#      tallies the same hits as cross-family contamination
#
# AND THE SCOPE CORRECTION, which is arm 2. `NEARBY_CLIP` was measured as moving "exactly one
# locus in the whole genome" -- INSIDE GATE [B2] (`probe_nearby_clip.py`). The same constant also
# filters the clip SAM parse at `x_transduction.py:865`, ahead of the unmapped and clip-qualified
# tests, and that is gate [A]'s DOMAIN. So the second arm re-parses at unbounded distance. It is
# a different population from the one the earlier sweep bounded, and this is the run that says
# whether the earlier claim needs its scope stated.
#
# CHECK 0 IS TWO CHECKS AND THE JOB FAILS ON EITHER. At the shipping (NEARBY_CLIP=50,
# n_clip_cutoff from the run's own ledger) pair the probe must reproduce the run's ledger `A_clip`
# reason split AND its `.unique_trsdct_half_clip` file as a SET of loci. Both PASSED on chr1
# (L1 1184/499 and 217 loci; SVA 1360/919 and 270). A mismatch here means this harness is
# misconfigured -- it is NOT a finding.
#
# WHAT IT CANNOT SAY: whether any locus becomes a CALL. Gate [A]'s survivors still face [B], [B2],
# [C], [D] and [E], and 306 of BASELINE-0's 322 false negatives never reach this decision at all.
# Expect a precision-shaped effect. Nothing here justifies moving a threshold on its own; a
# truth-set-scored two-arm run does.
# =============================================================================
set -euo pipefail

REPO=${REPO:-/n/data1/hms/dbmi/park/corinne/claude_sandbox/xTea}
# The post-flip WGS payoff run (job 49792012), which is the CURRENT default behaviour: the decoy
# test is structural there. Its L1 and Alu output is byte-identical to BASELINE-0, so for gate [A]
# -- which is upstream of the decoy classification entirely -- L1 is the same population either
# arm would have given. SVA is the family where the run choice could matter.
RUN=${RUN:-/n/scratch/users/c/cos689/xtea2_testdata/gate_decoy_structural_wgs_49792012/HG002.wgs}
OUT=${OUT:-/n/scratch/users/c/cos689/xtea2_testdata/probe_clip_gate_wgs_${SLURM_JOB_ID:-manual}}
SWEEP=${SWEEP:-1,2,3,4,5,6,8,10,15,20}
NEARBY=${NEARBY:-50,1000000000}

mkdir -p "$OUT"
echo "=== probe_clip_gate WGS ==="
echo "job     : ${SLURM_JOB_ID:-manual} on $(hostname)"
echo "repo    : $REPO ($(cd "$REPO" && git rev-parse --short HEAD))"
echo "run dir : $RUN"
echo "sweep   : $SWEEP"
echo "nearby  : $NEARBY"
echo "out     : $OUT"
echo

# The inputs this depends on, sized -- so an empty or truncated SAM is visible in the log rather
# than showing up as a family with no transduction population.
for fam in L1 SVA; do
  sam="$RUN/$fam/tmp/transduction/temp_transduction_clip.sam"
  if [[ ! -s "$sam" ]]; then
    echo "FATAL: missing or empty clip SAM: $sam"
    exit 2
  fi
  echo "  $fam clip SAM: $(stat -c %s "$sam") B"
done
echo

t0=$SECONDS
set +e
python3 "$REPO/docs/handoff/runs/probe_clip_gate.py" \
  --run-dir "$RUN" --sweep "$SWEEP" --nearby "$NEARBY" \
  2>&1 | tee "$OUT/probe_clip_gate.txt"
rc=${PIPESTATUS[0]}
set -e
echo
echo "elapsed: $((SECONDS - t0)) s"

# The probe returns non-zero iff a check 0 half failed, so the job's verdict is its verdict.
if [[ $rc -ne 0 ]]; then
  echo "RESULT: FAIL -- check 0 did not reproduce the run (probe exit $rc). This harness is"
  echo "        misconfigured; do not read the sweep as a measurement."
  exit 1
fi
echo "RESULT: PASS -- both halves of check 0 reproduced the run; the sweep is a measurement."
echo "artifacts: $OUT/probe_clip_gate.txt"
