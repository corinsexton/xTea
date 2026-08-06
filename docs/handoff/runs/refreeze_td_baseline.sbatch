#!/bin/bash
#SBATCH -J td_refreeze
#SBATCH -c 8
#SBATCH -t 3:00:00
#SBATCH --mem=32G
#SBATCH -p short
#SBATCH --account=park
#SBATCH -o /n/scratch/users/c/cos689/xtea2_testdata/refreeze_td_baseline_%j.out
#SBATCH --mail-type=NONE
# =============================================================================
# Re-freeze the 3-arm HG002 transduction baseline with the self-window guard ON (decision 6).
#
# WHY. `verify_transduction` now disqualifies an entry whose source window CONTAINS its own
# insertion point: reads are gathered AT the insertion and counted for aligning to the window, so
# a window spanning that point verifies itself off ordinary local coverage. The frozen 3-arm
# baseline contains **2 such entries of 115** -- and they are its **two highest-support entries
# by 3x**, `chr17:45585847` at support 143 and `chr6:31329020` at 76, against a #3 of 26 and a
# median of 1. Both were VERIFIED. So the frozen precision range 59.1%-81.0% includes them and
# the baseline has to move. Corinne's decision, 2026-08-06.
#
# **THE POOL DOES NOT CHANGE, ONLY THE ADJUDICATION.** `td_pool_hg002_3arm.tsv` is what the three
# callers found; the guard is a statement about what the READS can prove. Re-running the pool
# builder here would conflate "the baseline moved because the verifier got stricter" with "the
# baseline moved because the pool was rebuilt", and only the first is intended. The pool file is
# copied across and its sha256 asserted UNCHANGED.
#
# THE THREE CHECKS, and check 1 is the one that makes the other two trustworthy:
#
#   1. **CONTAINMENT: `--allow-self-window` must reproduce the OLD verdicts EXACTLY.** Same pool,
#      same BAM, same operating point, guard disabled => byte-identical to the frozen
#      `verdicts.tsv`. This is what proves the only behavioural change is the guard. Without it a
#      re-freeze cannot tell "the guard moved 2 entries" from "something else moved 2 entries and
#      the guard moved none" -- and the verifier's own history has exactly that failure mode
#      (a gate that PASSED while measuring the wrong pair, HANDOFF_2026-08-05 sec 4.4).
#   2. **The guard-ON run differs in EXACTLY the 2 known entries**, both VERIFIED -> UNDECIDABLE.
#      More than 2 means the guard is catching something unintended; fewer means it did not run.
#   3. **The new baseline round-trips**: re-adjudicating the new verdicts against themselves with
#      `--baseline` must report no regression, which is the reader agreeing with what was written.
#
# THE OLD BASELINE IS NOT DELETED. It moves to `td_gate_baseline_hg002_preguard/` and stays
# reproducible via `--allow-self-window`, because every gate measured against it (the frozen
# precision range, the decoy flip's no-regression arm) quoted numbers that only mean anything
# next to the baseline they were taken against.
#
# USAGE  sbatch refreeze_td_baseline.sbatch [worktree]
# =============================================================================
set -euo pipefail

SCRATCH=/n/scratch/users/c/cos689/xtea2_testdata
WT="${1:-/n/data1/hms/dbmi/park/corinne/claude_sandbox/xTea}"
OUT="${XTEA2_OUT:-$SCRATCH/refreeze_td_baseline_${SLURM_JOB_ID:-manual}}"
FROZEN="${XTEA2_TD_FROZEN:-$SCRATCH/td_gate_baseline_hg002}"
PREGUARD="${XTEA2_PREGUARD:-$SCRATCH/td_gate_baseline_hg002_preguard}"
TESTDATA=/n/data1/hms/dbmi/park/corinne/claude_sandbox/test_data/HG002
BAM="$TESTDATA/alignment/HG002.novaseq.pcr-free.30x.dedup.grch38.bam"
REF=/n/data1/hms/dbmi/park-smaht_dac/ref/GRCh38_no_alt/hg38_no_alt.fa
REPLIB=/home/cos689/xTea/rep_lib_annotation
CORES="${SLURM_CPUS_PER_TASK:-8}"
# Set to 1 to actually replace $FROZEN. Default DRY: the job measures and reports, and a human
# reads check 1 before anything overwrites a frozen artifact three gates cite.
COMMIT_REFREEZE="${XTEA2_COMMIT_REFREEZE:-0}"

source /n/app/conda/miniforge3/24.11.3-0/etc/profile.d/conda.sh
conda activate xtea2
export PATH=/home/cos689/software/bin:$PATH
test "$PYTHONNOUSERSITE" = 1

VERIFY="$WT/benchmark/verify_transduction.py"
POOL="$FROZEN/td_pool_hg002_3arm.tsv"

[ -f "$BAM" ]    || { echo "FATAL: $BAM missing"; exit 1; }
[ -f "$POOL" ]   || { echo "FATAL: $POOL missing"; exit 1; }
[ -f "$FROZEN/verdicts.tsv" ] || { echo "FATAL: no frozen verdicts to compare against"; exit 1; }
[ -f "$VERIFY" ] || { echo "FATAL: $VERIFY missing"; exit 1; }
grep -q 'allow-self-window' "$VERIFY" || {
    echo "FATAL: $VERIFY has no --allow-self-window; wrong commit"; exit 1; }
grep -q 'def self_window' "$VERIFY" || {
    echo "FATAL: $VERIFY has no self_window guard; wrong commit"; exit 1; }

mkdir -p "$OUT"
RC=0

echo "=============================================================="
echo "RE-FREEZE: the 3-arm HG002 td baseline, with the self-window guard ON"
echo "  job    : ${SLURM_JOB_ID:-<none>} on $(hostname)  cores=$CORES"
echo "  commit : $(git -C "$WT" rev-parse --short HEAD 2>/dev/null || echo '?')"
echo "  pool   : $POOL"
echo "  frozen : $FROZEN"
echo "  commit-refreeze: $COMMIT_REFREEZE (0 = measure and report only)"
echo "  started: $(date '+%F %T')"
echo "=============================================================="

echo
echo "--- check 0: the pool is the frozen one, unchanged ---"
# The guard is a statement about the READS, not about which entries exist. If this moved, the
# re-freeze would be measuring two changes at once.
sha_now=$(sha256sum "$POOL" | awk '{print $1}')
sha_rec=$(awk '/td_pool_hg002_3arm.tsv/ {print $1}' "$FROZEN/SHA256SUMS" 2>/dev/null || echo "")
echo "  pool sha256 : $sha_now"
echo "  recorded    : ${sha_rec:-<none recorded>}"
if [ -n "$sha_rec" ] && [ "$sha_now" != "$sha_rec" ]; then
    echo "  FATAL: the pool has changed since it was frozen. Re-freezing now would conflate a"
    echo "  stricter verifier with a rebuilt pool."; exit 1
fi
echo "  entries: $(( $(grep -vc '^#' "$POOL") - 1 ))"

echo
echo "--- check 1: CONTAINMENT -- guard OFF must reproduce the frozen verdicts exactly ---"
R_OFF="$OUT/guard_off"; mkdir -p "$R_OFF"
set +e
python "$VERIFY" --pool-tsv "$POOL" --bam "$BAM" --ref "$REF" --rep-lib-dir "$REPLIB" \
    --threads "$CORES" --out "$R_OFF" --allow-self-window \
    --baseline "$FROZEN/verdicts.tsv" > "$OUT/guard_off.report.txt" 2>"$OUT/guard_off.err.txt"
rc_off=$?
set -e
tail -25 "$OUT/guard_off.report.txt" | sed 's/^/    /'
if [ "$rc_off" -eq 3 ]; then
    echo "  FATAL: the operating point did not match -- the sweep midpoints moved."; exit 1
fi
if cmp -s "$R_OFF/verdicts.tsv" "$FROZEN/verdicts.tsv"; then
    echo "  PASS: byte-identical to the frozen verdicts. The ONLY behavioural change is the guard."
else
    echo "  FAIL: guard OFF does not reproduce the frozen baseline, so nothing below is"
    echo "  attributable to the guard. Differing lines:"
    diff "$FROZEN/verdicts.tsv" "$R_OFF/verdicts.tsv" | head -20 || true
    RC=1
fi

echo
echo "--- check 2: guard ON -- exactly 2 entries move, both VERIFIED -> UNDECIDABLE ---"
R_ON="$OUT/guard_on"; mkdir -p "$R_ON"
set +e
python "$VERIFY" --pool-tsv "$POOL" --bam "$BAM" --ref "$REF" --rep-lib-dir "$REPLIB" \
    --threads "$CORES" --out "$R_ON" > "$OUT/guard_on.report.txt" 2>"$OUT/guard_on.err.txt"
rc_on=$?
set -e
[ "$rc_on" -eq 0 ] || { echo "  FATAL: the guard-ON run exited $rc_on"; tail -20 "$OUT/guard_on.err.txt"; exit 1; }
tail -30 "$OUT/guard_on.report.txt" | sed 's/^/    /'

python - "$FROZEN/verdicts.tsv" "$R_ON/verdicts.tsv" <<'EOF' || RC=1
import sys
def read(p):
    out = {}
    for line in open(p):
        if line.startswith("#"):
            continue
        f = line.rstrip("\n").split("\t")
        if len(f) < 15 or f[0] == "family":
            continue
        out[(f[0], f[1], f[2], f[15] if len(f) > 15 else "")] = (f[13], f[14])
    return out
old, new = read(sys.argv[1]), read(sys.argv[2])
if set(old) != set(new):
    print("    FAIL: the entry SET moved, which the pool sha said it could not")
    raise SystemExit(1)
moved = [(k, old[k], new[k]) for k in old if old[k][0] != new[k][0]]
print("    entries whose STATE changed: {0} (expect exactly 2)".format(len(moved)))
for k, o, n in sorted(moved):
    print("      {0} {1}:{2}  {3} ({4})  ->  {5} ({6})".format(k[0], k[1], k[2], o[0], o[1], n[0], n[1]))
ok = len(moved) == 2
for k, o, n in moved:
    if o[0] != "VERIFIED" or n[0] != "UNDECIDABLE" or not n[1].startswith("self_window"):
        print("    FAIL: {0}:{1} did not move VERIFIED -> UNDECIDABLE(self_window)".format(k[1], k[2]))
        ok = False
expected = {("chr17", "45585847"), ("chr6", "31329020")}
got = set((k[1], k[2]) for k, _o, _n in moved)
if got != expected:
    print("    FAIL: the moved entries are {0}, expected {1}".format(sorted(got), sorted(expected)))
    ok = False
if not ok:
    raise SystemExit(1)
print("    PASS: exactly the two known circular entries moved, and only their state")
EOF

echo
echo "--- check 3: the new baseline round-trips through its own reader ---"
R_RT="$OUT/roundtrip"; mkdir -p "$R_RT"
set +e
python "$VERIFY" --pool-tsv "$POOL" --bam "$BAM" --ref "$REF" --rep-lib-dir "$REPLIB" \
    --threads "$CORES" --out "$R_RT" \
    --baseline "$R_ON/verdicts.tsv" > "$OUT/roundtrip.report.txt" 2>"$OUT/roundtrip.err.txt"
rc_rt=$?
set -e
grep -aE "NO-REGRESSION|REGRESSION|totals:" "$OUT/roundtrip.report.txt" | sed 's/^/    /' || true
[ "$rc_rt" -eq 0 ] || { echo "  FAIL: the new baseline does not validate against itself (rc=$rc_rt)"; RC=1; }
cmp -s "$R_RT/verdicts.tsv" "$R_ON/verdicts.tsv" &&
    echo "  PASS: deterministic -- two guard-ON runs agree byte-for-byte" || {
    echo "  FAIL: two guard-ON runs of the same pool disagree"; RC=1; }

echo
echo "--- the new stratified precision, which is what downstream gates will quote ---"
sed -n '/PRIMARY STRATUM/,/^$/p' "$OUT/guard_on.report.txt" | sed 's/^/  /' | head -30 || true

if [ "$RC" -ne 0 ]; then
    echo
    echo "RE-FREEZE: FAIL -- nothing was replaced."
    exit 1
fi

echo
if [ "$COMMIT_REFREEZE" = "1" ]; then
    echo "--- replacing the frozen baseline (XTEA2_COMMIT_REFREEZE=1) ---"
    # The old one is MOVED, not deleted: every number quoted against it (the 59.1%-81.0% range,
    # the decoy flip's no-regression arm) is meaningless without the baseline it was taken from.
    if [ -e "$PREGUARD" ]; then
        echo "  FATAL: $PREGUARD already exists -- refusing to overwrite a preserved baseline"; exit 1
    fi
    cp -a "$FROZEN" "$PREGUARD"
    cp "$R_ON/verdicts.tsv" "$FROZEN/verdicts.tsv"
    cp "$OUT/guard_on.report.txt" "$FROZEN/verify_report.txt"
    ( cd "$FROZEN" && sha256sum td_pool_hg002_3arm.tsv verdicts.tsv > SHA256SUMS )
    echo "  pre-guard baseline preserved at : $PREGUARD"
    echo "  new SHA256SUMS:"; sed 's/^/    /' "$FROZEN/SHA256SUMS"
    echo "  reproduce the pre-guard verdicts any time with --allow-self-window"
else
    echo "--- NOT replacing the frozen baseline (XTEA2_COMMIT_REFREEZE=0, the default) ---"
    echo "  The candidate is $R_ON/verdicts.tsv and every check above passed."
    echo "  Re-run with XTEA2_COMMIT_REFREEZE=1 to install it."
fi

echo
echo "RE-FREEZE: PASS"
echo "  finished: $(date '+%F %T')"
