#!/bin/bash
#SBATCH -J tdev_chr1
#SBATCH -A park_contrib
#SBATCH -p priopark
#SBATCH -t 2:00:00
#SBATCH -c 12
#SBATCH --mem 64G
#SBATCH -o /n/scratch/users/c/cos689/xtea2_testdata/gate_transduction_evidence_chr1_%j.out
# =============================================================================
# The transduction decision, moved into xtea/transduction_evidence.py -- gated on HG002 chr1
# as a PURE REFACTOR. Two arms, two commits, ONE command line: BEFORE (the parent, which
# still has the six in-place cascade methods) and AFTER (the five named gate groups).
#
# WHY THIS GATE'S SHAPE IS DIFFERENT FROM THE OTHER chr1 GATES. They compare two
# TOPOLOGIES at one commit, so they compare CALLS and accept that intermediates churn.
# This change is supposed to move nothing at all, so the gate is the strongest one
# available: **every file both arms produce must be byte-identical**, not just the VCFs.
# A dry run over an earlier gate's two arms (job 49599837) established that 98 files under
# the sample directory are already byte-identical between two arms of the SAME commit, and
# exactly two are not, for reasons that have nothing to do with any decision:
#
#   * cns_manifest.txt   -- contains each arm's own absolute output paths. Compared with
#                           the arm prefix stripped, so it is checked rather than skipped.
#   * scan_disc_pos.tsv  -- a PERMUTATION only (`sort` makes the two identical); the
#                           multithreaded scan writes its rows in completion order.
#                           Compared SORTED, which is stronger than excluding it.
#
# So this gate compares ~98 files per arm and tolerates two known, characterised
# differences. If a decision moved by one locus, an intermediate changes and this fails.
#
# THE ONE NEW FILE. AFTER writes `<out>.decision_ledger.tsv` per flanked family (L1, SVA --
# Alu has no transduction flank so the step never runs for it). Check 3 is about that file:
# a refactor that quietly produced an EMPTY ledger would pass every byte comparison here,
# because nothing reads it. Its rows are the only new information the change adds, so they
# get their own assertions -- present, non-empty, every reason in the module's vocabulary,
# and the per-phase counts echoed into the log for the handoff to quote.
#
# WHAT THIS GATE CANNOT SEE, and why a WGS arm follows it: chr1 yields 17 calls and no
# CREATED records at all, so the populations the decision actually differs on -- the 52 of
# 1,515 records that exist only through site creation, and the ~41.6K-locus universe -- are
# not represented. chr1 gates the WIRING; WGS gates the populations.
#
# USAGE  sbatch gate_transduction_evidence_chr1.sbatch <before-worktree> <after-worktree>
#   BEFORE may be any worktree whose xtea/ matches the parent commit --
#   /n/scratch/users/c/cos689/xtea2_wt/staged_barrier_d565445 qualifies (d565445..a26987e is
#   docs-only), which is also why its already-built xtea_scan serves BOTH arms: the C++ is
#   provably identical, so one binary removes it as a variable.
# =============================================================================
set -euo pipefail

SCRATCH=/n/scratch/users/c/cos689/xtea2_testdata
WT_BEFORE="${1:?usage: sbatch gate_transduction_evidence_chr1.sbatch <before-wt> <after-wt>}"
WT_AFTER="${2:?usage: sbatch gate_transduction_evidence_chr1.sbatch <before-wt> <after-wt>}"
OUT="${XTEA2_OUT:-$SCRATCH/gate_transduction_evidence_chr1_${SLURM_JOB_ID:-manual}}"
TESTDATA=/n/data1/hms/dbmi/park/corinne/claude_sandbox/test_data/HG002
BAM="$TESTDATA/chr1/HG002.chr1.bam"
SAMPLE=HG002.chr1
TYPES="${XTEA2_TYPES:-L1,Alu,SVA}"
REF=/n/data1/hms/dbmi/park-smaht_dac/ref/GRCh38_no_alt/hg38_no_alt.fa
CORES="${SLURM_CPUS_PER_TASK:-12}"
SCAN_BIN="${XTEA2_SCAN_BIN:-$WT_BEFORE/src/csrc/xtea_scan}"

source /n/app/conda/miniforge3/24.11.3-0/etc/profile.d/conda.sh
conda activate xtea2
export PATH=/home/cos689/software/bin:$PATH
export XTEA_BASIC_INFO_SEED=42
test "$PYTHONNOUSERSITE" = 1

for wt in "$WT_BEFORE" "$WT_AFTER"; do
    [ -d "$wt" ] || { echo "FATAL: worktree $wt missing"; exit 1; }
    [ -z "$(git -C "$wt" status --porcelain)" ] || { echo "FATAL: $wt is dirty"; exit 1; }
done
[ -x "$SCAN_BIN" ] || { echo "FATAL: $SCAN_BIN not built"; exit 1; }
[ -f "$BAM" ]      || { echo "FATAL: $BAM missing"; exit 1; }

mkdir -p "$OUT"

echo "=============================================================="
echo "TRANSDUCTION EVIDENCE GATE  --  HG002 chr1, types=$TYPES"
echo "  job     : ${SLURM_JOB_ID:-<none>} on $(hostname)  cores=$CORES"
echo "  BEFORE  : $(git -C "$WT_BEFORE" rev-parse --short HEAD)  $WT_BEFORE"
echo "  AFTER   : $(git -C "$WT_AFTER"  rev-parse --short HEAD)  $WT_AFTER"
echo "  scan bin: $SCAN_BIN (shared by both arms)"
echo "  out     : $OUT"
echo "  started : $(date '+%F %T')"
echo "=============================================================="

run_arm() {   # $1 = arm name, $2 = worktree
    # Separate statements, not one `local a=$1 b=...` line: under `set -u` the whole line
    # expands before `local` assigns (jobs 49452897 / 49516776).
    local arm=$1
    local wt=$2
    local t0 wall armout
    armout="$OUT/$arm"
    rm -rf "$armout"; mkdir -p "$armout"
    echo
    echo "=== arm $arm ($(git -C "$wt" rev-parse --short HEAD)) ==="
    t0=$(date +%s)
    ( cd "$wt" && python xtea/run_xtea.py \
      -i "$BAM" -s "$SAMPLE" -o "$armout" \
      --reference "$REF" \
      --rep-lib-dir /home/cos689/xTea/rep_lib_annotation \
      --gene-gff3 /home/cos689/data1/corinne/ref/gencode.v45.basic.annotation.gff3 \
      --blacklist /home/cos689/xTea/rep_lib_annotation/blacklist/hg38/centromere.bed \
      --types "$TYPES" \
      --scan-bin "$SCAN_BIN" \
      --bwa /home/cos689/software/bin/bwa \
      --scan-min-clip 2 \
      --mate-in-rep-families alu,sva \
      --cns-emit-scan \
      --td-gather-scan \
      --parallel-downstream \
      --keep-intermediates \
      -n "$CORES" ) > "$armout/driver.log" 2>&1 || {
          echo "FAIL: arm $arm exited $?"; tail -60 "$armout/driver.log"; return 1; }
    wall=$(( $(date +%s) - t0 ))
    printf '  wall %02d:%02d:%02d\n' $((wall/3600)) $((wall%3600/60)) $((wall%60))
}

RC=0
run_arm BEFORE "$WT_BEFORE" || RC=1
run_arm AFTER  "$WT_AFTER"  || RC=1
[ "$RC" -eq 0 ] || { echo; echo "TRANSDUCTION EVIDENCE GATE: FAIL (an arm did not complete)"; exit 1; }

echo
echo "--- check 0: the arms are the code they claim to be ---"
# The lesson from the barrier gate's check 0: an arm that silently ran the OTHER version
# would pass every comparison below. Prove the module is present on one side and absent on
# the other, and that the ledger line only appears where the module exists.
b_mod=$([ -f "$WT_BEFORE/xtea/transduction_evidence.py" ] && echo y || echo n)
a_mod=$([ -f "$WT_AFTER/xtea/transduction_evidence.py" ] && echo y || echo n)
echo "  transduction_evidence.py: BEFORE=$b_mod AFTER=$a_mod"
[ "$b_mod" = n ] && [ "$a_mod" = y ] || { echo "    FATAL: the arms are not before/after"; RC=1; }
n_led_before=$(grep -ah 'TD-LEDGER' "$OUT/BEFORE/$SAMPLE"/*.steps.log 2>/dev/null | wc -l || true)
n_led_after=$(grep -ah 'TD-LEDGER' "$OUT/AFTER/$SAMPLE"/*.steps.log 2>/dev/null | wc -l || true)
echo "  [TD-LEDGER] lines: BEFORE=$n_led_before AFTER=$n_led_after (expect 0 and 2)"
[ "$n_led_before" -eq 0 ] && [ "$n_led_after" -eq 2 ] || {
    echo "    FATAL: the ledger did not run exactly once per flanked family"; RC=1; }

echo
echo "--- check 1: EVERY file both arms produce is byte-identical ---"
A="$OUT/BEFORE/$SAMPLE"
B="$OUT/AFTER/$SAMPLE"
n_same=0; n_diff=0; n_missing=0; n_tol=0
while IFS= read -r f; do
    rel=${f#"$A/"}
    case "$rel" in
        tmp/*|*/tmp/*|tmp_basic/*|pub_clip/*) continue ;;   # scratch, not results
        *.log|step_timings.tsv) continue ;;                 # wall clocks and stderr
    esac
    if [ ! -f "$B/$rel" ]; then echo "  MISSING in AFTER: $rel"; n_missing=$((n_missing+1)); continue; fi
    if cmp -s "$f" "$B/$rel"; then n_same=$((n_same+1)); continue; fi
    case "$rel" in
      cns_manifest.txt)
        # Differs only by each arm's own output prefix; strip it and compare.
        if cmp -s <(sed "s|$OUT/BEFORE/||g" "$f") <(sed "s|$OUT/AFTER/||g" "$B/$rel"); then
            n_tol=$((n_tol+1)); echo "  $rel: identical after stripping the arm prefix"
        else
            echo "  DIFFERS (beyond the path prefix): $rel"; n_diff=$((n_diff+1))
        fi ;;
      scan_disc_pos.tsv)
        # A permutation only -- the scan writes rows in thread-completion order.
        if cmp -s <(sort "$f") <(sort "$B/$rel"); then
            n_tol=$((n_tol+1)); echo "  $rel: identical when sorted ($(wc -l < "$f") rows)"
        else
            echo "  DIFFERS (not a permutation): $rel"; n_diff=$((n_diff+1))
        fi ;;
      *)
        echo "  DIFFERS: $rel"; n_diff=$((n_diff+1))
        diff <(grep -v '^#' "$f" | cut -f1,2) <(grep -v '^#' "$B/$rel" | cut -f1,2) | head -20 || true ;;
    esac
done < <(find "$A" -type f | sort)
echo "  identical: $n_same    tolerated (characterised): $n_tol    DIFFERING: $n_diff    missing: $n_missing"
[ "$n_diff" -eq 0 ] && [ "$n_missing" -eq 0 ] || RC=1
# A run that produced almost nothing would otherwise "pass" with n_diff 0.
[ "$n_same" -ge 90 ] || { echo "    FATAL: only $n_same files compared, expected >=90"; RC=1; }

echo
echo "--- check 2: the only NEW files in AFTER are the ledgers ---"
n_new=0; n_new_led=0
while IFS= read -r f; do
    rel=${f#"$B/"}
    case "$rel" in
        tmp/*|*/tmp/*|tmp_basic/*|pub_clip/*|*.log|step_timings.tsv) continue ;;
    esac
    [ -f "$A/$rel" ] && continue
    case "$rel" in
        *.decision_ledger.tsv) n_new_led=$((n_new_led+1)) ;;
        *) echo "  UNEXPECTED new file: $rel"; n_new=$((n_new+1)) ;;
    esac
done < <(find "$B" -type f | sort)
echo "  new ledgers: $n_new_led    unexpected new files: $n_new"
[ "$n_new" -eq 0 ] && [ "$n_new_led" -eq 2 ] || RC=1

echo
echo "--- check 3: the ledger says something, and only things it is allowed to say ---"
# The reason vocabulary, from the module itself rather than retyped here -- a renamed reason
# must fail this, not silently widen it.
mapfile -t VOCAB < <(python - "$WT_AFTER" <<'EOF'
import sys, os
sys.path.insert(0, os.path.join(sys.argv[1], "xtea"))
import transduction_evidence as te
for r in sorted(te.PHASE_OF_REASON):
    print(r)
EOF
)
echo "  vocabulary (${#VOCAB[@]} reasons): ${VOCAB[*]}"
[ "${#VOCAB[@]}" -ge 14 ] || { echo "    FATAL: could not read the reason vocabulary"; RC=1; }
n_led=0
while IFS= read -r led; do
    n_led=$((n_led+1))
    rows=$(( $(wc -l < "$led") - 1 ))
    hdr=$(head -1 "$led")
    echo "  $(basename "$(dirname "$led")")/$(basename "$led"): $rows rows"
    [ "$hdr" = "#chrom	pos	phase	reason	detail" ] || {
        echo "    FATAL: unexpected header: $hdr"; RC=1; }
    [ "$rows" -gt 0 ] || { echo "    FATAL: empty ledger -- the gates recorded nothing"; RC=1; }
    # every reason in the file must be in the vocabulary, and no phase may be "?"
    bad=$(awk -F'\t' 'NR>1 {print $4}' "$led" | sort -u | grep -vxF -f <(printf '%s\n' "${VOCAB[@]}") || true)
    [ -z "$bad" ] || { echo "    FATAL: reasons outside the vocabulary: $bad"; RC=1; }
    n_q=$(awk -F'\t' 'NR>1 && $3=="?"' "$led" | wc -l)
    [ "$n_q" -eq 0 ] || { echo "    FATAL: $n_q rows with an unknown phase"; RC=1; }
    echo "    by phase:  $(awk -F'\t' 'NR>1 {c[$3]++} END {for (p in c) printf "%s=%d ", p, c[p]}' "$led" | tr -d '\n')"
    echo "    by reason: $(awk -F'\t' 'NR>1 {c[$4]++} END {for (r in c) printf "%s=%d ", r, c[r]}' "$led")"
done < <(find "$B" -name '*.decision_ledger.tsv' | sort)
[ "$n_led" -eq 2 ] || { echo "    FATAL: found $n_led ledgers, expected 2"; RC=1; }

echo
echo "--- check 4: call counts, for the record ---"
for arm in BEFORE AFTER; do
    l1=$(grep -vc '^#' "$OUT/$arm/$SAMPLE/merged/${SAMPLE}_LINE1.vcf" 2>/dev/null || echo 0)
    alu=$(grep -vc '^#' "$OUT/$arm/$SAMPLE/merged/${SAMPLE}_ALU.vcf" 2>/dev/null || echo 0)
    sva=$(grep -vc '^#' "$OUT/$arm/$SAMPLE/merged/${SAMPLE}_SVA.vcf" 2>/dev/null || echo 0)
    mrg=$(grep -vc '^#' "$OUT/$arm/$SAMPLE/${SAMPLE}.vcf" 2>/dev/null || echo 0)
    echo "  $arm: L1=$l1 Alu=$alu SVA=$sva merged=$mrg"
done

echo
echo "=============================================================="
if [ "$RC" -eq 0 ]; then
    echo "TRANSDUCTION EVIDENCE GATE: PASS"
else
    echo "TRANSDUCTION EVIDENCE GATE: FAIL"
fi
echo "  finished: $(date '+%F %T')"
echo "=============================================================="
exit "$RC"
