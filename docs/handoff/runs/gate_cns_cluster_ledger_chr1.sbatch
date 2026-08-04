#!/bin/bash
#SBATCH -J b2_chr1
#SBATCH -A park_contrib
#SBATCH -p priopark
#SBATCH -t 2:00:00
#SBATCH -c 12
#SBATCH --mem 64G
#SBATCH -o /n/scratch/users/c/cos689/xtea2_testdata/gate_cns_cluster_ledger_chr1_%j.out
# =============================================================================
# Gate [B2] -- `check_seprt_disc_consistency`'s decision moved into
# `transduction_evidence.gate_cns_disc_cluster`, with a caller-supplied ledger. Gated on HG002
# chr1 as a PURE INSTRUMENTATION CHANGE: two arms, two commits, one command line.
#
# WHAT IS DIFFERENT FROM `gate_transduction_evidence_chr1.sbatch`, WHICH THIS IS BUILT FROM.
# There, AFTER had a module BEFORE did not, and the ledgers were NEW FILES. Here both arms have
# the module and both write ledgers, so the ledger cannot be a new file and cannot be excluded:
# it is the one artifact this change is *supposed* to alter, by exactly the rows gate B2 adds.
# So check 1 makes that precise instead of tolerating it --
#
#   **AFTER's ledger with every `B2_cns_cluster` row deleted must be BYTE-IDENTICAL to
#   BEFORE's.**
#
# which asserts both halves at once: the change added rows, and it moved nothing else in the
# file -- not a row, not a reason, not an order. Everything else under the sample directory
# must be byte-identical outright, with the same two characterised exceptions the parent gate
# established over a dry run of two arms of the SAME commit (job 49599837): `cns_manifest.txt`
# embeds each arm's own output paths, and `scan_disc_pos.tsv` is a permutation because the
# multithreaded scan writes rows in completion order.
#
# WHAT chr1 CAN AND CANNOT SAY ABOUT GATE B2 -- read this before trusting a PASS.
# chr1's gate B survivors are 5 (L1) and 0 (SVA) and its cns-consistent set is empty for both
# (job 49634289's AFTER arm). So B2's DOMAIN at chr1 is 5 loci for L1 and **nothing at all for
# SVA**: this gate proves the wiring, the vocabulary and the file-level containment, at N=5, on
# one family. The population -- 135 of L1's 159 picked sources and 49 of SVA's 105 -- exists
# only at WGS, which is `gate_cns_cluster_ledger_wgs.sbatch`. A passing chr1 gate here says
# nothing about whether the WGS counts are right.
#
# CHECK 4 IS THE ONE THAT MATTERS MOST AND IT COSTS SECONDS: `probe_decision_funnel.py` on the
# AFTER arm reconciles the ledger's B2 rows against `picked sources` minus `cluster-consistent`,
# computed from two of the step's own site files. Two independent derivations of the same
# population from different files agreeing is what makes "the gate is where the ledger says it
# is" a measurement rather than a claim.
#
# USAGE  sbatch gate_cns_cluster_ledger_chr1.sbatch <before-worktree> <after-worktree>
#   BEFORE may be any worktree whose xtea/ is behaviourally the parent --
#   /n/scratch/users/c/cos689/xtea2_wt/td_evidence_48ec4c7 qualifies (48ec4c7..7921950 is docs
#   plus one comment block), which is also why its already-built xtea_scan serves BOTH arms:
#   src/ is byte-identical across the range, so one binary removes the C++ as a variable.
# =============================================================================
set -euo pipefail

SCRATCH=/n/scratch/users/c/cos689/xtea2_testdata
WT_BEFORE="${1:?usage: sbatch gate_cns_cluster_ledger_chr1.sbatch <before-wt> <after-wt>}"
WT_AFTER="${2:?usage: sbatch gate_cns_cluster_ledger_chr1.sbatch <before-wt> <after-wt>}"
OUT="${XTEA2_OUT:-$SCRATCH/gate_cns_cluster_ledger_chr1_${SLURM_JOB_ID:-manual}}"
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
# The code under test must be present BEFORE the compute, not after it. Written as `if`, not
# `grep && { ... }`: a failing AND-list is the whole statement's status, so under `set -e` the
# negative form aborts the script on the case it is meant to accept.
if ! grep -q 'def gate_cns_disc_cluster' "$WT_AFTER/xtea/transduction_evidence.py"; then
    echo "FATAL: AFTER has no gate_cns_disc_cluster; wrong commit"; exit 1
fi
if grep -q 'def gate_cns_disc_cluster' "$WT_BEFORE/xtea/transduction_evidence.py"; then
    echo "FATAL: BEFORE already has gate_cns_disc_cluster; wrong commit"; exit 1
fi

mkdir -p "$OUT"

echo "=============================================================="
echo "GATE B2 (cns-cluster ledger)  --  HG002 chr1, types=$TYPES"
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
[ "$RC" -eq 0 ] || { echo; echo "GATE B2: FAIL (an arm did not complete)"; exit 1; }

A="$OUT/BEFORE/$SAMPLE"
B="$OUT/AFTER/$SAMPLE"

echo
echo "--- check 0: the arms are the code they claim to be ---"
# An arm that silently ran the other version would pass every comparison below. Both arms have
# the module here, so the discriminator is the gate function and the wiring at the call site.
b_gate=$(grep -c 'def gate_cns_disc_cluster' "$WT_BEFORE/xtea/transduction_evidence.py" || true)
a_gate=$(grep -c 'def gate_cns_disc_cluster' "$WT_AFTER/xtea/transduction_evidence.py" || true)
# `ledger.record(` with the paren, so the sentence about it in a comment does not count.
b_wire=$(grep -c 'ledger\.record(' "$WT_BEFORE/xtea/x_clip_disc_filter.py" || true)
a_wire=$(grep -c 'ledger\.record(' "$WT_AFTER/xtea/x_clip_disc_filter.py" || true)
echo "  gate_cns_disc_cluster defined: BEFORE=$b_gate AFTER=$a_gate (expect 0 and 1)"
echo "  x_clip_disc_filter records:    BEFORE=$b_wire AFTER=$a_wire (expect 0 and 2)"
{ [ "$b_gate" -eq 0 ] && [ "$a_gate" -eq 1 ] && [ "$b_wire" -eq 0 ] && [ "$a_wire" -eq 2 ]; } || {
    echo "    FATAL: the arms are not before/after"; RC=1; }
# BOTH arms must write both ledgers -- unlike the parent gate, where only AFTER did.
for arm in BEFORE AFTER; do
    n=$(grep -ah 'TD-LEDGER' "$OUT/$arm/$SAMPLE"/*.steps.log 2>/dev/null | wc -l || true)
    echo "  [TD-LEDGER] lines in $arm: $n (expect 2)"
    [ "$n" -eq 2 ] || { echo "    FATAL: the ledger did not run once per flanked family"; RC=1; }
done

echo
echo "--- check 1: every file both arms produce is byte-identical, ledgers modulo B2 rows ---"
n_same=0; n_diff=0; n_missing=0; n_tol=0; n_led=0
while IFS= read -r f; do
    rel=${f#"$A/"}
    case "$rel" in
        tmp/*|*/tmp/*|tmp_basic/*|pub_clip/*) continue ;;   # scratch, not results
        *.log|step_timings.tsv) continue ;;                 # wall clocks and stderr
    esac
    if [ ! -f "$B/$rel" ]; then echo "  MISSING in AFTER: $rel"; n_missing=$((n_missing+1)); continue; fi
    if cmp -s "$f" "$B/$rel"; then n_same=$((n_same+1)); continue; fi
    case "$rel" in
      *.decision_ledger.tsv)
        # THE ONE ARTIFACT THIS CHANGE MAY ALTER, and only by adding B2 rows. Delete them and
        # the file must be byte-identical, which also pins row ORDER: B2's rows are appended
        # between gate B's and gate C's, so anything that reordered the rest fails here.
        n_b2=$(awk -F'\t' 'NR>1 && $3=="B2_cns_cluster"' "$B/$rel" | wc -l)
        if cmp -s "$f" <(awk -F'\t' 'NR==1 || $3!="B2_cns_cluster"' "$B/$rel"); then
            n_led=$((n_led+1))
            echo "  $rel: identical after removing $n_b2 B2 rows  (added: $n_b2)"
        else
            echo "  DIFFERS BEYOND THE B2 ROWS: $rel  (B2 rows present: $n_b2)"
            diff "$f" <(awk -F'\t' 'NR==1 || $3!="B2_cns_cluster"' "$B/$rel") | head -20 || true
            n_diff=$((n_diff+1))
        fi ;;
      cns_manifest.txt)
        if cmp -s <(sed "s|$OUT/BEFORE/||g" "$f") <(sed "s|$OUT/AFTER/||g" "$B/$rel"); then
            n_tol=$((n_tol+1)); echo "  $rel: identical after stripping the arm prefix"
        else
            echo "  DIFFERS (beyond the path prefix): $rel"; n_diff=$((n_diff+1))
        fi ;;
      scan_disc_pos.tsv)
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
echo "  identical: $n_same    ledgers reconciled: $n_led    tolerated: $n_tol    DIFFERING: $n_diff    missing: $n_missing"
[ "$n_diff" -eq 0 ] && [ "$n_missing" -eq 0 ] || RC=1
# A run that produced almost nothing would otherwise "pass" with n_diff 0. The parent gate
# measured 98 identical + 2 characterised; two ledgers move out of the identical count here.
[ "$n_same" -ge 88 ] || { echo "    FATAL: only $n_same files compared, expected >=88"; RC=1; }
# A ledger that came out byte-identical would land in n_same and never be reconciled. That is
# possible for SVA at chr1 (empty B2 domain) but not for both.
[ "$n_led" -ge 1 ] || {
    echo "    FATAL: neither ledger gained a B2 row -- the instrumentation did not run"; RC=1; }

echo
echo "--- check 2: NO new files, and none removed ---"
n_new=0
while IFS= read -r f; do
    rel=${f#"$B/"}
    case "$rel" in
        tmp/*|*/tmp/*|tmp_basic/*|pub_clip/*|*.log|step_timings.tsv) continue ;;
    esac
    [ -f "$A/$rel" ] && continue
    echo "  UNEXPECTED new file: $rel"; n_new=$((n_new+1))
done < <(find "$B" -type f | sort)
echo "  unexpected new files: $n_new (expect 0 -- both arms already write the ledgers)"
[ "$n_new" -eq 0 ] || RC=1

echo
echo "--- check 3: the B2 rows say only things the module allows ---"
mapfile -t VOCAB < <(python - "$WT_AFTER" <<'EOF'
import sys, os
sys.path.insert(0, os.path.join(sys.argv[1], "xtea"))
import transduction_evidence as te
for r in sorted(te.PHASE_OF_REASON):
    print(r)
EOF
)
echo "  vocabulary (${#VOCAB[@]} reasons): ${VOCAB[*]}"
# 14 before this change, 18 after. A widened vocabulary is expected; a shrunk one is not.
[ "${#VOCAB[@]}" -ge 18 ] || { echo "    FATAL: expected >=18 reasons, got ${#VOCAB[@]}"; RC=1; }
# The four B2 reasons must be the module's, read from the module rather than retyped.
mapfile -t B2VOCAB < <(python - "$WT_AFTER" <<'EOF'
import sys, os
sys.path.insert(0, os.path.join(sys.argv[1], "xtea"))
import transduction_evidence as te
for r, ph in sorted(te.PHASE_OF_REASON.items()):
    if ph == "B2_cns_cluster":
        print(r)
EOF
)
echo "  B2 reasons (${#B2VOCAB[@]}): ${B2VOCAB[*]}"
[ "${#B2VOCAB[@]}" -eq 4 ] || { echo "    FATAL: expected 4 B2 reasons"; RC=1; }
while IFS= read -r led; do
    fam=$(basename "$(dirname "$led")")
    tot=$(( $(wc -l < "$led") - 1 ))
    b2=$(awk -F'\t' 'NR>1 && $3=="B2_cns_cluster"' "$led" | wc -l)
    echo "  $fam: $tot rows, of which $b2 are B2"
    [ "$(head -1 "$led")" = "#chrom	pos	phase	reason	detail" ] || {
        echo "    FATAL: unexpected header"; RC=1; }
    bad=$(awk -F'\t' 'NR>1 {print $4}' "$led" | sort -u | grep -vxF -f <(printf '%s\n' "${VOCAB[@]}") || true)
    [ -z "$bad" ] || { echo "    FATAL: reasons outside the vocabulary: $bad"; RC=1; }
    n_q=$(awk -F'\t' 'NR>1 && $3=="?"' "$led" | wc -l)
    [ "$n_q" -eq 0 ] || { echo "    FATAL: $n_q rows with an unknown phase"; RC=1; }
    # Any B2-phase row must carry one of the four B2 reasons, and no other phase may.
    xr=$(awk -F'\t' 'NR>1 && $3=="B2_cns_cluster" {print $4}' "$led" | sort -u |
         grep -vxF -f <(printf '%s\n' "${B2VOCAB[@]}") || true)
    [ -z "$xr" ] || { echo "    FATAL: B2 rows with a non-B2 reason: $xr"; RC=1; }
    xp=$(awk -F'\t' -v OFS='\t' 'NR>1 {print $3, $4}' "$led" | sort -u |
         awk -F'\t' '$1!="B2_cns_cluster"' | cut -f2 | sort -u |
         grep -xF -f <(printf '%s\n' "${B2VOCAB[@]}") || true)
    [ -z "$xp" ] || { echo "    FATAL: a B2 reason under another phase: $xp"; RC=1; }
    echo "    by phase : $(awk -F'\t' 'NR>1 {c[$3]++} END {for (p in c) printf "%s=%d ", p, c[p]}' "$led")"
    echo "    by reason: $(awk -F'\t' 'NR>1 {c[$4]++} END {for (r in c) printf "%s=%d ", r, c[r]}' "$led")"
done < <(find "$B" -name '*.decision_ledger.tsv' | sort)

echo
echo "--- check 4: the funnel probe reconciles the B2 rows against the step's own files ---"
# The payoff check. `picked sources` minus `cluster-consistent` is derived from two site FILES;
# the B2 rows come from the ledger. They must be the same set of loci. Note this check is the
# one that can only SKIP at chr1 for SVA, whose B2 domain is empty -- the probe says so.
probe_out="$OUT/probe_after.txt"
if python "$WT_AFTER/docs/handoff/runs/probe_decision_funnel.py" --run-dir "$B" \
        > "$probe_out" 2>&1; then
    sed 's/^/  /' "$probe_out"
    if grep -q 'reconciliation checks: FAILED' "$probe_out"; then
        echo "    FATAL: the probe's reconciliation failed -- the ledger and the site files"
        echo "    disagree about which loci gate B2 dropped"; RC=1
    fi
    if grep -q 'SKIPPED, not passed' "$probe_out"; then
        echo "    FATAL: the probe found no B2 rows at all in a family that had picked sources"
        RC=1
    fi
else
    echo "  the probe exited nonzero:"; sed 's/^/  /' "$probe_out"; RC=1
fi

echo
echo "--- check 5: call counts, for the record (expect 1 / 14 / 2 / 17 in both arms) ---"
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
    echo "GATE B2 (chr1): PASS"
else
    echo "GATE B2 (chr1): FAIL"
fi
echo "  finished: $(date '+%F %T')"
echo "=============================================================="
exit "$RC"
