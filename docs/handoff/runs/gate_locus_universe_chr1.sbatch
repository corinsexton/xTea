#!/bin/bash
#SBATCH -J lu_chr1
#SBATCH -A park_contrib
#SBATCH -p priopark
#SBATCH -t 2:00:00
#SBATCH -c 12
#SBATCH --mem 64G
#SBATCH -o /n/scratch/users/c/cos689/xtea2_testdata/gate_locus_universe_chr1_%j.out
# =============================================================================
# --locus-universe, end to end on HG002 chr1. TWO arms at ONE commit.
#
# WHAT THIS GATES, and why the offline gate is not enough. `gate_locus_universe.py`
# proves the driver's enumeration reproduces a finished run's `.slct` byte for byte at
# WGS scale, but it runs the module directly against kept files. It cannot see the
# WIRING, and the wiring is the part that can move calls:
#
#   - `--locus-universe` opens the SAME cns_filter -> transduction barrier
#     `--family-assign-early` uses, so the per-family chain runs as TWO x_TEA_main
#     invocations instead of one (build_type_steps start_from/stop_after). A step that
#     opens an input it was never handed dies only on the real chain -- that is exactly
#     how job 49498881 failed (`SVA:transduction`, bare `KeyError`, because
#     x_TEA_main.py:551 derives `.high_confident` from its own `-i`).
#   - unlike `--family-assign-early`, this flag does NOT get the single-family
#     exemption, so the split can happen where it never has before.
#
# THE TWO CHECKS:
#   1. ON vs OFF -- every VCF, cns2 and sibling list byte-identical. The flag is
#      diagnostic; a moved call is a bug, not a finding.
#   2. per flanked family, the driver's replica vs the step's own file:
#        merged/locus_universe/<T>.slct  ==  <T>/candidate_list_from_disc.txt.clip_sites_raw_disc.txt.slct
#      Alu has no step `.slct` (no transduction flank), so the expected comparison
#      count is the FLANKED family count, and it is printed -- job 49516776's lesson:
#      a gate that only counts at comparison time cannot tell "wrong" from "compared
#      nothing".
#
# Native gather only (the default since job 49480961). This change touches no gather
# code, and the legacy leg is covered by refactor_gate_chr1.sbatch.
#
# USAGE  sbatch gate_locus_universe_chr1.sbatch <worktree>
# =============================================================================
set -euo pipefail

SCRATCH=/n/scratch/users/c/cos689/xtea2_testdata
WT="${1:?usage: sbatch gate_locus_universe_chr1.sbatch <worktree>}"
OUT="${XTEA2_OUT:-$SCRATCH/gate_locus_universe_chr1_${SLURM_JOB_ID:-manual}}"
TESTDATA=/n/data1/hms/dbmi/park/corinne/claude_sandbox/test_data/HG002
BAM="$TESTDATA/chr1/HG002.chr1.bam"
SAMPLE=HG002.chr1
TYPES="${XTEA2_TYPES:-L1,Alu,SVA}"
REF=/n/data1/hms/dbmi/park-smaht_dac/ref/GRCh38_no_alt/hg38_no_alt.fa
CORES="${SLURM_CPUS_PER_TASK:-12}"

source /n/app/conda/miniforge3/24.11.3-0/etc/profile.d/conda.sh
conda activate xtea2
export PATH=/home/cos689/software/bin:$PATH
export XTEA_BASIC_INFO_SEED=42
test "$PYTHONNOUSERSITE" = 1

[ -d "$WT" ]                    || { echo "FATAL: worktree $WT missing"; exit 1; }
[ -x "$WT/src/csrc/xtea_scan" ] || { echo "FATAL: $WT/src/csrc/xtea_scan not built"; exit 1; }
[ -f "$BAM" ]                   || { echo "FATAL: $BAM missing"; exit 1; }
[ -z "$(git -C "$WT" status --porcelain)" ] || { echo "FATAL: worktree is dirty"; exit 1; }

mkdir -p "$OUT"
cd "$WT"

echo "=============================================================="
echo "LOCUS UNIVERSE GATE  --  HG002 chr1, types=$TYPES"
echo "  job    : ${SLURM_JOB_ID:-<none>} on $(hostname)  cores=$CORES"
echo "  commit : $(git rev-parse --short HEAD) $(git log -1 --format=%s | cut -c1-60)"
echo "  out    : $OUT"
echo "  started: $(date '+%F %T')"
echo "=============================================================="

run_arm() {   # $1 = arm name, $2... = extra driver flags
    # Separate statements, not `local arm=$1 armout=...`: under `set -u` the whole line
    # expands before `local` assigns and $arm is still unset (job 49452897 / 49516776).
    local arm=$1; shift
    local t0 wall armout
    armout="$OUT/$arm"
    rm -rf "$armout"; mkdir -p "$armout"
    echo
    echo "=== arm $arm: extra flags: $* ==="
    t0=$(date +%s)
    python xtea/run_xtea.py \
      -i "$BAM" -s "$SAMPLE" -o "$armout" \
      --reference "$REF" \
      --rep-lib-dir /home/cos689/xTea/rep_lib_annotation \
      --gene-gff3 /home/cos689/data1/corinne/ref/gencode.v45.basic.annotation.gff3 \
      --blacklist /home/cos689/xTea/rep_lib_annotation/blacklist/hg38/centromere.bed \
      --types "$TYPES" \
      --scan-bin "$WT/src/csrc/xtea_scan" \
      --bwa /home/cos689/software/bin/bwa \
      --scan-min-clip 2 \
      --mate-in-rep-families alu,sva \
      --cns-emit-scan \
      --td-gather-scan \
      --parallel-downstream \
      --keep-intermediates \
      -n "$CORES" \
      "$@" \
      > "$armout/driver.log" 2>&1 || {
          echo "FAIL: arm $arm exited $?"; tail -60 "$armout/driver.log"; return 1; }
    wall=$(( $(date +%s) - t0 ))
    printf '  wall %02d:%02d:%02d\n' $((wall/3600)) $((wall%3600/60)) $((wall%60))
}

RC=0
run_arm OFF || RC=1
run_arm ON --locus-universe || RC=1
[ "$RC" -eq 0 ] || { echo; echo "LOCUS UNIVERSE GATE: FAIL (an arm did not complete)"; exit 1; }

echo
echo "--- the arms ran the topologies they claim ---"
# OFF is one wave: exactly one x_TEA_main --transduction per family. ON is two waves, so
# the chain is split -- cns_filter ends wave 1 and transduction starts wave 2. Checked from
# the driver log's own barrier line, which only the ON arm can print.
n_off=$(grep -c 'locus universe (cns_filter -> transduction' "$OUT/OFF/driver.log" || true)
n_on=$(grep -c 'locus universe (cns_filter -> transduction' "$OUT/ON/driver.log" || true)
echo "  barrier lines: OFF=$n_off  ON=$n_on"
[ "$n_off" -eq 0 ] && [ "$n_on" -eq 1 ] || {
    echo "  FATAL: arms did not run different topologies"; RC=1; }

echo
echo "--- check 1: ON vs OFF, no call may move ---"
n_same=0; n_diff=0; n_new=0
while IFS= read -r a; do
    rel=${a#"$OUT/OFF/"}
    b="$OUT/ON/$rel"
    if [ ! -f "$b" ]; then echo "  MISSING in ON: $rel"; n_diff=$((n_diff+1)); continue; fi
    if cmp -s "$a" "$b"; then n_same=$((n_same+1)); else
        echo "  DIFFERS: $rel"; n_diff=$((n_diff+1))
        diff <(grep -v '^#' "$a" | cut -f1,2) <(grep -v '^#' "$b" | cut -f1,2) | head -20 || true
    fi
done < <(find "$OUT/OFF" \( -name '*.vcf' -o -name 'candidate_disc_filtered_cns2.txt' \
                            -o -name 'candidate_sibling_transduction2.txt' \) -type f | sort)
# A file ON produces that OFF does not is expected ONLY under merged/locus_universe/.
while IFS= read -r b; do
    rel=${b#"$OUT/ON/"}
    case "$rel" in */merged/locus_universe/*) continue ;; esac
    [ -f "$OUT/OFF/$rel" ] || { echo "  NEW file: $rel"; n_new=$((n_new+1)); }
done < <(find "$OUT/ON" \( -name '*.vcf' -o -name 'candidate_disc_filtered_cns2.txt' \
                           -o -name 'candidate_sibling_transduction2.txt' \) -type f | sort)
echo "  identical: $n_same    differing: $n_diff    unexpected new: $n_new"
[ "$n_diff" -eq 0 ] && [ "$n_new" -eq 0 ] || RC=1
[ "$n_same" -gt 0 ] || { echo "  FATAL: nothing was compared"; RC=1; }

echo
echo "--- check 2: the driver's replica == the step's own site list ---"
S="$OUT/ON/$SAMPLE"
# Flanked families only: Alu has no transduction flank, so its step writes no .slct.
FLANKED=""
for t in ${TYPES//,/ }; do
    case "$t" in L1|SVA) FLANKED="$FLANKED $t" ;; esac
done
n_exp=$(echo $FLANKED | wc -w)
n_cmp=0; n_ident=0
for t in $FLANKED; do
    step="$S/$t/candidate_list_from_disc.txt.clip_sites_raw_disc.txt.slct"
    drv="$S/merged/locus_universe/$t.slct"
    if [ ! -f "$step" ]; then echo "  $t: FATAL: step .slct missing ($step)"; continue; fi
    if [ ! -f "$drv" ]; then echo "  $t: FATAL: driver replica missing ($drv)"; continue; fi
    n_cmp=$((n_cmp+1))
    if cmp -s "$step" "$drv"; then
        n_ident=$((n_ident+1))
        echo "  $t: IDENTICAL ($(wc -c < "$step") bytes, $(wc -l < "$step") rows)"
    else
        echo "  $t: DIFFERS"
        diff <(cat "$step") <(cat "$drv") | head -20 || true
    fi
done
echo "  products compared: $n_cmp (expected $n_exp)   identical: $n_ident"
[ "$n_cmp" -eq "$n_exp" ] && [ "$n_ident" -eq "$n_exp" ] && [ "$n_exp" -gt 0 ] || RC=1
# Alu must NOT have a step .slct -- if it grows one, the flankless assumption is wrong.
if [ -f "$S/Alu/candidate_list_from_disc.txt.clip_sites_raw_disc.txt.slct" ]; then
    echo "  FATAL: Alu grew a step .slct; the flankless assumption is wrong"; RC=1
fi

echo
echo "--- the universe, as the run reported it ---"
grep -A 8 'locus universe (cns_filter' "$OUT/ON/driver.log" | sed 's/^/  /' || true
echo "  table rows: $(grep -vc '^#' "$S/merged/locus_universe/locus_universe.tsv" || echo 0)"

echo
echo "--- call counts ---"
n_records() { [ -f "$1" ] || { echo 0; return 0; }; grep -vc '^#' "$1" || true; }
first_match() { local m; for m in $1; do [ -e "$m" ] && { printf '%s' "$m"; return 0; }; done; return 0; }
for arm in OFF ON; do
    d="$OUT/$arm/$SAMPLE"; line="  $arm:"
    for pair in "L1:LINE1" "Alu:ALU" "SVA:SVA"; do
        f=$(first_match "$d/merged/*_${pair##*:}.vcf")
        [ -n "$f" ] && line="$line ${pair%%:*}=$(n_records "$f")"
    done
    mf=$(first_match "$d/*.vcf"); [ -n "$mf" ] && line="$line merged=$(n_records "$mf")"
    echo "$line"
done

echo
if [ "$RC" -eq 0 ]; then echo "LOCUS UNIVERSE GATE: PASS"; else echo "LOCUS UNIVERSE GATE: FAIL"; fi
echo "finished $(date '+%F %T')"
exit "$RC"
