#!/bin/bash
#SBATCH -c 24
#SBATCH -t 8:00:00
#SBATCH --mem=64G
#SBATCH -p priopark
#SBATCH --account=park_contrib
#SBATCH -o /n/scratch/users/c/cos689/xtea2_testdata/b10_validate_%j.out
#SBATCH --mail-type=NONE

# =============================================================================
# B10 validation: does fixing the disc-mate double-count make the FAST cns path
# equivalent to the slow C++ fan, so the 31 minutes can be banked?
#
# Background. Two completed WGS runs differed by 23 call positions and nobody knew which
# was right:
#   48825976  step45_wgs_out  --cns-emit-scan  (C++ fan, whole-genome mate-gather)  1:31:21
#   48842713  fix_nofan_out   no flag          (Python indexed per-site fetch)        59:56
# Static + record-level comparison of those two runs (2026-07-28) showed they collect the
# SAME READS -- clip set-identical for all three types, disc key-identical -- and that the
# only substantive difference was that the indexed path emitted 4.8% of Alu disc records
# TWICE or THRICE (bug B10: one FASTA record per mate ALIGNMENT rather than per mate READ).
# Commit ee0742d fixes that. This run is the fix, on the FAST path.
#
# THE TWO THINGS THIS RUN DECIDES
#   1. Is the duplication gone, and does the fast path's disc record set now match the fan's?
#      That is asserted directly on temp_disc.sam, below, and it is the sharp test -- the key
#      excludes bwa's placement, so it is unaffected by the ordering noise in (2).
#   2. Does the 31-minute win survive the fix? (It should: the fix only adds two guards
#      inside a loop that already runs.)
#
# WHAT THIS RUN CANNOT DECIDE, so it is reported and never asserted: exact call-set equality
# with the fan run. The two paths feed bwa the same clip parts in a DIFFERENT ORDER, and
# bwa's tie-break among equal-scoring consensus placements is order-dependent (measured:
# AS=23 for three different placements of one 33bp query). That reshuffles which ALU_vN a
# multi-mapping clip part lands on and tips borderline candidates. Expect the churn vs the
# fan to SHRINK a lot from 23 positions, not to reach zero.
#
# Runs from a frozen git worktree, not the live checkout: run_xtea.py spawns a fresh
# x_TEA_main.py per step, so a job sharing the checkout would pick up mid-run edits.
# =============================================================================
set -euo pipefail

WT="${XTEA2_WORKTREE:-/n/scratch/users/c/cos689/xtea2_testdata/b10_wt}"
OUT="${XTEA2_OUT:-/n/scratch/users/c/cos689/xtea2_testdata/b10_out}"
SCRATCH=/n/scratch/users/c/cos689/xtea2_testdata
FAN="$SCRATCH/step45_wgs_out/NA12878.wgs"        # 48825976, --cns-emit-scan
OLDIDX="$SCRATCH/fix_nofan_out/NA12878.wgs"      # 48842713, indexed WITHOUT the B10 fix

source /n/app/conda/miniforge3/24.11.3-0/etc/profile.d/conda.sh
conda activate xtea2
export PATH=/home/cos689/software/bin:$PATH
export XTEA_BASIC_INFO_SEED=42   # pin the insert-size draw (matches the seeded golden)

[ -d "$WT" ] || { echo "FATAL: worktree $WT missing"; exit 1; }
[ -x "$WT/src/csrc/xtea_scan" ] || { echo "FATAL: $WT/src/csrc/xtea_scan missing or not executable"; exit 1; }

rm -rf "$OUT"; mkdir -p "$OUT"
cd "$WT"

echo "=============================================================="
echo "xTea2 WGS B10 VALIDATION: fast indexed cns path + disc-mate dedup fix"
echo "  job          : ${SLURM_JOB_ID:-<none>}  on $(hostname)"
echo "  worktree     : $WT"
echo "  commit       : $(git rev-parse --short HEAD) $(git log -1 --format=%s | cut -c1-70)"
echo "  worktree diff: $(git status --porcelain | wc -l) file(s) (expect 0)"
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
S="$OUT/NA12878.wgs"

echo
echo "##############################################################"
echo "# 1. TIMING -- does the 31-minute win survive the fix?"
echo "##############################################################"
printf 'this run: %02d:%02d:%02d wall (%ds)\n' $((WALL/3600)) $((WALL%3600/60)) $((WALL%60)) "$WALL"
grep -E "Percent of CPU|Elapsed \(wall clock\)|Maximum resident" "$OUT/time_v.txt" 2>/dev/null || true
cat <<'REF'

--- the three comparators (same BAM, same node class, one job, -c 24 --mem=64G) ----
  48825977  master 46520ab, serial -n 24 (the honest denominator)   50:20   40.6 GB
  48825976  xTea2 --cns-emit-scan (C++ fan)                       1:31:21   24.3 GB
  48842713  xTea2 indexed, NO B10 fix                               59:56   22.9 GB
  <this>    xTea2 indexed, WITH B10 fix -- expect ~= 59:56
If this lands near 1:00 the win holds and the fan's whole-genome mate-gather can go.
REF

echo
echo "--- per-phase / per-step wall time ---"
if [ -f "$S/step_timings.tsv" ]; then
  grep -E "^phase" "$S/step_timings.tsv" || true
  echo "  slowest steps:"
  awk -F'\t' '$1 ~ /^step/ {printf "    %-32s %8.2f min\n", $2, $3/60}' "$S/step_timings.tsv" \
    | sort -k2,2 -rn -t' ' | head -8 || true
else
  echo "  MISSING $S/step_timings.tsv"
fi

echo
echo "##############################################################"
echo "# 2. THE B10 ASSERTION -- disc records, one per mate read"
echo "##############################################################"
echo "Key = qname~is_first~lclip~rclip~anchor_rc~mate_rc~anchor_pos~ins_chrm~ins_pos over"
echo "PRIMARY records only. It excludes bwa's placement, so it isolates collection."
echo
B10_FAIL=0
printf '%-5s %12s %12s %12s %12s\n' type keys dup_keys dup_recs "vs_fan_keys"
for t in L1 Alu SVA; do
  f="$S/$t/tmp/cns/temp_disc.sam"
  g="$FAN/$t/tmp/cns/temp_disc.sam"
  if [ ! -f "$f" ] || [ ! -f "$g" ]; then
    printf '%-5s %12s\n' "$t" "MISSING"; continue
  fi
  key='!/^@/ && and($2,0x900)==0 {n=split($1,q,"~"); if(n>=10) print q[1]"~"q[2]"~"q[3]"~"q[4]"~"q[5]"~"q[6]"~"q[7]"~"q[8]"~"q[9]}'
  awk -F'\t' "$key" "$f" | sort > "$OUT/.k.this.$t"
  awk -F'\t' "$key" "$g" | sort > "$OUT/.k.fan.$t"
  nkeys=$(sort -u "$OUT/.k.this.$t" | wc -l)
  dupk=$(uniq -c "$OUT/.k.this.$t" | awk '$1>1' | wc -l)
  dupr=$(uniq -c "$OUT/.k.this.$t" | awk '$1>1{s+=$1-1} END{print s+0}')
  onlythis=$(comm -23 <(sort -u "$OUT/.k.this.$t") <(sort -u "$OUT/.k.fan.$t") | wc -l)
  onlyfan=$(comm -13 <(sort -u "$OUT/.k.this.$t") <(sort -u "$OUT/.k.fan.$t") | wc -l)
  printf '%-5s %12s %12s %12s %12s\n' "$t" "$nkeys" "$dupk" "$dupr" "+$onlythis/-$onlyfan"
  [ "$dupk" -eq 0 ] || B10_FAIL=1
  [ "$onlythis" -eq 0 ] && [ "$onlyfan" -eq 0 ] || B10_FAIL=1
  rm -f "$OUT/.k.this.$t" "$OUT/.k.fan.$t"
done
echo
if [ "$B10_FAIL" -eq 0 ]; then
  echo "  ==> PASS: no duplicate disc records, and the disc record set matches the C++ fan"
  echo "      exactly for every type. The two cns collection paths are now equivalent, so the"
  echo "      fast one is the one to keep."
else
  echo "  ==> FAIL: see the table. dup_keys>0 means the fix did not take; a nonzero vs_fan_keys"
  echo "      delta means the two paths still disagree about WHICH reads to collect, which the"
  echo "      static analysis said they do not -- that would be a new finding, not a regression."
fi
echo
echo "--- for contrast, the SAME table for 48842713 (indexed, pre-fix) ---"
for t in L1 Alu SVA; do
  f="$OLDIDX/$t/tmp/cns/temp_disc.sam"
  [ -f "$f" ] || { echo "  $t: MISSING (output may have been cleaned up)"; continue; }
  awk -F'\t' '!/^@/ && and($2,0x900)==0 {n=split($1,q,"~"); if(n>=10) print q[1]"~"q[2]"~"q[3]"~"q[4]"~"q[5]"~"q[6]"~"q[7]"~"q[8]"~"q[9]}' "$f" \
    | sort | uniq -c | awk -v T="$t" '{k++; if($1>1){d++; r+=$1-1}} END{printf "  %-5s keys=%d dup_keys=%d dup_recs=%d\n", T, k, d+0, r+0}'
done

echo
echo "##############################################################"
echo "# 3. CLIP records -- must be unchanged (B10 touches disc only)"
echo "##############################################################"
for t in L1 Alu SVA; do
  f="$S/$t/tmp/cns/temp_clip.sam"; g="$FAN/$t/tmp/cns/temp_clip.sam"
  if [ -f "$f" ] && [ -f "$g" ]; then
    a=$(grep -avc '^@' "$f" || true); b=$(grep -avc '^@' "$g" || true)
    echo "  $t: this=$a  fan=$b   (small delta expected: bwa emits supplementary lines whose"
    echo "        count follows its order-dependent tie-break; primary counts should match)"
  fi
done

echo
echo "##############################################################"
echo "# 4. CALL SET -- reported, NOT gated"
echo "##############################################################"
echo "vs the C++ fan run (48825976). Prediction: much less than the 23 positions that separated"
echo "the fan from the PRE-fix indexed run, with the residue being bwa tie-break noise at"
echo "borderline sites. A LARGER churn would mean the double-count was not the driver."
echo
for pair in "L1:LINE1" "Alu:ALU" "SVA:SVA"; do
  t="${pair%%:*}"; alt="${pair##*:}"
  a=$(ls "$S"/merged/*_${alt}.vcf 2>/dev/null | head -1)
  for cmp_name in fan oldidx; do
    root=$([ "$cmp_name" = fan ] && echo "$FAN" || echo "$OLDIDX")
    b=$(ls "$root"/merged/*_${alt}.vcf 2>/dev/null | head -1)
    if [ -n "$a" ] && [ -n "$b" ]; then
      na=$(grep -vc '^#' "$a"); nb=$(grep -vc '^#' "$b")
      lost=$(comm -13 <(grep -v '^#' "$a" | cut -f1,2 | sort -u) <(grep -v '^#' "$b" | cut -f1,2 | sort -u) | wc -l)
      gain=$(comm -23 <(grep -v '^#' "$a" | cut -f1,2 | sort -u) <(grep -v '^#' "$b" | cut -f1,2 | sort -u) | wc -l)
      printf '  %-4s vs %-7s this=%-6s that=%-6s  +%s / -%s\n' "$t" "$cmp_name" "$na" "$nb" "$gain" "$lost"
    fi
  done
done
echo
echo "--- merged VCF and per-type counts vs the seeded golden (context only) ---"
mf=$(ls "$S"/*.vcf 2>/dev/null | head -1)
[ -n "$mf" ] && echo "  merged: $(grep -vc '^#' "$mf") calls" || echo "  MISSING merged VCF"
printf '  %-6s %8s %8s %8s\n' type this golden delta
for pair in "L1:LINE1" "Alu:ALU" "SVA:SVA"; do
  t="${pair%%:*}"; alt="${pair##*:}"
  f=$(ls "$S"/merged/*_${alt}.vcf 2>/dev/null | head -1)
  g=$(ls "$SCRATCH/master_golden_wgs_seeded/NA12878/$t"/*.vcf 2>/dev/null | head -1)
  n=$([ -n "$f" ] && grep -vc '^#' "$f" || echo NA)
  m=$([ -n "$g" ] && grep -vc '^#' "$g" || echo NA)
  d=NA; [ "$n" != NA ] && [ "$m" != NA ] && d=$((n - m))
  printf '  %-6s %8s %8s %8s\n' "$t" "$n" "$m" "$d"
done
echo
echo "  The goldens were produced WITH the double count, so they are expected to sit ABOVE"
echo "  this run. That is the re-baselining debt, not a regression -- do not read these as a"
echo "  recall loss without checking the dropped sites' disc support first."

echo
echo "=============================================================="
echo "finished $(date '+%F %T')"
printf 'total wall %02d:%02d:%02d\n' $((WALL/3600)) $((WALL%3600/60)) $((WALL%60))
echo "B10 assertion: $([ "$B10_FAIL" -eq 0 ] && echo PASS || echo FAIL)"
echo "=============================================================="
exit "$B10_FAIL"
