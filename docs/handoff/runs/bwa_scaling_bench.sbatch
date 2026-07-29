#!/bin/bash
#SBATCH -c 24
#SBATCH -t 2:00:00
#SBATCH --mem=32G
#SBATCH -p priopark
#SBATCH --account=park_contrib
#SBATCH -o /n/scratch/users/c/cos689/xtea2_testdata/bwa_scaling_%j.out
#SBATCH --mail-type=NONE

# =============================================================================
# How does bwa mem scale on xTea2's realign workload, and is `threads/nfamilies`
# the right shape?
#
# WHY THIS EXISTS. realign.hpp splits the core budget across the three families
# (bwa_threads_for: 24 cores -> 8 threads each, 3 concurrent bwa processes). That shape
# assumes bwa mem scales well to 8 threads on THIS input, which is not ordinary bwa
# input: the reads are soft-clip fragments (13-70 bp), the reference is a tiny MEI
# consensus multi-FASTA, and the params are extremely permissive (-T 7 -k 7 -c 70), so
# each read produces ~4.5 alignment records. Nobody has measured whether that scales.
#
# The realign phase is ~18-20 min of a ~60 min WGS run and the newly added
# instrumentation shows the parent's own parse costs ~0.1-0.3 s per family, so
# essentially all of it is bwa. If bwa saturates below 8 threads, the fix is more bwa
# PROCESSES (shard the FASTQ), and part 3 measures exactly that -- including whether
# sharding is output-safe.
#
# This MUST run as a batch job: the interactive Claude session holds only 2 cores, which
# silently caps every measurement at ~1.9x and looks exactly like a serial bottleneck in
# bwa. (That is how the first pass at this was misread.)
# =============================================================================
set -euo pipefail

SCRATCH=/n/scratch/users/c/cos689/xtea2_testdata
WT="${XTEA2_WORKTREE:-/n/data1/hms/dbmi/park/corinne/claude_sandbox/xTea}"
R=/home/cos689/xTea/rep_lib_annotation/consensus
WORK=$(mktemp -d "$SCRATCH/bwabench.XXXXXX")
trap 'rm -rf "$WORK"' EXIT
export PATH=/home/cos689/software/bin:$PATH

CORES=${SLURM_CPUS_PER_TASK:-24}
echo "=============================================================="
echo "bwa mem scaling on the xTea2 realign workload"
echo "  job    : ${SLURM_JOB_ID:-<none>} on $(hostname)"
echo "  cores  : $CORES   (nproc=$(nproc), affinity=$(taskset -pc $$ | sed 's/.*: //'))"
echo "  bwa    : $(bwa 2>&1 | grep -i version | head -1)"
echo "  started: $(date '+%F %T')"
echo "=============================================================="

# ---------------------------------------------------------------------------------------
# Input: a real clip FASTQ from the scan pass. A 20 Mb chr20 window gives ~2M clip parts,
# big enough that startup is negligible and small enough to sweep many configurations.
# ---------------------------------------------------------------------------------------
echo
echo "--- generating the clip FASTQ from a real scan (chr20:20-40Mb) ---"
"$WT/src/csrc/xtea_scan" --bam "$SCRATCH/NA12878.wgs.bam" \
  --region chr20:20000000-40000000 --threads "$CORES" --min-clip 2 \
  --clip-fastq "$WORK/clip.fq" >/dev/null
NREADS=$(( $(wc -l < "$WORK/clip.fq") / 4 ))
echo "  $NREADS reads, $(du -h "$WORK/clip.fq" | cut -f1)"

run() {  # run <label> <threads> <ref> ; prints wall + cpu, output to /dev/null
  local T=$2 ref=$3
  /usr/bin/time -f "%e %U %S %M" bwa mem -t "$T" -K 1000000 -T 7 -k 7 -c 70 \
    "$ref" "$WORK/clip.fq" 2>"$WORK/e" >/dev/null
  tail -1 "$WORK/e"
}

echo
echo "##############################################################"
echo "# 1. THREAD SCALING, one bwa alone (per family)"
echo "##############################################################"
echo "If speedup plateaus well below the thread count, threads/nfamilies wastes cores."
for fam in SVA ALU LINE1; do
  echo
  echo "[$fam]"
  printf '  %-4s %9s %9s %9s %9s\n' -t wall cpu speedup cpu/wall
  base=""
  for T in 1 2 4 8 12 24; do
    [ "$T" -le "$CORES" ] || continue
    o=$(run "$fam" "$T" "$R/$fam.fa")
    w=$(echo "$o" | awk '{print $1}')
    [ -z "$base" ] && base=$w
    echo "$o" | awk -v T="$T" -v b="$base" \
      '{printf "  %-4s %9s %9.1f %9.2f %9.2f\n", T, $1, $2+$3, b/$1, ($2+$3)/$1}'
  done
done

echo
echo "##############################################################"
echo "# 2. THE CURRENT SHAPE: 3 families concurrent at CORES/3 threads each"
echo "##############################################################"
PER=$(( CORES / 3 ))
echo "bwa_threads_for($CORES, 3, i) = $PER threads per family"
S=$(date +%s)
for fam in LINE1 ALU SVA; do
  bwa mem -t "$PER" -K 1000000 -T 7 -k 7 -c 70 "$R/$fam.fa" "$WORK/clip.fq" \
    >/dev/null 2>"$WORK/e.$fam" &
done
wait
echo "  3 families concurrent, $PER threads each: $(( $(date +%s) - S ))s wall"

echo
echo "##############################################################"
echo "# 3. SHARDED: more bwa PROCESSES, fewer threads each"
echo "##############################################################"
echo "Only worth doing if part 1 shows a low plateau. Sharding splits the FASTQ, so it"
echo "changes each read's index within its bwa run -- and bwa's tie-break among EQUAL-"
echo "SCORING alignments is index-dependent (measured: AS=23 for three placements of one"
echo "33bp query). So it can change WHICH consensus copy a multi-mapper lands on."
echo "realign_sites_file only consumes the set of qnames whose PRIMARY alignment mapped,"
echo "which should be tie-break invariant -- that is asserted below, not assumed. It would"
echo "NOT be safe for the cns path (temp_clip.sam/temp_disc.sam), which reads RNAME/POS."
echo
mapped() { awk -F'\t' '!/^@/ && and($2,0x904)==0 {print $1}' "$@" | sort -u; }
bwa mem -t "$PER" -K 1000000 -T 7 -k 7 -c 70 "$R/SVA.fa" "$WORK/clip.fq" \
  2>/dev/null > "$WORK/whole.sam"
mapped "$WORK/whole.sam" > "$WORK/mapped.whole"
echo "  reference (unsharded, -t $PER): $(wc -l < "$WORK/mapped.whole") mapped qnames"
echo
printf '  %-8s %-10s %9s %9s  %s\n' shards thr/shard wall speedup "mapped-set"
for NS in 1 2 3 4 6 12; do
  TPS=$(( CORES / NS )); [ "$TPS" -lt 1 ] && TPS=1
  rm -f "$WORK"/s*.fq "$WORK"/s*.sam
  awk -v ns="$NS" -v w="$WORK" 'NR%4==1{r++} {print > (w "/s" (r%ns) ".fq")}' "$WORK/clip.fq"
  S=$(date +%s)
  for i in $(seq 0 $(( NS - 1 ))); do
    bwa mem -t "$TPS" -K 1000000 -T 7 -k 7 -c 70 "$R/SVA.fa" "$WORK/s$i.fq" \
      > "$WORK/s$i.sam" 2>/dev/null &
  done
  wait
  W=$(( $(date +%s) - S ))
  [ "$NS" -eq 1 ] && BASE=$W
  mapped "$WORK"/s*.sam > "$WORK/mapped.shard"
  same=$(cmp -s "$WORK/mapped.whole" "$WORK/mapped.shard" && echo IDENTICAL \
         || echo "DIFFERS($(comm -3 "$WORK/mapped.whole" "$WORK/mapped.shard" | wc -l))")
  awk -v n="$NS" -v t="$TPS" -v w="$W" -v b="$BASE" -v s="$same" \
    'BEGIN{printf "  %-8s %-10s %8ss %9.2f  %s\n", n, t, w, (w>0? b/w : 0), s}'
done
echo
echo "  Reading it: if speedup keeps climbing with shard count while part 1's plateaued,"
echo "  the win is process-level parallelism and realign_sites_file should shard rather"
echo "  than hand each family CORES/3 threads. IDENTICAL in the last column is the"
echo "  precondition -- without it, sharding is off the table for the detection realign too."

echo
echo "##############################################################"
echo "# 4. bwa -K: peak RSS (measured in-session already; confirmed here)"
echo "##############################################################"
echo "-K is the batch size in bases. It exists here for determinism (output independent of"
echo "thread count) and any FIXED value keeps that. Smaller = less resident, and the output"
echo "is byte-identical, so the only question is whether it costs time."
printf '  %-12s %9s %12s\n' -K wall maxRSS
for K in 10000000 1000000 200000; do
  o=$(/usr/bin/time -f "%e %U %S %M" bwa mem -t "$PER" -K "$K" -T 7 -k 7 -c 70 \
        "$R/SVA.fa" "$WORK/clip.fq" 2>"$WORK/e" >/dev/null; tail -1 "$WORK/e")
  echo "$o" | awk -v K="$K" '{printf "  %-12s %8ss %10.2f GB\n", K, $1, $4/1048576}'
done
echo
echo "  (in-session, 2 cores, 193k reads: 1.87 GB / 0.38 GB / 0.099 GB, output"
echo "   byte-identical across all three -- a 19x resident reduction for free)"

echo
echo "=============================================================="
echo "finished $(date '+%F %T')"
echo "=============================================================="
