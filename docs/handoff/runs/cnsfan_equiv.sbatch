#!/bin/bash
#SBATCH -c 24
#SBATCH -t 4:00:00
#SBATCH --mem=64G
#SBATCH -p priopark
#SBATCH --account=park_contrib
#SBATCH -o /n/scratch/users/c/cos689/xtea2_testdata/cnsfan_equiv_%j.out
#SBATCH --mail-type=NONE

# =============================================================================
# Does the reworked cns fan emit the SAME temp SAMs as the old one, and is it faster?
#
# The rework is three changes to emit_cns_from_spill:
#   1. spill mate_pos, so mates can be fetched by position instead of by a genome-wide
#      name scan;
#   2. drop clips/anchors that overlap no candidate of ANY type BEFORE gathering -- the
#      old code gathered every disc mate in the genome (~4.15M) to serve the ~0.5M any
#      type actually reads;
#   3. gather in parallel over coordinate-sorted merged regions instead of one
#      single-threaded whole-genome pass.
#
# (2) is the one that has to be argued rather than measured: the per-type fan only emits
# when J.cand.overlapping(...) is non-empty, so anything overlapping no candidate of any
# type cannot reach any output. That makes the filter output-neutral BY CONSTRUCTION --
# and this job is what checks the construction against reality.
#
# CORRECTNESS IS THE POINT, speed is secondary. If the temp SAMs differ at all, the
# rework is wrong and the timing is irrelevant. Both runs are full --cns-emit-scan runs
# over the same chr1 BAM, differing ONLY in the scan binary.
# =============================================================================
set -euo pipefail

SCRATCH=/n/scratch/users/c/cos689/xtea2_testdata
REPO=/n/data1/hms/dbmi/park/corinne/claude_sandbox/xTea
OLDBIN=${OLDBIN:?set OLDBIN to the pre-rework xtea_scan}
NEWBIN="$REPO/src/csrc/xtea_scan"
BAM="$SCRATCH/NA12878.chr1.bam"

source /n/app/conda/miniforge3/24.11.3-0/etc/profile.d/conda.sh
conda activate xtea2
export PATH=/home/cos689/software/bin:$PATH
export XTEA_BASIC_INFO_SEED=42

for b in "$OLDBIN" "$NEWBIN"; do [ -x "$b" ] || { echo "FATAL: $b not executable"; exit 1; }; done
[ -f "$BAM" ] || { echo "FATAL: $BAM missing"; exit 1; }

echo "=============================================================="
echo "cns fan equivalence: old vs reworked, --cns-emit-scan, chr1"
echo "  job : ${SLURM_JOB_ID:-<none>} on $(hostname)  cores=${SLURM_CPUS_PER_TASK:-24}"
echo "  old : $OLDBIN"
echo "  new : $NEWBIN ($(cd "$REPO" && git rev-parse --short HEAD), $(cd "$REPO" && git status --porcelain | wc -l) file(s) dirty)"
echo "=============================================================="

run_one() {  # run_one <tag> <binary>
  local tag=$1 bin=$2 out="$SCRATCH/cnsfan_${tag}_out"
  rm -rf "$out"; mkdir -p "$out"
  local t0=$(date +%s)
  /usr/bin/time -v python "$REPO/xtea/run_xtea.py" \
    -i "$BAM" -s NA12878.chr1 -o "$out" \
    --reference /home/cos689/hg38_no_alt.fa \
    --rep-lib-dir /home/cos689/xTea/rep_lib_annotation \
    --gene-gff3 /home/cos689/data1/corinne/ref/gencode.v45.basic.annotation.gff3 \
    --blacklist /home/cos689/xTea/rep_lib_annotation/blacklist/hg38/centromere.bed \
    --types L1,Alu,SVA --scan-bin "$bin" --bwa /home/cos689/software/bin/bwa \
    --scan-min-clip 2 --col5-families alu,sva --cns-emit-scan --parallel-downstream \
    -n "${SLURM_CPUS_PER_TASK:-24}" > "$out/log.txt" 2>&1 || true
  echo "$(( $(date +%s) - t0 ))" > "$out/wall.txt"
}

for pair in "old:$OLDBIN" "new:$NEWBIN"; do
  tag=${pair%%:*}; bin=${pair#*:}
  echo; echo "--- running $tag ---"
  run_one "$tag" "$bin"
  o="$SCRATCH/cnsfan_${tag}_out"
  printf '  wall %ss   peak RSS %s\n' "$(cat "$o/wall.txt")" \
    "$(grep -a 'Maximum resident' "$o/log.txt" | awk '{printf "%.2f GB", $NF/1048576}')"
  grep -a '^\[cns_fan\]' "$o/log.txt" | sed 's/^/  /' || echo "  (no [cns_fan] line -- expected for the old binary)"
  grep -a '^phase' "$o/NA12878.chr1/step_timings.tsv" 2>/dev/null | sed 's/^/  /' || true
done

O="$SCRATCH/cnsfan_old_out/NA12878.chr1"
N="$SCRATCH/cnsfan_new_out/NA12878.chr1"

echo
echo "##############################################################"
echo "# CORRECTNESS: the emitted temp SAMs must be identical"
echo "##############################################################"
echo "Compared on the records, not the bytes: bwa's \@PG CL: header differs whenever the"
echo "invocation differs, and record ORDER can differ because the fan's FASTA is built by"
echo "iterating maps. So: sorted record multisets, which is what the parsers consume."
FAIL=0
for t in L1 Alu SVA; do
  for k in clip disc; do
    a="$O/$t/tmp/cns/temp_$k.sam"; b="$N/$t/tmp/cns/temp_$k.sam"
    if [ ! -f "$a" ] || [ ! -f "$b" ]; then printf '  %-4s %-5s MISSING\n' "$t" "$k"; FAIL=1; continue; fi
    na=$(grep -avc '^@' "$a"); nb=$(grep -avc '^@' "$b")
    if cmp -s <(grep -av '^@' "$a" | sort) <(grep -av '^@' "$b" | sort); then
      printf '  %-4s %-5s IDENTICAL   (%s records)\n' "$t" "$k" "$na"
    else
      d=$(comm -3 <(grep -av '^@' "$a" | sort) <(grep -av '^@' "$b" | sort) | wc -l)
      printf '  %-4s %-5s DIFFERS     old=%s new=%s  %s differing records\n' "$t" "$k" "$na" "$nb" "$d"
      comm -3 <(grep -av '^@' "$a" | sort) <(grep -av '^@' "$b" | sort) | head -4 | cut -c1-140 | sed 's/^/        /'
      FAIL=1
    fi
  done
  # low_mapq_clip_cnt.txt is also emitted by the fan and read by the bypass
  a="$O/$t/tmp/cns/low_mapq_clip_cnt.txt"; b="$N/$t/tmp/cns/low_mapq_clip_cnt.txt"
  if [ -f "$a" ] && [ -f "$b" ]; then
    cmp -s <(sort "$a") <(sort "$b") && printf '  %-4s lowq  IDENTICAL\n' "$t" \
      || { printf '  %-4s lowq  DIFFERS\n' "$t"; FAIL=1; }
  fi
done

echo
echo "--- call sets (a consequence of the above; reported for completeness) ---"
for pair in "L1:LINE1" "Alu:ALU" "SVA:SVA"; do
  t=${pair%%:*}; alt=${pair##*:}
  a=$(ls "$O"/merged/*_${alt}.vcf 2>/dev/null | head -1); b=$(ls "$N"/merged/*_${alt}.vcf 2>/dev/null | head -1)
  if [ -n "$a" ] && [ -n "$b" ]; then
    printf '  %-4s old=%s new=%s  +%s/-%s\n' "$t" "$(grep -vc '^#' "$a")" "$(grep -vc '^#' "$b")" \
      "$(comm -13 <(grep -v '^#' "$a" | cut -f1,2 | sort -u) <(grep -v '^#' "$b" | cut -f1,2 | sort -u) | wc -l)" \
      "$(comm -23 <(grep -v '^#' "$a" | cut -f1,2 | sort -u) <(grep -v '^#' "$b" | cut -f1,2 | sort -u) | wc -l)"
  fi
done

echo
echo "##############################################################"
echo "# SPEED / MEMORY / DISK"
echo "##############################################################"
printf '  %-6s %10s %12s %14s\n' run wall peakRSS spill_on_disk
for tag in old new; do
  o="$SCRATCH/cnsfan_${tag}_out"; s="$o/NA12878.chr1"
  printf '  %-6s %9ss %12s %14s\n' "$tag" "$(cat "$o/wall.txt")" \
    "$(grep -a 'Maximum resident' "$o/log.txt" | awk '{printf "%.2f GB", $NF/1048576}')" \
    "$(du -cb "$s"/cns_spill_*.txt 2>/dev/null | tail -1 | awk '{printf "%.1f MB", $1/1048576}')"
done
echo
echo "  (spill size is expected to grow slightly -- mate_pos adds ~9 bytes per anchor line."
echo "   That is the disk cost of the speed fix and is reported, not hidden.)"

echo
echo "=============================================================="
echo "correctness: $([ "$FAIL" -eq 0 ] && echo 'PASS -- temp SAMs identical' || echo 'FAIL -- see above')"
echo "finished $(date '+%F %T')"
echo "=============================================================="
exit "$FAIL"
