#!/bin/bash
#SBATCH -c 24
#SBATCH -t 6:00:00
#SBATCH --mem=64G
#SBATCH -p priopark
#SBATCH --account=park_contrib
#SBATCH -o /n/scratch/users/c/cos689/xtea2_testdata/cnsport_%j.out
#SBATCH --mail-type=NONE

# =============================================================================
# The ported cns collection (C++ / htslib, per-candidate indexed fetch) vs the fixed
# Python step. chr1, everything else held fixed, so the ONLY difference is which code
# collects the cns reads.
#
# Two runs of the SAME commit:
#   python : no --cns-emit-scan  -> x_clip_disc_filter collects in Python (B10-fixed)
#   ported : --cns-emit-scan     -> xtea-scan --cns-manifest collects in C++
#
# CORRECTNESS FIRST. The port keeps master's access pattern (an indexed fetch per
# candidate window) precisely so the emitted temp SAMs can be compared record for record
# against the Python step. If they differ, the port is wrong and the timing is irrelevant.
#
# Then the three axes that decide whether the port is worth keeping (user's bar): speed,
# peak memory, and disk. The deleted fold-in path failed all three; this one has to win
# them or the Python step stays.
# =============================================================================
set -euo pipefail

SCRATCH=/n/scratch/users/c/cos689/xtea2_testdata
REPO=/n/data1/hms/dbmi/park/corinne/claude_sandbox/xTea
WT="${XTEA2_WORKTREE:-$REPO}"
BAM="$SCRATCH/NA12878.chr1.bam"

source /n/app/conda/miniforge3/24.11.3-0/etc/profile.d/conda.sh
conda activate xtea2
export PATH=/home/cos689/software/bin:$PATH
export XTEA_BASIC_INFO_SEED=42

[ -x "$WT/src/csrc/xtea_scan" ] || { echo "FATAL: no xtea_scan in $WT"; exit 1; }
[ -f "$BAM" ] || { echo "FATAL: $BAM missing"; exit 1; }

echo "=============================================================="
echo "ported C++ cns collection vs fixed Python -- chr1"
echo "  job : ${SLURM_JOB_ID:-<none>} on $(hostname)  cores=${SLURM_CPUS_PER_TASK:-24}"
echo "  tree: $WT ($(cd "$WT" && git rev-parse --short HEAD), $(cd "$WT" && git status --porcelain | wc -l) dirty)"
echo "=============================================================="

run_one() {  # run_one <tag> <extra args...>
  local tag=$1; shift
  local out="$SCRATCH/cnsport_${tag}_out"
  rm -rf "$out"; mkdir -p "$out"
  local t0=$(date +%s)
  /usr/bin/time -v python "$WT/xtea/run_xtea.py" \
    -i "$BAM" -s NA12878.chr1 -o "$out" \
    --reference /home/cos689/hg38_no_alt.fa \
    --rep-lib-dir /home/cos689/xTea/rep_lib_annotation \
    --gene-gff3 /home/cos689/data1/corinne/ref/gencode.v45.basic.annotation.gff3 \
    --blacklist /home/cos689/xTea/rep_lib_annotation/blacklist/hg38/centromere.bed \
    --types L1,Alu,SVA --scan-bin "$WT/src/csrc/xtea_scan" --bwa /home/cos689/software/bin/bwa \
    --scan-min-clip 2 --col5-families alu,sva --parallel-downstream \
    -n "${SLURM_CPUS_PER_TASK:-24}" "$@" > "$out/log.txt" 2>&1 || true
  echo "$(( $(date +%s) - t0 ))" > "$out/wall.txt"
}

echo; echo "--- python (B10-fixed, no --cns-emit-scan) ---"; run_one python
echo "--- ported (C++ --cns-emit-scan) ---";                run_one ported --cns-emit-scan

P="$SCRATCH/cnsport_python_out/NA12878.chr1"
C="$SCRATCH/cnsport_ported_out/NA12878.chr1"

echo
echo "##############################################################"
echo "# 1. CORRECTNESS -- emitted temp SAMs, record for record"
echo "##############################################################"
echo "Sorted record multisets, headers excluded: bwa's @PG CL: line differs whenever the"
echo "invocation does, and record order can differ, but the parsers consume records."
FAIL=0
for t in L1 Alu SVA; do
  for k in clip disc; do
    a="$P/$t/tmp/cns/temp_$k.sam"; b="$C/$t/tmp/cns/temp_$k.sam"
    if [ ! -f "$a" ] || [ ! -f "$b" ]; then printf '  %-4s %-5s MISSING (py=%s cpp=%s)\n' "$t" "$k" \
      "$([ -f "$a" ] && echo y || echo n)" "$([ -f "$b" ] && echo y || echo n)"; FAIL=1; continue; fi
    na=$(grep -avc '^@' "$a"); nb=$(grep -avc '^@' "$b")
    if cmp -s <(grep -av '^@' "$a" | sort) <(grep -av '^@' "$b" | sort); then
      printf '  %-4s %-5s IDENTICAL  (%s records)\n' "$t" "$k" "$na"
    else
      printf '  %-4s %-5s DIFFERS    python=%s cpp=%s\n' "$t" "$k" "$na" "$nb"
      echo "        python-only:"; comm -23 <(grep -av '^@' "$a" | sort) <(grep -av '^@' "$b" | sort) | head -3 | cut -c1-130 | sed 's/^/          /'
      echo "        cpp-only:";    comm -13 <(grep -av '^@' "$a" | sort) <(grep -av '^@' "$b" | sort) | head -3 | cut -c1-130 | sed 's/^/          /'
      FAIL=1
    fi
  done
  a="$P/$t/tmp/cns/low_mapq_clip_cnt.txt"; b="$C/$t/tmp/cns/low_mapq_clip_cnt.txt"
  if [ -f "$a" ] && [ -f "$b" ]; then
    if cmp -s <(sort "$a") <(sort "$b"); then printf '  %-4s lowq  IDENTICAL\n' "$t"
    else printf '  %-4s lowq  DIFFERS (%s lines)\n' "$t" "$(comm -3 <(sort "$a") <(sort "$b") | wc -l)"; FAIL=1; fi
  fi
done

echo
echo "--- call sets (a consequence; reported) ---"
for pair in "L1:LINE1" "Alu:ALU" "SVA:SVA"; do
  t=${pair%%:*}; alt=${pair##*:}
  a=$(ls "$P"/merged/*_${alt}.vcf 2>/dev/null | head -1); b=$(ls "$C"/merged/*_${alt}.vcf 2>/dev/null | head -1)
  [ -n "$a" ] && [ -n "$b" ] && printf '  %-4s python=%s cpp=%s  +%s/-%s\n' "$t" \
    "$(grep -vc '^#' "$a")" "$(grep -vc '^#' "$b")" \
    "$(comm -13 <(grep -v '^#' "$a"|cut -f1,2|sort -u) <(grep -v '^#' "$b"|cut -f1,2|sort -u)|wc -l)" \
    "$(comm -23 <(grep -v '^#' "$a"|cut -f1,2|sort -u) <(grep -v '^#' "$b"|cut -f1,2|sort -u)|wc -l)"
done

echo
echo "##############################################################"
echo "# 2. SPEED / MEMORY / DISK -- the bar the port has to clear"
echo "##############################################################"
printf '  %-8s %10s %12s %14s\n' run wall peakRSS extra_scratch
for tag in python ported; do
  o="$SCRATCH/cnsport_${tag}_out"; s="$o/NA12878.chr1"
  # extra scratch = cns working files that are NOT the temp SAMs both paths write
  extra=$(du -cb "$s"/*/tmp/cns/all_disc_pos.disc_pos "$s"/cns_*.txt 2>/dev/null | tail -1 | awk '{print $1}')
  printf '  %-8s %9ss %12s %14s\n' "$tag" "$(cat "$o/wall.txt")" \
    "$(grep -a 'Maximum resident' "$o/log.txt" | awk '{printf "%.2f GB", $NF/1048576}')" \
    "$(awk -v b="${extra:-0}" 'BEGIN{printf "%.1f MB", b/1048576}')"
done
echo
echo "--- phase breakdown (the cns step is cns_emit for ported, inside downstream for python) ---"
for tag in python ported; do
  echo "  [$tag]"; grep -a '^phase' "$SCRATCH/cnsport_${tag}_out/NA12878.chr1/step_timings.tsv" 2>/dev/null | sed 's/^/    /'
  grep -a '^\[cns_fan\]\|^\[bwa_sam\]' "$SCRATCH/cnsport_${tag}_out/log.txt" 2>/dev/null | head -6 | sed 's/^/    /' || true
done
echo
echo "--- python's per-type cns step wall (for the like-for-like collection comparison) ---"
for tag in python ported; do
  echo -n "  [$tag] cns steps: "
  awk -F'\t' '$1 ~ /^step/ && $2 ~ /cns_filter|:cns/ {printf "%s=%.1fs ", $2, $3}' \
    "$SCRATCH/cnsport_${tag}_out/NA12878.chr1/step_timings.tsv" 2>/dev/null; echo
done

echo
echo "=============================================================="
echo "correctness: $([ "$FAIL" -eq 0 ] && echo 'PASS -- temp SAMs identical to the Python step' || echo 'FAIL')"
echo "finished $(date '+%F %T')"
echo "=============================================================="
exit "$FAIL"
