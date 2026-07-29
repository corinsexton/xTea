#!/bin/bash
#SBATCH -c 24
#SBATCH -t 3:00:00
#SBATCH --mem=64G
#SBATCH -p priopark
#SBATCH --account=park_contrib
#SBATCH -o /n/scratch/users/c/cos689/xtea2_testdata/cnsport_rerun_%j.out
#SBATCH --mail-type=NONE

# Ported C++ cns collection, chr1, compared against the ALREADY-COMPLETED Python run
# (cnsport_python_out, 323 s, from job 49068968 -- its Python side succeeded; only the
# ported side failed, on a wrong option name). Rerunning just the one side that failed.
#
# Python baseline from that run, for reference:
#   wall 323 s | cns_filter steps L1 52.0s / Alu 57.4s / SVA 51.2s (concurrent, so ~57 s of
#   wall) | extra scratch 2.0 MB | candidates L1 89 / Alu 469 / SVA 134
#
# Note what cns_filter contains: collection + bwa + parse + filter. The C++ emit replaces
# only the first two; the step still runs afterwards, on the bypass, to parse and filter.
# So the honest comparison is TOTAL wall, not step-vs-step.
set -euo pipefail

SCRATCH=/n/scratch/users/c/cos689/xtea2_testdata
REPO=/n/data1/hms/dbmi/park/corinne/claude_sandbox/xTea
OUT="$SCRATCH/cnsport_ported_out"
P="$SCRATCH/cnsport_python_out/NA12878.chr1"

source /n/app/conda/miniforge3/24.11.3-0/etc/profile.d/conda.sh
conda activate xtea2
export PATH=/home/cos689/software/bin:$PATH
export XTEA_BASIC_INFO_SEED=42

[ -d "$P" ] || { echo "FATAL: python baseline $P missing"; exit 1; }
rm -rf "$OUT"; mkdir -p "$OUT"

echo "=== ported C++ cns collection, chr1 -- job ${SLURM_JOB_ID:-<none>} on $(hostname)"
echo "    tree $(cd "$REPO" && git rev-parse --short HEAD), $(cd "$REPO" && git status --porcelain | wc -l) dirty"
T0=$(date +%s)
/usr/bin/time -v python "$REPO/xtea/run_xtea.py" \
  -i "$SCRATCH/NA12878.chr1.bam" -s NA12878.chr1 -o "$OUT" \
  --reference /home/cos689/hg38_no_alt.fa \
  --rep-lib-dir /home/cos689/xTea/rep_lib_annotation \
  --gene-gff3 /home/cos689/data1/corinne/ref/gencode.v45.basic.annotation.gff3 \
  --blacklist /home/cos689/xTea/rep_lib_annotation/blacklist/hg38/centromere.bed \
  --types L1,Alu,SVA --scan-bin "$REPO/src/csrc/xtea_scan" --bwa /home/cos689/software/bin/bwa \
  --scan-min-clip 2 --col5-families alu,sva --parallel-downstream --cns-emit-scan \
  -n "${SLURM_CPUS_PER_TASK:-24}" > "$OUT/log.txt" 2>&1 || true
WALL=$(( $(date +%s) - T0 ))
C="$OUT/NA12878.chr1"

echo
echo "##############################################################"
echo "# 1. CORRECTNESS -- temp SAMs vs the Python step, record for record"
echo "##############################################################"
FAIL=0
for t in L1 Alu SVA; do
  for k in clip disc; do
    a="$P/$t/tmp/cns/temp_$k.sam"; b="$C/$t/tmp/cns/temp_$k.sam"
    if [ ! -f "$b" ]; then printf '  %-4s %-5s MISSING from the ported run\n' "$t" "$k"; FAIL=1; continue; fi
    na=$(grep -avc '^@' "$a"); nb=$(grep -avc '^@' "$b")
    if cmp -s <(grep -av '^@' "$a"|sort) <(grep -av '^@' "$b"|sort); then
      printf '  %-4s %-5s IDENTICAL  (%s records)\n' "$t" "$k" "$na"
    else
      printf '  %-4s %-5s DIFFERS    python=%s cpp=%s  (%s records differ)\n' "$t" "$k" "$na" "$nb" \
        "$(comm -3 <(grep -av '^@' "$a"|sort) <(grep -av '^@' "$b"|sort)|wc -l)"
      echo "        python-only:"; comm -23 <(grep -av '^@' "$a"|sort) <(grep -av '^@' "$b"|sort)|head -3|cut -c1-130|sed 's/^/          /'
      echo "        cpp-only:";    comm -13 <(grep -av '^@' "$a"|sort) <(grep -av '^@' "$b"|sort)|head -3|cut -c1-130|sed 's/^/          /'
      FAIL=1
    fi
  done
  a="$P/$t/tmp/cns/low_mapq_clip_cnt.txt"; b="$C/$t/tmp/cns/low_mapq_clip_cnt.txt"
  [ -f "$a" ] && [ -f "$b" ] && { cmp -s <(sort "$a") <(sort "$b") && printf '  %-4s lowq  IDENTICAL\n' "$t" \
     || { printf '  %-4s lowq  DIFFERS (%s lines)\n' "$t" "$(comm -3 <(sort "$a") <(sort "$b")|wc -l)"; FAIL=1; }; }
done

echo
echo "##############################################################"
echo "# 2. SPEED / MEMORY / DISK"
echo "##############################################################"
printf '  %-8s %8s %12s %14s\n' run wall peakRSS extra_scratch
printf '  %-8s %7ss %12s %14s\n' python 323 \
  "$(grep -a 'Maximum resident' "$SCRATCH/cnsport_python_out/log.txt"|awk '{printf "%.2f GB",$NF/1048576}')" \
  "$(du -cb "$P"/*/tmp/cns/all_disc_pos.disc_pos 2>/dev/null|tail -1|awk '{printf "%.1f MB",$1/1048576}')"
printf '  %-8s %7ss %12s %14s\n' ported "$WALL" \
  "$(grep -a 'Maximum resident' "$OUT/log.txt"|awk '{printf "%.2f GB",$NF/1048576}')" \
  "$(du -cb "$C"/*/tmp/cns/all_disc_pos.disc_pos "$C"/cns_*.txt 2>/dev/null|tail -1|awk '{printf "%.1f MB",$1/1048576}')"
echo
echo "--- phases ---"
for tag in python ported; do
  d=$([ $tag = python ] && echo "$P" || echo "$C")
  echo "  [$tag]"; awk -F'\t' '$1=="phase"{printf "    %-22s %7.1fs  childRSS %.2f GB\n",$2,$3,$4}' "$d/step_timings.tsv" 2>/dev/null
  awk -F'\t' '$1 ~ /^step/ && $2 ~ /cns/ {printf "    %-22s %7.1fs\n",$2,$3}' "$d/step_timings.tsv" 2>/dev/null
done
echo
echo "--- C++ collection internals ---"
grep -a '^\[cns_fan\]\|^\[bwa_sam\]' "$OUT/log.txt" 2>/dev/null | head -8 | sed 's/^/  /' || echo "  (none)"
echo
echo "--- call sets ---"
for pair in "L1:LINE1" "Alu:ALU" "SVA:SVA"; do
  t=${pair%%:*}; alt=${pair##*:}
  a=$(ls "$P"/merged/*_${alt}.vcf 2>/dev/null|head -1); b=$(ls "$C"/merged/*_${alt}.vcf 2>/dev/null|head -1)
  [ -n "$a" ] && [ -n "$b" ] && printf '  %-4s python=%s cpp=%s  +%s/-%s\n' "$t" \
    "$(grep -vc '^#' "$a")" "$(grep -vc '^#' "$b")" \
    "$(comm -13 <(grep -v '^#' "$a"|cut -f1,2|sort -u) <(grep -v '^#' "$b"|cut -f1,2|sort -u)|wc -l)" \
    "$(comm -23 <(grep -v '^#' "$a"|cut -f1,2|sort -u) <(grep -v '^#' "$b"|cut -f1,2|sort -u)|wc -l)"
done
echo
echo "correctness: $([ "$FAIL" -eq 0 ] && echo 'PASS' || echo 'FAIL')"
echo "finished $(date '+%F %T')"
exit "$FAIL"
