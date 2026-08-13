#!/bin/bash
# Sample a run tree's disk footprint until its slurm job ends.
#   sample_disk.sh <jobid> <dir> <out.tsv> [interval_s]
# Used for runs that started before mosaic_hapmap_wgs.sbatch grew its own sampler.
# Peak matters more than final: intermediates are removed as steps complete, so a
# final-size measurement can be a small fraction of the high-water mark.
JOB="${1:?jobid}"; DIR="${2:?dir}"; OUT="${3:?out.tsv}"; INT="${4:-60}"
printf 'elapsed_s\ttotal_bytes\ttotal_h\tn_files\n' > "$OUT"
T0=$SECONDS
while squeue -j "$JOB" -h -o %T 2>/dev/null | grep -qE 'RUNNING|PENDING|COMPLETING'; do
    B=$(du -sb "$DIR" 2>/dev/null | cut -f1); B=${B:-0}
    N=$(find "$DIR" -type f 2>/dev/null | wc -l)
    printf '%s\t%s\t%s\t%s\n' "$((SECONDS-T0))" "$B" "$(numfmt --to=iec "$B" 2>/dev/null)" "$N" >> "$OUT"
    sleep "$INT"
done
awk -F'\t' 'NR>1 && $2+0>p {p=$2; a=$1; n=$4}
     END {printf "PEAK %.2f GB at t+%ds (%d files)\n", p/1e9, a, n}' "$OUT"
