#!/bin/bash
# One-screen status for the mosaic run fleet. Safe to run any time.
#   docs/handoff/runs/mosaic_status.sh
#
# Reads the live slurm queue plus each run's own log, disk sampler and step timings, so it
# reports what the runs ACTUALLY did rather than what they were asked to do -- the cutoff bug
# (run 50427774) was visible only in the step log, while the launcher's own banner said 2/2.
SCR=/n/scratch/users/c/cos689/xtea2_testdata

printf '%-12s %-10s %-9s %-10s %-9s %s\n' JOB NAME STATE ELAPSED DISK STEP
squeue -u cos689 -h -o "%i %j %T %M" 2>/dev/null | grep -v claude_ses | while read -r id name state el; do
    case "$name" in
        chr1_nocns) d=$SCR/mosaic_hapmap_bcm_chr1_nocnsemit ;;
        mosaic_hap*) d=$SCR/mosaic_hapmap_bcm_chr1 ;;
        *) d=$SCR/mosaic_hapmap_wgs_$name ;;
    esac
    disk=$(tail -1 "$d/disk_usage.tsv" 2>/dev/null | cut -f3)
    log=$(ls -t $SCR/mosaic_wgs_${name}_*.out $SCR/mosaic_hapmap_chr1_${id}.out 2>/dev/null | head -1)
    step=$(grep -E "^\[xtea2 \+" "$log" 2>/dev/null | tail -1 | sed -E 's/^\[xtea2 \+ *([0-9.]+)s\] //; s|/n/scratch.*||' | cut -c1-46)
    printf '%-12s %-10s %-9s %-10s %-9s %s\n' "$id" "$name" "$state" "$el" "${disk:-–}" "${step:-scan pass}"
done

echo
echo "--- finished since 2026-08-12 ---"
sacct -u cos689 -S 2026-08-12T20:00 -o JobID,JobName%12,State,Elapsed,MaxRSS,ExitCode --units=G 2>/dev/null \
  | grep -vE "\.extern|\.batch|claude_sess|^JobID|^---" | grep -vE "RUNNING|PENDING" | head -12 | sed 's/^/  /'

echo
echo "--- mosaic call sets produced so far ---"
found=0
for d in $SCR/mosaic_hapmap_wgs_* $SCR/mosaic_hapmap_bcm_chr1*; do
    [ -d "$d" ] || continue
    for f in "$d"/*/*/candidate_disc_filtered_cns_post_filtering_mosaic.txt; do
        [ -f "$f" ] || continue
        found=1
        printf '  %-58s %s\n' "$(echo "$f" | sed "s|$SCR/||; s|/candidate.*||")" "$(wc -l < "$f") calls"
    done
done
[ "$found" = 1 ] || echo "  (none yet -- no run has reached postF_mosaic)"

echo
echo "--- the A/B: --cns-emit-scan on vs off, same chr1 slice ---"
for tag in "with:mosaic_hapmap_bcm_chr1" "without:mosaic_hapmap_bcm_chr1_nocnsemit"; do
    lab=${tag%%:*}; d=$SCR/${tag#*:}
    t=$(awk -F'\t' '$1=="phase" && $2=="scan_pass" {printf "scan %.0fs", $3}' "$d"/*/step_timings.tsv 2>/dev/null | head -1)
    pk=$(awk -F'\t' 'NR>1 && $2+0>p {p=$2; a=$1} END {if (p) printf "peak %.1f GB at t+%ds", p/1e9, a}' "$d/disk_usage.tsv" 2>/dev/null)
    printf '  %-8s %-18s %s\n' "$lab" "${t:-–}" "${pk:-–}"
done
