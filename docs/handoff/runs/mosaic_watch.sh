#!/bin/bash
# Block until the FIRST job in the mosaic fleet leaves the queue, then report it.
# First-exit rather than all-exit on purpose: the two chr1 canaries finish long before the
# WGS arms, and they carry the answers that matter earliest (does the mode produce calls at
# all, and the --cns-emit-scan A/B). Waiting for all seven would sit on that for many hours.
JOBS="${*:?usage: mosaic_watch.sh <jobid> [jobid ...]}"
SCR=/n/scratch/users/c/cos689/xtea2_testdata
alive () { squeue -h -j "$(echo $JOBS | tr ' ' ,)" -o %i 2>/dev/null | wc -l; }
N0=$(alive)
echo "watching $N0 jobs: $JOBS"
while [ "$(alive)" -ge "$N0" ] && [ "$N0" -gt 0 ]; do sleep 120; done
echo "=== a job finished at $(date '+%F %T') ==="
for j in $JOBS; do
    squeue -h -j "$j" -o %i 2>/dev/null | grep -q . && continue
    sacct -j "$j" -o JobID,JobName%14,State,Elapsed,MaxRSS,ExitCode --units=G -n 2>/dev/null \
      | grep -E "\.batch" | sed 's/^/  DONE /'
done
echo "--- full status ---"
bash /n/data1/hms/dbmi/park/corinne/claude_sandbox/xTea/docs/handoff/runs/mosaic_status.sh
