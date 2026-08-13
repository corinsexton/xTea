#!/bin/bash
# Block until a job in the mosaic fleet REALLY leaves the queue, then report it.
#
# First-exit rather than all-exit on purpose: the chr1 canaries finish long before the WGS
# arms and carry the answers that matter earliest.
#
# THE BUG THIS VERSION FIXES. The first version counted `squeue | wc -l` and exited when the
# count dropped. A transient squeue failure returns no lines, wc gives 0, 0 < N, and the
# watcher reported a completion that had not happened -- observed once, with all six jobs
# still RUNNING. Two changes make that impossible:
#   * a job counts as gone only when sacct confirms a TERMINAL state for it, so an empty or
#     failed squeue proves nothing;
#   * the absence must persist across two polls, so a single flaky call cannot trigger it.
JOBS="${*:?usage: mosaic_watch.sh <jobid> [jobid ...]}"
POLL=${POLL:-120}
declare -A SUSPECT

terminal () {   # sacct state for one job, or empty if not yet terminal
    sacct -j "$1" -o State -n 2>/dev/null | head -1 | awk '{print $1}' \
      | grep -Eo 'COMPLETED|FAILED|CANCELLED|TIMEOUT|OUT_OF_MEMORY|NODE_FAIL|PREEMPTED'
}

echo "watching: $JOBS  (poll ${POLL}s, terminal state confirmed via sacct)"
while :; do
    Q=$(squeue -h -j "$(echo $JOBS | tr ' ' ,)" -o %i 2>/dev/null)
    for j in $JOBS; do
        if echo "$Q" | grep -qw "$j"; then SUSPECT[$j]=0; continue; fi
        st=$(terminal "$j")
        [ -n "$st" ] || { SUSPECT[$j]=0; continue; }        # absent but not terminal -> ignore
        SUSPECT[$j]=$(( ${SUSPECT[$j]:-0} + 1 ))
        if [ "${SUSPECT[$j]}" -ge 2 ]; then
            echo "=== job $j finished ($st) at $(date '+%F %T') ==="
            sacct -j "$j" -o JobID,JobName%14,State,Elapsed,MaxRSS,ExitCode --units=G -n 2>/dev/null \
              | grep -E "\.batch" | sed 's/^/  /'
            echo "--- full status ---"
            bash "$(dirname "$0")/mosaic_status.sh"
            exit 0
        fi
    done
    sleep "$POLL"
done
