#!/bin/bash
# =============================================================================
# compare_runs.sh -- side-by-side report for the two WGS runs launched 2026-07-26.
#
#   xtea2  (trim steps 1-4 + instrumentation)  job 48824433
#   master (pristine 46520ab, SERIAL, 24-core) job 48825778   <- the true comparator
#
# Run it any time. It reports what is available and says plainly what is not, so it
# is useful while the jobs are still going, not only after.
#
#   ./compare_runs.sh                      # uses the default job ids/dirs below
#   ./compare_runs.sh 48824433 48825778    # explicit job ids
#
# Why the master run and not the golden. The frozen golden was THREE concurrent
# 8-core jobs, so its 38:31 is max(37:15, 29:06, 38:31) over a 24-core total
# allocation. xtea2 is one process with one 24-core allocation. Comparing the two
# charges xtea2 for concurrency SLURM gave master for free. Job 48825778 is master's
# own code run in the same shape as xtea2 -- serial, one job, 24 cores -- which is
# the only denominator that isolates code speed from job-shape luck.
# =============================================================================
set -uo pipefail

SCRATCH=/n/scratch/users/c/cos689/xtea2_testdata
X2_JOB="${1:-48824433}"
MA_JOB="${2:-48825778}"
X2_OUT="${XTEA2_OUT:-$SCRATCH/step45_wgs_out}"
MA_OUT="${MASTER_OUT:-$SCRATCH/master_serial_out}"
X2_LOG="$SCRATCH/step45_wgs_${X2_JOB}.out"
MA_LOG="$SCRATCH/master_serial_wgs_${MA_JOB}.out"
X2_S="$X2_OUT/NA12878.wgs"
MA_S="$MA_OUT/NA12878"

hr() { printf '%s\n' "----------------------------------------------------------------------"; }
hdr() { echo; printf '======================================================================\n'; echo "$1"; printf '======================================================================\n'; }

# elapsed seconds from sacct (works for running and finished jobs)
job_secs() {
  local e
  e=$(sacct -j "$1" --format=Elapsed -X -n 2>/dev/null | head -1 | tr -d ' ')
  [ -z "$e" ] && { echo ""; return; }
  # Elapsed is [DD-]HH:MM:SS
  python - "$e" <<'PY' 2>/dev/null
import sys
s = sys.argv[1]
d = 0
if "-" in s:
    d, s = s.split("-", 1)
    d = int(d)
p = [int(x) for x in s.split(":")]
while len(p) < 3:
    p.insert(0, 0)
print(d*86400 + p[0]*3600 + p[1]*60 + p[2])
PY
}
fmt() { [ -z "${1:-}" ] && { echo "--:--:--"; return; }; printf '%02d:%02d:%02d' $(($1/3600)) $(($1%3600/60)) $(($1%60)); }

hdr "0. JOB STATUS"
sacct -j "$X2_JOB,$MA_JOB" --format=JobID%14,JobName%26,State%12,Elapsed,MaxRSS%10,AllocCPUS 2>/dev/null \
  | grep -vE '\.extern|\.batch$' || echo "sacct unavailable"
echo
echo "  xtea2  log: $X2_LOG"
echo "  master log: $MA_LOG"
for f in "$X2_LOG" "$MA_LOG"; do
  [ -f "$f" ] || echo "  NOT YET WRITTEN: $f"
done

hdr "1. WALL CLOCK -- the headline"
X2_S_SEC=$(job_secs "$X2_JOB"); MA_S_SEC=$(job_secs "$MA_JOB")
printf '%-46s %10s\n' "master, 3 CONCURRENT 8-core jobs (frozen golden)" "38:31"
printf '%-46s %10s\n' "  (summed job time of those three)" "1:44:52"
printf '%-46s %10s\n' "master, SERIAL, one 24-core job [$MA_JOB]" "$(fmt "$MA_S_SEC")"
printf '%-46s %10s\n' "xtea2, one 24-core job, pre-trim (48110538)" "1:18:23"
printf '%-46s %10s\n' "xtea2, after trim steps 1-4 [$X2_JOB]" "$(fmt "$X2_S_SEC")"
echo
if [ -n "$X2_S_SEC" ] && [ -n "$MA_S_SEC" ] && [ "$MA_S_SEC" -gt 0 ]; then
  python - "$X2_S_SEC" "$MA_S_SEC" <<'PY'
import sys
x, m = int(sys.argv[1]), int(sys.argv[2])
r = m / x if x else 0
print(f"  xtea2 vs master-serial: {r:.2f}x  ({'xtea2 FASTER' if r > 1 else 'xtea2 SLOWER'})")
print(f"  xtea2 vs golden 38:31 : {2311/x:.2f}x  (golden had 3-way job concurrency)")
if r >= 1.0:
    print("  => on equal footing xtea2 is at parity or ahead. The remaining gap to the")
    print("     golden's 38:31 is a CONCURRENCY deficit, not slower code: recovering it")
    print("     means using more of the 2400% available, not rewriting algorithms.")
else:
    print("  => xtea2 is slower than master even serially. This is a real regression;")
    print("     use section 2 to find which phase/step owns it.")
PY
else
  echo "  (need both jobs to have elapsed time before the ratio is meaningful)"
fi
echo
echo "  NOTE: sacct Elapsed for a RUNNING job is time-so-far, not the final number."

hdr "2. WHERE THE TIME GOES"
echo "--- xtea2 phases (non-overlapping; these sum to the run) ---"
if [ -f "$X2_S/step_timings.tsv" ]; then
  python - "$X2_S/step_timings.tsv" <<'PY' || true
import sys, collections
rows=[l.split("\t") for l in open(sys.argv[1]).read().splitlines()[1:] if l.strip()]
ph=[(n,float(s),float(g)) for k,n,s,g in rows if k=="phase"]
st=[(n,float(s)) for k,n,s,g in rows if k.startswith("step")]
tot=sum(s for _,s,_ in ph)
for n,s,g in ph:
    print(f"  {n:<22}{s/60:8.2f} min {100*s/tot if tot else 0:5.1f}%  peak-so-far {g:5.2f} GB")
print(f"  {'TOTAL':<22}{tot/60:8.2f} min")
agg=collections.defaultdict(float)
for n,s in st:
    agg[n.split(':')[-1]] += s
print("\n  xtea2 per-step-label totals across types (types run CONCURRENTLY -- these")
print("  sum to more than the wall clock, by design; use them to rank, not to add):")
for k,v in sorted(agg.items(), key=lambda kv:-kv[1]):
    print(f"    {k:<16}{v/60:8.2f} min")
PY
else
  echo "  MISSING $X2_S/step_timings.tsv (written as the run progresses)"
fi
echo
echo "--- master per-step-label totals (types run SEQUENTIALLY -- these DO sum) ---"
if [ -f "$MA_OUT/step_timings.tsv" ]; then
  awk -F'\t' 'NR>1{s[$3]+=$4; n[$3]++} END{for(k in s) printf "    %-16s%8.2f min  (%d steps)\n", k, s[k]/60, n[k]}' \
    "$MA_OUT/step_timings.tsv" | sort -k2 -rn || true
  echo
  echo "  master per type:"
  awk -F'\t' 'NR>1{s[$1]+=$4} END{for(k in s) printf "    %-6s%8.2f min\n", k, s[k]/60}' \
    "$MA_OUT/step_timings.tsv" | sort || true
else
  echo "  MISSING $MA_OUT/step_timings.tsv"
fi
echo
echo "--- head-to-head on the steps BOTH run (master's clip+disc have no xtea2"
echo "    equivalent: xtea2 replaced them with the C++ scan pass, so they are"
echo "    compared against xtea2's scan_pass + cns_fan phases instead) ---"
if [ -f "$X2_S/step_timings.tsv" ] && [ -f "$MA_OUT/step_timings.tsv" ]; then
  python - "$X2_S/step_timings.tsv" "$MA_OUT/step_timings.tsv" <<'PY' || true
import sys, collections
x=collections.defaultdict(float); m=collections.defaultdict(float); xph={}
for l in open(sys.argv[1]).read().splitlines()[1:]:
    if not l.strip(): continue
    k,n,s,g=l.split("\t")
    if k=="phase": xph[n]=float(s)
    elif k.startswith("step"): x[n.split(':')[-1]]+=float(s)
for l in open(sys.argv[2]).read().splitlines()[1:]:
    if not l.strip(): continue
    t,i,lbl,s=l.split("\t"); m[lbl]+=float(s)
print(f"    {'step':<16}{'xtea2':>10}{'master':>10}{'ratio':>9}")
for k in ["cns_filter","transduction","sibling","postF_all","postF_hc","gene","genotype","gvcf"]:
    if k in x or k in m:
        xa, ma = x.get(k,0)/60, m.get(k,0)/60
        r = f"{ma/xa:.2f}x" if xa else "--"
        print(f"    {k:<16}{xa:9.2f}m{ma:9.2f}m{r:>9}")
cand_m = (m.get("clip",0)+m.get("disc",0))/60
cand_x = (xph.get("scan_pass",0)+xph.get("cns_fan",0)+xph.get("disc_pos_load",0))/60
if cand_m or cand_x:
    r = f"{cand_m/cand_x:.2f}x" if cand_x else "--"
    print(f"\n    candidate generation (the actual xTea2 rewrite):")
    print(f"    {'  scan vs clip+disc':<16}{cand_x:9.2f}m{cand_m:9.2f}m{r:>9}")
PY
else
  echo "    (need both timing files)"
fi

hdr "3. MEMORY"
printf '%-46s %12s\n' "golden L1 / Alu / SVA (8 cores each)" "9.17/11.30/13.68 GB"
for pair in "xtea2 [$X2_JOB]:$X2_JOB" "master-serial [$MA_JOB]:$MA_JOB"; do
  lbl="${pair%%:*}"; j="${pair##*:}"
  rss=$(sacct -j "$j" --format=MaxRSS -n 2>/dev/null | tr -d ' K' | grep -E '^[0-9]+$' | sort -rn | head -1)
  if [ -n "$rss" ]; then printf '%-46s %9.2f GB\n' "$lbl" "$(echo "scale=4; $rss/1048576" | bc)"
  else printf '%-46s %12s\n' "$lbl" "pending"; fi
done
echo
echo "  master-serial runs one step at a time, so its peak is the single heaviest step."
echo "  xtea2 holds three concurrent type chains, so its peak is inherently higher --"
echo "  that is a design consequence of --parallel-downstream, not a leak. The like-for-"
echo "  like read is xtea2 vs the LARGEST golden job (13.68 GB), not vs their sum."

hdr "4. INTERMEDIATE FILE SIZES"
printf '%-30s %12s %10s\n' "tree" "size" "files"
for pair in "xtea2:$X2_S" "master:$MA_S"; do
  lbl="${pair%%:*}"; d="${pair##*:}"
  if [ -d "$d" ]; then
    printf '%-30s %12s %10s\n' "$lbl total" "$(du -sh "$d" 2>/dev/null | cut -f1)" "$(find "$d" -type f 2>/dev/null | wc -l)"
  else
    printf '%-30s %12s %10s\n' "$lbl total" "absent" "-"
  fi
done
echo
for t in L1 Alu SVA; do
  xd="$X2_S/$t"; md="$MA_S/$t"
  xs=$([ -d "$xd" ] && du -sh "$xd" 2>/dev/null | cut -f1 || echo "-")
  ms=$([ -d "$md" ] && du -sh "$md" 2>/dev/null | cut -f1 || echo "-")
  printf '  %-6s xtea2 %10s   master %10s\n' "$t" "$xs" "$ms"
done
echo
echo "--- xtea2-only artefacts (the C++ scan path's disk cost; master has no analogue) ---"
if [ -d "$X2_S" ]; then
  find "$X2_S" -maxdepth 1 -type f -printf '%s\t%p\n' 2>/dev/null | sort -rn \
    | awk -v r="$X2_S/" '{s=$1;p=$2;sub(r,"",p); if (s>10485760) printf "    %9.1f MB  %s\n", s/1048576, p}' || true
else echo "    (absent)"; fi
echo
echo "--- master-only: shared pub_clip ---"
du -sh "$MA_S/pub_clip" 2>/dev/null || echo "    (absent)"
echo
echo "--- deletable scratch in each ---"
du -sh "$X2_S"/*/tmp 2>/dev/null | sed 's/^/    xtea2  /' || true
du -sh "$MA_S"/*/tmp 2>/dev/null | sed 's/^/    master /' || true

hdr "5. CALL SETS -- reported, never gated"
echo "The goldens predate the trim step 2 B8 fix, which deliberately drops slightly"
echo "more calls (min-divergence rule; strictly more conservative). xtea2 differing"
echo "from master here is EXPECTED for this pass -- the user accepted a shifted"
echo "baseline. master-serial, by contrast, should reproduce the golden closely; if it"
echo "does not, the comparator is not measuring what the golden measured and the"
echo "timing above needs a caveat."
echo
printf '%-6s %9s %9s %9s %9s\n' type xtea2 master golden x2-vs-ma
for pair in "L1:LINE1:203" "Alu:ALU:1262" "SVA:SVA:73"; do
  t="${pair%%:*}"; rest="${pair#*:}"; alt="${rest%%:*}"; g="${rest##*:}"
  xf=$(ls "$X2_S"/merged/*_${alt}.vcf 2>/dev/null | head -1)
  mf=$(ls "$MA_S/$t"/*_${alt}.vcf 2>/dev/null | head -1)
  xn=$([ -n "$xf" ] && grep -vc '^#' "$xf" || echo NA)
  mn=$([ -n "$mf" ] && grep -vc '^#' "$mf" || echo NA)
  d="NA"; { [ "$xn" != NA ] && [ "$mn" != NA ]; } && d=$((xn - mn))
  printf '%-6s %9s %9s %9s %9s\n' "$t" "$xn" "$mn" "$g" "$d"
done
echo
echo "--- positions gained/lost, xtea2 vs master-serial (same input, same day) ---"
for pair in "L1:LINE1" "Alu:ALU" "SVA:SVA"; do
  t="${pair%%:*}"; alt="${pair##*:}"
  xf=$(ls "$X2_S"/merged/*_${alt}.vcf 2>/dev/null | head -1)
  mf=$(ls "$MA_S/$t"/*_${alt}.vcf 2>/dev/null | head -1)
  if [ -n "$xf" ] && [ -n "$mf" ]; then
    lost=$(comm -23 <(grep -v '^#' "$mf" | cut -f1,2 | sort -u) <(grep -v '^#' "$xf" | cut -f1,2 | sort -u) | wc -l)
    gain=$(comm -13 <(grep -v '^#' "$mf" | cut -f1,2 | sort -u) <(grep -v '^#' "$xf" | cut -f1,2 | sort -u) | wc -l)
    echo "  $t: xtea2 -$lost / +$gain vs master (B8 predicts lost >> gained)"
  else
    echo "  $t: need both VCFs (xtea2=${xf:-missing} master=${mf:-missing})"
  fi
done

hdr "6. WHAT TO DO WITH THIS"
cat <<'EOF'
  * The one number that matters: section 1's "xtea2 vs master-serial" ratio. That is
    code-vs-code on equal footing. Everything else is diagnosis.
  * If xtea2 >= 1.0x there, the open work is CONCURRENCY (925% of 2400% used), and
    section 2's phase table says which phase is serial when it need not be.
  * Trim step 5 (dead locals) is still NOT done -- see
    docs/handoff/HANDOFF_2026-07-26_wgs_run.md. Neither is re-baselining the goldens,
    which is what would turn section 5 from "reported" into a real gate.
  * Both runs read frozen git worktrees, so the live checkout was safe to edit while
    they ran. Clean up when done:
      git worktree remove /n/scratch/users/c/cos689/xtea2_testdata/step45_wt
      git worktree remove /n/scratch/users/c/cos689/xtea2_testdata/master_cmp_wt
EOF
echo
