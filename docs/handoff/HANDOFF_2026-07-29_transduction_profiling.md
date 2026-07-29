# xTea2 handoff — 2026-07-29 evening (transduction profiling: the 300-candidate assumption was wrong)

**This is the current entry point.** Supersedes `HANDOFF_2026-07-29_downstream_concurrency.md` for
anything about what to do next; that handoff's WGS three-way comparison, downstream-concurrency
analysis, and Alu-skip finding are still accurate and still the reference for how we got to
"transduction is the lever." What's new: the live job it left pending (`49075238`) finished, a new
profiling job (`49078859`) was built and run, and together they overturn that handoff's central
per-candidate cost estimate and pin down (most of) where the real time — and the real port-vs-master
gap — actually goes.

**One-line status: `49075238` confirms port really is slower than master (not a `-n8` concurrency
artifact — it holds at `-n24` too, isolated per family). But the previous handoff's `1481s / 300 ≈
4.9s per candidate` framing is now known to be wrong: direct instrumentation of
`parse_td_sibling_by_site` shows it costs ~0.05-0.1s/candidate and sums to under 30s total — the
"300" was the wrong denominator. The real cost lives earlier, in `call_candidate_transduction_v3`'s
raw-disc-site processing (L1 has 53,953 raw sites, not 300), and file-timestamp bounding localizes
essentially the entire port-vs-master gap to two sub-phases there. The mechanism inside those
sub-phases (why port is ~1.2-2x slower specifically there) is not yet found — that's next.**

## 1. Job `49075238` finished: port-slower is real, not a `-n8` artifact

Full results (real elapsed times from each leg's `time.txt` — `results.tsv`'s `elapsed_s` column is
still literal `mm` for every row, the awk bug flagged in the last handoff was never fixed and
doesn't matter now that the numbers below were pulled directly from `time.txt`):

| impl | family | threads | wall time |
|---|---|---|---|
| master | L1 | 8 | 18:46 |
| master | SVA | 8 | 22:59 |
| **port** | **L1** | **8** | **24:51** |
| **port** | **SVA** | **8** | **29:37** |
| master | L1 | 24 | 13:13 |
| master | SVA | 24 | 17:28 |
| **port** | **L1** | **24** | **18:53** |
| **port** | **SVA** | **24** | **21:19** |

The `-n24` leg (each family alone, full 24-core allocation, no cross-family contention at all) shows
the same direction as `-n8`: port +42.9% on L1, +22.0% on SVA. **This rules out "it's just `-n8`
concurrent-family contention" as the explanation** — the previous handoff's open question ("is this
real or a node-contention artifact") is resolved in favor of real, though the mechanism is now known
to be almost entirely inside `call_candidate_transduction_v3`, not `parse_td_sibling_by_site` (see
§3). The confound job (`train_fixed_proportions_persample.py`, another user, pid `3120866`, ~217%
CPU) was still present on `compute-p-17-178` throughout — it can't be ruled out as a shared drag on
absolute wall-clock, but it can't explain a *consistent directional* bias against port specifically,
since both arms shared the same node and ran through the same multi-hour window.

**Output equivalence: master and port outputs differ at every leg** (sha256 mismatch + `diff` reports
differences for L1/SVA at both n8 and n24). This includes master's own n8 output differing from
master's own n24 output on the same inputs — some result content is sensitive to thread count/timing
(likely dict/scheduling-order noise, not incorrectness). **This is expected and already tracked
separately** under the 07-28 handoff's "re-baseline the goldens" item (still not started, still the
most overdue item across all handoffs) — do not treat this as a new finding or start chasing it from
this angle.

## 2. Added: opt-in per-candidate phase timing in `parse_td_sibling_by_site`

`xtea/x_orphan_transduction.py` — **uncommitted**, in the working tree now. Gated behind
`XTEA_TD_PROFILE=1` (checked once at import time into `_B_TD_PROFILE`; zero overhead when unset).
When set, each call to `parse_td_sibling_by_site` writes one line to stderr:

```
TD_PROFILE\t{pid}\t{chrm}\t{insertion_pos}\t{n_reads_in_window}\t{t_open_s}\t{t_fetch_s}\t{t_cluster_s}\t{t_total_s}
```

`t_open` = `BamInfo`+`is_chrm_contain_chr()`+second `AlignmentFile` open; `t_fetch` = the
`samfile.fetch()` iterate loop; `t_cluster` = `_form_unique_cluster` x2 +
`is_disc_cluster_consist_with_clip_position`. **Note the field layout has `TD_PROFILE` itself as
`$1`** — an awk script that forgets this and starts counting from pid as `$1` will misparse every
column (this bit the profiling sbatch script's own summary table the first time; see §4).

The same 4-hunk instrumentation was also applied, identically, to a copy of master's
`x_orphan_transduction.py` (see §4) to get a true side-by-side. Master's version of this function
differs from port's only in dead-local removal and the already-verified `fall_in_region` boolean
trim (confirmed again this session, unchanged from the 07-28 finding) — the instrumentation hunks
line up in both trees.

This instrumentation is safe to keep or strip; it's inert unless the env var is set. Not committed
because it's debug-only and the pending `--keep-intermediates` diff (§6) shouldn't get entangled
with it in one commit.

## 3. The central finding: "300 candidates" was the wrong denominator

The 07-29-morning handoff computed L1's per-candidate transduction cost as `step_timings.tsv`'s
`transduction: 1481s` divided by `candidate_disc_filtered_cns.txt`'s 300 lines ≈ 4.9s/candidate, and
flagged this as unverified. **It's now verified, and it's wrong.** Profiling job `49078859` (below)
measured `parse_td_sibling_by_site` directly: **mean total ~0.053-0.095s/candidate, summing to
14-29 seconds across the entire candidate set** (L1 n=302→305, SVA n=274) — see the table in §4.
That function accounts for essentially none of the 1481s.

The reason: `parse_td_sibling_by_site` is called from `call_sibling_TD_from_existing_list`, which
runs on the **300-ish post-filter candidate list**. But that's the *last* stage of the
`--transduction` step. Earlier in the same step, `call_candidate_transduction_v3` calls
`collect_clipped_reads`/`collect_disc_reads`, which operate on the **raw, pre-filter discordant-site
list** — confirmed by direct line count:

```
L1/candidate_list_from_disc.txt.clip_sites_raw_disc.txt:  53,953 lines
SVA/candidate_list_from_disc.txt.clip_sites_raw_disc.txt: 53,962 lines
```

~180x more sites than the 300 used in the old back-of-envelope math. That single wrong denominator
is what made "reduce L1/SVA's per-candidate transduction cost" look like it pointed at
`parse_td_sibling_by_site` — it doesn't. The real work, and the real port-vs-master gap, is upstream
in `call_candidate_transduction_v3`.

## 4. Profiling job `49078859`: phase-level localization of the port-vs-master gap

To avoid disturbing `49075238`'s live trees, two **separate, isolated** copies were made and
instrumented (identical 4-hunk timing patch from §2 applied to both):

- `/n/scratch/users/c/cos689/xtea2_testdata/master_wt_profile` (copied from `master_cmp_wt` @
  `46520ab`)
- `/n/scratch/users/c/cos689/xtea2_testdata/port_wt_profile` (copied from `port_wt` @ `b44a85f`)

Run via `docs/handoff/runs/transduction_profile.sbatch` (new, untracked) — same staged inputs as
`49075238`, `-n 8`, L1+SVA concurrent per impl, `XTEA_TD_PROFILE=1`. Completed in 59:03.
**The sbatch script's own built-in summary table has an off-by-one awk bug** (started field-counting
at pid instead of at the literal `TD_PROFILE` in `$1`) — its printed "mean_open_s ≈ 465s,
mean_n_reads ≈ 77M" numbers are garbage; ignore the job's own stdout summary and re-derive from the
raw `stderr.txt` files as shown below.

**Corrected `parse_td_sibling_by_site` phase means** (re-derived with correct field indices,
`$5`=n_reads, `$6`=open, `$7`=fetch, `$8`=cluster, `$9`=total):

| impl | family | n | mean_open_s | mean_fetch_s | mean_cluster_s | mean_total_s | sum_total_s | mean_n_reads |
|---|---|---|---|---|---|---|---|---|
| master | L1 | 302 | 0.0509 | 0.0024 | 0.0000 | 0.0533 | 16.10 | 465.1 |
| master | SVA | 274 | 0.0512 | 0.0025 | 0.0000 | 0.0537 | 14.71 | 479.4 |
| port | L1 | 305 | 0.0928 | 0.0026 | 0.0000 | 0.0955 | 29.11 | 462.4 |
| port | SVA | 274 | 0.0914 | 0.0026 | 0.0000 | 0.0940 | 25.77 | 481.0 |

Confirms §3: totals are seconds, not minutes. Side note: port's `open` phase (the double
`BamInfo.is_chrm_contain_chr()` + second `AlignmentFile` open — invariant per-BAM, redone every
call, present identically in both trees) is ~1.8x master's here too, but the absolute gap is only
~13-15s across the whole candidate set — real, but not where the minutes are.

**File-mtime bounding of `call_candidate_transduction_v3`'s intermediates** (`ls --time-style=full-iso`
on each run's `tmp/` dir) localizes essentially the entire gap:

| phase | master L1 | port L1 | gap L1 | master SVA | port SVA | gap SVA |
|---|---|---|---|---|---|---|
| bwa-index-done → `all_disc_pos.disc_pos` (`collect_clipped_reads`, over the ~54K raw sites) | 17:04 | 19:53 | +2:49 | 18:13 | 21:32 | +3:19 |
| clip.sam → disc.sam (`parse_transduction_clip_algnmt` + `collect_disc_reads`/`extract_mate_reads_by_name` + bwa disc realign) | 2:57 | 5:19 | **+2:22 (~1.8x)** | 4:50 | 7:11 | **+2:21 (~1.5x)** |
| disc.sam → cns realign + wrapup | ~5s | ~7s | negligible | ~7s | ~10s | negligible |

Sums: L1 gap ≈ 5:11 (total gap was 27:38−22:32=5:06); SVA gap ≈ 5:40 (total gap was
32:24−26:34=5:50). **These two phases account for essentially all of the port-vs-master difference**
in both families — nothing outside them needs explaining.

The second phase is the more interesting one: it's not just slower, it's **relatively** slower
(~1.5-1.8x) than the first phase's ~1.15-1.2x. That phase does per-read-name BAM mate lookups
(`extract_mate_reads_by_name` in `x_alignments.py`) and parses a ~2.4GB SAM file
(`parse_transduction_clip_algnmt`) — plausible candidates for an environment-level explanation
(pysam/htslib build differences between the `xTea` and `xtea2` conda envs) rather than an algorithmic
one, especially given master's and port's *code* for the earlier-verified functions is essentially
identical (§2). **Not yet confirmed — this is the next thing to instrument, not a conclusion.**

## 5. NEXT SESSION — start here

### 1. Instrument phase 2 specifically: SAM-parse vs mate-lookup vs bwa-subprocess time

Add the same opt-in timing pattern (§2) around three things inside the clip.sam→disc.sam window:
`parse_transduction_clip_algnmt` (SAM parsing), `collect_disc_reads`/`BamInfo.extract_mate_reads_by_name`
(`x_alignments.py:410`, per-read-name BAM mate lookup), and the `bwa_align.realign_disc_reads` bwa
subprocess call itself (just wrap with `time.time()` around the `cmd_runner` call, or diff bwa's own
reported time from stderr if it prints one). Whichever of the three carries the ~1.5-1.8x will tell
you whether to look at Python-side code, htslib/pysam version, or the bwa binary/version in each
conda env.

### 2. Compare pysam/htslib/bwa versions between the two conda envs

Given master's code for phase 2 is essentially the same as port's (only dead-local removal, per the
07-28 line-by-line diff), and phase 2 (not phase 1) carries the *bigger relative* slowdown, check
whether `xTea` and `xtea2` envs differ in pysam/htslib/bwa versions:
```
conda list -n xTea 'pysam|htslib|bwa'
conda list -n xtea2 'pysam|htslib|bwa'
```
This is the most direct way to test the "environmental, not algorithmic" hypothesis from §4 without
more profiling runs.

### 3. Re-baseline the goldens (carried over from 07-28 and 07-29-morning, still the most overdue item)

Still untouched. Now additionally supported by this session's finding that master/port outputs
differ at *every* leg of `49075238`, including master-vs-itself across thread counts — consistent
with, not contradicting, the existing "tens of sites expected to move" pre-work. See the 07-28
handoff for the full blast-radius numbers.

### 4. Commit the pending work

Still nothing committed. Three independent pieces sitting in the tree:
- `--keep-intermediates` (`xtea/run_xtea.py`, `tests/test_output_cleanup.py`, `README.md`) — finished,
  tested, from before this session even started. Safe to commit alone.
- The `XTEA_TD_PROFILE` instrumentation (`xtea/x_orphan_transduction.py`, §2) — debug-only, inert by
  default, this session's addition. Could go in its own commit, or be stripped if not wanted
  long-term now that its diagnostic job has already run.
- `docs/handoff/runs/transduction_attribution.sbatch` and the new
  `docs/handoff/runs/transduction_profile.sbatch` — both scripts, safe to commit together.

## What NOT to redo

- **Do not re-derive the per-candidate transduction cost from `step_timings.tsv` divided by the
  300-line filtered candidate list.** That arithmetic is now known to be wrong (§3) — the real
  denominator for `call_candidate_transduction_v3`'s read-collection passes is the ~54K-line raw
  disc-site list, not the post-filter list `parse_td_sibling_by_site` runs on.
- **Do not re-profile `parse_td_sibling_by_site` looking for the port-vs-master gap.** Directly
  measured this session (§4): ~14-29s total across the whole candidate set for either impl. It's not
  where the minutes are. The `open`-phase 1.8x there is real but too small in absolute terms (~13-15s)
  to matter.
- **Do not re-run `49075238`'s `-n8` vs `-n24` comparison to check whether port-slower is a
  concurrency artifact.** Answered this session: it holds at `-n24` (isolated per family), so it
  isn't just `-n8` cross-family contention. (Carried forward from the prior handoff's own
  "what NOT to redo" list, now confirmed rather than just argued from arithmetic.)
- **Do not chase the master-vs-port output differences from this angle.** Expected, tracked under
  "re-baseline the goldens" (07-28 handoff), not a new bug surfaced by this session.
- **Do not trust the profiling sbatch script's own printed summary table** (in
  `transduction_profile_49078859.out`) — off-by-one awk field bug (§4). Re-derive from raw
  `stderr.txt` `TD_PROFILE` lines using `$5`=n_reads, `$6`=open, `$7`=fetch, `$8`=cluster, `$9`=total
  (not `$4`/`$5`/`$6`/`$7`).

## Where things live

| | |
|---|---|
| the completed attribution job | `49075238`, results/logs at `/n/scratch/users/c/cos689/xtea2_testdata/transduction_attribution_49075238/` (22G) |
| the completed profiling job | `49078859`, `docs/handoff/runs/transduction_profile.sbatch`, results/logs at `/n/scratch/users/c/cos689/xtea2_testdata/transduction_profile_49078859/` (11G) |
| isolated instrumented trees (do not confuse with `49075238`'s live `master_cmp_wt`/`port_wt`) | `/n/scratch/users/c/cos689/xtea2_testdata/master_wt_profile` (86M), `/n/scratch/users/c/cos689/xtea2_testdata/port_wt_profile` (7.0M) |
| the `XTEA_TD_PROFILE` instrumentation | `xtea/x_orphan_transduction.py:parse_td_sibling_by_site` (uncommitted) |
| the raw disc-site counts that broke the old "300 candidates" assumption | `port_wgs_out/NA12878.wgs/{L1,SVA}/candidate_list_from_disc.txt.clip_sites_raw_disc.txt` (53,953 / 53,962 lines) |
| the read-collection code now implicated | `xtea/x_transduction.py:call_candidate_transduction_v3` (line 173), `xtea/x_clip_disc_filter.py:collect_clipped_reads`/`collect_disc_reads` (lines 3307/3333), `xtea/x_alignments.py:extract_mate_reads_by_name` (line 410) |
| the cns/perf findings that still stand | `PERF_2026-07-26_why_xtea2_is_slower.md` |
| the WGS three-way comparison, cns decision, C++ port | `HANDOFF_2026-07-28_cns_decision.md` |
| downstream-concurrency analysis, Alu-skip finding (both still accurate) | `HANDOFF_2026-07-29_downstream_concurrency.md` |
| `--keep-intermediates` implementation (uncommitted, unrelated to this session) | `xtea/run_xtea.py:cleanup_successful_run`, `tests/test_output_cleanup.py` |
