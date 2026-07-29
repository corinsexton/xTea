# xTea2 handoff — 2026-07-29 (downstream concurrency, the Alu skip, transduction as the next lever)

**This is the current entry point.** Supersedes `HANDOFF_2026-07-28_cns_decision.md` for anything
about what to do next; that handoff's cns-port history and WGS run table are still accurate and
still the reference for how xTea2 got here. What's new: the WGS three-way comparison it left
pending has been read, a downstream-concurrency fix was explored and **not** built (arithmetic
killed it before code did), and a real, previously-unexamined lever surfaced instead —
transduction's own per-candidate cost for L1/SVA, which is now ~93% of the whole downstream phase.

One-line status: **the C++ cns port clears the bar (memory/speed/disk) the user set against the
Python step — that question is closed. The "run downstream serially at -n 24" idea was correctly
declined (already refuted 07-28). A "shared core pool across types" idea was designed but abandoned
before being written, once real step timings showed its ceiling is ~15-20 seconds against a 9.4
minute gap. Alu's near-zero transduction cost turned out to be a deliberate, correct design choice,
not a bug or an imbalance — which reframes the open question as "is L1/SVA's transduction cost
itself reducible," unexplored. A live job (`49075238`) is mid-flight with an unexplained result:
port's transduction came out slower than master's at `-n 8`, which contradicts the "near-identical
code" finding below and needs to be resolved, not assumed away.**

## 1. The WGS three-way comparison is read, and the C++ port wins

From `compare_wgs_three_way.py --quick` (jobs: master `48825977`, python-cns `49065573`, cpp-cns
`49070030`):

| arm | elapsed | MaxRSS | output |
|---|---|---|---|
| master (pre-B10) | 50:20 | 40.6 GB | 6.0 G |
| python cns | 60:37 | 21.99 GB | 7.1 G |
| **cpp cns** | **54:39** | **20.56 GB** | 7.0 G |

Cpp beats the Python step on all three axes the user named (memory, speed, disk): −6 min wall
(−9.5%), −1.4 GB peak, disk flat. **This closes item 2 of the 07-28 handoff's NEXT SESSION list.**
The 07-28 handoff's own framing applies: "if it clears the bar, the Python collection path becomes
the thing to delete." It cleared the bar. Calls move a little (python→cpp: L1 +1/−0, Alu +14/−1,
SVA +1/−0), consistent with the known bwa tie-break order noise, though the +14 Alu is larger than
chr1 suggested and has not been individually checked — the read-level sections of
`compare_wgs_three_way.py` (run without `--quick`) haven't been run yet. **Item 1 of the 07-28
handoff (re-baselining the goldens) is still untouched and is still the most overdue item** — see
that handoff's pre-worked blast-radius numbers, unchanged by anything this session found.

## 2. Uncommitted: `--keep-intermediates` retention policy

Present in the tree at session start, finished and tested this session, **not yet committed**:

- `xtea/run_xtea.py`: `cleanup_successful_run()` — after a successful run, deletes everything in
  `sample_root` except the final merged VCF, unless `--keep-intermediates` is passed. All-or-nothing:
  refuses to touch anything if the merged VCF is missing, so failed/incomplete runs are never
  cleaned regardless of the flag.
- `tests/test_output_cleanup.py` — 2 tests (retains-only-final-vcf, refuses-when-vcf-missing), both
  passing.
- `README.md` updated with the flag and a size comparison (the completed WGS C++-scan tree was
  ~7 GB before cleanup; the merged VCF is ~1 MB).

Motivation: the 7 GB vs 1 MB gap seen in this session's own WGS runs. Nothing here depends on
anything else in this handoff — safe to commit independently whenever convenient.

## 3. Downstream concurrency: one idea declined, one designed then killed by arithmetic

The user asked to test running the three types' downstream chains **serially at `-n 24`** instead
of the current concurrent `--parallel-downstream` (`cores // ntypes` each). **This was declined
without running it**: it is exactly the test already run and refuted on 07-28
(`48842713` concurrent 34.90 min vs `48842714` serial-24 39.23 min — serial was *slower*, see
`PERF_2026-07-26_why_xtea2_is_slower.md`'s RESULTS section), and it's already on that handoff's own
"what NOT to redo" list. Downstream code (transduction/sibling/postF) is unchanged since then, so
there was no reason to expect a different answer on the C++-port pipeline.

Given three options (cost-weighted static split / re-run serial-24 anyway / shared worker pool),
the user picked **shared worker pool**: instead of a fixed `cores // ntypes` slice per type for the
whole chain, redistribute a finished type's cores to whichever types are still running.

**Design sketched (not written to any file):** a `SharedDownstreamBudget` object shared across the
three concurrent per-type threads — `cores()` returns `total // active_types`, recomputed fresh
every time any thread asks; `finished()` decrements `active_types` once a type's whole chain (through
`postF_hc`) completes. Each type would re-read the budget **before every step** (via
`build_type_steps(..., n_cores=dcores, stop_after=name, start_from=name)` to get one step's argv at
a time), rather than locking in one thread count for the whole chain the way `run_type_downstream`
does today. Bounded safely by construction — `active * (total // active) <= total`, no
oversubscription possible regardless of timing races.

**Killed before being written**, by reading `port_wgs_out/NA12878.wgs/step_timings.tsv` — the actual
per-step costs from the completed cpp WGS run (`49070030`), which nobody had looked at yet:

| type | cns_filter | transduction | sibling | postF_all | postF_hc | **type total** |
|---|---|---|---|---|---|---|
| Alu | 116s | 1.5s | 38s | 7s | 6s | **169s** |
| L1 | 55s | 1481s | 55s | 5s | 6s | 1602s |
| SVA | 57s | 1772s | 52s | 1s | 1s | **1884s** (= the whole phase) |

Alu really does finish at t≈169s while its 8 cores sit idle for the rest of the 31.4-minute phase —
the original diagnosis (07-26, root cause 2) was right about *that*. But **the redistribution can
only take effect at a step boundary**, and each type's `transduction` step is one subprocess spawned
once with a fixed `-n` — it cannot pick up more threads mid-flight. L1's transduction subprocess is
already running from t≈55s to t≈1536s; Alu's cores free up at t≈169s, deep inside that window, and
do nothing for it. The freed cores would only reach L1/SVA's next steps (`sibling`+`postF*`, ~60-110s
each) — arithmetic ceiling **~15-20 seconds**, against a 9.4-minute deficit. Not worth building.
**The only design that could reach the real deficit is merging all three types' candidate records
into one shared `pool.map` call inside `x_orphan_transduction.py` itself** (so idle worker capacity
gets picked up mid-transduction, not just between steps) — this touches correctness-sensitive
per-candidate code (type-specific `--fflank`/`--rtype`/`-r`/`-a` args would need reconciling inside
one call), not just driver orchestration, and was not attempted this session.

## 4. Alu's near-zero transduction cost is by design, not a bug

Before deciding whether to build the cross-type merge, the user asked whether Alu is *supposed* to
skip transduction. **Yes, confirmed in code:**

- `xtea/run_xtea.py:74` — `REP_TYPES["Alu"]["flank"] = None  # no transduction flank for Alu`
- `resolve_annot()` (`run_xtea.py:97`) turns `None` into the literal string `"null"` for `--fflank`.
- `x_TEA_main.py`'s `--transduction` branch gates on
  `os.path.isfile(sf_flank) == True  # for Alu and many others, there is no transduction` — `"null"`
  is never a real file, so Alu takes the empty branch.

This is biologically correct: Alu is a non-autonomous SINE with no internal Pol II promoter, so it
has no mechanism for 3′ read-through past its own poly-A signal the way L1 (internal promoter, weak
poly-A) and SVA (weak poly-A, hijacks L1 machinery) do. **Alu's 1.5s is doing the right amount of
work, not idling.** This reframes root cause 2 entirely: it was never a load-balancing problem across
three similar workloads — Alu has essentially none, and L1/SVA carry unavoidable, legitimate cost.
Cross-type pool-sharing (§3) would still help *if* built, since it's about reusing Alu's idle cores
for L1/SVA regardless of why Alu finishes early — but the real lever now looks like reducing L1/SVA's
own per-candidate transduction cost.

## 5. NEXT SESSION — start here

### 1. Resolve the live job first: `49075238`

Purpose: a controlled master-vs-port comparison, feeding **identical staged cns-filter inputs**
(from the completed cpp WGS run) to L1/SVA transduction only, at `-n 8` (both concurrent) and
`-n 24` (each serial), to isolate whether xTea2's transduction code itself is faster/slower/equal to
master's — independent of anything upstream. Script: `docs/handoff/runs/transduction_attribution.sbatch`
(fixed since the first attempt — `49073948` died on an awk `system` reserved-word collision in the
results-formatting `awk`, which the `time.txt` files bypass; **`results.tsv`'s `elapsed_s` column
still reads literal `mm` for every row**, a cosmetic awk bug — read `time.txt` per cell instead, or
fix the awk, before trusting that column).

State at handoff (job at ~1:05 elapsed of a 6h limit, ~4:55 remaining, healthy):

| impl | family | wall (`-n 8`) |
|---|---|---|
| master | L1 | 18:46 |
| master | SVA | 22:59 |
| **port** | **L1** | **24:51** |
| **port** | **SVA** | **29:37** |

**This is the opposite of what should happen, and it is not yet explained.** Section 4's earlier
finding in this same session (before the job had these numbers) traced this job's exact code path —
`--transduction` dispatches to `call_sibling_TD_from_existing_list` → `parse_td_sibling_by_site`, and
diffing that function line-for-line between the master and port trees shows **only** dead-local
removal (`working_folder`, `b_first`, `query_name`/`query_seq`/`query_quality`, `icnt`) plus the
already-verified-equivalent `fall_in_region` boolean-only trim (spec §4 B8c). There is no algorithmic
difference that should make the port slower, yet it was 6-6.5 minutes slower on both families.
**Do not assume this is real without ruling out confounds first**: an unrelated job
(`train_fixed_proportions_persample.py`, pid 3120866, a different user's workload judging by the
path, running the entire time this session at 200%+ CPU) was observed on the same node
(`compute-p-17-178`) throughout monitoring. Whether it caused the slowdown, or whether it's
reproducible without contention, is unknown — check node occupancy/`seff` for this job, and consider
whether the `-n 24` block (now running) reproduces the gap or narrows it before concluding anything.
The job is also mid-flight into the `-n 24` leg as of this handoff (master L1 just started, 24
workers observed at ~65% CPU each) — the `-n 24` numbers and the script's own sha256sum/diff
output-equivalence check are both still to come.

### 2. If the port-slower result holds up: profile transduction's per-candidate cost

L1 has 300 cns candidates (`port_wgs_out/NA12878.wgs/L1/candidate_disc_filtered_cns.txt`), SVA has
269, Alu has 2220 — Alu's candidate count is the *largest* of the three, yet (per §4) it does
effectively no transduction work while L1 and SVA do real work per site: **L1 ≈ 1481s / 300 ≈ 4.9s
per candidate; SVA ≈ 1772s / 269 ≈ 6.6s per candidate**, at `-n 8`. Together L1+SVA transduction is
~93% of the entire 31.4-minute `downstream_parallel` phase in the cpp WGS run — this is now the
single largest concentration of time in the pipeline outside of `scan_pass`'s bwa realignment.

Nothing has profiled what that per-candidate cost actually is. Two live possibilities, and the
difference matters for whether there's anything to fix:
- A fixable inefficiency analogous to the per-record RepeatMasker-tree rebuild the trim found and
  deleted (see 07-28 handoff, "best single find") — i.e. something redone per-candidate that doesn't
  need to be.
- Genuine bwa realignment against the 3kb flank library, done per site — in which case, like
  `scan_pass`'s realign (07-28 handoff, root cause 3), it may be near-irreducible compute, and the
  only lever would be fewer candidates or less permissive params, not a code fix.

**This was flagged as the recommended next step but not started.** Suggested first move: instrument
`parse_td_sibling_by_site` (or whatever it calls) with the same spawn/first-record/total/getline-vs-
body split that `realign.hpp` got in the 07-28 session, since that same instrumentation is what
turned root cause 3 from a hypothesis into a measurement there.

### 3. Re-baseline the goldens (carried over, still the most overdue item)

Unchanged from the 07-28 handoff — nothing this session touched the cns/B10 question. See that
handoff's "NEXT SESSION item 1" for the full pre-worked blast-radius numbers (tens of sites
expected to move, 95-site Alu inspection list). Still not done.

### 4. Commit the pending work

`--keep-intermediates` (§2) is finished, tested, and independent of everything else in this handoff
— safe to commit whenever. `docs/handoff/runs/transduction_attribution.sbatch` should be committed
alongside it (or once job `49075238` is fully resolved and the script's final form is known, in case
the `results.tsv` awk bug gets fixed first).

## What NOT to redo

- **Do not run downstream types serially at `-n 24` as a standalone test.** Refuted twice now:
  07-28 (`48842713` vs `48842714`, 34.90 vs 39.23 min) and by inference from job `49075238`'s
  `-n 24` leg (transduction-only, same code, currently running) — if that leg's serial numbers don't
  beat the concurrent `-n 8` numbers either, that's a third confirmation, not a new question.
- **Do not build the step-boundary "shared core pool" as designed in §3.** Arithmetic ceiling is
  ~15-20 seconds against a 9.4-minute deficit, derived from `49070030`'s actual `step_timings.tsv`
  (Alu's whole chain finishes in 169s; SVA's alone is 1884s = the entire phase). It was never
  written to any file — `run_xtea.py`'s diff is unchanged except for `--keep-intermediates`.
- **Do not treat Alu's near-instant transduction as a bug, an imbalance, or something to load-balance
  away.** `REP_TYPES["Alu"]["flank"] = None` is deliberate and biologically correct (§4). If a future
  cross-type pool merge (the one design that *could* reach the 9.4-min gap) gets built, it should
  still route through this — Alu contributing ~0 candidates to a merged pool is the correct
  behavior, not something to special-case around.
- **Do not assume the port's transduction code differs from master's.** Traced line-for-line for the
  exact function this session's live job exercises (`parse_td_sibling_by_site`) — only dead-local
  removal and two already-verified-equivalent optimizations (B8c, P1-EVALUATED). If job `49075238`'s
  port-slower result persists after ruling out node contention, the cause is environmental
  (interpreter/package versions between the two conda envs used, `xTea` vs `xtea2`) or something not
  yet found — not a known code difference.

## Where things live

| | |
|---|---|
| this session's live job | `49075238`, `docs/handoff/runs/transduction_attribution.sbatch`, output at `/n/scratch/users/c/cos689/xtea2_testdata/transduction_attribution_49075238/` |
| the WGS three-way comparison | `compare_wgs_three_way.py` (no arguments; `--quick` for calls+timing) — read this session, cpp wins |
| the cpp WGS run's real per-step costs | `/n/scratch/users/c/cos689/xtea2_testdata/port_wgs_out/NA12878.wgs/step_timings.tsv` |
| the cns/perf findings that still stand | `PERF_2026-07-26_why_xtea2_is_slower.md` |
| the cns decision, C++ port, re-baselining pre-work | `HANDOFF_2026-07-28_cns_decision.md` |
| the bug/optimization log (B8c, P1-EVALUATED cited this session) | `../clipdisc_contract_spec.md` §4 |
| `--keep-intermediates` implementation | `xtea/run_xtea.py:cleanup_successful_run`, `tests/test_output_cleanup.py` |
| Alu's transduction skip | `xtea/run_xtea.py:74` (`REP_TYPES`), `resolve_annot` (`:97`), master/port `x_TEA_main.py`'s `--transduction` branch |
