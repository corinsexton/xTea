# xTea2 handoff — 2026-07-26 (WGS run session)

**This is the current entry point.** Supersedes `HANDOFF_2026-07-25_late_night.md` for state.

One-line status: **Trim step 4 is DONE (scan-only, −2,573). All four WGS jobs are COMPLETE.
The performance question is answered and one 31-minute fix is confirmed — see
[`PERF_2026-07-26_why_xtea2_is_slower.md`](PERF_2026-07-26_why_xtea2_is_slower.md), which is the
substantive document from this work. Step 5 (dead locals) is NOT done. One open question blocks
adopting the fix: it moves 15 calls.**

## Check this first

```bash
/n/scratch/users/c/cos689/xtea2_testdata/compare_runs.sh   # baseline vs master; no arguments
```

**Four completed jobs.** All NA12878 whole genome, one job, `-c 24 --mem=64G -p priopark`.

| job | config | wall | MaxRSS |
|---|---|---|---|
| `48825977` | **master `46520ab`, serial `-n 24` — the true comparator** | **50:20** | 40.6 GB |
| `48825976` | xTea2 baseline: `--cns-emit-scan` + `--parallel-downstream` | 1:31:21 | 24.3 GB |
| `48842713` | xTea2 **without `--cns-emit-scan`** (best known) | **59:56** | 22.9 GB |
| `48842714` | xTea2 without `--cns-emit-scan`, serial downstream | 1:02:19 | 22.1 GB |

The 64G ceiling was fine for everything: master peaked at 40.6 GB (it was unmeasured at `-n 24`
beforehand, hence the caution in earlier versions of this file) and every xTea2 config at ~23 GB.
No OOM.

`compare_runs.sh` still defaults to the `48825976` / `48825977` pair. For the fix runs, pass them
explicitly — `./compare_runs.sh 48842713 48825977` — or read their phase tables directly:

```bash
grep ^phase /n/scratch/users/c/cos689/xtea2_testdata/fix_nofan_out/NA12878.wgs/step_timings.tsv
grep ^phase /n/scratch/users/c/cos689/xtea2_testdata/fix_nofan_serial_out/NA12878.wgs/step_timings.tsv
```

**`sacct` "FAILED" on `48825976` means nothing** — the last thing that script runs is the
concordance gate, which asserts `>= 1.0` against a baseline the B8 fix deliberately moved.

## The headline, and a correction to what this file used to say

**xTea2 was 1.81× slower than master on equal footing. Dropping `--cns-emit-scan` recovers 31.4
minutes, taking it to 1.19×.**

Earlier versions of this file said the gap was "likely shape, not algorithm" — a concurrency
deficit, because xTea2 used 925% of an available 2400%. **That was wrong, and it was wrong in a
way that would have sent the next person after the wrong target.** Both pipelines ran at ~1000%
CPU (xTea2 910%, master 1036%); xTea2 was not idling, it was doing 1.60× the CPU *work*. The
cause was algorithmic: `--cns-emit-scan`'s mate-gather does
`sam_itr_queryi(idx, HTS_IDX_START, 0, 0)` — a whole-genome BAM pass — to service ~2,800 candidate
sites, where master fetches an indexed window per site. 34.96 min vs 5.47 min for the same result.

Also corrected: the trim did not cause this. xTea2's best pre-trim run was 1:18:23, already 1.56×
slower than master-serial. **xTea2 had been slower than the code it replaces for its entire
existence**, and the only comparator on file — the golden's three concurrent jobs — could not show
it. The `-n 24` serial comparator is what made it visible.

A second proposed fix (run types serially at `-n 24`) was **refuted by measurement**: 1:02:19,
2.4 min *slower*. Transduction does not scale linearly to 24 threads, so `--parallel-downstream`
should stay on. Details and the surviving weaker version of that idea are in the PERF doc.

## ⚠ The open question that blocks adopting the fix

Dropping `--cns-emit-scan` is **not call-neutral**:

| | L1 | Alu | SVA | merged |
|---|---|---|---|---|
| with `--cns-emit-scan` (`48825976`) | 199 | 1278 | 68 | 1545 |
| without (`48842713`, `48842714`) | 199 | 1264 | 67 | 1530 |

**15 calls move.** The two fix runs agree exactly with each other, so the threading change is
call-neutral and the cns emit path is what differs. Which of the two is *correct* is unknown.
Do not bank the 31 minutes until that is settled: compare `temp_clip.sam` / `temp_disc.sam`
between the runs at the 15 differing positions. This is independent of the B8 baseline shift,
which both runs carry equally.

## NEXT SESSION — start here

Work in this order. Item 1 gates a confirmed 31-minute speedup; item 2 is the only remaining
large perf target; items 3-5 are the trim's unfinished business. Nothing here is blocked on
anything outside the repo.

### 1. Decide which cns collection path is CORRECT (blocks a 31-min win)

Dropping `--cns-emit-scan` makes xTea2 31.4 min faster (1:31:21 → 59:56) but changes calls. The
churn is **23 positions, not 15** — 15 was the net. Full list committed at
[`runs/cnsemit_diff_positions.tsv`](runs/cnsemit_diff_positions.tsv):

- **19 present only with `--cns-emit-scan`** (lost when the flag is dropped) — 17 Alu, 1 L1, 1 SVA
- **4 present only without it** (gained) — 3 Alu, 1 L1
- **but one of those is a type flip, not a call change**: `chr12:99387860` is the SVA "loss" *and*
  the Alu "gain". Same locus, reassigned between types. So the real churn is 18 lost / 3 gained /
  1 retyped. Do not count it twice, and note it exercises move A.2's type-resolution rule rather
  than the cns path.

The two candidate outputs are both still on disk (`temp_clip.sam` ~105 MB each for Alu):

```bash
A=/n/scratch/users/c/cos689/xtea2_testdata/step45_wgs_out/NA12878.wgs      # with --cns-emit-scan
B=/n/scratch/users/c/cos689/xtea2_testdata/fix_nofan_out/NA12878.wgs       # without (indexed)
# for a differing Alu position, compare the reads each path collected:
for f in Alu/tmp/cns/temp_clip.sam Alu/tmp/cns/temp_disc.sam; do
  diff <(grep -c . $A/$f) <(grep -c . $B/$f)
done
```

**The question to answer:** at those positions, does the C++ fan collect reads master's indexed
fetch misses, or does it collect reads it should not? The bias is informative — 18 of 19 losses
means the fan is *more permissive*, so the likely story is that the whole-genome mate-gather finds
mates outside the `cns_extnd` window that master's per-site fetch never sees. If those mates are
legitimate evidence, the fan is more sensitive and the 31 min is a real trade, not free. If they
are noise, drop the flag and bank it.

Either way the decision is a **judgment call about sensitivity, not a bug hunt** — frame it that
way, and record which was chosen and why, because it changes the call set permanently.

### 2. Instrument `scan_pass`, then decide whether to optimise it

Now the largest single cost: **20-22 min of a ~60 min run**, of which only ~2 min is detection
(the C++ scan's spill files complete ~2 min in). The remaining ~18-20 min is bwa realignment.

It cannot be attacked yet, and that is the actual task: `realign.hpp` sends bwa's stderr to
`/dev/null` and has no internal timers, so **nothing inside that phase is observable.** Add, per
family: time to bwa spawn, time to first SAM record, total parse time, records parsed.

Only then test the hypothesis in the PERF doc — that `realign_one_file`'s parent parses bwa's SAM
single-threaded via `fgets` on a pipe, one reader per family, so three readers must consume what
24 bwa threads emit. **That is unverified.** If it holds, the fix is bwa → file + parallel parse,
or more bwa processes with fewer threads each.

Getting master's clip-step realign time separately would make this comparable; master's 15.63 min
clip step bundles its BAM scan and its realign together, so the two are not yet directly
comparable.

### 3. Trim step 5 — dead locals (the trim's last step)

See outstanding item 2 below for the counts and the RHS hazard. Re-derive the list against the
post-step-4 tree; step 4 deleted modules the old count included.

### 4. Re-baseline the goldens — but AFTER item 1, not before

This is owed from step 2 (checklist item 8) and has been deferred twice. **Sequence it after item
1 deliberately:** the cns-path decision changes the call set, so re-baselining first means doing it
twice. Re-baseline once, with the final cns config chosen.

When it happens, the goldens to regenerate are the WGS seeded set; `tests/test_downstream_from_cns.py`
and `tests/test_e2e_concordance.py` both compare against pre-B8-fix goldens today, which is why
every gate currently reports expected differences.

### 5. Housekeeping

- `pip uninstall pyflakes` from the `xtea2` env, or add it to `envs/xtea2.yml` — it was installed
  ad hoc this session and the env is meant to be a pinned spec.
- Trim step 3 (attribute B7) remains skipped by user decision. If picked up, note `TRIM_PLAN`'s
  master path is wrong: it is `/home/cos689/xtea_exploration/xTea/`.
- Both worktrees can be removed once items 1-2 are done (commands at the bottom of this file).
  Keep them until then — the fix runs used `step45_wt` and a different commit would not be
  comparable.

### What NOT to redo

- **Do not re-litigate the concurrency hypothesis.** Both pipelines run at ~1000% CPU. It was
  measured and refuted; see the correction below.
- **Do not run types serially at `-n 24`.** Measured 2.4 min *slower* than `--parallel-downstream`.
  Transduction does not scale linearly to 24 threads.
- **Do not re-establish the master baseline.** It is `48825977`, 50:20, and
  `master_serial_wgs.sbatch` reproduces it.

## What this session did

| commit | what |
|---|---|
| `254d0a5` | trim step 4: scan-only, **−2,573 net** |
| `57f433a` | per-phase / per-step timing + peak-RSS instrumentation in `run_xtea` |

`xtea/` is now **20,273 lines / 42 modules**, down from 32,963 / 58 (**−38.5%**).

### Step 4, and where it exceeded the checklist

Deleted: `x_TEI_locator.py` (1,315), `clip_read.ClipReadInfo` (1,062), `OneClipRead` + 5 pool
trampolines (44), the `-C`/`-D`/`--mosaic`/`-B`/`-Q` branches of `x_TEA_main` (178), 14 dead
`parse_option` entries (36), the clip/disc steps in `build_type_steps` (13).

The plan said "the `-C`/`-D` branches". Deleting only those would have left `x_TEI_locator`
alive through three other branches (`--mosaic` and `-Q` use `TE_Multi_Locator`, `-B` uses
`TELocator`), so the module could not die — which is what the plan's own −2,700 estimate
assumed. `--mosaic` does **not** block `MOSAIC_PLAN.md`: its M.-1 reference run is against
master's own checkout, and M.1 replaces this path with a `mosaic_adapter.py` over
xtea-scan's polyA tally rather than calling master's mosaic clip function at all.

`--scan-clip`/`--scan-disc` are **kept as accepted no-ops** (forced True, nothing branches on
them) because ~15 scripts under `/n/scratch/users/c/cos689/xtea2_testdata` pass them.
`--nclip`/`--cr`/`--nd` are **removed** — no caller passed them, and a silently-ignored cutoff
flag is worse than one that fails loud.

### The one verification worth reusing

Deleting a `from X import *` line cannot be checked with pyflakes: it sees the star-import,
prints "unable to detect undefined names", and reports nothing further — 19 such lines in
`x_TEA_main.py`, zero findings. So a `NameError` in a branch no test runs would surface only
when that branch ran. `docs/handoff/check_star_import_safety.py` builds the real namespaces
instead: of the 281 names `x_TEI_locator` exported, **6 were unique to it**, and exactly **2**
were read by `x_TEA_main` — `TE_Multi_Locator` and `TELocator`, at the 5 sites in the 5 deleted
branches. Reports SAFE after the deletions.

Other gates: step plan **byte-identical** (`--dry-run` produces exactly the pre-step-4 18-step
plan; and *without* the scan flags it now produces that same 18-step plan instead of the old
24-step one, which is the real assertion that scan-only is the only path) · import sweep 42/42
· fast suite 83 passed / 1 pre-existing skip, unchanged.

### The timing baseline, which had never been established

The comparison the run reports against — **the naive original golden WGS run** — was not
recorded anywhere before this session. It is jobs `48013142`/`48013143`/`48013144`, Jul 17
16:36, master's own per-type SLURM pipeline, `-c 8 --mem=192G` each:

| | elapsed | MaxRSS |
|---|---|---|
| L1 | 37:15 | 9.17 GB |
| Alu | 29:06 | 11.30 GB |
| SVA | 38:31 | 13.68 GB |

They ran **concurrently**, so the golden's wall clock is **38:31** across 24 total cores.
xTea2's best WGS to date is **1:18:23 at 925% CPU / 18.68 GB** (`48110538`, speed3).

> **⚠ SUPERSEDED.** This section used to conclude here that the gap was "shape, not algorithm"
> — a concurrency deficit. Measurement refuted that; see the headline section above and the PERF
> doc. Both pipelines run at ~1000% CPU. The cause was `--cns-emit-scan` doing a whole-genome BAM
> pass for ~2,800 point queries. The per-phase instrumentation added this session is what
> settled it, which is the one thing this section got right.

## Why there is a second job

The golden is **not** a fair comparator for xTea2's runtime, and this had been going
unnoticed. It was three per-type jobs the scheduler ran **concurrently** at 8 cores each, so
its "38:31" is `max(37:15, 29:06, 38:31)` over a 24-core *total* allocation. xTea2 is one
process with one 24-core allocation. Comparing the two charges xTea2 for concurrency SLURM
handed master for free — it is not a code-speed measurement at all.

Job `48825977` removes the asymmetry: master's own code, same input, same allocation, three
types **sequentially** in one job. That is the honest denominator.

Its command lines are **not retyped** — they are the golden's own three
`run_xTEA_pipeline.sh` files transformed by `sed` (paths, `XTEA_PATH`, `-n`). Retyping 30
invocations would have been a silent-divergence source, and a comparator is only worth
anything if exactly one thing differs. Exhaustively, what differs from the golden:
sequential instead of 3 jobs; `-n 8` → `-n 24` so each step gets the whole allocation;
`XTEA_PATH` points at a worktree of this repo's `master` (the golden's
`/home/cos689/src/xtea/xtea/` **no longer exists**, which is what made TRIM_PLAN's path stale);
fresh output dir. Cutoffs, annotations, reference, the py3.8 `xTea` env with deep-forest,
step order and `--sva` placement are byte-identical.

The faithful `-n 8` serial figure needs no run: it is 37:15 + 29:06 + 38:31 = **1:44:52** by
summing the golden's own jobs. `-n 24` measures the thing that is not already derivable.

**The result:** `48825977` landed at **50:20**, well below xTea2's 1:18:23 — so the second
branch of that test is what happened. It was not a concurrency deficit. The per-step table is
what located it; see the PERF doc.

## Outstanding

0. **Settle the 15-call difference** from dropping `--cns-emit-scan` (see the warning above).
   This gates a confirmed 31-minute win, so it is the highest-value open item, ahead of step 5.
1. **Instrument the clip realignment** — `scan_pass` is now the largest remaining cost (20-22 min
   of a 60-min run, of which only ~2 min is detection; the rest is bwa realign). It cannot be
   attacked yet because `realign.hpp` sends bwa's stderr to `/dev/null` and has no internal
   timers, so nothing inside that phase is observable. Add the timers before touching it. A
   hypothesis is recorded in the PERF doc (single-threaded parent parsing bwa's SAM through a
   pipe, one reader per family) and is explicitly **unverified**.
2. **Step 5 (dead locals) is not done.** 184 locals, clustered in `x_orphan_transduction` (55),
   `x_post_filter` (24), `x_clip_disc_filter` (19), `x_transduction` (10), `x_somatic_calling`
   (10). Each RHS must be checked by hand for side effects — deleting `x = f()` where `f`
   writes a file changes behavior. Re-derive the list against the post-step-4 tree.
   Also fold in the **11 `parse_option` dests that were already dead before step 4**
   (`af1`, `af2`, `analysis`, `contig_realign`, `copies`, `duplication`, `iniclip`, `spectrum`,
   `sv`, `trace`, `visualization`) — trivially safe, deliberately left out of step 4 so its
   diff stayed attributable.
3. **Re-baselining the goldens is still owed** (step-2 checklist item 8). The user accepted a
   shifted baseline *for this pass only*.
4. **Step 3 (attribute B7) was skipped**, per the same decision. If it is picked up: `TRIM_PLAN`
   says master lives at `/home/cos689/src/xtea` — **that path does not exist**. The real
   untouched checkout is `/home/cos689/xtea_exploration/xTea/` (what
   `cc3_master*_isolation.sh` already runs against). Note it is on branch `cli` with one
   modified file.
5. **`pyflakes` was pip-installed into the `xtea2` env** for this session's analysis. It is a
   dev-only tool with no runtime effect, but the env is pinned as `envs/xtea2.yml`, so either
   `pip uninstall pyflakes` or add it to the spec.

## Where things live

| | |
|---|---|
| **the performance findings** | **`docs/handoff/PERF_2026-07-26_why_xtea2_is_slower.md`** |
| the runs | `step45_wgs.sbatch` (`48825976`), `master_serial_wgs.sbatch` (`48825977`), `fix_nofan.sbatch` (`48842713`), `fix_nofan_serial.sbatch` (`48842714`) — all under `/n/scratch/users/c/cos689/xtea2_testdata/`, copies in `docs/handoff/runs/` |
| its worktree (frozen at `57f433a`) | `/n/scratch/users/c/cos689/xtea2_testdata/step45_wt` |
| its output | `/n/scratch/users/c/cos689/xtea2_testdata/step45_wgs_out` |
| per-step timings | `<out>/NA12878.wgs/step_timings.tsv` |
| star-import checker | `docs/handoff/check_star_import_safety.py` |
| the plan | `docs/handoff/TRIM_PLAN.md` |

All four jobs are finished, so both worktrees are now safe to remove:

```bash
git worktree remove /n/scratch/users/c/cos689/xtea2_testdata/step45_wt
git worktree remove /n/scratch/users/c/cos689/xtea2_testdata/master_cmp_wt
```

Keep them if you intend to act on outstanding item 0 or 1 — the fix runs used `step45_wt`, and
re-running against a different commit would not be comparable.
