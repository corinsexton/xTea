# xTea2 handoff — 2026-07-26 (WGS run session)

**This is the current entry point.** Supersedes `HANDOFF_2026-07-25_late_night.md` for state.

One-line status: **Trim step 4 is DONE (scan-only, −2,573). Timing/memory instrumentation added.
The full-genome NA12878 run is submitted as job `48824433` from a frozen git worktree.
Step 5 (dead locals) is NOT done — the session ran out of budget and the run was launched
without it, deliberately.**

## Check this first

```bash
/n/scratch/users/c/cos689/xtea2_testdata/compare_runs.sh     # <-- one command, both runs
```

That prints the whole side-by-side: wall clock, per-phase/per-step breakdown, memory,
intermediate file sizes, and call sets. It works while the jobs are still running and says
plainly what is not available yet, so there is nothing to time correctly.

**Two jobs, not one.**

| job | what |
|---|---|
| `48824433` | xTea2 after trim steps 1-4, one 24-core job |
| `48825778` | **master `46520ab`, SERIAL, one 24-core job — the true comparator** |

Raw logs if needed:
```bash
sacct -j 48824433,48825778 --format=JobID%14,JobName%26,State%12,Elapsed,MaxRSS%10,AllocCPUS
tail -120 /n/scratch/users/c/cos689/xtea2_testdata/step45_wgs_48824433.out
tail -120 /n/scratch/users/c/cos689/xtea2_testdata/master_serial_wgs_48825778.out
```

The job prints its own report at the end: timing vs the golden, memory, intermediate file
sizes, and the call-set diff. **A non-zero exit / "FAILED" in `sacct` is expected and means
nothing here** — the last thing the script runs is the concordance gate, which asserts
`>= 1.0` against a baseline the B8 fix deliberately moved. Read the report, not the status.
Expect the whole run to take roughly 1h20m based on history.

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

**So xTea2 is currently ~2× slower in wall clock than the code it replaced, at ~1.4× the peak
memory of the largest golden job.** Do not let the per-component speedups obscure this. The
likely reason is shape, not algorithm: master got 3 processes × 8 cores *free* from the
scheduler, while xTea2 is one process that has to manufacture the same concurrency
(`--parallel-downstream`), and 925% of an available 2400% means most of the allocation sits
idle. That is the number to attack next, and this run's per-phase breakdown is the first data
that can point at which phase to attack.

## Why there is a second job

The golden is **not** a fair comparator for xTea2's runtime, and this had been going
unnoticed. It was three per-type jobs the scheduler ran **concurrently** at 8 cores each, so
its "38:31" is `max(37:15, 29:06, 38:31)` over a 24-core *total* allocation. xTea2 is one
process with one 24-core allocation. Comparing the two charges xTea2 for concurrency SLURM
handed master for free — it is not a code-speed measurement at all.

Job `48825778` removes the asymmetry: master's own code, same input, same allocation, three
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

**How to read the result.** If `48825778` lands near or above xTea2's 1:18:23, xTea2 is at
parity or ahead on equal footing, and the remaining gap to 38:31 is a **concurrency** deficit
(xTea2 has peaked at 925% of an available 2400%) rather than slower code — which points at a
completely different fix. If it lands well below, the trim has a real regression and the
per-step-label table in section 2 says where.

## Outstanding

1. **Step 5 (dead locals) is not done.** 184 locals, clustered in `x_orphan_transduction` (55),
   `x_post_filter` (24), `x_clip_disc_filter` (19), `x_transduction` (10), `x_somatic_calling`
   (10). Each RHS must be checked by hand for side effects — deleting `x = f()` where `f`
   writes a file changes behavior. Re-derive the list against the post-step-4 tree.
   Also fold in the **11 `parse_option` dests that were already dead before step 4**
   (`af1`, `af2`, `analysis`, `contig_realign`, `copies`, `duplication`, `iniclip`, `spectrum`,
   `sv`, `trace`, `visualization`) — trivially safe, deliberately left out of step 4 so its
   diff stayed attributable.
2. **Re-baselining the goldens is still owed** (step-2 checklist item 8). The user accepted a
   shifted baseline *for this pass only*.
3. **Step 3 (attribute B7) was skipped**, per the same decision. If it is picked up: `TRIM_PLAN`
   says master lives at `/home/cos689/src/xtea` — **that path does not exist**. The real
   untouched checkout is `/home/cos689/xtea_exploration/xTea/` (what
   `cc3_master*_isolation.sh` already runs against). Note it is on branch `cli` with one
   modified file.
4. **`pyflakes` was pip-installed into the `xtea2` env** for this session's analysis. It is a
   dev-only tool with no runtime effect, but the env is pinned as `envs/xtea2.yml`, so either
   `pip uninstall pyflakes` or add it to the spec.

## Where things live

| | |
|---|---|
| the run | `/n/scratch/users/c/cos689/xtea2_testdata/step45_wgs.sbatch`, job `48824433` |
| its worktree (frozen at `57f433a`) | `/n/scratch/users/c/cos689/xtea2_testdata/step45_wt` |
| its output | `/n/scratch/users/c/cos689/xtea2_testdata/step45_wgs_out` |
| per-step timings | `<out>/NA12878.wgs/step_timings.tsv` |
| star-import checker | `docs/handoff/check_star_import_safety.py` |
| the plan | `docs/handoff/TRIM_PLAN.md` |

**The live checkout is safe to edit while `48824433` runs** — it reads the worktree, not this
tree. Do not `git worktree remove` it until the job finishes.
