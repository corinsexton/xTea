# xTea2 handoff — 2026-07-29 late (the transduction gap is root-caused and fixed; validation in flight)

**This is the current entry point.** Supersedes `HANDOFF_2026-07-29_transduction_profiling.md` for
anything about what to do next. That handoff's phase-level localization is still correct and is what
made this session short — it had narrowed the gap to two sub-phases and named the right next two
experiments. Both were run. Both answered.

**One-line status: the port-vs-master transduction gap was never algorithmic. It is a per-BAM-open
constant that pysam 0.24/py3.10 pays at ~2.05x pysam 0.17/py3.8, multiplied by ~108K opens that the
code was doing once per candidate SITE for a value invariant per BAM. Fixed in `fb1d766` by a
per-process handle/meta cache, which removes ~5,112s of worker CPU per family on the port (~10.7 min
of wall at -n 8) and ~2,845s on master — i.e. it more than closes the ~5-6 min gap, and speeds up
both arms. A second, independent fix in the same commit stops splitting cores three ways for the
transduction step, since Alu has no transduction. chr1 runs clean and the core split verifiably takes
effect. The WGS timing/memory/disk run has NOT been submitted yet — it is gated on an A/B that proves
the cache is output-neutral.**

## 1. Root cause: a per-open constant times a per-site loop

Two measurements, and the agreement between them is the whole argument.

**(a) The env difference is real and it is in pysam, not xTea.** Single-threaded microbenchmark, no
xTea code involved, on the real 39 GB `NA12878.wgs.bam` (3,366 contigs):

| env | python | pysam | `AlignmentFile` open+close |
|---|---|---|---|
| `xTea` (master) | 3.8.15 | 0.17.0 | 0.0185s |
| `xtea2` (port) | 3.10.20 | 0.24.0 | **0.0380s (2.05x)** |

Newer htslib builds a full header dictionary on open, which scales with contig count — so this BAM's
3,366 contigs make the constant unusually large. Script:
`docs/handoff/scratchpad/bench_pysam_open.py` (reproduce with the two `conda activate`s).

**(b) That constant was being paid per candidate site.**
`XClipDisc.collect_clipped_disc_reads_by_region` (`x_clip_disc_filter.py:79`) is a `pool.map` task
run **once per raw disc site**, and the transduction step feeds it the RAW list — 53,953 sites for L1
on WGS, not the ~300 post-filter candidates (this is the denominator the previous handoff corrected).
Each call did:

- `BamInfo(...).is_chrm_contain_chr()` → **BAM open #1** (header only)
- `pysam.AlignmentFile(...)` → **BAM open #2**
- `_get_chrm_id_name(samfile)` → `get_tid()` for all 3,366 contigs

All three invariant per BAM. That is **~108K opens and ~181M `get_tid` calls** of pure overhead.
Measured whole-preamble cost per site: **0.0527s (master) vs 0.0948s (port)**.

**The agreement that confirms the attribution:** the previous handoff's in-situ profiling job
(`49078859`) independently measured the same preamble inside `parse_td_sibling_by_site` at **0.0509s
(master) / 0.0928s (port)** — within 4% of the microbenchmark, under real 8-way concurrency. And in
that function the preamble was **97% of total time** (0.0928 of 0.0955). Same preamble shape, same
cost, in the 54K-site loop.

Arithmetic against the isolated transduction baselines from `49075238`:

| family | port before | − cache saving | port projected | master (no fix) | outcome |
|---|---|---|---|---|---|
| L1 | 24:51 (1491s) | −639s | **~14:12** | 18:46 (1126s) | **274s faster than master** |
| SVA | 29:37 (1777s) | −639s | **~18:58** | 22:59 (1379s) | **241s faster than master** |

(639s = 53,953 sites × 0.0948s ÷ 8 workers.) The gap goes from +365s to −274s on L1 — the fix is
larger than the gap because it deletes work master was also doing.

## 2. What `fb1d766` changes

**Fix 1 — per-process BAM handle/meta cache.** `x_alignments.get_cached_bam` /
`get_cached_bam_meta`, applied at the two per-task sites: the clip/disc collection above, and
`BamInfo.extract_mate_reads_of_region` (per-bin, phase 2). Neither closes its handle any more —
both call sites carry a comment saying so, because closing a shared handle is the obvious way to
break this.

**Keyed on `os.getpid()`, and that is not decoration.** A pysam handle owns a file offset; one
inherited across `fork()` and used by both parent and child interleaves seeks and returns garbage.
The pid in the key means a forked child misses the cache and opens its own. Verified on the real BAM
(`docs/handoff/scratchpad/verify_bam_cache.py`, all three claims pass):
1. cached `b_with_chr` and tid map are **identical** to what the old code computed;
2. reused-handle fetches match fresh-open fetches across interleaved repeated regions;
3. forked children match a fresh open **even with the parent's cache already warm** — which is the
   dangerous case, and is why the test warms it deliberately before forking.

**Fix 2 — transduction gets its own core split.** `--parallel-downstream` gave every type
`cores // ntypes`. That is right for every step except transduction: `REP_TYPES["Alu"]["flank"] is
None`, so Alu's transduction is a ~1.5s no-op while L1's and SVA's are ~25 and ~30 min, and
transduction is ~93% of the whole downstream phase. Splitting 24 cores three ways parked a third of
the machine on nothing for ~25 minutes. Now: transduction alone splits over the **flanked** types
(12 each at -c 24), flankless types get 1. New param `n_cores_transduction` on `build_type_steps` /
`run_type_downstream`; 6 tests in `tests/test_transduction_core_split.py`.

**This is not the design that was rejected.** The rejected one redistributed cores at *step
boundaries* and was costed at a ~15-20s ceiling because a running subprocess cannot widen its own
`-n` mid-flight. This one is static — which types have a flank is known from `REP_TYPES` before
anything runs — which is exactly why deciding up front works where redistributing later did not.

**Known, deliberate wart:** the sum can briefly exceed `opt.cores` while a flankless type is still
in an earlier step (Alu's `cns_filter` runs to t≈116s) and L1/SVA are already in transduction on the
wider split. Bounded by the flankless types' pre-transduction chain, and much cheaper than idling
those cores for the ~25 minutes that follow. Documented at the call site.

## 3. Validation status — READ THIS BEFORE CLAIMING THE FIX WORKS

**chr1 (`49099643`) ran clean at `fb1d766`:** wall 4:50, MaxRSS 1.33 GB, disk 496 M with
intermediates / 24 K merged VCF. L1 transduction 101.8s, SVA 111.7s, Alu 0.8s. The core split
verifiably took effect (`transduction n=12` for L1/SVA, `n=1` for Alu, in the run log).

**chr1 cannot measure the win, and the reason is worse than "fewer sites" — MEASURED, not assumed.**
Two effects compound:
- 4,118 L1 raw disc sites vs WGS's 53,953 (13.1x fewer), and the saving is strictly per-site;
- **the chr1 BAM opens 7.8x cheaper than the WGS BAM** — 0.0050s vs 0.0386s, at the *same* 3,366
  contigs. So open cost tracks BAM/index size, not contig count, which was the earlier assumption.

Together that is ~100x less to save. chr1's theoretical cache saving is **~4s** (4,118 sites x
~0.012s preamble / 12 workers) — below run-to-run noise. And indeed the measured cache-only arms
came out 5% *slower* with the cache (L1 97.4 -> 102.4s, SVA 107.7 -> 112.5s), which is variance, not
signal. **Do not treat that as evidence the cache is slow, and do not treat chr1 as a perf gate for
this fix at all.** The 0.0948s/site figure underpinning the WGS projection was measured on the WGS
BAM, so it is unaffected — but it stays a projection until job `49101214` reports.

The chr1 speedup that IS real came from fix 2, not fix 1: parent at transduction `-n 8` vs
`fb1d766` at `-n 12` gave L1 115.8 -> 101.8s (-12.1%) and SVA 126.6 -> 111.7s (-11.8%).

**Output equivalence: ANSWERED, and it is clean.** Two comparisons, both on chr1:
- **Cache isolated** (arms B vs C, both transduction `-n 12`, L1+SVA): all 3 VCFs byte-identical and
  **all 70 per-type intermediates byte-identical, zero differences.** The cache is output-neutral.
- **Both fixes together** (parent 3-type `-n 8` vs `fb1d766` 3-type `-n 12`): **all 4 VCFs
  byte-identical**, call counts unchanged (L1 0 / Alu 26 / SVA 3 / merged 29). So even the
  thread-count change did not move chr1's calls, despite the known thread-count sensitivity.

chr1 at `fb1d766` gives
L1 0 / Alu 26 / SVA 3 / merged 29. Two things to keep straight:
- **`L1 = 0` on chr1 is correct** — master's own chr1 golden is also L1 0 (Alu 27, SVA 5). Not a
  regression, do not chase it.
- The last xtea2 chr1 run on file (`a2_chr1_cnsemit_out`) gave L1 2 / Alu 29 / SVA 4 / merged 35.
  **That delta is not attributable to this session** — dozens of known call-moving commits sit
  between (B10, the C++ cns port, the candidate-major clip FASTQ). It is also not *cleared* by this
  session. Hence the A/B below.

**Why an A/B is required rather than optional:** fix 2 changes thread counts, and this pipeline is
known to move calls with thread count — `49075238` found master's own `-n 8` output differing from
master's own `-n 24` output on identical inputs. So a plain parent-vs-`fb1d766` 3-type comparison
confounds the cache (must be output-neutral) with the core split (known ordering noise).

Three arms, `49100173` / `49100174` / `49100175` — **all three have now run; result above.**
(B and C exited non-zero, but only in the report section, *after* the pipeline had completed and
written every output: `set -o pipefail` plus `f=$(ls glob|head -1)` makes a non-matching glob exit 2,
and a failing command substitution in an assignment trips `set -e`. The Alu VCF does not exist in a
2-type run. Fixed; the outputs were valid and were compared as-is.)

| arm | commit | types | transduction `-n` | purpose |
|---|---|---|---|---|
| A | `6fb6247` (parent) | L1,Alu,SVA | 8 | baseline for the combined effect |
| B | `6fb6247` (parent) | L1,SVA | 12 | — |
| C | `fb1d766` (mine) | L1,SVA | 12 | **B vs C isolates the cache** |

B and C run transduction at the *same* thread count (parent: 24÷2 types; mine: 24÷2 flanked types),
so **B vs C must be byte-identical**. If it is not, the cache has a real bug and the WGS run must
wait. A vs the 3-type `fb1d766` run shows the combined effect, where residual diffs are the known
thread-count noise.

**WGS is submitted: job `49101214`** (gate B == C passed). Command:
```
sbatch --export=ALL,XTEA2_SCALE=wgs,XTEA2_WORKTREE=/n/scratch/users/c/cos689/xtea2_testdata/perf_wt \
       docs/handoff/runs/perf_validate.sbatch
```
Baselines it must be read against (three-way WGS table): master (pre-B10) 50:20 / 40.6 GB / 6.0 G;
python cns 60:37 / 21.99 GB / 7.1 G; cpp cns 54:39 / 20.56 GB / 7.0 G. Projected with both fixes:
downstream phase 31.4 → ~19.5 min, whole run **~42-43 min**.

### Two harness bugs that cost a round of jobs — do not reintroduce

1. **`sbatch --export=ALL,FOO=a,b,c` splits on every comma.**
   `--export=...,XTEA2_TYPES=L1,Alu,SVA` set `XTEA2_TYPES=L1` and tried to export `Alu` and `SVA` as
   variable *names*. All three arms silently ran **L1-only while printing success-shaped output** —
   the failure was only visible in the script's own `types :` echo line. The type list is now a
   **positional arg**. Any future comma-bearing value must be too.
2. **The parent worktree predates `--keep-intermediates`** (committed in `fb1d766`), so A/B arms
   against it exit 2 on an unrecognized flag. The script now probes `--help` rather than assuming.
   Generally: an A/B arm on an older commit cannot be assumed to accept the current flag set.

Also fixed: `set -e` plus a trailing `[ -n "$f" ] && printf` made a missing per-type VCF kill the
script *after* the report printed, so a fine run looked failed.

## 4. NEXT SESSION — start here

### 1. Read the WGS result (`49101214`) — the only remaining open item on this work
Output equivalence is already settled (§3). What WGS decides is whether the projected timing win is
real: L1 transduction ~14:12 and SVA ~18:58 against the 1481s/1772s baselines, downstream phase
31.4 → ~19.5 min, whole run ~42-43 min against master's 50:20 and the current 54:39. Also collect
MaxRSS and disk against 40.6 GB / 6.0 G (master), 20.56 GB / 7.0 G (cpp cns).

If the win does NOT materialise, the first thing to check is whether the per-site preamble is
actually as large in `collect_clipped_disc_reads_by_region` as it is in `parse_td_sibling_by_site`
(where it was 97% of the function). The chr1 result showed the projection is sensitive to the BAM's
open cost, which varies ~8x with BAM size — so verify against the WGS BAM's 0.0386s, not a general
assumption. If it does not hold, the structural items in §3 below are the fallback, and item 1 there
(filter earlier) is worth more than every constant-factor fix combined.

### 2. Re-baseline the goldens (carried over from 07-28, 07-29 morning, 07-29 evening — still the most overdue item)
Untouched again. This session made it more pressing, not less: the chr1 call counts moved relative to
`a2_chr1_cnsemit_out` and there is no current baseline to say which of those moves are intended. See
the 07-28 handoff for the pre-worked blast-radius numbers (tens of sites expected to move, 95-site
Alu inspection list).

### 3. The structural refactor of transduction (discussed, designed, not built)
The fixes in `fb1d766` are constant-factor. The step's shape is the real cost. Purpose, so the
design has a reason: 3′ transduction is when an L1's weak poly-A is skipped, transcription reads
through into downstream flanking sequence, and that flank is co-mobilized to the new insertion site
— so detection means asking, per candidate, whether its reads align to a *source element's 3′ flank*
rather than to the element consensus. Orphan transduction is the flank arriving with little or no L1
body. Alu is excluded because a non-autonomous SINE with no internal Pol II promoter cannot read
through — correct, not a gap.

Ranked, highest value first:

1. **Filter earlier — potentially worth more than everything else combined.** The step grinds 53,953
   raw disc sites to report on ~300 candidates. If any later filtering is computable from the
   disc-site file alone, moving it ahead of read collection is a >100x work reduction. **Needs a
   correctness answer first** — which filters qualify is not known. Answer this before building
   anything below it.
2. **Range-partition instead of site-partition.** Today: 54K pool tasks, `chunksize=1`, one random
   `fetch()` each, over heavily-overlapping ±extnd windows (so reads get re-read). Instead: coalesce
   the windows, shard the genome into N contiguous chunks of roughly equal *read volume*, one
   sequential sweep per shard, testing reads against a local interval index. `GenomeIntervalIndex`
   is already in the tree and is the right structure. Also load-balances by data, not site count.
3. **Stop round-tripping 108K temp files.** Each site writes two tiny files the parent re-reads and
   unlinks — on a shared filesystem this is the next biggest cost after the opens. One file per
   shard, or return records through the pool.
4. **One `Pool` for the step.** `collect_clipped_reads_of_given_list` creates a fresh `Pool` **per
   chromosome** (~25 fork storms). A single pool with an `initializer` that opens the BAM once per
   worker is the from-scratch version of the cache added here.
5. **Fuse the three sweeps** (clip, disc, mate gather) into one pass with three output streams — ~3x
   less BAM I/O.
6. **Drop the 2.4 GB SAM text intermediate** in `parse_transduction_clip_algnmt` (chr1 already shows
   `temp_transduction_clip.sam` at 171 M as the single largest file).

Precedent: this is the same shape as the C++ `cns_emit.hpp` port — candidate-major sorted emission
replacing per-site random gathers — so the pattern already exists in the codebase.

## What NOT to redo

- **Do not look for an algorithmic port-vs-master difference in transduction.** There isn't one. The
  code is master's (verified line-by-line across three handoffs); the gap was a pysam-version
  constant times a per-site loop (§1). Both measurements are on file and agree to within 4%.
- **Do not re-derive the per-candidate cost from the 300-line filtered candidate list.** Wrong
  denominator, corrected in the previous handoff and confirmed here: 53,953 raw disc sites.
- **Do not re-profile `parse_td_sibling_by_site` for the gap.** ~14-29s total across the whole
  candidate set either way. Its 97%-preamble share was useful only as corroboration for §1.
- **Do not treat chr1 timings as evidence about the cache.** ~100x less to save there: 13.1x fewer
  sites AND a 7.8x cheaper BAM open (0.0050s chr1 vs 0.0386s WGS, same 3,366 contigs — open cost
  tracks BAM/index size, not contig count). chr1's theoretical saving is ~4s, under noise; the
  measured cache-only arms landed 5% slower, which is variance. chr1 IS a valid gate for output
  equivalence (and it passed, byte-for-byte) and for fix 2's core split (-12%).
- **Do not re-run the chr1 A/B for output equivalence.** Done and clean: 3 VCFs + 70 per-type
  intermediates byte-identical with the cache isolated at matched `-n 12`, and all 4 VCFs identical
  with both fixes together.
- **Do not treat chr1's `L1 = 0` as a bug.** Master's own chr1 golden is L1 0.
- **Do not run downstream types serially at `-n 24`.** Refuted three times now (07-28 `48842713` vs
  `48842714`; `49075238`'s `-n 24` leg).
- **Do not build the step-boundary shared core pool.** ~15-20s ceiling. Fix 2 gets the same win
  statically, which is why it works.
- **Do not put comma-bearing values in `sbatch --export`.** §3.

## Where things live

| | |
|---|---|
| the fix | `fb1d766`; `xtea/x_alignments.py:get_cached_bam`/`get_cached_bam_meta`, `xtea/x_clip_disc_filter.py:collect_clipped_disc_reads_by_region`, `xtea/run_xtea.py` (core split) |
| the core-split tests | `tests/test_transduction_core_split.py` (6 tests) |
| the pysam microbenchmark | `docs/handoff/scratchpad/bench_pysam_open.py`, `docs/handoff/scratchpad/bench_preamble.py` |
| the cache equivalence + fork-safety proof | `docs/handoff/scratchpad/verify_bam_cache.py` (3 claims, all pass on the real BAM) |
| the validation harness | `docs/handoff/runs/perf_validate.sbatch` (`XTEA2_SCALE=chr1\|wgs`, types as `$1`) |
| chr1 run at the fix | `49099643`, `/n/scratch/users/c/cos689/xtea2_testdata/perf_chr1_out` |
| the A/B arms (complete, clean) | `49100173` (parent 3type), `49100174` (parent 2type), `49100175` (mine 2type); outputs `ab_parent_3type`, `ab_parent_2type`, `ab_mine_2type` |
| the WGS timing run (in flight at handoff) | `49101214`, output `/n/scratch/users/c/cos689/xtea2_testdata/perf_wgs_out` |
| frozen worktrees | `perf_wt` (`fb1d766`), `parent_wt` (`6fb6247`) — both with `src/csrc/xtea_scan` built via `PKG_CONFIG_PATH=/home/cos689/software/lib/pkgconfig make` |
| the phase-level localization that pointed here | `HANDOFF_2026-07-29_transduction_profiling.md` |
| downstream-concurrency analysis, the Alu-skip finding | `HANDOFF_2026-07-29_downstream_concurrency.md` |
| the cns decision, C++ port, re-baselining pre-work | `HANDOFF_2026-07-28_cns_decision.md` |
| the cns/perf findings that still stand | `PERF_2026-07-26_why_xtea2_is_slower.md` |
