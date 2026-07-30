# xTea2 handoff — 2026-07-29 late-morning (WGS timing run cancelled; transduction's INPUT SET is the new finding)

> **SUPERSEDED** by `HANDOFF_2026-07-30_hg002_acquired.md`, and more importantly by the plan
> `docs/handoff/DOWNSTREAM_OVERHAUL_PLAN.md` (written 11:52 the same day, an hour after this
> doc). **Every measurement below still stands and is still worth trusting.** What is obsolete is §5's
> agenda: the plan commits to a full downstream rewrite with HG002/SMHT005 benchmarking and **retires
> NA12878**, which makes the "if you want the cache measured" resubmission at the end of this doc
> pointless — it is a constant-factor fix to code slated for rewrite, measured against a retired
> comparator. Read this doc for the transduction *findings* (§3 and §4 feed the plan's Stage 1); read
> the plan for what to do.

**~~This is the current entry point.~~** Supersedes `HANDOFF_2026-07-29_transduction_fixed.md` for what to
do next. Everything that handoff established about the *fixes* in `fb1d766` still stands — root cause,
output equivalence, the core split — and none of it needs redoing. What changed is the agenda: the
WGS timing run was **deliberately cancelled mid-transduction at the user's request**, so the cache's
timing win is still unmeasured, and a separate finding surfaced that is probably worth more than the
cache: **transduction's input site set is family-agnostic, and the two flanked types grind almost the
same 50K sites independently.**

**One-line status: the two fixes in `fb1d766` are correct and output-neutral (proven on chr1), fix 2's
core split is measured (-12% chr1) and verified live at WGS scale, but fix 1's cache remains a
PROJECTION — the run that would have measured it was cancelled at 39:25 with L1/SVA transduction ~15
minutes in and no step timing emitted. Separately: transduction is fed 53,953 raw disc sites of which
53,945 are shared across all three types, because the raw-disc gate keys on family-independent
columns. L1 and SVA then collect reads over `.slct` sets of 49,042 and 53,569 sites that overlap in
48,290 — so one shared collection pass over the 54,321-site union would replace 102,611 per-site BAM
gathers with 54,321, a 1.89x cut in the dominant loop with NO filtering-correctness question attached.
That is the thing to discuss before writing more transduction code.**

## 0. State at handoff — read this first

| | |
|---|---|
| branch / HEAD | `xTea2` @ `ca1390d`, working tree clean before this doc |
| the fixes | `fb1d766` (unchanged, still correct — see §1) |
| tests | 82 passed, 18 skipped (skips are the unbuilt C scan bin, pre-existing) |
| **nothing in flight** | WGS job `49101214` **CANCELLED by user request** at 39:25 elapsed |
| its partial output | `/n/scratch/users/c/cos689/xtea2_testdata/perf_wgs_out` (L1/SVA transduction incomplete), log `perf_validate_49101214.out` |

The cancellation was intentional, not a failure. The job was healthy: correct commit, clean tree, all
three types, and reproducing the baseline to within 0.5% on every phase it completed.

## 1. What `49101214` did and did not establish before it was cancelled

**Established — the harness and fix 2 are sound.** Every unchanged phase reproduced the pre-fix
baseline almost exactly, which is what makes this a valid comparison rather than noise:

| step | baseline `49070030` | run `49101214` | Δ |
|---|---|---|---|
| `scan_pass` | 1203.15s | 1197.07s | −6.1s (−0.5%) |
| `disc_pos_load` | 4.36s | 4.50s | +0.1s |
| `cns_emit` | 40.62s | 41.15s | +0.5s |
| `L1:cns_filter` | 55.08s | 61.60s | +6.5s |
| `SVA:cns_filter` | 57.37s | 63.63s | +6.3s |
| `Alu:cns_filter` | 116.22s | 129.64s | +13.4s |
| `Alu:transduction` | 1.50s | 2.26s | +0.8s |
| `Alu:sibling` | 38.47s | 38.37s | −0.1s |
| children peak RSS | 7.19 GB | 7.19 GB | identical |

The `cns_filter` steps being 6–13s slower is **fix 2 working as designed, not a regression**: Alu's
`cns_filter` now runs on fewer cores because L1/SVA took 12 each for transduction. Trading 13s of Alu
`cns_filter` for 7 freed cores across ~20 minutes of transduction is the entire point.

Fix 2 verified live at WGS scale, from the run log:
```
Phase 2: 3 per-type downstreams CONCURRENT (n=8 each; transduction n=12 for 2 flanked type(s), n=1 for the rest)
  L1  (rtype 1): downstream (n=8, transduction n=12)
  Alu (rtype 2): downstream (n=8, transduction n=1)
  SVA (rtype 4): downstream (n=8, transduction n=12)
```
And `Alu:transduction` at 2.26s on 1 core confirms the premise fix 2 was built on — it is a no-op, so
giving it 1 core instead of 8 cost 0.8s and freed 7 cores for the ~20 minutes that follow.

**NOT established — the cache's timing win.** L1 and SVA transduction started at 10:24:20 and were
~15 minutes in when the job was cancelled at 10:39:36. Neither emitted a step timing. **Fix 1 remains
a projection.** The pass/fail table in §4.1 of the previous handoff is still the right gate; it just
has not been run.

**One arithmetic correction to that projection.** The previous handoff computed the saving as
53,953 sites ÷ 8 workers. Both inputs were wrong in the same direction: the loop runs on the
**49,042 selected (`.slct`) sites, not 53,953**, and transduction now runs at **n=12, not 8**. So the
cache's projected wall saving is 49,042 × 0.0948s ÷ 12 ≈ **387s**, not 639s. This does not change any
verdict threshold, but the projection was optimistic on both counts.

**Benign noise to expect, do not chase:** `Error happen at merge clip and disc feature step: chrM not
exist`, thousands of times. Pre-existing on the WGS path — the baseline `49070030` emits 4,395 of them
and ran to completion. chr1 emits none (no chrM in its disc list). Despite the word "Error" it aborts
nothing. Any monitor filter on `Error` must exclude this line or it will drown.

## 2. Which "master" is which — this cost real confusion, fix it in your head first

There are three reference points and they are **not** interchangeable. The previous handoff's
transduction table said "master 1126s", which is true but is not master's real-world time.

| label | job | what it actually is |
|---|---|---|
| "baseline" in `perf_validate.sbatch` output | `49070030` | **xtea2/port before this session's fixes** (C++ cns). Isolates *this session's* two fixes. NOT master. |
| **master, whole pipeline** | `48825977` | master @ `46520ab`, py3.8/pysam 0.17, **serial, one 24-core job**, `-n 24`/step. **The true end-to-end comparator: 50:18 (3018s).** |
| "master" in the transduction table | `49075238` | master's transduction code **fed the port's staged inputs** at `-n 8` → L1 1126s, SVA 1379s. A code-vs-code test at equal work. NOT a real master run. |

Master's own per-step numbers from `48825977`, for whole-pipeline comparison:

| | L1 | Alu | SVA |
|---|---|---|---|
| clip | 414s | 236s | 288s |
| disc | 135s | 168s | 131s |
| cns_filter | 87s | 156s | 85s |
| **transduction** | **448s** | 3s | **587s** |

Master's clip+disc totals 1372s (serial, `-n 24`); the port's C++ `scan_pass` replaces all of it at
1197s, so **the scan already beats master's clip+disc by ~13%**. Also on file: master's frozen golden
(3 concurrent jobs × 8 cores) was 38:31 wall for 1:44:52 of summed job time — faster in wall but
using the same 24 cores less efficiently. The fair single-job comparator is **50:18**.

Master's transduction (448s on 36,687 sites at n=24) vs the port's (1481s on 53,953 at n=8) is not a
like-for-like speed comparison — see §3, the port is fed ~1.5x more sites. Normalized, master spent
0.293 CPU-s/site; the port at the projected ~900s × 12 threads over 49,042 sites would be ~0.220.

## 3. THE NEW FINDING: transduction's input set is family-agnostic

### 3.1 The funnel, from this run's own log

```
L1 :disc  324797 peaks -> 2375 pass disc cutoff, 53953 raw-disc rows
Alu:disc  324813 peaks -> 9295 pass disc cutoff, 53960 raw-disc rows
SVA:disc  324831 peaks -> 1568 pass disc cutoff, 53962 raw-disc rows
```

- **Clip peaks are family-agnostic** (324,797 / 324,813 / 324,831 — within 34). One shared C++ scan
  pass produces them. Master's per-family clip step produced *different* lists per type
  (624,700 / 661,513 / 695,379 candidate clip lines).
- **The disc cutoff is the family filter** — 2,375 / 9,295 / 1,568, sharply type-specific. It is the
  mate-in-repeat test against the type's repeat annotation (`disc_adapter.py:162-168`, with
  `SVA_ANNOTATION_EXTND=200` for SVA). Its output is `candidate_list_from_disc.txt`, which
  `cns_filter` consumes.
- **So the family is decided BEFORE cns, not by cns.** cns narrows *within* a family; it does not
  choose one. (This corrects a natural assumption — worth stating because it changes where you would
  intervene.)

### 3.2 The raw-disc list bypasses that filter by construction — MEASURED

| raw disc **positions** | count |
|---|---|
| L1 | 53,953 |
| Alu | 53,960 |
| SVA | 53,962 |
| **shared by all three** | **53,945** |
| L1 positions not in Alu | **8** |

99.98% identical. The three files are not byte-identical and their md5s differ, but only in the
`nL`/`nR` count columns — the *site set* is the same.

**Why, in code.** `disc_adapter.py:200` appends to `raw_rows` **outside** the `if keep:` disc-cutoff
gate, and the raw-disc write gate at `disc_adapter.py:227-229` is:
```python
b_l = (nrawL >= i_half_cutoff) and l_ok
b_r = (nrawR >= i_half_cutoff) and r_ok
```
`nraw*` are pre-annotation raw counts and `*_ok` is cluster consistency. **The family-specific
`nL`/`nR` are written to the file but never gate it.** So raw_disc is pre-cutoff *and* pre-family by
construction.

**This is deliberate and must not be "fixed".** A 3′ transduction's inserted sequence is the *source
element's 3′ flank* — unique genomic sequence, not the repeat consensus — and an orphan transduction
carries little or no element body at all. The candidate list is built by asking the opposite question
("do the mates land in the family's repeat annotation?"), so **the sites transduction exists to find
are precisely the ones that fail the test used to build the candidate list.** Being rejected there is
a transduction signature, not evidence against one.

**This is measured, not just reasoned — the obvious optimization would delete 20% of the calls.**
On WGS L1 (`port_wgs_out`), the natural-looking shortcut is "just use the 300-site candidate list
instead of 50K raw sites." What that would cost:

| L1, WGS | sites |
|---|---|
| `candidate_disc_filtered_cns.txt` (the candidate list) | 300 |
| `candidate_disc_filtered_cns2.txt` (transduction output) | 322 |
| **cns2 positions NOT present in the candidate list** | **64 (20% of output)** |
| `candidate_sibling_transduction2.txt` (orphan/sibling) | 59 |
| **orphan/sibling positions NOT in the candidate list** | **59 — all of them** |

Those 64 sites come from the raw 50K list and nowhere else. Every orphan transduction call is absent
from the candidate list, which is exactly what the biology predicts.

**So: filtering `raw_disc` by family, or by the disc cutoff, or down to the candidate list, does NOT
qualify** as a safe "filter earlier" — it answers part of the §4.3-item-1 gate ("which filters qualify
— needs a correctness answer first") in the negative. That was the highest-ranked lever in the
previous handoff; this narrows it considerably and redirects the effort to §4.

### 3.3 Why the port does ~1.5x more transduction work than master

| raw disc sites | master `48825977` | port | ratio |
|---|---|---|---|
| L1 | 36,687 | 53,953 | 1.47x |
| Alu | 40,681 | 53,960 | 1.33x |
| SVA | 44,373 | 53,962 | 1.22x |

Master's raw-disc list inherited family specificity from its **per-family clip step**. The port's
single shared C++ scan pass removes that, so the port's list is the family-agnostic superset. **This
is traceable entirely to the scan replacing master's clip step — not to anything in the transduction
code.** It is consistent with the previous handoff's conclusion that transduction has no algorithmic
port-vs-master difference; the difference is upstream, in what transduction is *handed*.

Open question for the discussion, not answered here: **is the superset intended?** A larger input is
more sensitive (won't miss transductions whose sites master's clip step dropped) but costs 1.5x and
may admit false positives. Nobody has decided this deliberately; it fell out of the scan port.

## 4. THE LEVER: one shared collection pass across the flanked types (1.89x, no correctness question)

L1 and SVA each run `collect_clipped_reads_of_given_list` → `collect_clipped_disc_reads_by_region`
over their own `.slct` site set. Those sets nearly coincide:

| | sites |
|---|---|
| L1 `.slct` | 49,042 |
| SVA `.slct` | 53,569 |
| shared positions | 48,290 |
| **union** | **54,321** |
| **today's total gathers** | **102,611** |

So a single collection pass over the union would do **54,321 per-site BAM gathers instead of 102,611
— a 1.89x cut in per-site gather COUNT.** (Not exactly 2x: the `.slct` re-selection *is* type-specific,
because its repetitive-region filter uses the type's annotation, which is why 49,042 ≠ 53,569.)

### 4.1 Two things found while designing this that change how it must be built

**(a) The redundant gathers are currently absorbed by parallelism, so a lazy shared cache wins
nothing.** L1 and SVA transduction are two *concurrent subprocesses* on 12 cores each, walking nearly
the same site list in nearly the same order. Doing a gather twice therefore costs CPU, not wall time.
Wall today ≈ 53,569·t/12 ≈ 4,464·t. A "first type to reach a site computes it, the other waits on a
lock" design leaves 12 cores blocked and lands at ≈ 53,967·t/12 — **no improvement**. The win exists
only if the shared pass runs **once at full 24-core width**: 53,967·t/24 ≈ 2,249·t, i.e. **~1.98x off
the collection**. That forces it *outside* the per-type steps, into the parent before
`downstream_parallel`. Do not build the lazy-cache version.

**(b) `.slct` is NOT a subset of `raw_disc`, so the parent cannot know the exact union up front.**

| | L1 | SVA |
|---|---|---|
| `.slct` positions | 49,042 | 53,569 |
| `raw_disc` positions | 53,953 | 53,962 |
| **`.slct` positions absent from `raw_disc`** | **206** | **151** |
| union of `raw_disc` (L1∪SVA) | 53,967 | |
| union of `.slct` (L1∪SVA) | 54,321 | |

The re-selection (`re_slct_with_clip_raw_disc_sites`) can introduce positions, and it runs *inside*
each type's transduction step using that type's annotation. So either the step gets split into
re-select → collect → rest (bigger change), or the parent collects over `union(raw_disc)` = 53,967
sites (known before any type runs — the parent writes those files) and each type computes its own
~150-200 residual sites itself. **The latter needs no step split and still gets ~99.6% of the win.**

### 4.2 The design that works, if this is built

1. **Run the shared pass as a subprocess of `x_TEA_main.py`**, not in-process in `run_xtea.py`. The
   collector reads runtime globals (`READ_LENGTH`, `DISC_THRESHOLD`, the clip/mapq cutoffs) and
   `extnd = max(400, mean_is + 3·std_var)` from `automatic_gnrt_parameters`. Reproducing that setup in
   the parent invites silent divergence, and a divergent shared temp would be *used*, not detected —
   the file exists, its content is just wrong for the window. A subprocess through the same branch
   gets the identical setup for free.
2. **Write a manifest** recording `extnd`, the BAM, and the site count. Each consumer verifies the
   manifest's `extnd` equals its own and **ignores the shared folder entirely on mismatch**. This is
   the guard that makes (1) safe rather than merely likely-correct.
3. **Consumers merge, never compute-and-delete.** In `collect_clipped_reads_of_given_list`, partition
   the site list into *have* (temp present in the shared folder) and *need*; pool over *need* only,
   into the type's own folder; then merge in the **unchanged original per-chrm/per-site order**,
   reading each site from whichever folder holds it and deleting only the type's own temps. Preserving
   merge order is what makes the output byte-identical — and this pipeline is known to move calls with
   ordering/thread count, so it is not optional.
4. **The parent removes the shared folder after `downstream_parallel`.**
5. Files touched: `x_clip_disc_filter.py` (collector + `collect_clipped_reads*`),
   `x_transduction.py:196` (thread the folder), `x_TEA_main.py` (new shared-collect branch, new
   `--td-shared-folder` option on the transduction branch), `run_xtea.py` (`build_type_steps` argv +
   the new phase before `downstream_parallel`), plus tests.

Precedent for the pattern already exists in the tree: the cns-emit fan unions all types' candidate
windows "so overlapping ones are read once rather than once per type" (`run_xtea.py`, cns_emit block).

### 4.3 BLOCKER — item 1's payoff is unquantified, measure before building

**Transduction has no internal phase breakdown.** `XTEA_TD_PROFILE=1` instruments only
`parse_td_sibling_by_site` in `x_orphan_transduction.py` (~305 calls) — **not** the ~49K-site
collection loop where this saving lives. So of L1's 1481s we do not know the split between:
- per-site clip/disc collection (what item 1 halves) — bounded below by ~390s (49,042 × 0.0948s ÷ 12),
- the bwa realignments against the flank library,
- parsing the 2.4 GB `temp_transduction_clip.sam`,
- disc mate extraction.

Item 1 is a ~6-file change to the most fragile module in the pipeline. **Cutting 2x off 26% of the
step is a different decision from 2x off 80%, and right now nobody knows which it is.** Add
`time.time()` phase timers around the four blocks in `call_candidate_transduction_v3` (cheap, inert)
and run `docs/handoff/runs/transduction_attribution.sbatch` (~25 min, stages pre-computed cns inputs)
before committing to the refactor.

### 4.4 "Run the whole thing once and split by family afterwards?" — partly yes

Asked in-session; worth recording because it is the right instinct and item 1 is its safe subset.

- **Collection: yes.** Family-independent, byte-identical achievable. That is item 1.
- **Realignment + parsing: not free.** Each family aligns to *its own* flank library
  (`hg38_FL_L1_flanks_3k.fa` vs `hg38_FL_SVA_flanks_3k.fa`), each extended by polymorphic full-length
  insertions drawn from that family's own cns candidates via `construct_novel_flanks`, with
  `ileftmost_cns` = 500 for L1 but **5000 for SVA**. Concatenating the libraries and splitting hits by
  which family's flank they landed on is implementable, but bwa reports best hits — so adding SVA
  flanks changes which reads L1 retains. **That is a semantic change to the calls, not a refactor.**
  It may well be an improvement (cross-family ambiguity becomes explicit), but it requires an
  equivalence study and should ride on the overdue golden re-baseline, not precede it.
- **`.slct` re-selection: stays per-family** (type-specific annotation), and is cheap anyway.

**Why this is safe, verified in code:**
- `collect_clipped_disc_reads_by_region` (`x_clip_disc_filter.py:80`) takes only
  `(chrm, pos, extnd, sf_bam, working_folder)`. Nothing type-specific enters the gather.
- `extnd` is `iextnd = max(400, mean_is + 3*std_var)` (`x_TEA_main.py:465,490-491`), derived from the
  BAM's insert-size distribution by `automatic_gnrt_parameters`. **There is no `IS_CALL_SVA` or
  `rep_type` branch — same BAM gives L1 and SVA the same value.** This was the one thing that could
  have killed the idea; it does not.
- The type-specific work is entirely downstream: the realignment target
  (`all_with_polymerphic_flanks.fa`, built per type from that type's flank reference) and the
  per-type working folder.

This composes with, and is independent of, the `fb1d766` cache — the cache removes per-call constant
overhead, this removes ~half the calls.

**Unresolved before building:** the gather writes two temp files per site into a **per-type** working
folder, keyed by `chrm_pos`. A shared pass needs either one shared folder that both types then read,
or records returned through the pool (§5 item 3). That interacts with the existing 108K-temp-file
problem, so design them together rather than separately.

## 5. NEXT SESSION — the transduction module discussion

**The user explicitly asked to discuss this module again before more code lands.** Agenda, in the
order the answers unblock each other:

1. **Is the family-agnostic superset (§3.3) intended?** Sensitivity vs 1.5x cost vs false positives.
   This is a *correctness/design* call, and it changes the value of everything below. Nobody decided
   it; it fell out of the scan port.
2. **The shared collection pass (§4).** 1.89x, no correctness question, `extnd` already verified
   type-independent. Highest value-per-risk item currently known.
3. **Then re-rank the previous handoff's §4.3 list**, which is still substantially right but whose
   item 1 has now been narrowed by §3.2:
   - ~~Filter earlier by family / disc cutoff~~ — **ruled out on correctness grounds (§3.2)**. If any
     filter still qualifies it is not one of these two.
   - Range-partition instead of site-partition (`GenomeIntervalIndex` is already in the tree).
   - Stop round-tripping 108K temp files — design jointly with §4.
   - One `Pool` for the step instead of one per chromosome (~25 fork storms).
   - Fuse the three sweeps (clip, disc, mate gather) into one pass.
   - Drop the 2.4 GB SAM text intermediate in `parse_transduction_clip_algnmt`.
4. **Still overdue, fourth handoff running: re-baseline the goldens.** Untouched. chr1 call counts
   moved relative to `a2_chr1_cnsemit_out` and there is no current baseline to say which moves are
   intended. Pre-worked blast-radius numbers are in the 07-28 handoff.

**If you want the cache measured** (optional now — it is a constant-factor fix whose ceiling is ~387s,
and §4 is worth ~2x more), resubmit:
```
sbatch --export=ALL,XTEA2_SCALE=wgs,\
XTEA2_WORKTREE=/n/scratch/users/c/cos689/xtea2_testdata/perf_wt,\
XTEA2_OUT=/n/scratch/users/c/cos689/xtea2_testdata/perf_wgs_out \
  docs/handoff/runs/perf_validate.sbatch L1,Alu,SVA
```
The type list is a POSITIONAL arg — never `--export` (commas split). Worktree is frozen at `fb1d766`
with `src/csrc/xtea_scan` already built. Cheaper alternative that costs ~25 min instead of a full
pipeline: `docs/handoff/runs/transduction_attribution.sbatch`, which stages pre-computed cns inputs.

The one gap in the cache's evidence, if it is measured and underdelivers: the per-site preamble's
dominance was verified in `parse_td_sibling_by_site` (97% of it, but only ~305 calls), **not** in
`collect_clipped_disc_reads_by_region` (the ~49K-site loop where the saving is claimed). Set
`XTEA_TD_PROFILE=1` (already in the tree, inert by default) to check.

## What NOT to redo

Carried forward from `HANDOFF_2026-07-29_transduction_fixed.md` — all still true:

- **Do not look for an algorithmic port-vs-master difference in transduction.** There isn't one. The
  gap was a pysam-version constant times a per-site loop. §3.3 *reinforces* this: the real difference
  is in the input volume handed to it, upstream.
- **Do not re-derive the per-candidate cost from the ~300-line filtered candidate list.** Wrong
  denominator. It is 53,953 raw / 49,042 selected.
- **Do not re-profile `parse_td_sibling_by_site` for the gap.** ~14-29s total either way.
- **Do not treat chr1 timings as evidence about the cache.** ~100x less to save: 13.1x fewer sites AND
  a 7.8x cheaper BAM open (0.0050s chr1 vs 0.0386s WGS, same 3,366 contigs — open cost tracks
  BAM/index size, not contig count). chr1's theoretical saving is ~4s, under noise.
- **Do not re-run the chr1 A/B for output equivalence.** Done and clean: 3 VCFs + 70 per-type
  intermediates byte-identical with the cache isolated at matched `-n 12`; all 4 VCFs identical with
  both fixes together.
- **Do not treat chr1's `L1 = 0` as a bug.** Master's own chr1 golden is L1 0.
- **Do not run downstream types serially at `-n 24`.** Refuted three times.
- **Do not build the step-boundary shared core pool.** ~15-20s ceiling; fix 2 gets it statically.
- **Do not put comma-bearing values in `sbatch --export`.**

New this session:

- **Do not "fix" the family-agnostic raw-disc list by adding a family or disc-cutoff filter, and do
  not restrict transduction to the few-thousand-site candidate list.** It is pre-family and pre-cutoff
  by construction, deliberately: on WGS L1, 64 of 322 output sites (20%) and **all 59** orphan/sibling
  calls are absent from the 300-site candidate list. Filtering first deletes the answers (§3.2).
- **Do not compare raw-disc files across types by whole line.** The `nL`/`nR` columns are
  family-specific, so line-wise `comm` reports ~66K "differences" and hides the fact that the site
  *positions* are 99.98% identical. **Compare `cut -f1,2` only.** This mistake was made and corrected
  in-session; it inverts the conclusion.
- **Do not treat the previous handoff's "master 1126s / 1379s" as master's real transduction time.**
  That is master's code on the port's inputs (`49075238`). Master's own run did 448s / 587s (§2).
- **Do not filter monitor output on bare `Error` for WGS runs** — `chrM not exist` fires thousands of
  times benignly and will drown the signal (§1).
- **Do not assume `extnd` differs between L1 and SVA.** Checked: `max(400, mean_is + 3*std_var)`,
  no type branch (§4).

## Where things live

| | |
|---|---|
| the fixes | `fb1d766`; `xtea/x_alignments.py:get_cached_bam`/`get_cached_bam_meta`, `xtea/x_clip_disc_filter.py:collect_clipped_disc_reads_by_region`, `xtea/run_xtea.py` (core split) |
| the core-split tests | `tests/test_transduction_core_split.py` (6 tests) |
| the family filter (disc cutoff) | `xtea/disc_adapter.py:162-168` (annotation), `:196-200` (`keep` gate, and `raw_rows.append` outside it) |
| the family-agnostic raw-disc gate | `xtea/disc_adapter.py:218-232` |
| the per-site collector | `xtea/x_clip_disc_filter.py:80`; callers `:336`, `:3309`; transduction caller `xtea/x_transduction.py:196` |
| `extnd` derivation | `xtea/x_TEA_main.py:465`, `:490-491` |
| the pysam microbenchmarks | `docs/handoff/scratchpad/bench_pysam_open.py`, `bench_preamble.py` |
| the cache equivalence + fork-safety proof | `docs/handoff/scratchpad/verify_bam_cache.py` (3 claims, all pass) |
| the validation harness | `docs/handoff/runs/perf_validate.sbatch` (`XTEA2_SCALE=chr1\|wgs`, types as `$1`) |
| chr1 run at the fix | `49099643`, `perf_chr1_out` |
| the chr1 A/B arms (complete, clean) | `49100173` / `49100174` / `49100175`; outputs `ab_parent_3type`, `ab_parent_2type`, `ab_mine_2type` |
| **the cancelled WGS run** | `49101214`, partial output `perf_wgs_out` |
| master whole-pipeline comparator | `48825977`, `master_serial_out`, log `master_serial_wgs_48825977.out` |
| master-vs-port transduction at equal inputs | `49075238`, `transduction_attribution_49075238/` |
| the port baseline these were diffed against | `49070030`, `port_wgs_out` |
| frozen worktrees | `perf_wt` (`fb1d766`), `parent_wt` (`6fb6247`) — both with `src/csrc/xtea_scan` built via `PKG_CONFIG_PATH=/home/cos689/software/lib/pkgconfig make` |
| the root-cause handoff this supersedes | `HANDOFF_2026-07-29_transduction_fixed.md` |
| phase-level localization | `HANDOFF_2026-07-29_transduction_profiling.md` |
| downstream-concurrency analysis, the Alu-skip finding | `HANDOFF_2026-07-29_downstream_concurrency.md` |
| the cns decision, C++ port, re-baselining pre-work | `HANDOFF_2026-07-28_cns_decision.md` |
| the cns/perf findings that still stand | `PERF_2026-07-26_why_xtea2_is_slower.md` |
