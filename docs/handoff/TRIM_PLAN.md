# xTea2 — plan for trimming and optimizing the Python

Written 2026-07-25. Companion to `docs/clipdisc_contract_spec.md` §4 (the bug/dead-code log), which
holds the measured evidence behind every claim here as entries **B6, B7, P1, D2**.

> ## STATUS: steps 1, 2 and 4 are DONE. Step 3 SKIPPED; step 5 is the only work left.
>
> **See [`HANDOFF_2026-07-26_wgs_run.md`](HANDOFF_2026-07-26_wgs_run.md) first.** Step 4 landed
> (`254d0a5`, −2,573) and the full-genome NA12878 run is job `48825976`. `xtea/` is
> **20,273 lines / 42 modules** (−38.5% from 32,963). Step 3 was skipped and re-baselining
> deferred by user decision, for this pass only. **Step 5 (dead locals) is NOT done.**
>
> Correction: this file says master lives at `/home/cos689/src/xtea`. That path does not
> exist. The untouched checkout is `/home/cos689/xtea_exploration/xTea/`.
>
> | step | | lines | commit |
> |---|---|---|---|
> | 1 | pure deletions | **−8,389 actual** (est. −8,679) | `120e207`, `bceff69`, `52bb453` |
> | 2 | interval tree | **−1,739 actual** (est. −1,990) | `4a9e824` (evidence), `65fca02` (code) |
> | 3 | attribute B7 | — | **SKIPPED** (user decision) |
> | 4 | scan-only | **−2,573 actual** (est. −2,700) | `254d0a5` |
> | 5 | dead locals | small | **← the only step left** |
>
> `xtea/` is now **22,845 lines / 43 modules**, down from 32,963 / 58 (−31%).
>
> **Step 2 outcome.** The AVL tree is gone: −1,996 deleted, +147 for `xtea/interval_index.py`,
> +110 net across the 11 edited modules (mostly explanatory comments) = **−1,739 net**.
> Construction dropped ~100× (23.6 s → 0.22 s on Alu) and is now cached, so the ~12× rebuild per
> 3-type run is gone too. The design changed from the checklist: the API is **split** into a
> boolean-only entry point (seven hot sites) and a position+divergence one (three cold sites),
> which makes B8 and B8c non-bugs by construction rather than by convention. B8's min-divergence
> rule and B8b's ordering fix are in.
>
> **The differential passed** over all 2,188,624 interval starts, run against the tree before
> deleting it: zero boolean differences, zero hit-set differences, and every filter-verdict flip
> in the newly-dropped direction (zero newly-kept), confirming the rule is monotone.
>
> **⚠ THE OUTSTANDING OBLIGATION IS RE-BASELINING** — checklist item 8. Calls have deliberately
> moved; B6 and B7 were measured against the old baseline and must be re-measured. See
> "Step 2's carry-forward" at the bottom of this file.
>
> **From step 1 (still current):**
> Two of D2's four counts were overstated and were corrected during execution — imports 26→**23**,
> commented-out code 739→**438**; modules 13 / 7,940 confirmed exactly. See §4 **D2-EXECUTED** for
> the measurements and, more importantly, for the four comment-classifier false-positive classes
> that had to be excluded — one of them is `interval.py`'s half-open-interval documentation, **which
> step 2 depends on**.
>
> Also cleared: `src/csrc/xtea_scan` has been **rebuilt** (it was stale vs `xtea_scan.cpp`), and job
> `48774776` has been **read — it passes**, with one remaining unrecovered Alu target
> (`chr7:64580547`, lost inside cns filtering) now logged in §4.

## The situation in numbers

`xtea/` is **32,963 lines** across 58 modules. Of that:

| | lines | |
|---|---|---|
| new in xTea2 | 1,511 | `run_xtea.py` 953, `disc_adapter` 237, `scan_adapter` 164, `locus_table` 157 |
| edited from upstream | ~320 | 8 files |
| **upstream, untouched** | **~31,000** | everything else |

So xTea2 today is a thin orchestration layer plus a C++ scan engine (2,042 lines) bolted onto an
essentially unmodified xTea. The trimming work is about the ~31,000.

**The measured headline (P1): the vendored pure-Python interval tree is simultaneously the biggest
deletion target and the biggest speedup in the codebase.** It is 96% of every postF invocation.

## What makes this tractable now: the fast downstream gate

`tests/test_downstream_from_cns.py` (added this session) stages the golden's own post-cns outputs
and runs only the tail — transduction → sibling → postF ×2 → gene → genotype → gVCF — against the
golden's own downstream outputs.

- **59 files byte-identical** across L1 (31), Alu (7), SVA (21) — reported, not asserted (see 5c).
- **The gate is VCF columns 1-8.** Intermediate differences print as diagnostics, labelled when
  they are a recognised class (B6 genotype column, B7 SVA transduction intermediates).
- Cost: **~11 min** for all three types; a single type via `-k Alu` is a few minutes.
- Argv comes from `run_xtea.build_type_steps(..., start_from="transduction")` — the same function
  the real pipeline uses — so the harness cannot drift from the real chain.

This is the gate for everything below. It is not fast in absolute terms *yet*, and the reason is
P1: fix the interval tree and this gate drops to ~2-3 minutes, which then makes every subsequent
trim cheap to verify.

> **⚠ That last sentence turned out to be wrong** (step 2, 2026-07-25). The gate went ~11 min →
> **8m25s**, not 2-3 min: the tree was 96% of *postF*, not of the whole downstream chain, so the
> saving was always bounded by postF's share of it. Construction itself did improve ~100×. The
> remaining time is other work — profile it before optimizing, which is how P1 was found.

> **P1-P5 are identifiers, not execution order.** They are numbered by when they were identified so
> that the bug log and this plan agree. The order to actually work in is the sequence table at the
> bottom — pure deletions first (risk-free, and they shrink what you have to read), interval tree
> second (it makes the gate fast), then B7, scan-only, dead locals.

---

## P1 — Replace the interval tree  ·  *execution step 2*

**Deletes ~1,990 lines and is the largest single speedup available.**

`intervaltree.py` (1,142) + `interval.py` (305) + `node.py` (593) = **2,040 lines** of vendored
pure-Python AVL tree. Profile of one chr1 Alu `--postF` (evidence in P1):

```
104.8 s total
 100.9 s  construct_interval_tree                     (96%)
  77.8 s    index_rmsk_annotation_interval_tree       955,824 single inserts
                                                      13.4M node.add, 14.3M rotate
  23.1 s    load_rmsk_annotation_...                  parsing the annotation
   ~4 s   the actual post-filtering
```

Rebuilt in **~12 places per 3-type germline run** (postF ×2/type, transduction ×1/type,
disc_adapter ×1/type, plus gene and blacklist). The cost is **fixed** — the annotation is
genome-wide, so it is the same ~25 s on WGS as on chr1.

**Why a tree is the wrong structure here:** the annotations are static. Built once from the repeat
`.out` files, then queried read-only — 17 `is_within_repeat_region_interval_tree` call sites, no
mutation after construction. A tree buys insert/delete that nothing uses.

**The change, in two parts:**
1. Replace the tree with sorted arrays + `np.searchsorted` (~50 lines). Intervals per chromosome,
   sorted by start, with a running max-end for correct overlap queries.
2. **Build once and cache.** Even a fast build shouldn't run 12×. Key the cache on
   (annotation path, boundary_extnd, min_copy_len) — those are the only inputs.

**Risk: medium, and concentrated in semantics, not performance.** The replacement must reproduce
`is_within_repeat_region_interval_tree` exactly: boundary extension (`boundary_extnd`, and SVA's
200 vs Alu/L1's 0), open/closed interval convention, and the tie behavior at exact boundaries. Get
this wrong and calls move. Mitigation: write a differential test that runs both implementations
over the real annotation and asserts identical answers for a large sample of query positions,
*before* deleting the old one.

> ### ⚠ UPDATED 2026-07-25 evening — step 2 as written above will NOT work. See §4 **B8**.
>
> "Sorted arrays + `np.searchsorted`" returns the **leftmost** overlapping interval. The original
> returns **an arbitrary** one — `for rcd in set_rslt: return True, rcd[0]`, the first element of a
> *set* iteration. Measured on the real hg38 annotations at the real parameters, those two answers
> differ on **9.2% of Alu queries and 24.3% of L1 queries**. The "tie behavior at exact boundaries"
> note above understated this by an order of magnitude: boundary extension (Alu/SVA ±100, L1 ±400)
> makes neighbouring copies overlap constantly, so 17.7% (Alu) / 45.2% (L1) of queries are
> multi-overlap in the first place.
>
> **And it is results-determining.** The returned position feeds `get_div_subfamily`, whose
> divergence rate gates a call-dropping filter (`x_post_filter.py:90,97,153,160`, cutoff 15%).
> Overlapping copies land on opposite sides of 15% often enough that the arbitrary pick flips the
> verdict at **2.95% of Alu and 5.79% of L1 positions**. The *boolean* is unaffected, so the call
> sites that discard the position are safe — only the position-consuming ones are exposed.
>
> **Revised design (still gets the speedup, and is bit-exact):** the 78 s is interval *construction*,
> not querying (~4 s of 104.8 s). So keep sorted arrays + `bisect` to find candidates fast, then when
> there is more than one candidate reproduce the original choice by putting the survivors in a real
> `set()` and taking `next(iter(...))` — cheap at ≤5 elements, and bit-exact because
> `Interval.__hash__` is `hash((begin, end))` (`interval.py:128`), a plain int tuple, stable across
> runs and independent of `PYTHONHASHSEED`. Residual risk to test: set insertion order can matter on
> a hash collision, and the original inserts in tree-traversal order vs. the replacement's sorted
> order. Run the differential test over all ~1.45M positions, not a sample.
>
> **DECISION (user, 2026-07-25 evening): FIX the bug.** Step 2 is therefore no longer a pure
> optimization — it is an optimization **plus a deliberate behavior change**, and it carries a
> re-baselining obligation. The rule is the **minimum-divergence** overlapping copy, tie-broken by
> leftmost start; see §4 B8 for why that is the only rule the filter's purpose justifies, and note
> that it is strictly more conservative, so xTea2 will drop somewhat more calls than master here.
>
> **Revised step-2 checklist:**
> 1. Sorted arrays + `bisect` (prefix-max-end pruning) replacing the AVL tree — the actual speedup.
> 2. Cache the built annotation on (path, extnd, min_copy_len) — kills the ~12× rebuild.
> 3. Apply the min-divergence rule in `is_within_repeat_region_interval_tree`.
> 4. Sort the hits in `is_within_repeat_region_interval_tree2` and in **B8b**
>    `x_gene_annotation.query_by_position` (arbitrary `GENE_INFO` ordering).
> 5. Neutralize **B8c** `x_black_list.fall_in_region` so the trap cannot be re-armed.
> 6. Differential test over all ~1.45M positions: the **boolean must be identical everywhere**
>    (it is order-independent, so any boolean difference is a bug in the new code, not the fix).
>    Positions are *expected* to differ — characterize and record where, do not assert equality.
> 7. Delete `intervaltree.py` + `interval.py` + `node.py` (2,040 lines) once nothing imports them.
>    **Note:** `x_black_list.py` and `x_gene_annotation.py` also `from intervaltree import
>    IntervalTree` — they must be migrated too, or the deletion breaks them.
> 8. Re-baseline the goldens; re-measure B6 and B7 against the new baseline.

> ### ⚠ REVISED AGAIN 2026-07-25 late night, after measuring. See §4 **P1-EVALUATED** and **B9**.
>
> The direction above is right; three things about it are not.
>
> - **The cost figures were cProfile-inflated ~3×.** Real: AVL build 23.6 s (Alu) / 18.0 s (L1),
>   parse 2.75 s / 2.13 s — not 77.8 s / 23.1 s. Sorted arrays build in **0.22 s (105×)** and query
>   **12× faster**. Run-level: ~270 s of construction per 3-type run → ~33 s with sorted arrays →
>   ~14 s with the cache. **The structure swap dominates; the cache is a ~19 s second**, which is
>   the reverse of how items 1 and 2 above are weighted.
> - **Item 7's line count is stale**: the three modules are **1,996 lines**, not 2,040 — step 1's
>   comment stripping already took the difference.
> - **Items 3 and 5 are better solved structurally.** Rather than one `(bool, position)` function
>   with a selection rule bolted on, **split the API**: a boolean-only entry point for the seven
>   hot sites that bind the position and never read it, and a position+divergence entry point for
>   the three that do. B8 and B8c then cannot be re-armed *by construction* — there is no position
>   for a later edit to misuse. **Approved by the user, 2026-07-25 late night.**
>
> Two further facts that change surrounding work:
>
> - **B8's L1 row was measured at the wrong parameter.** postF's L1 branch uses
>   `L1_min_ref_copy_len=125`, not `REP_LINE_MIN_LEN=300`. At the real value the exposure is
>   **55.9% multi-overlap and 7.00% verdict flips**, not 45.2% / 5.79%. Item 8's re-baselining will
>   move more than recorded.
> - **B9: the master pysam clip/disc fallback is already py3-broken** (`x_annotation.py:582`,
>   float list index, `TypeError`). P5b's "accepted cost" of keeping it as the attribution control
>   is illusory — it cannot run. This weakens the case for sequencing P4 before P5b.

**Expected after:** postF ~29 s → ~1 s; downstream gate ~11 min → ~2-3 min; the full chr1 gate
should drop substantially from 18 min.

---

## P2 — Free deletions  ·  *execution step 1 (with 5a)*

**~2,800 lines here, plus 5,879 from 5a = 8,679 total in this step.** No behavior risk, no decision
outstanding. Do it first: it makes everything after easier to read.

| what | lines | notes |
|---|---|---|
| orphaned modules | 2,061 | `x_vcf_merger` 1006, `x_cmp` 520, `x_insert_size` 262, `x_sv` 100, `x_population_AF` 72, `union_find_set` 65, `x_ref_TE_deletion` 18, `x_somatic_panel_normal` 18 |
| commented-out code | 739 | verified by regex as code, not prose |
| unused imports | 26 | |

All eight modules are unreachable from `run_xtea` + `x_TEA_main`, verified by AST import-graph
analysis **and** by hand — `x_vcf_merger`, `x_insert_size` and `union_find_set` appear only in
comments and docstrings. `x_vcf_merger.py` is what `locus_table.py` replaced (see B5).

Gate: import sweep over all remaining modules + `run_xtea.py --dry-run` + the fast suite. No
compute needed.

---

## P3 — Dead local variables in the downstream  ·  *execution step 5 (last)*

**184 dead locals**, clustered exactly where the downstream lives:

```
x_orphan_transduction.py   55      x_transduction.py     10
x_post_filter.py           24      x_somatic_calling.py  10
x_clip_disc_filter.py      19      (56 more across 20 modules)
```

Gate: the downstream harness (call-set), and **read its intermediate-diff diagnostics** — dead-local
removal is exactly the kind of edit that should move nothing at all, so any `NOTE — intermediate
moved` line here means the removal was not as dead as it looked.

**The one real hazard:** a local whose right-hand side is a *call with side effects* is not dead —
deleting `x = f()` where `f` writes a file changes behavior. Check every RHS; never bulk-delete on
the name alone. This is why it is P3 and not P2.

---

## P4 — Attribute B7 before trimming transduction  ·  *execution step 3*

The SVA transduction step already differs from the golden at the intermediate level (extra site
`chr1:182305180`; HC gains it and loses `chr1:50381610` + `chr1:182510697`). It is deterministic and
washes out — the final SVA VCF is byte-identical, 5 calls.

**Do this before touching transduction**, otherwise a trim there gets blamed or credited for a
pre-existing difference. Likely causes: the 84 changed lines in `x_clip_disc_filter.py`, or
py3.8→3.10. Cheapest attribution: run the same staged SVA chain against master's untouched
`/home/cos689/src/xtea/xtea/x_TEA_main.py` under the same interpreter — that isolates our edits from
the interpreter in one run.

---

## P5 — Decisions: ALL MADE (user, 2026-07-25)  ·  *5a in step 1, 5b in step 4*

**5a. Drop the legacy entry points — YES, delete all 5,879 lines.**
`gnrt_pipeline_local_chm13.py` 1349, `gnrt_pipeline_local_v38.py` 1347, `gnrt_pipeline_local.py`
1251, `gnrt_pipeline_germline_cloud_v38.py` 828, `x_TEA_main_germline_cloud.py` 1104. Upstream's
"generate a shell script, then submit it" workflow — precisely what `run_xtea.py` replaced, and
unreachable from xTea2's entry points. xTea2 no longer supports the upstream invocation style or
the cloud path. Recoverable from git history.

**5b. Commit to scan-only — YES (~2,700 lines).**
`--scan-clip`/`--scan-disc` become always-on; the master pysam clip/disc fallback goes.
`x_TEI_locator.py` (1,356) and much of `clip_read.py` (1,316) become dead, along with the `-C`/`-D`
branches of `x_TEA_main.py`. **Accepted cost:** the fallback was the controlled variable in
attribution runs (`a2_validate_chr1_cnsemit.sh` varies the scan binary while holding everything
else fixed). Do this **after** P4, because B7's attribution wants to compare against master's
untouched path — once it is deleted, that comparison needs a git checkout instead.

**5c. Gate scope — call-set concordance only (VCF columns 1-8).**
Intermediate files may churn during the trim without failing the build. Implemented in
`tests/test_downstream_from_cns.py::test_downstream_call_set_matches_golden`.

**Deliberate hedge, worth keeping:** the byte-level comparison of every intermediate still *runs*
and prints what moved — as diagnostics, not a verdict. B7 was found only because intermediates were
compared, and a trim that silently changes an intermediate is exactly what should be visible in the
log before it becomes a call-set bug. Cost of keeping it: zero (the files are already on disk).
**If an intermediate starts moving, read it — do not assume a green gate means nothing changed.**

**5d. Package layout — LATER.** Keep `xtea/` for now; move to `src/xtea/` + pyproject after the
trim, on a smaller codebase. Avoids churning the paths in the scripts and docs just converted.

---

## Sequence and expected end state

Ordering changed now that 5a/5b are approved: the pure deletions are risk-free and no longer need a
decision, so they go first and make everything after easier to read. P1 still precedes the
finer-grained work because it is what makes the gate fast.

```
step                                   lines      gate
-----------------------------------------------------------------------------------
1. P2 + P5a  pure deletions           -8,389     import sweep + --dry-run + fast suite
   [DONE]    (8 orphaned modules, 4 gnrt generators, cloud main,     ** ACTUAL **
              438 commented-out lines, 23 unused imports)
2. P1        interval tree            -1,739     differential (2,188,624 positions) +
   [DONE]    (-1,996 tree, +147 index, +110 docs) import sweep + --dry-run + fast suite
                                       ** ACTUAL **
3. P4        attribute B7                 --     no deletion; unblocks transduction
   [NEXT]
4. P5b       scan-only                -2,700     downstream gate + full chr1 gate once
5. P3        dead locals              small      downstream gate (check every RHS)
-----------------------------------------------------------------------------------
    32,963 -> 24,584 -> 22,845 -> ~20,100         (steps 1-2 banked; ~39% smaller at the end)
```

## Step 2's carry-forward — READ BEFORE STEP 3

**0. What the downstream gate actually said after step 2 — it PASSED, in 8m25s (was ~11 min).**
No final call moved on the chr1 golden, despite the B8 fix changing the selected copy on 8.9%
(Alu) / 30.6% (L1) of *annotation positions* and flipping the filter verdict on 2.88% / 7.00%.
The gap is the point: those percentages are over annotation positions genome-wide, not over the
golden's ~30 actual candidate sites, and on this sample none of them landed on a straddling copy.
**Do not read that as "the fix is a no-op"** — it is one sample on one chromosome, and the
re-baselining below is still owed. Four L1 intermediates did move, all in transduction:

```
candidate_disc_filtered_cns2.txt.all_non_sibling_td.txt.unique_trsdct_half_clip
candidate_list_from_disc.txt.clip_sites_raw_disc.txt.slct
candidate_list_from_disc.txt.clip_sites_raw_disc.txt.slct2
candidate_list_from_disc.txt.clip_sites_raw_disc.txt.slct2.tmp
```

That these surfaced at all is the 5c hedge working as designed — the call-set gate could not see
them, and the intermediate diagnostics could. **They land in exactly the code path step 3 is about
to attribute (B7, SVA transduction), so read them before starting step 3**, or a step-2 effect will
be mistaken for a B7 finding.

**1. Re-baselining is owed and not yet done (checklist item 8).** The B8 fix deliberately moves
calls, always in the more-conservative direction. Until the goldens are re-baselined:

- `tests/test_downstream_from_cns.py` compares against goldens generated *before* the fix, so
  any call-set difference it reports is expected, not a regression. Read the diff, do not just
  read the exit code.
- **B6 and B7 are measured against the old baseline** and must be re-measured. B7 especially —
  step 3 is about attributing it, and attributing it against a stale baseline wastes the run.

**2. The gate did not get as fast as predicted, and that is fine.** P1 predicted "~11 min → ~2-3
min". The interval tree was 96% of *postF*, not 96% of the downstream chain, so the saving is
bounded by postF's share of the whole. Do not treat the remaining time as a regression; profile
before optimizing it further, exactly as P1 was found by profiling.

**3. B9 changes step 3's cheapest path.** The master pysam clip/disc fallback cannot run under
py3.10 (`x_annotation.py:582`, float list index). P4's "cheapest attribution" — run the staged
SVA chain against master's untouched `x_TEA_main.py` under the same interpreter — is still valid
*provided* the chain does not enter the clip/disc fallback. If it does, expect a `TypeError` and
budget for fixing `is_within_repeat_region` in the master checkout first.

**Carried into step 2 from step 1's execution:** `interval.py:46,48` document the half-open
convention (`begin <= C < end`, `self.begin <= C < self.end`) and were deliberately preserved when
the commented-out code was stripped. They are the written form of the semantics the
`np.searchsorted` replacement has to match, so write the differential test against them.

**Do P4 before P5b.** Attributing B7 wants a run against master's untouched clip/disc path; deleting
it first makes that a git-checkout exercise instead of a flag flip.

## What is NOT in this plan, deliberately

- **Rewriting the big downstream modules** (`x_clip_disc_filter` 3,442, `x_transduction` 1,562,
  `x_orphan_transduction` 1,421, `x_post_filter` 1,247). Trim them first, profile them second, and
  only then decide whether any deserves a rewrite. Nobody has profiled them yet; P1 was found by
  profiling and is worth more than a speculative rewrite of all four.
- **Porting more to C++.** The measured bottleneck in the downstream is a data-structure choice, not
  a language boundary. Fix P1 and re-profile before assuming C++ is the answer.

## Loose ends carried in — both now CLOSED

- ~~`src/csrc/xtea_scan` is stale relative to `xtea_scan.cpp`~~ — **rebuilt 2026-07-25 20:24**
  (`PKG_CONFIG_PATH=/home/cos689/software/lib/pkgconfig make` in `src/csrc`, exit 0, warnings only).
- ~~Job `48774776` result not yet read~~ — **read, and it passes.** Assertion 1b clears (raw-disc
  rows differ from pre-fix, so `b_tumor` reached `disc_adapter`), Alu `chr12:121357375` reaches the
  final VCF, germline Alu recall 96.6% → 99.4%. Full result and the two remaining gaps in §4.
- ~~The working tree has uncommitted changes~~ — committed as `4185655`; step 1 is `120e207`,
  `bceff69`, `52bb453`. Tree clean.
