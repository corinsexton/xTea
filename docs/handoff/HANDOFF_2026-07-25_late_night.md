# xTea2 handoff — 2026-07-25 (late night session)

**This is the current entry point.** The two earlier 2026-07-25 handoffs (afternoon, evening) are
history; where they describe state, this file supersedes them.

One-line status: **Trim step 2 is DONE — the vendored AVL interval tree is deleted, annotation
construction is ~100× cheaper and cached, and the B8 bug family is fixed. The tree is clean, all
gates pass, and step 3 (attribute B7) is next. One obligation is outstanding and deliberate:
re-baselining the goldens.**

Read [`TRIM_PLAN.md`](TRIM_PLAN.md) for the work — its STATUS block and its
**"Step 2's carry-forward"** section are authoritative. This file is what happened and what state
things are in.

## Start here

```bash
conda activate xtea2        # the entire env setup; nothing else to export
cd /n/data1/hms/dbmi/park/corinne/claude_sandbox/xTea

XTEA2_SCAN_BIN=$PWD/src/csrc/xtea_scan pytest -m "not slow"   # 83 passed, 1 skipped
pytest tests/test_interval_index.py -q                        # 14 tests, <1s, no external data
pytest tests/test_downstream_from_cns.py -q -s                # the call-set gate, 8m25s
```

The one skip (`test_e2e_concordance.py`, `XTEA2_OUTPUT_DIR not set`) is pre-existing and expected.

## 1. What happened this session

The session opened on "step 2: replace the interval tree" as written in `TRIM_PLAN.md`. Before
implementing, the plan's premise was re-examined on request — and that changed the design, so the
work split into three commits, evidence deliberately committed *before* the code that relies on it:

| commit | what |
|---|---|
| `4a9e824` | measurements, B9, corrected B8 numbers, two reproducible benchmark scripts |
| `65fca02` | the code change |
| `d25d41a` | outcome: differential results, gate result, step 3's inheritance |

**`xtea/` is now 22,845 lines / 43 modules** (was 24,584 / 45; 32,963 / 58 before the trim — **−31%**).
Net −1,739 = **−1,996** (the three vendored modules: `intervaltree.py` 1,137 + `interval.py` 305
+ `node.py` 554) **+147** (`xtea/interval_index.py`, stdlib-only) **+110** net across the 11 edited
modules, most of it the comments explaining the change.

## 2. Why the design differs from the checklist

**The tree was never needed for its data-structure properties.** `intervaltree.py` was a vendored
copy of the PyPI package; xTea used exactly one operation — `tree[point]` — after a build-once
phase, never inserting or deleting afterwards. Both files that adopted it say why:
`x_annotation.py:9-10` ("Use interval tree to replace binary search") and `x_gene_annotation.py:9`
("replace the query module from binary tree (with bug) to interval tree"). The bug is real and the
search is still in the tree: binary search on start positions returns at the first midpoint hit and
can walk past a genuine hit when copies overlap. So the tree was a **correctness** fix — and sorted
arrays + a prefix-max-end scan is correct under overlap *and* ~105× cheaper to build.

Three corrections to the plan, all now recorded in `clipdisc_contract_spec.md` §4 **P1-EVALUATED**:

- **The cost figures were cProfile-inflated ~3×.** AVL build is 23.6 s (Alu) / 18.0 s (L1) and
  parse 2.75 s / 2.13 s — not 77.8 s / 23.1 s. The 96%-of-postF ratio survives; the magnitudes don't.
- **The API is split**, rather than one `(bool, position)` function with a selection rule bolted on.
  Seven call sites bound the position and never read it; three consume it. The boolean-only entry
  point makes **B8 and B8c non-bugs by construction** — no position exists for a later edit to
  misuse. `x_black_list.fall_in_region` returns only a boolean for the same reason.
- **pandas parsing is slower** than the Python loop (3.5 s vs 2.75 s). Recorded with four other
  rejected options (numpy batch, NCLS/cgranges, on-disk index, bitset) so they aren't re-proposed.

## 3. The two new findings

**B8's L1 numbers were measured at the wrong parameter.** postF's L1 branch uses
`L1_min_ref_copy_len=125` (`x_post_filter.py:415`), not `REP_LINE_MIN_LEN=300`. Real exposure is
**55.9% multi-overlap and 7.00% verdict flips**, not 45.2% / 5.79%. The re-measuring harness
reproduces the original row exactly at 300, which is what validates it.

**B9 — the master pysam clip/disc fallback cannot run under py3.10 at all.**
`x_annotation.py:582` does `mid = (lo + hi) / 2`, a float list index; every call raises
`TypeError` (confirmed by execution). All three call sites are in the fallback that step 4 deletes,
so xTea2's default C++-scanner flow never trips it. **But P5b's "accepted cost" — keeping the
fallback as the controlled variable for attribution runs — is illusory, because it cannot run.**
That weakens the argument for sequencing P4 before P5b. Left unfixed on purpose: step 4 deletes it,
and fixing it would mean fixing a wrong algorithm.

## 4. Verification — what was actually checked

**Differential against the AVL tree, run before deleting it**, over every interval start in the
real hg38 annotations at the real parameters (Alu 100/200, L1 400/125, SVA params 200/300) —
**2,188,624 positions**:

| | boolean diffs | hit-set diffs | position changed | verdict flipped | **newly kept** |
|---|---|---|---|---|---|
| Alu | **0** | **0** | 86,997 (8.9%) | 28,238 (2.88%) | **0** |
| L1 | **0** | **0** | 235,048 (30.6%) | 53,856 (7.00%) | **0** |
| SVA params | **0** | **0** | 28,988 (6.6%) | 7,381 (1.67%) | **0** |

Boolean and hit-set identity is the correctness result. **Zero newly-kept** is the one worth
keeping: every verdict flip is in the newly-dropped direction, which *measures* rather than assumes
that min-divergence is monotone.

Other gates: import sweep 43/43 · `--dry-run` 24 steps exit 0, plan unchanged · fast suite 83
passed / 1 skip · downstream gate **passes in 8m25s** (was ~11 min) · cache verified (warm hit
free, differing parameters correctly rebuild).

## 5. Two results that invite the wrong conclusion

**The downstream gate passing is weaker evidence than it looks.** No final call moved on the chr1
golden — but the 8.9% / 30.6% rates above are over *genome-wide annotation positions*, while the
golden has ~30 candidate sites, none of which happened to sit on a straddling copy. One sample, one
chromosome. It is not evidence the fix is inert.

**Four L1 transduction intermediates DID move:**

```
candidate_disc_filtered_cns2.txt.all_non_sibling_td.txt.unique_trsdct_half_clip
candidate_list_from_disc.txt.clip_sites_raw_disc.txt.slct
candidate_list_from_disc.txt.clip_sites_raw_disc.txt.slct2
candidate_list_from_disc.txt.clip_sites_raw_disc.txt.slct2.tmp
```

They were visible only because decision 5c kept intermediate comparison running as *diagnostics*
rather than verdict — the call-set gate could not see them. **They sit in the code path step 3
attributes (B7, transduction), so read them before starting step 3**, or a step-2 effect will be
mistaken for a B7 finding.

**The gate did not reach the predicted "~2-3 min"** because the tree was 96% of *postF*, not of the
whole downstream chain, so the saving was bounded by postF's share. Construction itself did improve
~100×. Noted in the gate's own docstring so the next person profiles rather than assumes.

## 6. Outstanding — one item, deliberately not done

**Re-baseline the goldens (step-2 checklist item 8), and re-measure B6 and B7 against the new
baseline.** Not done for two reasons: it is a separate compute job, and freezing a new baseline
immediately before step 3 touches transduction seemed like the wrong moment — that is a judgment
call to confirm, not a blocker discovered late.

Until it happens: **the goldens predate the B8 fix**, so any call-set difference the downstream gate
reports is expected rather than a regression. Read the diff, not just the exit code. Flagged in
`TRIM_PLAN.md` STATUS, its "Step 2's carry-forward" section, and
`tests/test_downstream_from_cns.py`'s docstring.

## 7. Next: step 3 — attribute B7

No deletion; it unblocks trimming transduction. What it inherits from step 2, in order of how
easily it could waste a run:

1. **The four moved L1 intermediates above** — read them first.
2. **B9** — P4's "cheapest attribution" (run the staged SVA chain against master's untouched
   `x_TEA_main.py` under the same interpreter) is still valid *provided* the chain does not enter
   the clip/disc fallback. If it does, expect a `TypeError` and budget for fixing
   `is_within_repeat_region` in the master checkout first.
3. **B7 was measured against the old baseline** and may have moved.

After step 3: step 4 (scan-only, −2,700) then step 5 (dead locals). Expected end state ~20,100
lines, ~39% smaller than the 32,963 the trim started from.

## 8. Where things live

| | |
|---|---|
| the plan / next actions | `docs/handoff/TRIM_PLAN.md` (STATUS + "Step 2's carry-forward") |
| bug & evidence log | `docs/clipdisc_contract_spec.md` §4 — B8, B8b, B8c, **B9**, **P1-EVALUATED**, **P1-EXECUTED** |
| the new index | `xtea/interval_index.py` (module header explains why it is not a tree) |
| its tests | `tests/test_interval_index.py` (14, brute-force oracle, no external data) |
| the benchmarks | `docs/handoff/bench_interval_index.py`, `bench_b8_impact.py` — **historical**, they need the deleted tree: `git checkout 4a9e824 -- xtea/intervaltree.py xtea/interval.py xtea/node.py` |
| env | `docs/environments.md` — `conda activate xtea2`, one lever |
