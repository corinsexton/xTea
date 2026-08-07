# xTea2 handoff — 2026-08-07, the objective is met and measured, and the agenda is rebased

**This is the current entry point for state.** The agenda entry point is still
`docs/handoff/DOWNSTREAM_OVERHAUL_PLAN.md`; the gate's design entry point is still
`HANDOFF_2026-08-04_stage3_td_gate_spec.md`.

> Supersedes `HANDOFF_2026-08-06_two_decisions_built.md`. **Its §6 NEXT item 1 is DONE.** Items 2–6
> are PARKED, with their evidence, in §6 below.
>
> **§10's four decisions (2026-08-04) remain the standing agreement. Nothing here reopens one.**

**The one-sentence version: the project was started to make the downstream faster and leaner, that
was achieved and never written down, and the last five days went somewhere else.**

---

## 0. State at handoff

| | |
|---|---|
| branch | `xTea2` |
| commits this session | `eca38d2` the concat default flip; `308e1a0` check 2b out of the heredoc; `a091b8a` both gates test the default; `f551432` the objective probe |
| **`xtea/` DID change** | `run_xtea.py` (`resolve_decoy_concat`, the banner, two help texts), `global_values.py` and `td_decoy_library.py` (comments only). **The default of `--td-decoy-concat` is now ON** |
| `benchmark/` changed | **new** `adjudicate_dropped.py` — check 2b, lifted out of an sbatch heredoc |
| tests | **324 passed / 30 skipped**, was 313/30. 11 new in `tests/test_adjudicate_dropped.py` |
| jobs | `50038654` chr1 4-arm concat gate; `50038678` WGS concat gate on the default path |
| mutations | 5 against the adjudicator, 1 against the flipped resolver. **All 6 bite** |
| what is NOT done | NEXT items 2–6, all parked in §6; Stage 2c's per-locus loop; Stage 4 |

## 1. THE HEADLINE, and it needed no compute

Plan §1's intended outcome is a downstream that is faster and leaner and is judged on
sensitivity/precision against orthogonal truth sets. **Nobody had ever checked.** Every number
needed to answer it had been on disk since 2026-08-05, in `sacct` and in a `step_timings.tsv`
column nothing has ever read.

Same HG002 WGS BAM, same `-c 24`, same `--mem=64G`, same partition, whole pipeline on both sides:

| | wall | peak RSS |
|---|---|---|
| master (`46520ab`), job `49763841` | 1:10:15 | **38.03 GB** |
| xTea2 concat-OFF, job `49792012` | **38:32** | **19.50 GB** |
| xTea2 concat-ON — the new default, job `49879677` | 42:15 | 19.91 GB |

**1.82× faster on 51% of the peak memory**, and sensitivity never moved: L1 188 / Alu 1223 / SVA 82
byte-identical Stage 0 → Stage 2, SVA TP 77 / FN 13 unchanged across both decoy flips.

`docs/handoff/runs/probe_objective_metrics.py` reproduces this offline, in seconds, and prints
every caveat next to the number rather than in a footnote. Three of them matter:

- **master loops L1/Alu/SVA SERIALLY**, one family resident at a time; xTea2 runs all three
  concurrently. So the memory result is *stronger* than the ratio: **xTea2 holds three families in
  half the memory master needs for one.**
- **Conversely, the wall-clock ratio is not pure per-family speed.** Part of it is a parallelism
  master does not have. That *is* an xTea2 change and not a benchmark artifact, but it is not the
  C++ gather, which is separately measured at 15.25× (L1) / 17.82× (SVA) on the gather itself.
- `sacct` MaxRSS is polled and can miss a spike; neither job was near its 64 GB ceiling.

### 1.1 The per-step memory attribution came from an unread column

`record_timing` (`run_xtea.py:271`) has been writing `children_peak_rss_gb` into every run's
`step_timings.tsv` since the instrument was added. On the WGS OFF arm **exactly two rows of the
entire run raise the high-water**: `scan_pass` to 11.64 GB, then `SVA:transduction` to 18.18 GB.
That is an attribution `sacct` cannot give, and it cost nothing.

It is `getrusage(RUSAGE_CHILDREN).ru_maxrss` — the high-water of the **largest single child**,
monotone. Never a per-step delta, never a sum over concurrent children, so it is a different
quantity from `sacct`'s job-level number and is not expected to match it. That 18.18 and 19.50 land
close says one child dominates the job's footprint.

## 2. THE AGENDA IS REBASED, and the reason is a measurement

**Since 2026-08-03 essentially all work has been transduction correctness** — the decoy defect, the
Stage 3 gate, the verifier, the polyA decoupling, the re-freeze, the 2026-08-07 rule amendment. It
found a real bug (SVA ran with **zero** decoy competition) and it is good work. But it is
constraint-defense, and **all six items on the previous handoff's NEXT list are more of it. None
makes anything faster or lighter.**

**And there is no large speed item left to pivot to.** The obvious candidate — the ~1.9×
cross-family sharing that Stage 2c's single per-locus loop would cash in — **does not exist as a
speed win, and this was already measured.** `run_xtea.py:937`'s docstring records jobs `49550925` /
`49548949`: the hoist is **−13.7 s**, the barrier it needs costs **+52.9 s**, and
`downstream_parallel` went 462.7 → **505.3 s**, i.e. slower. `StagedBarrier` takes the +52.9 back,
so at best it is neutral. The docstring says it outright: *"The reason to do the hoist is not
speed."*

> **So: the speed story is COMPLETE, and the per-locus loop must not be sequenced as a performance
> item by a future session.** It remains worth doing on structural grounds — it is the front edge of
> Stage 3 and its dependencies are built (`site_universe.py`, `--locus-universe`,
> `transduction_evidence.py`) — but the case for it is architecture, not throughput.

## 3. `--td-decoy-concat` IS ON BY DEFAULT (`eca38d2`)

hg38's SVA flank library ships **0** of the 45 decoy records against both L1 libraries' 45, so SVA
has had no decoy channel on either the clip or the disc side. Appending is BY NAME, so the flip is a
structural no-op on L1 and the fix on SVA.

The case, as it stood: the recipe is verified as **bytes**; L1 and Alu unmoved at both scales; SVA
`decoy_cross_family` 0 → 44; exactly one call moves genome-wide and the truth set counts it a
**false positive** (SVA FP 8 → 7, TP 77 / FN 13 unchanged, position precision 90.6% → **91.7%**);
Stage 3 PASSES with identical 54.0%–76.2% ranges in both arms; cost ~55 s, ~12% of the transduction
phase.

**The subprocess seam needed no change**, and that is worth recording so nobody re-derives it:
`run_xtea` writes `ENV_TD_DECOY_CONCAT` explicitly as `"1"`/`"0"` on every path and is its single
writer. The env var's **own** rule stays "absent means OFF", exactly as `ENV_TD_DECOY_STRUCTURAL`'s
does, so a direct `x_TEA_main` invocation keeps pre-flip behaviour — which is what the frozen gate
harnesses depend on.

### 3.1 The flip surfaced a live defect its sibling left behind

**22 driver-invoking harnesses under `docs/handoff/runs/` pass neither decoy flag** (a 23rd was
`baseline0_hg002.sbatch`, fixed below; 3 more name `--td-decoy-structural` only). They inherit
whatever the defaults are, and when `--td-decoy-structural` flipped on 2026-08-05 nobody
retrofitted them. The consequence is not hypothetical:

> **`baseline0_hg002.sbatch` — whose entire contract is that re-running it reproduces the frozen
> BASELINE-0 VCFs — has since 2026-08-05 been silently producing a DIFFERENT baseline under the same
> name, exiting 0 while doing it.**

It now pins both `--no-td-decoy-*` explicitly, probed via `--help` rather than hardcoded so it still
runs against pre-flip worktrees, reusing the detect-rather-than-hardcode idiom already there for
`--mate-in-rep-families`. **The other 22 are NOT retrofitted** — most are historical gates whose
baselines are frozen and which nobody re-runs, and touching them would change arms in gates that
would then need re-verifying. Anyone re-running one of those for a *fresh* comparison must decide
which decoy configuration it means and say so on the command line.

### 3.2 Both gates now test the DEFAULT (`a091b8a`)

The flip's original evidence was three chr1 arms and one WGS arm, and **every one of them names a
flag**. All four would pass unchanged if the new default never reached the step — a resolver
returning the wrong branch, an option parsed and dropped, a banner printed from one variable while
the env is written from another. This is the argument `gate_polya_cutoff_chr1.sbatch` already had to
make for its own arm 4.

- **chr1 gains `AFTER_DEFAULT`**, an arm that passes nothing, plus check **1b** (`AFTER_DEFAULT ==
  AFTER_ON` byte-for-byte over the whole output tree) and check **1c** (`AFTER_DEFAULT !=
  AFTER_OFF`, without which 1b would also pass for a wholly inert flag).
- **WGS stops passing the flag entirely.** Its check 0 now asserts two things: that the
  concatenation is on, and that the **DEFAULT** turned it on — the driver's banner names its own
  cause, printing `(DEFAULT)` or `(--td-decoy-concat)`.
- check 1's tree comparison is now `compare_trees`, one definition called three times, rather than a
  fourth private copy of a shipping rule (plan §10).

<!-- GATE RESULTS -->

## 4. Check 2b left the heredoc (`308e1a0`)

Handoff §5.2 ends by telling its reader to **verify, rather than take**, that check 2b would still
FAIL a drop of a VERIFIED entry — and §5.3 concedes the verification was impossible: the
adjudication was a `python - <<EOF` block inside `td_gate_no_regression.sbatch`, the one surface in
this project with no tests and no mutation coverage, where both of 2026-08-06's harness defects
lived while the unit suite stayed green.

`benchmark/adjudicate_dropped.py` is that logic, importable. Behaviour is unchanged: run against job
`49995743`'s own `lost_keys.txt` and `verify_off/verdicts.tsv` it reproduces that job's check 2b
output **line for line**. `tests/test_adjudicate_dropped.py`, 11 tests:

| branch | expected |
|---|---|
| a dropped entry the OFF arm **VERIFIED** | **FAIL** — the hazard the old rule protected |
| UNDECIDABLE / REFUTED / SUB_THRESHOLD | PASS, count still printed — the amendment itself |
| absent from the verdicts entirely | **FAIL** — silence is not a pass |

Plus §4.1's lesson as a test: the frozen 16-column layout and the current 18-column one must give
the **same** answer, because `verdict` moved from field 13 to 15 and an index reader would compare
`n_clip_support` against it.

**Two things changed rather than moved**, both about not confusing a broken harness with a failed
gate. A verdicts file missing a required column, or empty, now exits **2** — distinct from the
gate's own 1 — instead of raising a `KeyError`; note that `SystemExit("msg")` exits **1**, which is
exactly the collision being avoided. The sbatch reads exit 2 as *"HARNESS error, do not read this as
either a PASS or a FAIL"*.

## 5. Numbers that are now stale, and what they should say

| stale | current | why |
|---|---|---|
| `59.1%–81.0%` | **57.8%–80.4%** | the old one is the **`_preguard`** baseline (50 VERIFIED); current is 48, after decision 6's self-window guard. The two removed entries were the pool's **highest-support** ones, so the range does not simply shift down |
| "`0 → 241`" unqualified | **227** is now the shipping number | 241 is concat-OFF, 227 is concat-ON, and concat-ON is the default as of `eca38d2`. Plan lines 56, 77, 88 and `HANDOFF_2026-08-04_sva_decoy_fix.md:364` all still say 241 |
| "SVA position precision 90.6%" | **91.7%** | 90.6% is the shipping number only before `eca38d2` |
| `--td-decoy-concat` "OFF BY DEFAULT" | ON | fixed in `run_xtea.py`, `global_values.py`, `td_decoy_library.py`; **the plan's §4 and §5 prose still says off** |

## 6. PARKED — the previous NEXT list, with its evidence

None of these advances the objective in §1. Each is recorded so it is not re-derived.

1. **~~`--td-decoy-concat`'s default~~ — DONE**, §3.
2. **Gate [A]'s polyA-membership veto, as a MEASUREMENT only** (Corinne, 2026-08-06). Drops 68% /
   78% of the clip-qualified population including 594 / 810 loci with `n_clip >= 10`. Expect
   *precision*: 306 of BASELINE-0's 322 FNs never reach this gate. **Not started.**
3. **`NEARBY_CLIP` at gate [A]** — the largest lever measured in the step, and **not** a one-line env
   change: **11 read sites across 3 modules, 6 of them in `x_orphan_transduction.py`**, which is
   neither gate [A] nor [B2].
4. **§11's bare-shape writer fix** — scoped, not built, per §10.3.
5. **§3.4's repeat cutoff by family** — explicitly NOT the route to justifying anything about the
   decoy fix (§10.1).
6. **A third arm that is neither master nor xTea2** — still the only thing that can take the
   no-regression margin out of PROVISIONAL (§10.4). Worth restating plainly: the gate compares
   xTea2-flag-off against xTea2-flag-on, so it can say *"turning this on lost no verified
   evidence"*, but it has never been fed a change it **should** fail. The PASS is real; its power is
   unmeasured.

And two structural items, neither of them a speed item (§2):

7. **Stage 2c's single per-locus loop** — the front edge of Stage 3; dependencies built; **no
   measured speed win**, see §2.
8. **Stage 4** — somatic, mosaic, tail, legacy deletion. Untouched.

## 7. Where things live

| | |
|---|---|
| **the agenda** | `docs/handoff/DOWNSTREAM_OVERHAUL_PLAN.md` |
| **the standing decisions** | `HANDOFF_2026-08-04_td_gate_built.md` §10 — four, unchanged |
| the previous state handoff | `docs/handoff/HANDOFF_2026-08-06_two_decisions_built.md` |
| **the objective's numbers** | `docs/handoff/runs/probe_objective_metrics.py` — offline, seconds, no allocation |
| **the concat default** | `run_xtea.resolve_decoy_concat`; `ENV_TD_DECOY_CONCAT`; `--no-td-decoy-concat` is the rollback |
| **check 2b** | `benchmark/adjudicate_dropped.py` + `tests/test_adjudicate_dropped.py` |
| **the concat gates** | `gate_decoy_concat_chr1.sbatch` (4 arms), `gate_decoy_concat_wgs.sbatch` (the default path) |
| **why the per-locus loop is not a speed item** | `run_xtea.py:937`'s docstring — jobs `49550925`, `49548949` |
| **a trap in the neighbouring harness** | `docs/handoff/runs/master_serial_wgs.sbatch` opens `rm -rf "$OUT"` and its default `OUT` is the NA12878 tree the gate spec cites. Never invoke it without `MASTER_OUT` |
| running the tests | `python3 -m pytest -m "not slow"` on the system 3.9.25 — **324 passed / 30 skipped** |
