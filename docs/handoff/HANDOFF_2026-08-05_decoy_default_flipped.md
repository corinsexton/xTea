# xTea2 handoff — 2026-08-05, the decoy default is FLIPPED, and `-1` means something narrower than we thought

**This is the current entry point for state.** The agenda entry point is still
`docs/handoff/DOWNSTREAM_OVERHAUL_PLAN.md`; the gate's design entry point is still
`HANDOFF_2026-08-04_stage3_td_gate_spec.md`.

> Supersedes `HANDOFF_2026-08-04_td_gate_built.md`. **Its §8 NEXT items 1, 2, 3 and 5 are DONE**;
> item 4 stays deliberately unbuilt per §10.3; item 6 stays explicitly-not-the-route per §10.1.
> Item 7 is what remains, and one of its three is now measured (§6).
>
> It also discharges `HANDOFF_2026-08-04_sva_decoy_fix.md`'s §0 agenda item **1** — "decide the
> default for `--td-decoy-structural`" — which had been open through two handoffs.

**Four things, and the second is the one that changes how to read an existing number.**

1. **`--td-decoy-structural` is ON by default**, measured first and flipped second, in the order
   §10.1 made an obligation rather than a preference. §2, §3.
2. **`start == -1` in a bare `TD_SRC` is a repeat-EXCUSED missing side, not merely a missing one** —
   pinned to `x_orphan_transduction.py:1120`. It turns the one-sided verification deficit from a
   correlation into a prediction. §5.
3. **The no-regression gate exists, and it refuses two comparisons rather than making them badly**:
   a baseline at a different operating point, and a single-arm pool against the 3-arm frozen
   baseline. §3.
4. **Shape × sidedness is live and it is reporting-only** — byte-identical verdicts, proven rather
   than asserted. §4.

---

## 0. State at handoff

| | |
|---|---|
| branch | `xTea2` |
| commits this session | `f10daa7` the restratification + the no-regression gate + its harness; `3d1437a` **the flip** (the first change to `xtea/` since `55055ba`); `f3dd6d0` the `-1` pinning |
| tests | **232 passed / 30 skipped** in the fast subset. 14 new: 13 in `tests/test_verify_transduction.py` (new file) and 1 in `tests/test_transduction_evidence.py` pinning the flipped driver default |
| **the environment changed under us** | the python that ran the suite last session — `/n/app/python/3.10.11.conda` — **no longer exists on the system**. The default `python3` is now 3.9.25, and `pytest` had to be reinstalled (`pip install --user pytest`). Measured: the pristine tree at this commit is **218 passed / 30 skipped** on 3.9.25, so the handoff's "31 skipped" was that python and not a regression. One env-gated skip became a pass |
| in flight | *(see §7 — filled at the end of the session)* |
| the flip's evidence | job **49791196**, `docs/handoff/runs/td_gate_no_regression.sbatch`, **PASS in 19 s**. Artifacts `/n/scratch/users/c/cos689/xtea2_testdata/td_gate_no_regression_49791196/` |
| the frozen baseline | unchanged and re-verified by sha256 at the top of that job: `/n/scratch/users/c/cos689/xtea2_testdata/td_gate_baseline_hg002/` |
| what is NOT done | §6's three: plan §5's semantic half on gate **[A]**, the decoy-consensus concatenation, and `ileftmost_cns` 500-vs-5000 (now measured — §6.1) |

## 1. What the restratification cost, which is nothing, and why that is a result

§10.2 promoted `family × call class × shape × sidedness` over the source-class axis. The obvious
risk in a restratification is that it quietly moves the thing it is reporting on, so the harness
checks it the strongest available way: **the per-entry verdicts file must be byte-identical to the
frozen one once the two new columns are stripped.** It is.

So the primary table is new *reporting* over unchanged *measurement* — and the table earns its
promotion on the first run:

| stratum | VERIFIED | REFUTED | SUB_THR | UNDEC | precision |
|---|---|---|---|---|---|
| L1 regular bare **two_sided** | 23 | 4 | 2 | 1 | **85.2%** |
| L1 regular bare **one_sided** | 9 | 7 | 1 | 3 | **56.2%** |
| L1 sibling sibling_point | 2 | 2 | 2 | 20 | 50.0% |
| SVA orphan tailed two_sided | 5 | 0 | 0 | 1 | 100.0% |
| SVA regular tailed two_sided | 2 | 0 | 0 | 1 | 100.0% |
| SVA sibling sibling_point | 9 | 7 | 5 | 5 | 56.2% |
| L1 regular tailed two_sided | 0 | 1 | 0 | 2 | 0.0% (n=1) |
| L1 orphan tailed two_sided | 0 | 0 | 0 | 1 | *empty* |

**The 29-point sidedness gap reproduces on HG002** (85.2% vs 56.2%) after being measured on
master's NA12878 (83.9% vs 55.0%). Two samples, same direction, same size — and §5 now supplies the
mechanism.

**One decision inside the axis that is not cosmetic.** A `sibling_point` gets sidedness `-`, not
`two_sided`. Its `src_start` is the pool builder's own `[POS-100, POS+100]` bound, so reading
sidedness off it would report an invented number — **and it would read as the CONFIDENT value**,
which is the direction that does damage. Pinned by a test.

## 2. The flip, and what it does NOT rest on

`--td-decoy-structural` is on by default as of `3d1437a`. Per §10.1 the justification is
**structural correctness** — a substring test that misclassifies 100% of hg38's 5,454 SVA flank
records as decoy consensus is wrong regardless of how many calls it moves — **plus the step-level
recovery 0 → 241 reference copies**, measured on `unique_trsdct_disc_only` column 3, the
flank-derived field.

**It does not rest on the +2 VERIFIED entries the gate found.** Those are orphan-class, and §10.1
explicitly refuses to rest a shipping default on the class plan §1 says is not an MEI insertion at
all. They are corroboration. Stated here because the temptation to promote them will recur: the
gate's job was to show no regression, and that is all it showed.

**And the structural blocker survived the flip rather than being retired by it.** `SVA regular
reference_copy` is still **0 VERIFIED / 0 REFUTED / 1 UNDECIDABLE** — SVA_D at 2.7% divergence,
disqualified by §3.4's repeat test by construction, exactly as §4 predicted. The one entry that
would have made the fix's *regular* class scorable is the one entry that cannot be scored.

### 2.1 Three things about the flip that are changes of KIND, not of default

These are the parts a future reader will trip on, so they are in the code as well as here.

1. **`run_xtea` is now the SINGLE WRITER of the policy** and sets `XTEA_TD_DECOY_STRUCTURAL` on
   every path, `"1"` or `"0"`. Before the flip it wrote only on opt-in, so that an inherited `"0"`
   could not fake an opt-in. With the fix on by default **the risk reverses**: an inherited value
   could fake an opt-OUT and silently restore the defect.
2. **The env var's own rule is UNCHANGED — absent still means OFF** (`x_transduction.py:78`, and
   `tests/test_transduction_evidence.py`'s existing pin is untouched, which is evidence rather than
   assertion). The default is a **driver** policy with a blast radius of one entry point. A dozen
   frozen gate harnesses under `docs/handoff/runs/` invoke `x_TEA_main` **directly**, with baselines
   measured under the substring test; if absent meant ON, this one flip would silently move all of
   them at once. **Do not "fix" this asymmetry without re-measuring those baselines** — both files
   say so at the line.
3. **`--(no-)td-decoy-structural` keeps THREE states** via `default=None`, and the resolution is a
   module-level `run_xtea.resolve_decoy_structural` rather than two lines inside `main` — a shipping
   default that can only be exercised by running the whole driver is a default with no test. The
   banner names the channel (`DEFAULT since 2026-08-05` vs the flag) because two gate harnesses grep
   it as their only evidence of which configuration ran.

### 2.2 The off arm now has to say so, and one gate needed changing for it

`gate_decoy_structural_chr1.sbatch`'s **AFTER_OFF arm passed no flag**, which was correct when
nothing was off. Left alone it would have become AFTER_ON under another name, and the gate would
have reported a byte-identity FAIL **it had itself caused**. It now passes
`--no-td-decoy-structural`.

Its check 0 also gained a **positive** check on a new SUBSTRING banner. "No STRUCTURAL banner"
cannot be told from "the driver died before it got there" — an off-arm gate whose only evidence is
an *absent* line passes when the run fails early. And a precondition grep refuses the pre-flip
commit *before* three arms of compute rather than after.

`gate_decoy_structural_wgs.sbatch` needed **no change**: it passes the flag explicitly, still an
accepted opt-in, and the banner prefix it greps is preserved verbatim.

## 3. The no-regression gate, and the two comparisons it refuses

`verify_transduction.py --baseline <verdicts.tsv>` gates on **absolute VERIFIED counts per
stratum** (§10.4). Two refusals are the load-bearing parts:

- **A baseline at a different operating point is REFUSED (exit 3)**, not warned about. Every
  VERIFIED count is a function of all five thresholds — that is decision 5 and the entire reason
  the sweep exists — so a baseline at N=3 against a run at N=2 would report growth and regression
  that are artifacts of N. The operating point is parsed from the header this program itself
  writes, and `format_op_point`/`_OP_POINT_RE` round-trip is a test.
- **The strata are RE-DERIVED from `td_src_raw` and `src_start`, never read from columns.** That is
  what lets the frozen baseline — written before the shape axis existed — be compared with one
  implementation of the rule and no migration. Verified against the real frozen file: it reproduces
  all 50 VERIFIED entries, and both axes partition the same set.

**And the third refusal is in the harness, not the tool.** `td_gate_no_regression.sbatch` does NOT
compare a post-flip single-arm pool against the frozen 3-arm baseline. That pool holds ~90 entries
against 115, of which master contributes 21 uniquely, so it would "lose" verified entries in almost
every stratum **from arm scope and not from the flag** — a large FAIL that is entirely an artifact
of which arms are in the pool. This is the same error the plan records three times. The frozen
3-arm baseline is the right reference for the gate's **precision range** and the wrong reference for
a **one-flag delta**.

So the harness runs two different questions as two checks:

| check | what it compares | result |
|---|---|---|
| 0 | the frozen baseline's sha256s | PASS |
| 1 | the frozen 3-arm pool, re-adjudicated | **PASS**, all strata +0, per-entry file byte-identical, range 59.1–81.0% unchanged |
| **2a** | is the flip **additive** at the pool level | **PASS — strict superset**, 90 → 94 |
| 2 | BASELINE-0 (OFF) vs the decoy arm (ON), per stratum | **PASS**, every pre-existing stratum +0, `SVA orphan tailed two_sided` 0 → **+2** |

**Check 2a is separate from check 2 for a structural reason.** A *dropped* hypothesis is invisible
to a per-stratum VERIFIED count — it is absent from both sides of it. So "no stratum lost a verified
entry" and "no hypothesis was removed" are different claims needing different checks, and only one
of them is what the flag could plausibly break.

Measured, and it is sharper than the previous handoff had it: the ON pool is a **strict superset**
of the OFF pool. All 90 OFF hypotheses appear verbatim — same locus, same source, same `TD_SRC` —
plus exactly the 4 SVA reference-copy entries §4 names. So on this pair of arms a VERIFIED loss
could not have come from the flag at all; it could only have come from the same entry, in the same
BAM, at the same operating point, getting a different verdict. **Check 2 is therefore also a
determinism check**, and a future FAIL should be read that way first.

Precision ranges, reported and not gated: OFF **54.2–76.7%**, ON **54.9–76.7%**. The low end rises;
the high end is identical.

**The margin prints PROVISIONAL and the tool says why.** Two of the three arms differ by one flag,
so 82-of-115 agreement is a floor on independence, not a measure of it. And PASS means *no stratum
lost a verified entry*, which is weaker and different from *precision did not fall* — a stratum can
hold its VERIFIED count while gaining REFUTED entries, and this gate says nothing about that by
construction.

## 4. What the gate cost to run, which is the surprising part

**19 seconds**, 8 cores on `short`, for three full adjudications (115-entry pool + 90 + 94) against
a 43 GB WGS BAM. It reads ~2,400 candidate read parts over 115 entries and runs two bwa passes.

This matters for how the gate gets used: it was implicitly treated as expensive enough to freeze a
baseline and compare against it. At 19 seconds it can be a **per-change gate**, and the reason to
freeze a baseline is reproducibility rather than cost. No `priopark` reservation is needed —
`-c 8 -p short -A park` is enough.

## 5. `start == -1` is a repeat-EXCUSED missing side — §8 step 5, and the answer is narrower than the question

The previous handoff had `-1` as "one of the two discordant-pair midpoints was never established",
with the mechanism inferred from `min(src_pos_start, src_pos_end)` and **explicitly flagged as not
pinned to a line**. It is now pinned, and the extra specificity changes what the stratum means.

`-1` is `mid_pos`'s initialiser in `_form_unique_cluster` (`x_orphan_transduction.py:1064`),
returned unchanged when a side has too few discordant pairs (`:1072`) or no dominant source
chromosome (`:1082`). **But a failed side alone does not produce a surviving record.** `b_consistent`
(`:596-619`) requires one side to have formed, and for the other it calls
`_check_form_repeat_or_same_cluster`, which returns:

| return | line | meaning |
|---|---|---|
| `True, -1` | `:1120` | the missing side's disc mates are **predominantly INSIDE REPEATS** (`n_in_rep > n_half_disc` and the ratio clears `f_rep_ratio`) |
| `True, <pos>` | `:1126` | they cluster with the other side instead — a **real** position |
| `False, -1` | `:1128` | no explanation; the record does not survive |

So in surviving output **`start == -1` means the second side never clustered AND its absence was
excused because its discordant mates fell inside repeat regions.**

**That predicts the verification deficit instead of accompanying it.** Those entries are precisely
the ones whose second side is repeat-driven — and repeat-driven discordant support is exactly what
`verify_transduction`'s decision 4 exists to distrust. The ~30-point gap (56.2% vs 85.2% here;
55.0% vs 83.9% on NA12878) now has a mechanism, which also means it should not be read as "the
caller is 30 points worse when it only has one side": it is "the caller accepted a repeat-region
excuse for the second side, and the reads mostly do not back that".

Byte-identical in master (`:1093`, `:1148`, `:1155`), so it is inherited, not an xTea2 behaviour.

**And the parser's assertion is now sourced from the writer.** `END` cannot be `-1`: the writer
takes `i_end = max(...)` (`:359-363`) and at least one side is a formed cluster, whose position
comes from `_position_of_same_chrm_form_cluster` (`:1086-1096`) as the median of a **non-empty**
list. So `parse_td_src`'s `100--1` branch is **unreachable for this writer** rather than merely
unobserved. Kept anyway, as the tripwire for a third writer.

## 6. NEXT

**In this order. Nothing below is blocked on compute.**

1. **The previous handoff's items 2 and 3, unchanged.** Plan §5's semantic half on gate **[A]** —
   97% of the rejecting, expect *precision*; and **the decoy-consensus concatenation**, which is what
   makes §5's SVA/Alu hazard testable at all. The concatenation was sequenced "after the default
   decision in item 1" so the two changes' effects on SVA stay separable — **item 1 is now done, so
   that sequencing constraint is discharged.**
2. **`ileftmost_cns` 500-vs-5000** — see §6.1, now measured rather than open.
3. **§11's bare-shape writer fix** — still scoped, not built, per §10.3, and §5 above strengthens the
   case for it: the source-class axis is empty for 108 of 115 entries *because* the bare writer emits
   disc midpoints. Read §11's trap first — the midpoints are what the verifier anchors on and must
   stay available.
4. **Whether §3.4's repeat cutoff should differ by family** remains explicitly NOT the route to
   justifying anything about the decoy fix (§10.1). §5 above raises its stakes rather than lowering
   them: if `-1` entries are repeat-excused by construction, then the repeat cutoff and the
   sidedness stratum are measuring overlapping things, and that is worth understanding *before*
   anyone touches the cutoff.
5. **A third arm that is neither master nor xTea2**, which is the only thing that can take the
   no-regression margin out of PROVISIONAL (§10.4). Nothing else on this list unblocks it.

### 6.1 `ileftmost_cns` — *(filled below)*

## 7. The two gate scales, re-run after the flip — *(filled below)*

## 8. Where things live

| | |
|---|---|
| **the agenda** | `docs/handoff/DOWNSTREAM_OVERHAUL_PLAN.md` |
| **the gate's design** | `docs/handoff/HANDOFF_2026-08-04_stage3_td_gate_spec.md` (corrections inline, `[CORRECTED]`) |
| the previous state handoff | `docs/handoff/HANDOFF_2026-08-04_td_gate_built.md` — **§10's four decisions are the standing agreement** |
| the pool builder | `docs/handoff/runs/probe_td_pool.py` |
| the verifier + the no-regression gate | `benchmark/verify_transduction.py` (`--baseline`) |
| **the no-regression harness** | `docs/handoff/runs/td_gate_no_regression.sbatch` — job 49791196, PASS, 19 s |
| the gate's frozen baseline | `/n/scratch/users/c/cos689/xtea2_testdata/td_gate_baseline_hg002/` (+ `SHA256SUMS`) |
| the flip | `xtea/run_xtea.py` (`resolve_decoy_structural`, the two flags, the two banners), `xtea/global_values.py:178-190`, `xtea/x_transduction.py:72-86`, `xtea/transduction_evidence.py`'s `decoy_tags` docstring |
| the flip's worktree | `/n/scratch/users/c/cos689/xtea2_wt/td_decoy_flip_3d1437a` (`3d1437a`); `src/` is byte-identical to `365f352` so one prebuilt `xtea_scan` serves every arm |
| the two decoy gates | `docs/handoff/runs/gate_decoy_structural_{chr1,wgs}.sbatch` — **chr1's AFTER_OFF now needs `--no-td-decoy-structural`** |
| **a trap in the neighbouring harness** | `docs/handoff/runs/master_serial_wgs.sbatch` opens `rm -rf "$OUT"` and its default `OUT` is the NA12878 tree the gate spec cites. Never invoke it without `MASTER_OUT` |
| the truth-set harness | `benchmark/benchmark_truth_set.py` — `parse_td_src`, `td_src_shape` (**now carries the `-1` mechanism**), `td_src_call_class`, `td_src_source_class` |
| running the tests | `python3 -m pytest -m "not slow"` on the system 3.9.25 (see §0 — the old 3.10 conda python is gone) |
