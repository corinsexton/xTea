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
| commits this session | `f10daa7` the restratification + the no-regression gate + its harness; `3d1437a` **the flip** (the first change to `xtea/` since `55055ba`); `f3dd6d0` the `-1` pinning; `53ca400` a self-review fix to the flip (§2.3) |
| tests | **233 passed / 30 skipped** in the fast subset. 15 new: 13 in `tests/test_verify_transduction.py` (new file) and 2 in `tests/test_transduction_evidence.py` pinning the flipped driver default and its parser contract |
| **the environment changed under us** | the python that ran the suite last session — `/n/app/python/3.10.11.conda` — **no longer exists on the system**. The default `python3` is now 3.9.25, and `pytest` had to be reinstalled (`pip install --user pytest`). Measured: the pristine tree at this commit is **218 passed / 30 skipped** on 3.9.25, so the handoff's "31 skipped" was that python and not a regression. One env-gated skip became a pass |
| in flight | **nothing.** All three jobs reported: `49791196` no-regression **PASS** 19 s, `49792011` chr1 containment **PASS** 18:19, `49792012` WGS payoff **PASS** 38:32 |
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

**The flip's cost, stated plainly.** The WGS gate's truth-set score (reported, never gated) has SVA
position precision going **93.9% → 90.6%**: three new calls, **all off-truth**, TP and FN unchanged.
L1 and Alu are byte-identical in every figure. §7.2 has the numbers and the reason that is *not* an
explanation for them — and it also records an inference I made and then refuted, which reaches
§10.1's own wording. Anyone re-opening the flip should start there rather than here.

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
6. **NEW — the two "orphan"s, and it is small and cheap.** `SUBTYPE == orphan_or_sibling_transduction`
   and `td_src_call_class == "orphan"` are disjoint populations (§7.2), and both get called "orphan"
   in prose, including in plan §1 and in §10.1's reasoning. Decide which one each existing claim
   meant, and give them different words. This is documentation plus possibly one scoring decision —
   **not a cleanup**, because re-keying the scorer's carve-out would move a reported precision
   number.
7. **NEW — why the three recovered SVA calls are all off-truth** (§7.2). The honest options are that
   they are correct transductions the MEI-position truth set cannot represent, or that they are wrong.
   The verifier can adjudicate two of them from the BAM and did (both VERIFIED, `chr20:45799810` on 15
   clip + 3 disc), which is evidence but not the same question. Worth resolving before the flip is
   quoted as free.

### 6.1 `ileftmost_cns` — measured, and the live 10× difference is ZERO

The item was framed as "an unexplained 10× difference" between L1's 500 and SVA's 5000. **That
framing does not survive the code, and the code already said so.** `family_params.py:18`:
"`leftmost_cns` was therefore **always 500** — SVA's 5000 branch". `x_TEA_main.py:311` resolves the
bundle from `from_sva_flag(options.sva)`, and `run_xtea.build_type_steps` passes `--sva` to
**cns_filter alone** (matching the golden `run_xTEA_pipeline.sh`), so the transduction step resolves
`GENERIC` = 500. **Today both families run 500.** The question was never "why do they differ" — it
is "should SVA's dormant 5000 be woken", which is a documented one-word change.

Measured through the shipping `pick_full_length_insertion`, identical across three independent
finished runs (`docs/handoff/runs/probe_leftmost_cns.py`):

| `ileftmost_cns` | L1 (of 300) | SVA (of 292) |
|---|---|---|
| 100 | 53 · 17.7% | 82 · 28.1% |
| **500** | **58 · 19.3%** ← LIVE | **162 · 55.5%** ← the value actually in use |
| 1000 | 60 · 20.0% | 251 · 86.0% |
| 2000 | 69 · 23.0% | 289 · 99.0% |
| **5000** | 167 · 55.7% | **292 · 100.0%** ← SVA canonical, DORMANT |
| 10000 | 300 · 100.0% | 292 · 100.0% |

**`family_params`' claim is CONFIRMED and is now a measurement rather than an inference from a
sequence length.** Waking SVA's 5000 takes it 162 → 292, **+130 candidates (+80%)** into the
polymorphic-flank library and every realignment built from it. The consensus-length mechanism
checks out independently: `consensus/SVA.fa`'s longest record is **3183 bp** (< 5000, so SVA
saturates) and `consensus/LINE1.fa`'s is **6342 bp** (> 5000, so L1 does not, reaching 100% only at
10000). The probe asserts that **asymmetry** rather than a round number, because the asymmetry is
what breaks if the libraries or the column layout move.

**NEW, and not in the docstring: at the value actually in use the two families are nowhere near
equivalent — L1 19.3% against SVA 55.5%.** So "both run 500" is not "both are filtered the same";
SVA's cns clusters sit much closer to the 5′ end. **That is the thing to explain before deciding
anything**, and it is arguably more interesting than the dormant 5000.

**What would settle the wake.** A truth-set-scored two-arm gate — SVA at 500 vs SVA at 5000 — on the
`gate_decoy_structural_wgs.sbatch` pattern. The known obstacle is the one §7.1 already records: the
HG002 truth set's derived transduction denominator is **0 (L1) and 2 (SVA)**, so it cannot score
transduction calls *as* transductions; what it can score is SVA position sensitivity/precision. Do
not flip it blind — `family_params`' own words.

> **A methodological finding worth more than the numbers.** This probe's first version carried its
> own copy of `_reach_5_prime_end` *and* its own candidate parser while its docstring claimed to
> import the shipping predicate. The copies were not equivalent, and not subtly:
> `load_in_candidate_list_str_version` stores `fields[2:]`, dropping chrom and pos, so the shipping
> `rcd[17..20]` are the line's **columns 20-23** — the four cluster fields. Reading the whole line
> put `rcd[17]` two columns early and fed plain **counts** into the cluster test (`0 < 4 < 500`
> passes), so L1 at 500 read **50.3% instead of 19.3%**. Plausible numbers, 2.6× wrong, and they
> were on their way into this handoff. This is the concrete case for the rule the other probes
> follow, and it is now in the probe's docstring so nobody re-copies the predicate.

## 7. The two gate scales, re-run after the flip — §8 step 3's third leg

Both re-run at `3d1437a` against the worktree `/n/scratch/users/c/cos689/xtea2_wt/td_decoy_flip_3d1437a`,
sharing `td_b2_365f352`'s prebuilt `xtea_scan` (`src/` is byte-identical across the range, so the
C++ is removed as a variable rather than introduced as one).

> **Neither gate exercises `53ca400`'s fix**, and that is checkable rather than assumed: every arm
> passes an explicit flag — AFTER_OFF `--no-td-decoy-structural`, AFTER_ON and the WGS arm
> `--td-decoy-structural` — and BEFORE is the parent commit, which has no banner at all. The
> not-asked path is not on any arm's route, so they were not resubmitted.

### 7.1 chr1 containment — job **49792011**, **PASS**, 18:19

**The claim that mattered: the containment claim survives the change of channel.** AFTER_OFF now
expresses "off" with a flag instead of with silence, and check 1 says that reproduces BEFORE
**byte-for-byte over the whole output tree** — 100 identical, 2 tolerated (the two characterised
exceptions), **0 differing, 0 missing**. So `--no-td-decoy-structural` is a real rollback and not
merely a different code path that happens to agree.

| check | result |
|---|---|
| 0 · the arms are the code and the config they claim | BEFORE STRUCTURAL banner **0**; AFTER_OFF STRUCTURAL **0** *and* SUBSTRING **1**; AFTER_ON STRUCTURAL **1** |
| 1 · AFTER_OFF ≡ BEFORE, whole tree | 100 identical / 2 tolerated / **0 differing** |
| 2 · no new or removed files in either AFTER arm | 0 / 0 |
| 3 · the prediction: AFTER_ON moves SVA and nothing else | **L1=0 Alu=0 SVA=1**, other=0 |
| 5 · the ledgers | AFTER_ON/SVA **2281 → 2315 rows**, with `no_dominant_source=1` and `low_disc_support` **2 → 35** appearing |
| 6 · call counts | identical across all three arms |

**The new positive banner check earned its place immediately**: AFTER_OFF's `SUBSTRING banner: 1` is
what distinguishes "the opt-out reached the driver" from "the driver never got there", and without it
check 1's PASS would have been consistent with an arm that died early.

**Check 5 is where the flag proves it is WIRED at chr1, and check 6 is where chr1's limits show.**
Call counts are identical in all three arms — but the ledger moved: 34 more rows, and
`low_disc_support` going 2 → 35 is exactly hits that used to be discarded at the decoy branch now
reaching the source gate and being rejected there instead. Per the gate's own header, **a zero call
delta at chr1 is not evidence the fix does nothing** — chr1 has 0 SVA gate-B survivors, so it has no
population to move. That is what the WGS leg is for.

### 7.2 WGS payoff — job **49792012**, **PASS**, 38:32

**The payoff check is unambiguous: SVA reference-copy picked sources 0 → 241**, out of 346, against a
5,454-record library built to be exactly that source set. The library measurement is on this arm's
own run: 5,454 of 5,454 SVA flank records match the substring test, 0 of them are genuine decoy
records — 100% collision, which is the defect the flip fixes.

| check | result |
|---|---|
| 0 · the arm ran with the fix ON | banner 1 |
| 1 · **the prediction** | L1 **IDENTICAL** to BASELINE-0, Alu **IDENTICAL**, SVA and merged DIFFER — 2 and 2, exactly as predicted |
| 2 · gated call counts | **L1 = 188**, **Alu = 1223** (both exact); SVA 82 → **85**, merged 1493 → 1496 |
| 3 · ledgers | vocabulary-clean, 18 reasons, no new reason, 2 ledgers, L1 = 21158 rows exactly |
| 4 · `downstream_parallel` | 471.5 s (bound ≤ 900 s), children peak RSS 18.18 GB |
| **5 · THE PAYOFF** | **SVA reference-copy picked sources 241** (BASELINE-0: 0 of 105) — **PASS** |

The three gained SVA positions are `chr12:39235233`, `chr20:45799810`, `chrY:10120435` — three of the
four entries §4 named, and the same three the no-regression pool saw as additions.

#### The truth-set score moved, it is REPORTED not gated, and it is not good news

| | BASELINE-0 | this arm |
|---|---|---|
| SVA calls | 82 | **85** |
| SVA TP / FN | 77 / 13 | **77 / 13** (unchanged) |
| SVA FP | 5 | **8** |
| SVA position precision | 93.9% | **90.6%** |
| all-families precision | 88.4% | 88.2% |

**All three new SVA calls are off-truth.** TP does not move, FN does not move, so the fix bought
three calls and none of them is a truth record. L1 and Alu are byte-identical in every figure.

Read it the way the gate's own header says to: the HG002 truth set enumerates MEI insertion
**positions**, and its derived transduction denominator is 0 (L1) and 2 (SVA), so it **cannot score
transduction calls as transductions**. What it can say is whether new SVA calls land on truth. They
do not. That is a real cost to weigh against the structural argument, and it should not be waved away
— but neither is it a verdict, because a transduction call landing off the MEI-insertion truth set is
what §7.1 predicts for a *correct* orphan-class transduction as much as for a wrong one.

> **AN INFERENCE I MADE AND THEN REFUTED, recorded because it is a trap.** Two of the three new calls
> carry `~orphan` in `TD_SRC`, so I expected them to fall inside the scorer's existing orphan
> carve-out — the one that excludes L1's 37 "orphan transduction SITES" with the reasoning "these are
> not MEI insertions and the truth set does not enumerate them, so they cannot match". **They do not,
> and the reason is that there are two different "orphan"s.** The carve-out keys on
> `SUBTYPE == orphan_or_sibling_transduction`. Cross-tabulated over every transduction-asserting call
> in this run:
>
> | family | SUBTYPE | `td_src_call_class` | n |
> |---|---|---|---|
> | L1 | `orphan_or_sibling_transduction` | **regular** | **37** |
> | SVA | `one_side_and_half_transduction` | orphan | 2 |
> | SVA | `one_half_side_tprt` | orphan | 1 |
>
> **The two fields do not agree on a single call.** All 37 of the carve-out's L1 calls are `regular`
> by `TD_SRC`; SVA has **zero** `orphan_or_sibling_transduction` calls, so the carve-out never fires
> for it at all. So SVA's precision drop cannot be explained away as "the same category L1 excludes",
> and my first reading of it was wrong.
>
> **The consequence reaches §10.1's own argument.** It rejected resting the flip on the 2 orphan-class
> VERIFIED entries because "plan §1 says [orphan] is not an MEI insertion at all (BASELINE-0's L1
> orphans are 0/37 on-truth)" — but the 0/37 is about **SUBTYPE** orphans and the 2 VERIFIED entries
> are **TD_SRC** orphans. Different populations, so the cited evidence was not evidence about them.
> **The decision is unaffected** — not resting on them is still the conservative choice, and the flip
> rests on structural correctness plus 0 → 241 — but the reason as written conflates two fields.
> Recorded at the carve-out in `benchmark_truth_set.py` so the next reader does not repeat it, and
> the carve-out is deliberately **left keyed on SUBTYPE**: that is the population 0/37 was measured
> on, and re-keying it would move a reported precision number under cover of a comment.

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
