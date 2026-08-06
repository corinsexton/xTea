# xTea2 handoff — 2026-08-06 (second session), two decisions are built and the flip's last blocker is gone

**This is the current entry point for state.** The agenda entry point is still
`docs/handoff/DOWNSTREAM_OVERHAUL_PLAN.md`; the gate's design entry point is still
`HANDOFF_2026-08-04_stage3_td_gate_spec.md`.

> Supersedes `HANDOFF_2026-08-06_timing_resolved.md`. **Its §5 NEXT items 1 and 2 are DECIDED
> (Corinne, 2026-08-06) and the first half of each is BUILT.** Items 3-6 are unchanged and
> re-listed in §6.
>
> **§10's four decisions (2026-08-04) remain the standing agreement. Nothing here reopens one.**

**The one-sentence version: the guard Corinne approved found a defect two entries wide instead of
one, and the mechanism is not the one the brief named — which is why the approved fix would not
have worked.**

---

## 0. State at handoff

| | |
|---|---|
| branch | `xTea2` |
| commits this session | `90a7c47` the polyA decoupling; `b8b2ac3` its 4-arm gate; `9769a54` the verifier's self-window guard; `f5e311a` the re-freeze harness; `460e13d` the plan's decision-6 record; `d608a30` the compare-by-name fixes |
| **`xtea/` DID change** | `transduction_evidence.py`, `x_transduction.py`, `run_xtea.py`, `global_values.py` — all for the polyA flag, which is **absent-means-DERIVED on both entry points**, so BASELINE-0 and every frozen harness baseline still reproduce |
| `benchmark/` changed | `verify_transduction.py` — decision 6, ON by default, `--allow-self-window` to reproduce a pre-guard baseline |
| tests | **313 passed / 30 skipped** (`python3 -m pytest -m "not slow"`), was 254/30. 48 new in `tests/test_td_polya_cutoff.py`, 11 added to `tests/test_verify_transduction.py` |
| jobs | `49966543` the polyA chr1 4-arm gate, **PASS**, 4 arms at ~6 min. `49968034` the re-freeze dry run, 16 s — **its check 1 reported a spurious FAIL, §4.1**. `49969684` the corrected re-freeze |
| mutations | 6 run against the polyA tests, 4 against the guard's. **All 10 bite** |
| what is NOT done | the polyA veto arm (§6 item 3), which Corinne scoped as *measure, do not change* |

## 1. The three decisions, as taken

Asked directly and answered 2026-08-06:

1. **The `chr17:45585847` instrument conflict → add the source-distance guard**, argued and gated
   as its own change, accepting that it moves the frozen 3-arm baseline. §3.
2. **Build the polyA decoupling next.** §2.
3. **The polyA-membership veto → measure it, do not change it yet.** Not started; §6 item 3.

## 2. Gate [A]'s two thresholds are separable, and the default is unchanged byte-for-byte

`n_polya_cutoff` was derived as `n_clip_cutoff / 2 - 1` inside `TransductionParams`, so a
clip-support sweep was single-variable only up to 4. It is now
`run_xtea --td-polya-cutoff` / `XTEA_TD_POLYA_CUTOFF`, **absent-means-DERIVED on both entry
points** — deliberately not absent-means-zero, because 0 is a threshold a caller may legitimately
want and must stay distinguishable from "no opinion". The env channel therefore carries the
literal `"derived"` rather than `"0"`; the step is its own `x_TEA_main` subprocess, so that
string is the whole contract across the seam. A malformed value **raises** rather than reverting
to the derivation.

### 2.1 The gate: four arms, and the fourth is why the other three mean anything

`docs/handoff/runs/gate_polya_cutoff_chr1.sbatch`, job `49966543`, **PASS**.

| check | result |
|---|---|
| **1. containment** — BEFORE vs AFTER_DEFAULT | **100 identical, 0 differing, 0 missing** (+2 characterised) |
| **2. the no-op arm** — AFTER_DEFAULT vs `--td-polya-cutoff 0.5` | **100 identical, 0 differing** |
| **3. the positive control** — `--td-polya-cutoff 2` | `low_polya_clip` **0 → 67 (L1), 0 → 88 (SVA)** |
| 4. Alu, which never runs the step | identical to BEFORE in all three after arms |

**Arms 1-3 all predict byte-identity, and three byte-identical arms are equally consistent with
"the flag works and is inert here" and with "the flag never reached the gate at all."** A typo in
the env var, a resolver returning `None`, a params object built before the value was read — checks
1 and 2 pass just as happily. So arm 4 demands the threshold BITE, and check 3 says explicitly
that a zero there makes checks 1 and 2 **vacuous** rather than passing.

**And what it bites on is a dead parameter coming alive.** `low_polya_clip` is *unreachable* at
the shipping operating point — every `m_polyA` record has `sum(rcd[:6]) >= 1` by construction, 0
zero-sum entries of 4,604 (L1) / 5,261 (SVA) — so its count must be 0 in arms 1-3 and nonzero in
arm 4. **155 rejections appeared.** That is a cleaner control than a VCF delta: chr1's whole call
set is 17 calls and stayed 17 in all four arms, far too small to register a gate [A] change, which
is exactly why the ledger was the instrument.

The step also prints what it compared against:
`[td-polya-cutoff] gate [A] polyA threshold DECOUPLED: 2.0 (derived would be 0.5, from
n_clip_cutoff=3)` — the "armed vs changed anything" split, and independent confirmation that
HG002's operating point is `n_clip_cutoff = 3`.

## 3. The verifier's guard — and the approved fix would not have worked

### 3.1 The brief said proximity. It is geometry, and that changes the remedy

The standing diagnosis (previous handoff §4.5, plan §10) was *"a source close enough to the
insertion can be verified by proximity"*, with the remedy *"mirror the caller's
`MIN_POLYMORPHIC_SOURCE_DIST = 1000`, which the verifier lacks"*. Right that the verifier lacked a
guard the caller has. **Wrong about the mechanism — and 1000 would have caught nothing**, because
the entry sits **2,602 bp** from its source anchor.

What is actually wrong is that with `W = flklen = 3000` the downstream window
`[src_end, src_end + W] = 45583245–45586245` **CONTAINS the insertion at 45585847**.
`gather_candidates` collects reads AT the insertion site and `n_support` counts those aligning to
the source window — so a window spanning that site is satisfied by ordinary local coverage. The
reads are real; the inference is circular.

**So the guard is a function of the window the verifier already builds, and adds no new
threshold.** A test pins the 1000-bp fact so that a later "simplification" into a distance
constant fails and says what it breaks.

### 3.2 It is two entries, and they are the pool's two highest-support entries

| | support | window | verdict |
|---|---|---|---|
| `chr17:45585847` | **143** | downstream, contains pos at **+2602** | VERIFIED → UNDECIDABLE |
| `chr6:31329020` | **76** | **upstream mirror**, +1624 | VERIFIED → UNDECIDABLE |
| next highest (#3) | 26 | — | unaffected |
| pool median | 1 | | |

**Ranks #1 and #2 of 115, by 3×.** That ranking is the circularity's signature rather than a
coincidence: coverage at an insertion site is deep, so a window covering the site collects all of
it. **The second entry is reachable only through the UPSTREAM mirror**, so `--one-sided` hides it
— worth carrying into the spec's §8 orientation question, which is otherwise argued on where
support lands.

Shipped as **decision 6**, ON by default, checked **before** the repeat-window test (circularity
disqualifies even a pristine unique region), with `--allow-self-window` to reproduce a pre-guard
baseline.

### 3.3 This retires the flip's last blocker

`chr17:45585847` is the call the truth set counts as a **false positive** and the verifier counted
as **VERIFIED** — the one call on which the two instruments disagreed, and the whole of the Stage
3 FAIL that blocked `--td-decoy-concat`'s default. **It is now UNDECIDABLE, so they no longer
disagree.** With the timing (previous handoff §2) and the 241 → 227 (§4.1 there) already closed,
**§5's three blockers are all discharged** — see §6 item 1.

## 4. The re-freeze, and a lesson about frozen artifacts

`docs/handoff/runs/refreeze_td_baseline.sbatch`. The pool is **not** rebuilt — check 0 asserts its
sha256 against the recorded one, because the guard is a statement about what the READS prove and
rebuilding the pool would conflate a stricter verifier with a different pool.

Measured (job `49968034` dry, then `49969684` **PASS, installed**, 16 s each):

| | VERIFIED |
|---|---|
| frozen baseline | 50 |
| guard OFF (`--allow-self-window`) | **50** — 0 of 115 verdicts differ |
| guard ON | **48** — exactly the 2 entries above |

**The re-freeze is INSTALLED.** All three checks passed: containment 0 of 115 verdicts differing,
exactly the 2 known entries moving, and two guard-ON runs agreeing byte-for-byte. The pool's
sha256 is `7e797af4…` in both directories — unchanged, as check 0 required. The pre-guard baseline
is preserved at `td_gate_baseline_hg002_preguard/` with its original `verdicts.tsv`
(`f01efc40…`), and stays reproducible from the new code with `--allow-self-window`.

**Anything quoting the frozen precision range 59.1%–81.0% must now say which baseline it means**
— that range was measured over the 50-VERIFIED pre-guard verdicts, and the current baseline has
48. The two entries removed were its two highest-support ones, so the range does not simply
shift down.

### 4.1 A FROZEN ARTIFACT PINS CONTENT, NOT SERIALISATION — and the index bug underneath it

**Check 1 reported a spurious FAIL**, and the cause is worth carrying: the frozen `verdicts.tsv`
has **16 columns** and the current writer emits **18**. The restratification added `shape` and
`sidedness` *after* that file was frozen, so **every data line differs while every verdict is
identical**. `cmp` was the wrong instrument.

**The second-order bug was worse and was live in check 2: reading those columns by INDEX.**
`verdict` is field 13 in the frozen file and field 15 in the current one, so an index-based reader
silently compares `n_clip_support` against `verdict` and can report any answer at all.

Both now parse by **header name**; check 1 reports the layout change instead of failing on it, and
check 2 compares guard-OFF against guard-ON so both sides come from the same writer version.

> **Note which instrument was immune.** The `--baseline` reader passed throughout — because it
> re-derives the strata from `td_src_raw` rather than reading columns, exactly as its docstring
> says. The generalised harness was right and the bespoke check bolted next to it was wrong, which
> is the third time this project has paid for a private copy of a shipping rule (plan §10).

## 5. A gate that printed "(none found)" where its evidence should be

The polyA gate's step-banner glob was one directory too shallow — `<arm>/<sample>_*.steps.log`
instead of `<arm>/<sample>/<sample>_*.steps.log`, the mirror image of `26504d2`. Harmless to the
verdict: check 0's driver banners and check 3's 155 rejections both stand on their own. **But a
silent "(none found)" in a check whose job is to display evidence is the failure mode that gate
exists to prevent**, so it is now an assertion — 4 step-level lines across the two decoupled arms,
0 in the derived one. Verified against `49966543`'s own artifacts rather than by re-running.

## 6. NEXT

1. **`--td-decoy-concat`'s default — all three blockers are now discharged**, so this is a clean
   decision rather than a blocked one. Re-run the Stage 3 no-regression gate on the same pair
   (`49792012` OFF vs `49879677` ON) **against the re-frozen baseline**; the call it failed on is
   now UNDECIDABLE in both arms. Cost ~20 s. If it passes, the flip is justified on: byte-verified
   recipe, total containment, L1/Alu unmoved at both scales, SVA's `decoy_cross_family` 0 → 44,
   one FP removed, and ~55 s (12%) of measured cost.
2. **Gate [A]'s second semantic decision — the polyA-membership veto — as a MEASUREMENT only**
   (Corinne, 2026-08-06). It drops 68% / 78% of the clip-qualified population including 594 / 810
   loci with `n_clip >= 10`. Build the arm behind a flag, produce truth-set-scored numbers, and
   make no recommendation about the default. Expect *precision*: 306 of BASELINE-0's 322 FNs never
   reach this gate.
3. **`NEARBY_CLIP` at gate [A]** — still the largest lever measured in the step, and still **not a
   one-line env change**: 11 read sites across 3 modules, **6 of them in
   `x_orphan_transduction.py`**, which is neither gate [A] nor [B2] (previous handoff §5.1).
4. **§11's bare-shape writer fix** — still scoped, not built, per §10.3.
5. **§3.4's repeat cutoff by family** — still explicitly NOT the route to justifying anything
   about the decoy fix (§10.1).
6. **A third arm that is neither master nor xTea2** — still the only thing that can take the
   no-regression margin out of PROVISIONAL (§10.4), and the re-freeze does not change that.

## 7. Where things live

| | |
|---|---|
| **the agenda** | `docs/handoff/DOWNSTREAM_OVERHAUL_PLAN.md` |
| **the standing decisions** | `HANDOFF_2026-08-04_td_gate_built.md` §10 — four, unchanged |
| the previous state handoff | `docs/handoff/HANDOFF_2026-08-06_timing_resolved.md` |
| **the polyA flag** | `transduction_evidence.TransductionParams(n_polya_cutoff=...)`; `x_transduction.resolve_env_polya_cutoff`; `run_xtea.resolve_polya_cutoff`; tests `tests/test_td_polya_cutoff.py` |
| **its gate** | `docs/handoff/runs/gate_polya_cutoff_chr1.sbatch` — 4 arms, job `49966543` PASS |
| **the verifier's guard** | `benchmark/verify_transduction.py` decision 6 — `PoolRow.self_window`, `--allow-self-window` |
| **the re-freeze** | `docs/handoff/runs/refreeze_td_baseline.sbatch`; the pre-guard baseline is preserved at `td_gate_baseline_hg002_preguard/` and stays reproducible with `--allow-self-window` |
| **the timing instrument** | `docs/handoff/runs/probe_step_timing.py` — reports null arms, refuses to attribute without one |
| **a trap in the neighbouring harness** | `docs/handoff/runs/master_serial_wgs.sbatch` opens `rm -rf "$OUT"` and its default `OUT` is the NA12878 tree the gate spec cites. Never invoke it without `MASTER_OUT` |
| running the tests | `python3 -m pytest -m "not slow"` on the system 3.9.25 — **313 passed / 30 skipped** |
