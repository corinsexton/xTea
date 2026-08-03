# xTea2 handoff — 2026-08-03 evening (the loop's front half is real; its middle box is not what §3 drew)

**This is the current entry point for state. It is NOT the entry point for agenda — the plan is.**

> **Read `docs/handoff/DOWNSTREAM_OVERHAUL_PLAN.md` first.** It still outranks this chain on agenda.
> Everything factual below is written into the plan at the section that owns it; this is the index.
> Supersedes `HANDOFF_2026-08-03_site_enumeration.md`, whose §6 NEXT item 1 is the per-locus loop.
> This session built the **first two boxes** of that loop and found that its **third** box cannot be
> a loop in the sense §3's diagram implies.

**One-line status: the per-locus loop's iteration domain is now enumerated by the driver
(`xtea/downstream_unified.py`, `--locus-universe`) — 41,604 loci at WGS, 3,271 on chr1 — and the
`--td-manifest` gather that has had no production caller since 2026-08-02 now has one
(`--td-gather-driver`). Both are OFF by default and answer-preserving: the universe pass is gated
PASS end to end on chr1 (job `49547574`, 9/9 identical, 2/2 replicas byte-identical) and offline at
WGS scale (`gate_locus_universe.py`, 4/4 arguments + 2/2 files). The barrier they run at is *forced*,
not chosen — the site-creation rule excludes against `cns_filter`'s own output, so "up front" means
before **transduction**, not before `cns_filter`, and two docstrings said otherwise. The 7-of-22
ceiling is unchanged, and there is now a number that says why it is structural: **38,815 of the
41,604 loci are reachable only by creation.**

## 0. State at handoff

| | |
|---|---|
| branch / HEAD | `xTea2` @ `1c31801` (universe `440e846`, gather caller `7dfae51`) |
| working tree | clean; all work committed (`440e846` universe, `7dfae51` gather caller) |
| tests | **187 passed, 1 skipped** with `XTEA2_SCAN_BIN` set, excluding `test_downstream_from_cns.py`. Was 173 — **+14** in the new `tests/test_downstream_unified.py`. The skip is `test_e2e_concordance` (needs `XTEA2_OUTPUT_DIR`) |
| in flight | **the WGS gate** (`gate_locus_universe_wgs.sbatch`) — see §5.4. Both chr1 gates are done and **PASSED**: job `49547574` (two arms) and job `49548949` (three arms, incl. the driver gather) |
| the agenda | the loop's **third** box: the per-locus decision. §4 says what shape it can have |

## 1. What was built, and what each piece is worth

### 1.1 `xtea/downstream_unified.py` — the iteration domain, driver-side (`--locus-universe`)

Runs `site_universe.write_slct_file` for every family in the driver, at the `cns_filter →
transduction` barrier, and unions the results into `merged/locus_universe/locus_universe.tsv` with
per-family provenance (did this family reach the locus by a **cns row** or by **creation**).

**Diagnostic, deliberately.** Each family's transduction step still derives its own `.slct` through
the same one copy of the rule, so the driver writes a *replica* and no call can move. What that buys
is a gate that is a `cmp` rather than an argument — and the replica is what `--td-gather-driver`
consumes.

### 1.2 `--td-gather-driver` — the caller `--td-manifest` was waiting for

`xtea_scan --td-manifest` has been gated since 2026-08-02 (job `49516920`, 6/6 byte-identical at WGS)
with **no production caller**, because hoisting the gather meant re-implementing the step's site
selection. It does not any more, and the last named gap (the step's `iextnd`) turned out to be
something **the driver already had** — see §2. So this gathers every flanked family's transduction
reads at the barrier out of ONE BAM decode, onto the exact two paths the step would have written, and
the step skips its own gather via `XTEA_TD_GATHER_PRECOMPUTED`.

**What it is worth, MEASURED and less flattering than expected: it returns CPU and costs wall
clock.** Job `49548949`, three arms at one commit, chr1:

| | OFF | universe only | GATHER |
|---|---|---|---|
| step `native_gather` L1 / SVA | 21.97s / 22.48s | 21.97s / 22.48s | **none / none** |
| driver gather (serialized) | — | — | **19.6s** |
| transduction step L1 / SVA | 30.50s / 41.18s | 30.82s / 41.73s | **9.27s / 19.82s** |
| `downstream_parallel` | 102.9s | 110.2s | **114.2s** |

So **~25s of CPU returned** (44.5 → 19.6) and **~12s of wall clock lost**. Two structural reasons:
the wave split alone costs ~7s at chr1 in extra `x_TEA_main` launches, and removing 22s from a ~100s
concurrent phase whose critical path is sibling/postF does not shorten it. The WGS ratio is different
(53-57s of gather per family against a 457s downstream) and is **unmeasured end to end** — §5.4.
**The reason to do it is not speed**: it is the per-locus loop's batched evidence stage, built once
instead of twice, and the only place one BAM decode can serve every family.

## 2. The four things that were wrong or unknown, and are now neither

1. **"Up front" meant the wrong place.** `site_universe.py` and `x_transduction.py` both said the
   rule's inputs are "available before `cns_filter` runs". raw_disc and the annotation are; **the cns
   candidate list is `cns_filter`'s own output**, and the rule excludes against it. So the universe is
   computable before **transduction**, not before `cns_filter` — the same point the Stage 2c barrier
   occupies. Both docstrings are corrected. Placing the pre-pass on that sentence would have it
   excluding against a file that does not exist yet.
2. **The `iextnd` the previous handoff named as the `--td-manifest` blocker was never missing.** The
   step's is `max(400, int(mean_is + 3*std_var))` (`x_TEA_main.py:473`, `:497-499`) and the driver's
   `cns_extnd` is that second term already; the 400 floor never binds at WGS insert sizes (954).
   `downstream_unified.transduction_extnd()` is the one copy.
3. **The `annotation1` trap.** Transduction is handed `{g}_L1.fa.out`; `cns_filter` is handed
   `{g}_L1_larger_500_with_all_L1HS.out`. Both exist, both look like an L1 annotation, and the wrong
   one changes the rule's third rejection. This is why the new offline gate **checks the arguments,
   not only the bytes** — see §5.1.
4. **§3's middle box cannot be a loop.** `call_candidate_transduction_v3` builds one flank reference
   and runs **four batched bwa passes over the whole `.slct` list** before any per-site record exists.
   One bwa per locus is a different program, not a refactor. §4.

## 3. The number that reframes the ceiling: 38,815 of 41,604

The universe carries the provenance, so the barrier's own limit can be measured on the *domain*
instead of on the final call set:

| | loci |
|---|---|
| union universe (WGS, exact position) | **41,604** — 7,910 / 33,139 / 555 by one / two / three families |
| of those, **no family has a cns row** | **38,815** |
| chr1 equivalents | 3,271 universe, 3,230 created-only |

**Quote it as a fact about the iteration domain, never as a call count.** The transduction evaluation
rejects almost all of them; only 52 become final records. What it says is that the created population
**dominates** the universe rather than being a tail — so the 7-of-22 assignment ceiling is
*structural*, not a small-numbers accident. `--family-assign-early` still stays OFF.

## 4. What the loop's remaining box can be

The corrected shape, recorded at the §3 diagram itself:

```
enumerate loci (DONE)  →  batched per-family evidence production (gather DONE)  →  per-locus decision  →  assign family  →  filter
```

Only the third box is a loop. **And its record type already exists** — this is the previous handoff's
open question ("what does a per-family evidence record look like at a locus only `.slct` knows
about?"), answered by reading rather than by a probe:

- `collect_ins_info` assembles `(n_l_disc, lcns_start, lcns_end, n_r_disc, rcns_start, rcns_end,
  n_cns_clip, n_polyA_clip, i_lcov2, i_rcov2, s_src)`;
- `_collect_gntp_features` returns the 10-field genotype tuple;
- `add_new_rcd_with_cutoff` writes a **created** site out of exactly those, in the same format as an
  updated one.

So a created locus's evidence has the **same shape as any other locus's**. Do not invent a record
type for the loop — this is it.

What is genuinely still unknown is only the *ordering inside the tail*: where family assignment sits
relative to the sibling/orphan pass, and whether the four bwa passes can be shared across families
the way the gather now is (the flank reference is per-family and built from that family's own
full-length copies, so probably not — but it has not been measured).

## 5. The gates

### 5.1 `docs/handoff/runs/gate_locus_universe.py` — offline, WGS scale, ~5s. **PASS**

It differs from `gate_site_universe.py` in the one way that matters, and both are needed. That gate
feeds the rule **the arguments a run's own step log records**, deliberately, so it cannot drift into
testing the driver. This one builds every argument **the driver's way** — `REP_TYPES` layout,
coverage-derived cutoff, the step's own 450 boundary extend, 225 min copy — and **fails on any
mismatch with the run's**, then byte-compares.

Both halves matter: agreeing bytes with a disagreeing annotation path would mean the rule is
insensitive to it *on this data*, not that the driver is right.

| | result |
|---|---|
| arguments matched (raw_disc, cns list, annotation, rtype) | **4/4 for all three families** |
| L1 `.slct` | **IDENTICAL**, 34,197 rows / 793,811 B |
| SVA `.slct` | **IDENTICAL**, 39,452 rows / 915,815 B |
| Alu | no step `.slct`, correctly (no transduction flank) — a pass, not a skip |
| universe | 41,604 loci from 75,853 claims, 1.82× |

### 5.2 `gate_locus_universe_chr1.sbatch` at `440e846`, job `49547574` — **PASS**

Two arms, one commit. This is the leg the offline gate cannot reach: `--locus-universe` opens the
two-wave split, and a step that opens an input it was never handed dies only on the real chain (job
`49498881` died exactly that way).

| | OFF | ON |
|---|---|---|
| VCF / cns2 / sibling artifacts | — | **9/9 identical** |
| calls (L1 / Alu / SVA / merged) | 1 / 14 / 2 / 17 | **1 / 14 / 2 / 17** |
| driver replica vs step's own `.slct` | — | **2/2 identical** (L1 2,815 rows, SVA 3,205) |
| universe pass cost | — | **3.0s** |
| wall | 6:04 | 5:59 |

### 5.3 `gate_locus_universe_chr1.sbatch` at `7dfae51`, job `49548949` — three arms

OFF / ON / **GATHER**, extended with a third check that reads the **step logs**: the GATHER arm must
show the bypass once and **zero** of its own `--td-sites` calls per flanked family, the ON arm the
reverse. That check exists because both failure modes are silent — both gathering is correct and buys
nothing, neither gathering is empty.

**RESULT: PASS.**

| | OFF | ON | GATHER |
|---|---|---|---|
| universe barrier / driver gather lines | 0 / 0 | 1 / 0 | **1 / 1** |
| vs OFF: VCF / cns2 / sibling | — | 9/9 identical | **9/9 identical** |
| calls (L1 / Alu / SVA / merged) | 1 / 14 / 2 / 17 | 1 / 14 / 2 / 17 | **1 / 14 / 2 / 17** |
| driver replica vs step `.slct` | — | 2/2 identical | 2/2 identical |
| bypass checks (per flanked family, both directions) | — | — | **4/4** |
| wall | 6:02 | 6:02 | 6:02 |

The bypass check is the one worth keeping: `precomputed=1 own-gather=0` in GATHER and the reverse in
ON, per flanked family. Both gathering would be silently correct and buy nothing; neither gathering
would be silently empty.

**Two bugs in this gate script itself, found by reading its own output and fixed:** the phase-timing
block grepped `<sample>_xtea2.log` when `run_xtea`'s `log()` writes to **stdout**, i.e. `driver.log`
(the identical mistake job `49502115` made), and the `[TD-PHASE]` reader used `$3`/`$4` when the
label is `$4` and the seconds `$5`. Both printed *empty* rather than wrong, which is how a timing
section can silently report nothing — the numbers in §1.2 were recovered by hand afterwards.

### 5.4 `gate_locus_universe_wgs.sbatch` — WGS, submitted

Owed even though chr1 passed, because **the last bug in this wave split was invisible on chr1**: job
`49498881` died at `SVA:transduction` on a step opening an input it was never handed, and chr1's
smaller candidate sets never reached it. One arm covers both flags (`--td-gather-driver` implies
`--locus-universe`), and it adds two legs chr1 cannot have: **OFF == frozen BASELINE-0** (isolates
commit drift — without it a clean GATHER-vs-OFF diff could be two equally-wrong arms) and the WGS
answer to the wall-clock question in §1.2.

## 6. What NOT to redo

Plan §10 is the authority and gained seven entries. The ones most likely to cost time:

- **The universe is computable before TRANSDUCTION, not before `cns_filter`** (§2.1). The two
  docstrings that said otherwise are fixed; do not re-derive from an old copy.
- **`transduction_extnd` is `max(400, cns_extnd)`** and the driver already had it. Do not go looking
  for a missing `iextnd` in `automatic_gnrt_parameters`.
- **Take `REP_TYPES[t]["annotation1"]`, not `["annotation"]`**, for anything reproducing transduction.
- **`--td-gather-driver` returns ~25s of CPU and LOSES ~12s of wall clock at chr1 scale** (§5.3).
  Measured, not argued. Do not describe it as a speedup in either direction until the WGS leg runs.
- **The step does NOT read the driver's cached basic-info file.** `X_BasicInfo` keys that cache on its
  own working folder, so every step writes its own copy; they agree because the sampling is **seeded
  by default** (all seven copies in the 2c gate run are identical). With `XTEA_BASIC_INFO_SEED=none`
  the gather window and the evaluation could come from different estimates.
- **`downstream_unified` and `x_transduction` are a deliberate import cycle, broken one way.**
  `x_transduction` imports the module for one env-var constant; `downstream_unified` defers its
  `from x_transduction import XTransduction` into the one function that needs it. Do not "tidy" either
  to a top-level import — both orders are verified by test.
- **The driver and the step agree on two gather file names with no shared constant**, because the step
  hard-codes them. `tests/test_downstream_unified.py` pins the literals **in the step's source**. If a
  rename breaks that test, change both sides.
- Everything in the previous handoff's §5 still stands, in particular the half-open `(-100, +100]`
  proximity window, `is_matched_rep_type` having no Alu arm, the transduction/sibling steps ignoring
  `--nd`, and that L1 precision must always be quoted with its stratum (74.5% overall, 92.7%
  excluding orphan sites).

## 7. NEXT

1. **The per-locus decision loop** — the third box, and now the only one left. §4 has its record type
   and the two things still unknown about its tail. This is the front edge of Stage 3, so its gate
   stops being agreement and becomes the truth set.
2. **Flip `--locus-universe` / `--td-gather-driver` on** when the loop consumes them; both are proven
   answer-preserving and buy nothing standing alone, so flipping them earlier costs a barrier for
   nothing (the same reasoning that keeps `--family-assign-early` off).
3. **`--emit-per-type`**, still deliberately last: flipping per-type VCFs to opt-in breaks every gate
   script that `cmp`s them.
4. **Stage 3**, whose gate is fully specified (`HANDOFF_2026-08-02_stage2c_barrier.md` §1):
   position-level per family plus the orphan / sibling / no-`TD_SRC` strata, no transduction pass/fail.

## 8. Where things live

| | |
|---|---|
| **the agenda** | `docs/handoff/DOWNSTREAM_OVERHAUL_PLAN.md` (§3 diagram + two correction blocks, §4 Stage 2c's last two blocks, §8, §10) |
| the loop's front half | `xtea/downstream_unified.py` — `enumerate_family_sites`, `LocusUniverse`, `transduction_extnd`, `write_td_manifest`, `td_gather_paths` |
| its driver entry points | `xtea/run_xtea.py` — `run_locus_universe`, `run_td_gather_driver`, `downstream_waves` |
| the step-side bypass | `xtea/x_transduction.py` — `call_candidate_transduction_v3`'s `b_td_precomputed` |
| unit gate | `tests/test_downstream_unified.py` (17 tests; `test_driver_enumeration_reproduces_the_step` is the load-bearing one) |
| offline WGS gate | `docs/handoff/runs/gate_locus_universe.py` |
| end-to-end gate | `docs/handoff/runs/gate_locus_universe_chr1.sbatch` (3 arms) |
| the site-creation rule, one copy | `xtea/site_universe.py` + `tests/test_site_universe.py` + `gate_locus_universe.py`'s sibling `gate_site_universe.py` |
| the multi-entry gather | `src/csrc/td_gather.hpp`; CLI `xtea_scan --td-manifest` |
| the run both offline gates read | `/n/scratch/users/c/cos689/xtea2_testdata/stage2c_family_assign_49502115/OFF/HG002.wgs` |
| worktrees for this session | `/n/scratch/users/c/cos689/xtea2_wt/locus_universe_440e846`, `…_7dfae51` |
| the family rule, one copy | `xtea/locus_table.py` |
| BASELINE-0 | `/n/scratch/users/c/cos689/xtea2_testdata/baseline0_hg002_L1-Alu-SVA/` |
| truth sets | `/n/data1/hms/dbmi/park/corinne/claude_sandbox/test_data/HG002/truth/` |
