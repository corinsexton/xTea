# xTea2 handoff — 2026-08-02 (one BAM decode for many families; the Stage 3 blocker is sized)

**This is the current entry point for state. It is NOT the entry point for agenda — the plan is.**

> **Read `docs/handoff/DOWNSTREAM_OVERHAUL_PLAN.md` first.** It still outranks this chain on agenda.
> Everything factual below is written into the plan at the section that owns it; this is the index.
> Supersedes `HANDOFF_2026-08-02_stage2c_barrier.md`, whose §4 NEXT items are what this session
> addressed: **item 1 is decided and closed** (§1), **2b's C++ half is built and gated** (§2),
> and the per-locus loop's blocker is now measured rather than named (§3).

**One-line status: `--family-assign-early` is decided — it stays OFF until unified site creation
makes it load-bearing; `td_gather.hpp` now serves several families from ONE decode of the BAM
(`xtea_scan --td-manifest`), with the single-family form reduced to a one-entry call into it so
byte-identity is structural; and the Stage 3 dependency §2.4 named without sizing now has a number —
52 of 1,515 final per-family records are created after `cns_filter`, all by transduction, which puts
a hard ceiling of 7 of 22 on any family decision made before transduction. That 7 is the previous
session's gate number, re-derived from the opposite direction.**

## 0. State at handoff

| | |
|---|---|
| branch / HEAD | `xTea2` @ this session's commits on top of `e161c38` |
| working tree | clean; all work committed |
| tests | **162 passed, 1 skipped** with `XTEA2_SCAN_BIN` set, excluding `test_downstream_from_cns.py` (the slow chr1 golden). Was 155 — +7 in `test_td_gather.py`. The skip is `test_e2e_concordance` (needs `XTEA2_OUTPUT_DIR`) |
| in flight | see §2.3 — the chr1 agreement gate for the gather change |
| the agenda | the per-locus loop (§3 says why it is Stage 3's front edge), then Stage 3 |

## 1. NEXT item 1 is CLOSED: `--family-assign-early` stays OFF

Not "still unproven" — **proven answer-preserving and still not worth turning on.** It is a *second*
decision point, not a relocated one: the barrier does not replace the postF_hc merge, both run. The
only case for flipping it is that it is a prerequisite for deleting the merge's resolution, so it
flips **in the same change that makes it load-bearing** (unified site creation), not before.

Recorded at the flag itself in `run_xtea.py` and in plan §4 Stage 2c, so it cannot be re-litigated
from the help text alone. Two places still asserted the simulation's disproven prediction — *"4 of
1493 loci move, all toward truth"* — and both now carry the measured zero **and the reason**, which
is the part that generalises.

## 2. 2b's C++ half: several families out of one BAM decode

### 2.1 What was built

`td_gather.hpp` took one `TdSites`; it now takes a `std::vector<TdGatherEntry>`, builds the fetch
plan over the **union** of the entries' windows at the **largest** entry's `extnd`, and fans each
entry out of the shared collection in its own site-file order. The CLI mirrors
`--cns-emit-manifest`: `xtea_scan --td-manifest F`, one line per family,
`tag sites extnd clip_fastq disc_fa [disc_pos]`, `-` for an output the entry does not want.

**`run_td_gather` — the single-family form `XClipDiscFilter.native_td_gather` still calls — is now a
one-entry vector into `run_td_gather_multi` and nothing else.** There is no separate single-entry
path left to drift, which is what makes "the multi-entry version cannot change the one-entry output"
a property rather than a claim.

### 2.2 Three invariants that make sharing invisible, each with a test

- **Collection is a superset.** The union index at the widest `extnd` holds reads a narrower entry
  must not emit. That is safe *only* because of the rule already in the header — **holding is
  permissive, emission is exact** — whose exact half is the fan's two `ref_end <= wbeg` tests.
  `test_entries_may_have_different_extnd` fails if either side of that is lost.
- **Order survives the wider plan.** A fetch region owns reads starting at or after the previous
  region's end, so `ev` comes out ascending for *any* merge boundaries; the `stable_sort` is a no-op
  and equal-start reads keep the BAM order **bwa's tie-break rides on**.
- **The mate gather** merges its region plan over the union of wanted mate positions — different
  boundaries, same ownership rule. Each stripe now writes one part file **per entry**, and a mate's
  sequence is decoded once however many entries claim it.

The manifest is **stricter than the cns one on purpose**: a short line is a hard error, not a silent
skip. A dropped field there means a family quietly does not run, which on a WGS gate reads as "that
family found nothing". Two entries sharing an output path are refused up front.

### 2.3 How it is gated

Unit: `tests/test_td_gather.py` (+7). Every one is the same shape — run N entries together, run each
alone, `cmp` all three products. Transitivity does the rest: the pre-existing tests already pin
*solo native == the Python oracle* on real chr1 data, so *multi == solo* makes *multi == the oracle*
without a second oracle run.

End-to-end: `docs/handoff/runs/refactor_gate_chr1.sbatch` — **job `49516538`**, HG002 chr1, both
gather paths, against the Stage 1 gate's arms (`stage1_gate_chr1_49453032/{A,B}`). RESULT: see §5.

> **The validity condition that script cannot check itself.** It requires
> `git diff <ref_commit>..HEAD -- xtea/ src/` to be nothing but the change under test, and that is
> **not literally true** here — 2a, the conflict sidecar and 2c all landed in between. It is true in
> the sense that matters: each of those was itself gated byte-neutral (2a's 9/9 chr1, 2c's 4/4 vs
> BASELINE-0 on the OFF arm), so a diff in this run points at this change. Stated rather than
> glossed, because the next person to reuse that reference run inherits the same caveat.

### 2.4 No production caller, deliberately

Hoisting the gather to the driver is the throwaway half plan §4 2b describes: the driver would have
to compute the step's `iextnd` and run `re_slct_with_clip_raw_disc_sites` itself, re-implementing
step logic the per-locus loop deletes. `--td-manifest` is what that loop will call.

**Anyone measuring this today will find nothing, and that is expected.** L1's and SVA's transduction
steps are already separate concurrent subprocesses, so their two native gathers (53.3s and 56.8s at
WGS scale) overlap; one shared gather is still ~57s of wall time. The saving is CPU, and it only
becomes wall time when the downstream is one job at full width.

## 3. The Stage 3 blocker is no longer just named — it is 7 of 22

Plan §2.4 said "Stage 2c cannot end with one decision point until site creation is unified" and left
it there. `docs/handoff/runs/probe_site_creation.py` sizes it offline from the 2c gate's kept
intermediates — no BAM, no run. Where each family's final `postF_hc` record **first exists** (±50 bp,
one-to-one, because the chain moves coordinates and exact matching would invent movement):

| | postF_hc | cns | cns2 (transduction) | sibling |
|---|---|---|---|---|
| L1 | 200 | 150 | **50** | 0 |
| Alu | 1230 | 1230 | 0 | 0 |
| SVA | 85 | 83 | **2** | 0 |

**52 of 1,515 per-family final records do not exist when the barrier runs** — all L1 or SVA, all from
the **transduction** step, **none from sibling**. And on the conflict population: **exactly 7 of the
22 have both sides present at cns** (12 loser-invented — the `Alu over L1 tprt_override` block — plus
2 winner-invented and 1 both).

**That 7 is the previous session's gate number, derived from the opposite direction.** The gate
counted which conflicts the barrier *actually* resolved (7 both, 0 disagreements); this counts which
ones it *could have*. Two independent derivations of the same ceiling, which is what makes it a fact
about **site creation** rather than about the rule: widening `MERGE_WINDOW`, moving the barrier
earlier, or reordering the priority cannot raise it above 7 of 22.

## 4. NEXT

1. **The per-locus loop, and read it as Stage 3's front edge rather than 2c leftovers.** Plan §8
   lists it under *Rewrite* next to `downstream_unified.py` / `transduction_evidence.py`, which are
   Stage 3 modules; §3 above adds unified site creation to its dependencies. Sequencing it as a
   Stage-2-sized change left over is the reading that will cost a week. It has two callers waiting
   for it: `--td-manifest` (§2) and the collapse of the two decision points (§3).
2. **`--emit-per-type`**, still deliberately last: flipping per-type VCFs to opt-in breaks every gate
   script that `cmp`s them.
3. **Stage 3**, whose gate is fully specified (previous handoff §1): position-level per family plus
   the orphan / sibling / no-`TD_SRC` strata, no transduction pass/fail.

## 5. The chr1 gate — job `49516538`

RESULT PENDING AT WRITE TIME. See §6 for how to read it.

## 6. What NOT to redo

Plan §10 is the authority and gained four entries this session. The ones most likely to cost a day:

- **"Holding is permissive, emission is exact" is now load-bearing, not a nicety.** Tightening the
  hold test in `collect_stripe`, or loosening the fan's `ref_end <= wbeg` tests, corrupts multi-entry
  output while leaving the single-entry path correct — the exact failure a one-entry test cannot see.
- **A shared decode cannot reorder the gather because of the region ownership rule, not the sort.**
  Do not "simplify" the `stable_sort` into a sort by anything but `ref_start`.
- **Sharing the td gather buys CPU, not wall clock, until the downstream is one job** (§2.4).
- **52 of 1,515 final records are created after `cns_filter`, all by transduction, none by sibling**
  — hence the 7-of-22 ceiling (§3). Do not re-derive it by re-running the pipeline; the probe reads
  intermediates.
- **`--family-assign-early` is decided, not pending** (§1). The "4 moved loci, all toward truth"
  prediction is dead in the tree as well as in the docs.
- Everything in the previous handoff's §5 still stands, in particular that `scan_clip_realign.tsv`
  carries no per-family signal, that `two_side_tprt_both` is written by `cns_filter`, that Stage 2c
  is not a speedup, that the barrier and the merge BOTH run, and that L1 precision must always be
  quoted with its stratum (74.5% overall, 92.7% excluding orphan sites).

## 7. Where things live

| | |
|---|---|
| **the agenda** | `docs/handoff/DOWNSTREAM_OVERHAUL_PLAN.md` (§4 Stage 2c and §10 updated) |
| the multi-entry gather | `src/csrc/td_gather.hpp` — `TdGatherEntry`, `union_sites`, `fan_entry`, `run_td_gather_multi`; CLI in `src/csrc/xtea_scan.cpp` (`--td-manifest`) |
| its tests | `tests/test_td_gather.py` (13 total, +7 this session) |
| the site-creation probe | `docs/handoff/runs/probe_site_creation.py` — no BAM; reads a run's intermediates |
| the family rule, one copy | `xtea/locus_table.py` — `load_rows`, `assign_rows`, `assign_and_filter` |
| the barrier + waves | `xtea/run_xtea.py` — `downstream_waves`, `run_family_assign_barrier` |
| the chr1 agreement gate | `docs/handoff/runs/refactor_gate_chr1.sbatch` (job `49516538`) |
| the 2c gate's arms | `/n/scratch/users/c/cos689/xtea2_testdata/stage2c_family_assign_49502115/{OFF,ON}` |
| BASELINE-0 | `/n/scratch/users/c/cos689/xtea2_testdata/baseline0_hg002_L1-Alu-SVA/` (intermediates kept) |
| truth sets | `/n/data1/hms/dbmi/park/corinne/claude_sandbox/test_data/HG002/truth/` |
