# xTea2 handoff — 2026-08-04 (the decision is a record and named gates; the reasons are data)

**This is the current entry point for state. It is NOT the entry point for agenda — the plan is.**

> **Read `docs/handoff/DOWNSTREAM_OVERHAUL_PLAN.md` first.** It still outranks this chain on agenda.
> Everything factual below is written into the plan at the section that owns it (§4 Stage 3's first
> block, §8, §10); this is the index.
> Supersedes `HANDOFF_2026-08-03_staged_barrier.md`, whose §5 NEXT item **1** — the per-locus decision
> — is what this session started. Its items 2, 3, 4 and 5 are untouched.

**One-line status: `xtea/transduction_evidence.py` exists, the transduction step's whole decision
cascade runs through it, and every drop now carries a named reason in a joinable file.** This is plan
§5's first sentence built as a *structural* change: the thresholds did not move, only their location,
so it is agreement-gated rather than truth-set gated. The truth-set half — changing the numbers §5
flags — is what this unblocks.

**And the shape §5 drew is wrong in one way that matters: `decide()` is five gate groups, not one
function, because two of the five sit behind batched BAM fetches over sets the earlier groups
define.** That is a property of the data flow, not a compromise, and §4 Stage 3's first block has
the diagram.

## 0. State at handoff

| | |
|---|---|
| branch | `xTea2` |
| the feature commits | `e2da5bd` (the module + the rewiring) and `48ec4c7` (a wiring bug and the tests that close its class). **No HEAD hash recorded on purpose** — an earlier handoff named one and an `--amend` invalidated it the same day. |
| working tree | clean |
| tests | **24 in `tests/test_transduction_evidence.py`**, whole fast suite green. |
| in flight | see §4 |
| the agenda | plan §5's *semantic* half: the four "interrogate" thresholds, now that each is a named parameter with a reason attached. |

## 1. What was built

### 1.1 `xtea/transduction_evidence.py` — thresholds as parameters, drops as data

Five gate groups replace the same decisions spread across six methods of `x_transduction.py`, all six
now deleted:

| gate group | replaces | thresholds it owns |
|---|---|---|
| **[A] clip** | `filter_by_clip_algnmt` | `n_clip_cutoff`, `n_polya_cutoff` |
| **[B] source** | `_chk_decoy_cns_cnt`, `_pick_trsdct_rescued`, `_slct_orphan_transduction` | decoy cutoff, `n_half_disc`, the 0.45 dominance ratio, the 0.5 uniqueness ratio |
| **[C] candidate** | `call_out_candidate_transduct` | none — three memberships and a cluster |
| **[D] sanity** | `filter_transduction`'s body | `i_max_cov`, `n_half_clip`, the 5% repeat divergence |
| **[E] genotype** | `add_new_rcd_with_cutoff`'s one `continue` | `ave_cov * max_indel_ratio` |

`TransductionParams` computes every *derived* cutoff with the step's own expressions, so there is one
copy of each: `n_polya_cutoff = n_clip_cutoff/2 - 1` (1.0 at `--nclip 4`), `n_half_disc =
n_disc_cutoff/2` with the `<= 1` fallback (2.5 at `--nd 5`), the decoy cutoff as `n_half_disc*2`, and
`n_half_clip = n_clip_cutoff/2` — a *different* quantity from `n_polya_cutoff`, which is easy to miss
when they are inline. **The floats are not rounded**, because `>= 2.5` and `>= 3` admit different
loci.

**A threshold a caller does not have is `None`, not 0.** Two of the three constructors are partial,
because `filter_transduction` and `add_new_rcd_with_cutoff` are shared with the sibling step
(`x_orphan_transduction`), which reaches those gates from its own per-site parse and never runs
[A]-[C]. A gate reached with a None threshold raises `TypeError` instead of silently admitting or
rejecting everything, which is what a 0 default would have done.

### 1.2 Why it is five groups, and why that is not fixable

```
clip SAM parse -> [A] ----.
disc SAM parse -> [B] --> [C]
                           |   ReadDepth over the SURVIVING candidates only
                           v  [D]
                           |   the genotyper over the SURVIVING sites only
                           v  [E]
```

`filter_transduction`'s `calc_coverage_of_two_regions` and `_collect_gntp_features`' genotyper each
open the BAM for the set handed to them. Flattening this into one `decide()` means computing local
coverage and genotype features for the whole `.slct` list — **34-40K loci per family instead of a few
thousand**. So "one decision function" means **one place per threshold**, not one call. Anyone
designing the per-locus loop off §5's picture should read this first; it is in plan §10 as a
do-not-re-derive.

### 1.3 The reason ledger — what this buys today

`<out>.decision_ledger.tsv` per flanked family: `chrom, pos, phase, reason, detail`, one row per
dropped locus, drops only. 14 reasons, five phases. §5: "`reason` is not decoration — it is what
makes the benchmark interpretable, because every missed truth call gets attributed to the specific
gate that dropped it." **BASELINE-0's 322 false negatives have no such attribution today**, and the
information existed only as `print()` calls in a step's stdout (nothing greps them — checked).

`benchmark/ledger_attribution.py` is the join: truth set → false negatives → the gate that determined
each one. **Read its docstring before quoting it**, because two honest caveats decide how the table
reads:

- **The ledger covers the transduction step only, and that step is not where most calls are made.**
  `cns_filter` produces the candidate list; transduction *adds* the 52-of-1,515 created records and
  *updates* the rest. An FN with **no ledger row never reached these gates** — it is reported as its
  own line (`no_ledger_row`), never folded into "unattributed", because the two statements are
  completely different. That population is expected to dominate.
- **A ledger row near a true positive is not a contradiction**: the step can reject its own candidate
  at a position `cns_filter` already called. Printed as `tp_with_ledger_row` so the FN table is not
  read as causal on its own.
- **Alu can never be attributed**: no transduction flank, so the step never runs for it. Biology, not
  a gap.

The ledger is diagnostic — nothing in the pipeline reads it. **The 14 reason strings are the only
public contract in the module**; rename one and a join against an older run's ledger silently misses,
which is why `ledger_attribution.py` validates every row against the module's own vocabulary and
raises rather than counting.

## 2. The two defects preserved on purpose, and the one comment fixed

`DEFECT-PRESERVED` in the module, and a test fails if either is "fixed":

1. **The cluster test sits inside the per-read loop** (`dominant_side_cluster`). It runs once per read
   instead of once per locus, and only the last read's evaluation survives for the side it picked
   while the other side keeps whatever an earlier iteration left there — so a locus can carry
   `b_l_trsdct` from read 3 next to boundaries from read 9, and **"both sides true" is reachable only
   through that stale flag**. `collect_ins_info` reads exactly that pair to decide which side's cns
   interval to blank, so hoisting the test out of the loop **moves calls**. Stage 3 with a truth-set
   measurement behind it, not a cleanup.
2. **`_pick_trsdct_rescued`'s second return value was both discarded and wrong** — its accumulation
   sat inside `if ins_chrm not in m_n_source_picked`, so it recorded only the first locus per
   chromosome. Not rebuilt; the ledger's `no_dominant_source` rows are that population, correctly.

**And a documentation defect worth knowing, because both windows are gated on.**
`filter_transduction` passes `focal_win = LOCAL_COV_WIN` (**200**) and `focal_win2 = COV_ISD_CHK_WIN`
(**900**) and its comments described them the other way round. So `i_lcov1/i_rcov1` are the 200 bp
window, `i_lcov2/i_rcov2` the 900 bp one, all four are tested against `i_max_cov`, and **cov2 is what
`collect_ins_info` reports as the site's local coverage**. Only the comment changed — swapping the
constants would move calls.

## 3. The unit gate, and the bug it did NOT catch

`tests/test_transduction_evidence.py`, 24 tests. Its centre is six `_legacy_*` functions that are
**verbatim copies of the deleted bodies**, driven against the library over randomized inputs — the
same pattern `test_site_universe.py` uses, and for the same reason: the safety argument is "moved, not
rewritten", and that has to be checked rather than asserted.

**Every differential compares `_ordered()` — nested (key, value) lists — not dicts.** Key order is
part of the answer: the step writes its output files by iterating these maps, and job `49596361`
failed on exactly that class of bug in the locus-universe table. Three mutations confirm the
differentials bite: a strict clip cutoff, hoisting the cluster test out of its per-read loop, and a
`sorted()` that changes **only** the order.

### AND THEN AN END-TO-END RUN FOUND A BUG ALL 21 OF THEM WERE BLIND TO

`filter_transduction` was left returning **three** values where both callers unpack two. Job
`49633492` (the NA12878 golden gate) died with `ValueError: too many values to unpack` inside
`SVA:transduction`; the fast suite was green.

**Why: the differentials test the LIBRARY, and nothing in the fast suite called the two methods the
sibling step shares.** Every differential drives `transduction_evidence.apply_site_filters` directly,
so a broken return arity in `x_transduction.filter_transduction` is invisible to all of them. The gap
was the *wiring*, and only a real chain exercised it.

`48ec4c7` closes the class, not just the instance: three tests now run both shared methods **in
process** with the BAM work stubbed — a fake `ReadDepth` returning fixed coverage, a fake annotation
with no repeat hits — asserting the return arity its callers unpack, that the depth map still carries
cov2 for `collect_ins_info`, that the coverage gate fires with the right reason, and that
`add_new_rcd_with_cutoff` still writes one line per surviving site with the large-indel rejection
intact. Reverting the fix fails two of them.

**Worth generalising: when a refactor moves logic OUT of a method, a differential on the new library
proves the logic and says nothing about the seam.** Test the seam separately, in process, with the
I/O stubbed — it is cheap, and the alternative is finding it 10 minutes into a cluster job.

## 4. The gates

### 4.1 `docs/handoff/runs/gate_transduction_evidence_chr1.sbatch` — the strongest agreement gate here

Two arms, two commits, one command line: BEFORE (the parent, six in-place cascade methods) and AFTER.
Because this change should move *nothing*, the gate is **every file both arms produce must be
byte-identical** — ~98 of them, not just the VCFs.

That threshold is available because of a cheap dry run over an earlier gate's two arms (job
`49599837`): **98 files under the sample directory are already byte-identical between two arms of the
same commit**, and the only two that are not are characterised —

- `cns_manifest.txt` — holds each arm's own absolute output paths; compared with the arm prefix
  stripped, so it is checked rather than skipped.
- `scan_disc_pos.tsv` — a **permutation** only (`sort` makes the two identical); the multithreaded
  scan writes rows in completion order. Compared sorted, which is stronger than excluding it.

Checks: 0 the arms are the code they claim (the module present on one side, absent on the other, and
`[TD-LEDGER]` appearing exactly twice in AFTER and never in BEFORE); 1 the full-tree byte compare,
with a floor of 90 files so a run that produced almost nothing cannot pass; 2 the only new files are
the two ledgers; 3 the ledgers are non-empty, header-correct, and every reason is in the module's own
vocabulary with no unknown phase — **a refactor that produced an EMPTY ledger passes checks 1 and 2 in
silence**, so that one has its own assertions; 4 call counts for the record.

**RESULT at `48ec4c7`, job `49634289`: PASS.** (Job `49633211` was cancelled after the golden gate
exposed the return-arity bug — its AFTER arm was going to fail the same way.)

| check | result |
|---|---|
| 0. the arms are the code they claim | module BEFORE=n AFTER=y; `[TD-LEDGER]` lines 0 and **2** |
| 1. every file both arms produce | **98 identical, 0 differing, 0 missing**, 2 characterised (`cns_manifest.txt` identical after stripping the arm prefix; `scan_disc_pos.tsv` identical sorted, 302,912 rows) |
| 2. new files in AFTER | **2 ledgers, 0 unexpected** |
| 3. the ledgers | L1 **1,690 rows**, SVA **2,281 rows**, header correct, every reason in the 14-word vocabulary, no unknown phase |
| 4. call counts | **1 / 14 / 2 / 17 in both arms** |

**This also settles the golden gate's four L1 notes**: AFTER is byte-identical to BEFORE, so nothing
this change did could have produced them.

#### AND IT EXPOSED WHAT chr1 CANNOT GATE — read this before trusting it alone

The ledgers carry **only phases A and B**. That is not a wiring hole; it is the slice. Measured on the
AFTER arm:

| | L1 | SVA |
|---|---|---|
| `.slct` rows (what gate A sees, via the clip parse) | 2,815 | 3,205 |
| clip-selected — gate A survivors | 217 | 270 |
| disc sources picked — gate B survivors | 5 | 0 |
| cns-consistent — **gate C's domain** | **0** | **0** |
| candidates / after sanity / new records | 0 / 0 / 0 | 0 / 0 / 0 |

**Gate C's domain is empty at chr1, so gates C, D and E never execute and the transduction step writes
zero records.** chr1 therefore gates [A] and [B] end to end (3,971 recorded drops between them) and
proves the step does not crash — and nothing more. [D] and [E] are covered in process by §3's wiring
tests; **[C] and [D] are exercised end to end only at WGS**, where the step writes the ~50 created
records. So the WGS leg is not a magnitude check that could be skipped, it is the only end-to-end
exercise of half the gates. If a future change touches [C]/[D]/[E], a passing chr1 gate says nothing
about it.

### 4.2 `docs/handoff/runs/gate_transduction_evidence_wgs.sbatch` — one arm, against frozen BASELINE-0

**One arm, not three, and the reason is worth keeping.** The barrier gates ran three arms because they
compared *topologies* and needed a same-node control for timings. This change's only correctness claim
is "no artifact moves", which §4.1 already establishes over ~98 files. What chr1 cannot supply is the
**population**: 17 calls and no created records, against WGS's 188/1223/82 (merged 1493) including the
52 that exist only through site creation. BASELINE-0 is frozen with sha256s and is the canonical
statement of those numbers, so one arm against it answers the population question at half the compute.

Checks: 4/4 BASELINE-0 VCFs byte-identical; the ledgers at WGS population (counts printed per phase
and per reason, for the record); calls exactly 188/1223/82/1493; and `downstream_parallel` inside a
**loose** bound (≤700s) — not a measurement, since one arm on one node has no control, but a guard
against the decision loop having become structurally expensive (an evidence object per locus over
~40K loci per family). Submitted `--dependency=afterok` on the chr1 job.

**RESULT at `48ec4c7`, job `49634290`: _see §4.4_.**

### 4.3 The NA12878 golden downstream gate — a second sample and a frozen external baseline

`tests/test_downstream_from_cns.py` at the refactor commit, native gather forced on. Independent of
the HG002 A/B in both sample and baseline: it stages the golden's own post-cns outputs and byte-
compares every downstream intermediate, gating on VCF columns 1-8. **This is the gate that found the
return-arity bug**, ~10 minutes in, which is the argument for running it first on any downstream
change: it is the cheapest end-to-end signal in the repo.

**RESULT at `48ec4c7`, job `49634291`: PASS**, all three families, gate contract (VCF cols 1-8) held.

Four L1 intermediates are printed as "moved (not gated)" and they are the **pre-existing B8
carry-forward** the test's own docstring names, not this change: all four are downstream of
`…clip_sites_raw_disc.txt.slct` being 3,272 rows where the stale golden has 3,302, a file this change
never touches (`.slct`, `.slct2`, `.slct2.tmp`, and `unique_trsdct_half_clip` which is 256 vs 260
because it has 30 fewer input sites). §4.1's byte-identity against the parent commit is the independent
confirmation. Alu's three genotype flips are the documented DF21-golden-vs-ONNX class.

### 4.4 Results

| gate | job | result |
|---|---|---|
| NA12878 golden downstream (`test_downstream_from_cns.py`) | `49634291` | **PASS** — 3/3 families |
| HG002 chr1, before vs after, full tree | `49634289` | **PASS** — 98/98 byte-identical, calls 1/14/2/17 |
| HG002 WGS vs frozen BASELINE-0 | `49634290` | _queued behind the chr1 job; see below_ |

**The fast suite** is green at 24 tests in this file, whole-suite green.

### 4.5 The attribution tool, run for the first time (chr1, job `49634289`'s AFTER arm)

`benchmark/ledger_attribution.py --chrom chr1` against the HG002 truth set. It works on real ledgers,
and the answer is a clean negative:

| | truth (chr1) | calls | TP | FN | FN attributed to a gate |
|---|---|---|---|---|---|
| L1 | 9 | 1 | 1 | 8 | **0** — all 8 `no_ledger_row` |
| Alu | 111 | 14 | 14 | 97 | **0** — no ledger exists (no flank) |
| SVA | 9 | 2 | 2 | 7 | **0** — all 7 `no_ledger_row` |

**3,971 loci were dropped by gates [A] and [B] at chr1, and not one of them sits within 200 bp of a
missed truth call.** Read it with both caveats: gates C/D/E did not execute here at all (§4.1), and
Alu can never be attributed. But as far as chr1 can say anything, it says the misses are **upstream of
the transduction decision** — which is the hypothesis the WGS run exists to test at a population where
it can actually fail. `--chrom` was added for exactly this: the truth sets are whole-genome, so
without it every off-contig record counts as a miss.

Do **not** read the chr1 sensitivity here as a quality number — a 3.4 GB chr1 slice is not comparable
to BASELINE-0's WGS 1223-of-1355 Alu.

## 5. What NOT to redo

Plan §10 is the authority and gained a block for this change. Everything in the previous handoff's §4
still stands.

- **The transduction decision cannot be one function.** Two batched BAM fetches interleave it (§1.2).
  Do not re-derive this as a design question.
- **Do not tidy the two preserved defects** (§2). The first moves calls; the second no longer exists.
- **Do not swap `LOCAL_COV_WIN` and `COV_ISD_CHK_WIN`** to match the old comments. The values are
  right and the comments were wrong.
- **Do not rename a ledger reason** without knowing that it is the module's public contract.
- **A differential against a verbatim copy proves the logic and says nothing about the seam** (§3).
- **`scan_disc_pos.tsv` differs between two runs of the same commit as a permutation only**, and
  `cns_manifest.txt` because it embeds output paths. Everything else is byte-identical, which is what
  makes a full-tree comparison usable as a gate.
- **`_extract_transduct_info` is a sixth dead function** (zero callers). Left in place — deleting it
  belongs to the §8 sweep, not to a decision refactor.

## 6. NEXT

1. **The semantic half of Stage 3** — now unblocked and now cheap to reason about, because each
   threshold is a named parameter and each rejection is a counted reason. §5's four "interrogate"
   items, in the order the ledger's own counts justify: `NEARBY_CLIP=50`, the `ileftmost_cns`
   500-vs-5000 asymmetry (dormant today — `family_params`' docstring says why waking it is not a flag
   flip), the two-sided polyA/polyT veto (make conditional on source strength), and the three-way
   decoy tally. **Each is now a one-line change plus a truth-set measurement**, which is the shape
   Stage 3's gate wants.
1b. **One known ledger gap, deliberately left open** (marked at the code, `call_candidates`). A locus
   with clip evidence and disc polyA but **no entry in the picked-source map** appears in no ledger at
   all when the reason is that it had *zero* disc reads realigning to any source — gate B's domain is
   the disc SAM's own map, so it never sees such a locus, and gate C skips it silently. Adding a
   `no_disc_evidence` reason would attribute it, and it is worth doing: "no disc evidence at all" and
   "disc evidence but no dominant source" are different findings for a missed truth call. It was not
   done in this session on purpose — it changes the ledger's contents after the two gates that
   measured them (1,690 / 2,281 rows at chr1), and the WGS numbers were being produced at `48ec4c7`.
   Do it as its own commit and re-run the chr1 gate, which is 14 minutes.

2. **Run `benchmark/ledger_attribution.py` on the WGS gate's output** and read the FN table with its
   caveats. Its first job is to say how much of the 322 is even reachable from these gates; the answer
   may well be "little", and that is a finding rather than a disappointment — it would say the
   sensitivity ceiling is upstream, in `cns_filter`.
3. **The per-locus decision loop's remaining box** — the *batched evidence stage* (the `--td-manifest`
   driver caller, which has no missing input left) is still the piece that would let one process serve
   all three families. The decision is now a library, so the loop's third box is a call rather than a
   rewrite.
4. **The staged-barrier handoff's items 2, 3, 4, 5**, unchanged: reconsider flipping
   `--td-gather-driver` on; the last ~18s (give the driver gather fewer threads while a family is
   still in `cns_filter` — **measure before building**, the +16.3s it is inferred from is inside
   ±18s of per-family variance); `--emit-per-type`; and Stage 3's own gate spec.

## 7. Where things live

| | |
|---|---|
| **the agenda** | `docs/handoff/DOWNSTREAM_OVERHAUL_PLAN.md` (§4 Stage 3's first block, §8, §10) |
| the decision | `xtea/transduction_evidence.py` — `TransductionParams`, `SiteEvidence`, `gate_*`, the five drivers, `DecisionLedger` |
| its unit + wiring gate | `tests/test_transduction_evidence.py` |
| its end-to-end gates | `docs/handoff/runs/gate_transduction_evidence_chr1.sbatch`, `…_wgs.sbatch`, and `tests/test_downstream_from_cns.py` |
| the FN attribution | `benchmark/ledger_attribution.py` (joins the ledger to `benchmark/benchmark_truth_set.py`'s misses) |
| the step | `xtea/x_transduction.py` — `call_candidate_transduction_v3`, `filter_transduction`, `add_new_rcd_with_cutoff`, `parse_transduction_disc_algnmt` |
| worktrees for this session | `/n/scratch/users/c/cos689/xtea2_wt/td_evidence_48ec4c7` (AFTER; its `xtea_scan` is **copied** from the BEFORE worktree, legitimately — `src/` is identical between the two commits, so one binary removes the C++ as a variable), BEFORE is `…/staged_barrier_d565445` |
| everything else | see `HANDOFF_2026-08-03_staged_barrier.md` §6 |
