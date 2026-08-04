# xTea2 handoff — 2026-08-04, later (the last big rejection has a reason; and SVA's decoy zero is a bug)

**This is the current entry point for state. It is NOT the entry point for agenda — the plan is.**

> **Read `docs/handoff/DOWNSTREAM_OVERHAUL_PLAN.md` first.** It still outranks this chain on agenda.
> Everything factual below is written into the plan at the section that owns it (§4 Stage 3's second
> block, §5, §8, §10); this is the index.
> Supersedes `HANDOFF_2026-08-04_transduction_evidence.md`, whose §6 NEXT item **1** — instrumenting
> `check_seprt_disc_consistency` — is what this session did. Its items 2, 3 and 4 stand, item 2 with a
> correction below.

**Two things happened, and the second was not on the agenda.**

1. **Gate `[B2]` exists: `check_seprt_disc_consistency`'s decision is now a named gate with four
   reasons in the ledger.** It drops 135 of L1's 159 picked sources and 49 of SVA's 105 at WGS — the
   step's largest threshold-driven rejection after gate A — and it was the last big population the
   attribution was blind to.
2. **The agenda said "look at `decoy_cross_family` before touching it" (144 for L1, 0 for SVA). Looked.
   The zero is not biology: the SVA arm of that gate cannot fire, and behind the same line of code,
   *every* reference SVA source flank is misclassified as a decoy and skipped.** Measured: **0 of
   SVA's 105 picked sources are reference copies; all 105 are polymorphic.** Not fixed — it is
   call-moving and BASELINE-0 will move.

## 0. State at handoff

| | |
|---|---|
| branch | `xTea2` |
| the commits that touch code | `ff7cde1` gate [B2] — the module's four reasons, `PHASE_ORDER`/`phase_rank`, the wiring in `x_clip_disc_filter.py` and `x_transduction.py`, and `ledger_attribution.py`'s switch off the lexical phase sort. Everything after it is docs, gates, tests and one probe. **No HEAD hash recorded on purpose** — an earlier handoff named one and an `--amend` invalidated it the same day. |
| working tree | clean |
| tests | **38 in `tests/test_transduction_evidence.py`** (was 25), +1 in `tests/test_benchmark_truth_set.py`; whole fast suite green at 232. |
| in flight | **nothing — all three gates have reported** (§4). |
| the agenda | the SVA decoy fix (§3) is now the biggest known *sensitivity* item in the transduction step and needs its own arm. Then plan §5's semantic half, still reprioritised onto gate A. |

## 1. Gate [B2] — what was built, and the one constraint that shaped it

`check_seprt_disc_consistency` is **shared**: `cns_filter` calls it twice (`x_clip_disc_filter.py:1480`,
`:1652`) and the transduction step once (`x_transduction.py:326`). The previous session named that as
the obstacle and it is the reason this is not shaped like the other five gates:

- **its body did not move.** It stays where both callers reach it; only the predicate moved, into
  `transduction_evidence.gate_cns_disc_cluster`.
- **the ledger is caller-supplied.** The transduction caller passes one; `cns_filter`'s two pass
  nothing; a `None` ledger records nothing and the gate's answer is identical either way. That is what
  keeps the instrumentation off a path the decision refactor deliberately did not touch.
- **its differential drives the SHIPPING method, not a library driver** — because here the method *is*
  the seam, and §3 of the previous handoff is the reason: a differential against a verbatim copy proves
  the logic and says nothing about the seam. It needs no BAM; `m_repsnt_pos` and `m_disc` are plain
  dicts from an earlier SAM parse.

### 1.1 Four reasons, and the first is not a threshold

| reason | what it is |
|---|---|
| `no_cns_disc_reads` | **an ABSENCE** — no disc read of this locus survived the third bwa pass to the consensus. Never to be added to the three below. |
| `no_cns_disc_cluster` | neither side's consensus hits form a cluster |
| `low_two_sided_cns_disc` | both sides cluster, total support below the caller's cutoff |
| `low_one_sided_cns_disc` | the `b_transduct=True` arm — **unreachable today**, implemented anyway |

The absence/threshold split is the funnel probe's rule made structural: §4.7 of the previous handoff
insisted "no evidence was produced" and "a threshold rejected it" are different statements. The
set-difference derivation could not tell them apart; the ledger can, and
`test_absence_of_cns_disc_evidence_is_not_a_count_of_zero` pins it.

`low_one_sided_cns_disc` is implemented and switched off (`b_one_sided_cns_min`, False in all three
callers) so that turning the arm on later — a **call-moving** change, and for the transduction caller a
*quarter* of `--nd`, since the caller already halved it — arrives with its own attribution instead of
needing this work done twice.

### 1.2 The cutoff is the CALLER's number, not `n_disc_cutoff`

`TransductionParams.for_cns_cluster(n_cns_cluster_cutoff, b_one_sided_cns_min=False)` stores what it is
given and derives nothing, because the two callers hand it **different numbers**: transduction passes
`n_half_disc` (2.5 at `--nd 5`), `cns_filter` passes the full `ndisc_cutoff`. Storing it under a name
that promised "the disc cutoff" would make every reader of a transduction ledger row halve it twice.
`test_the_cns_cluster_cutoff_is_the_caller_s_own_number` pins that, and that a partial params object
leaves every threshold it does not own as `None` so a wrong gate fails loud.

## 2. AND THE DECISION IS SIX GATE GROUPS, NOT FIVE — a third batched pass the diagram omitted

```
clip SAM parse -> [A] clip ---------------------------------.
disc SAM parse -> [B] source                                |
                       |  (B's SURVIVORS' disc reads         |
                       |   re-selected and realigned to      |
                       |   the CONSENSUS -- a 3rd bwa pass)  |
                       v [B2] cns cluster -----> [C] candidate
                                                     |  (ReadDepth over survivors)
                                                     v [D] sanity
                                                     |  (the genotyper over survivors)
                                                     v [E] genotype
```

`re_select_disc_reads` writes only the disc reads of gate [B]'s survivors, and they are realigned to
the consensus before [B2] can see them. That is the same barrier argument as the two BAM fetches, one
stage earlier — so "one decision function" means one place per threshold, over **six** groups and
**three** batched passes. The previous handoff's §1.2 and plan §10 both said two fetches and five
groups; both are corrected.

## 3. THE FINDING THAT WAS NOT ON THE AGENDA: every reference SVA source is skipped

The agenda said to look at `decoy_cross_family` firing 144 times for L1 and **0** for SVA, which §10
read as evidence *against* §5's SVA/Alu-competition hazard. **That reading was wrong.** The decoy test
in `parse_transduction_disc_algnmt` (`x_transduction.py:1160`) is a **substring** match on the flank
record's name, with `TD_DECOY_SVA = "SVA"`. Two independent consequences, both measured offline from
job `49634290`'s artifacts and the two flank libraries
(`docs/handoff/runs/probe_td_source_classes.py`, seconds, no BAM):

| library | records | called decoy by the substring test | of those, records with a FLANK's structure |
|---|---|---|---|
| `hg38_FL_L1_flanks_3k.fa` | 707 | 45 | **0** |
| `hg38_FL_SVA_flanks_3k.fa` | 5,454 | **5,454** | **5,454 (100%)** |

1. **The SVA flank library has no decoy records at all.** Only the L1 library has the consensus
   concatenated in — `>LINE1`, `>SVA`, `>ALU`, `>HERV-K`, 35 `>ALU_Alu*` and 6 `>SVA_[A-F]`, and all
   45 of its substring matches are genuine consensus records, so the test is exactly right *there*.
   `gate_decoy_cross_family`'s SVA arm reads `n_decoy_alu` and `n_decoy_l1` and **no code path can
   increment either on an SVA run**. So §5's hazard is **untested, not refuted**.
2. **Worse: every reference SVA source flank is misread as the SVA consensus.** The SVA library's
   records are named `chr13~51077120~51078923~SVA_D~0L`, so `"SVA" in hit_flank_id` is True for all
   5,454; line 1160's branch counts the hit into `m_SVA` and **`continue`s before the source is
   extracted**. So SVA transduction can only ever be called from **polymorphic** flanks
   (`~polymerphic~`, written with `start == end`). Confirmed end to end on the same run:

| | picked sources | reference copies | polymorphic |
|---|---|---|---|
| L1 | 159 | **114** | 45 |
| SVA | 105 | **0** | **105** |

A 5,454-record library built to be the reference SVA source set, and every hit to it discarded.

**It arrived from master byte-identical** (`46520ab`: same four constants, same substring test), so
xTea2 did not break it — the same provenance as `b_transduct`'s deadness.

**NOT FIXED HERE, deliberately.** The right test is on the record name's **structure** — a decoy
consensus record has no `~`-delimited coordinate fields, a flank record has at least three — and making
it structural **moves calls and BASELINE-0 with them**. How much it unlocks is **not derivable from
existing artifacts**: gate B's SVA domain is 210 loci today, all reached through polymorphic flanks,
and only a run with the test fixed can say how much larger it becomes. So: its own arm, with the probe
run on both sides of it.

## 4. The gates

| gate | job | result |
|---|---|---|
| NA12878 golden downstream (`gate_golden_downstream.sbatch`) | `49703320` | **PASS** — 3/3 families, 4:18 wall, exactly the documented notes: 4 L1 carry-forward, 3 Alu genotype flips |
| HG002 chr1, before vs after, whole tree + ledgers | `49703321` | **PASS** — 99 identical, 1 ledger reconciled after removing its 5 B2 rows, 0 differing, 0 new files, calls 1/14/2/17 in both arms |
| HG002 WGS vs frozen BASELINE-0, with a prediction | `49703322` | *(see below)* |

### 4.1 `gate_cns_cluster_ledger_chr1.sbatch` — and why its ledger check is stronger than its parent's

The parent gate (`gate_transduction_evidence_chr1.sbatch`) could treat the ledgers as **new files**,
because only its AFTER arm had the module. Here **both arms write ledgers**, so the ledger cannot be
excluded — it is the one artifact this change is *supposed* to alter, by exactly the rows [B2] adds. So
the check is exact:

> **AFTER's ledger with every `B2_cns_cluster` row deleted must be BYTE-IDENTICAL to BEFORE's.**

That asserts both halves in one comparison: rows were added, and nothing else moved — not a row, not a
reason, not an order. (Order is a real property here: [B2]'s rows are appended between gate B's and
gate C's, so anything that reordered the rest fails.) Everything else under the sample directory must
be byte-identical outright, with the same two characterised exceptions the parent gate established over
a dry run of two arms of the *same* commit (job `49599837`): `cns_manifest.txt` embeds each arm's own
output paths, `scan_disc_pos.tsv` is a permutation because the multithreaded scan writes rows in
completion order.

**WHAT chr1 CANNOT SAY ABOUT [B2] — read this before trusting a PASS.** chr1's gate B survivors are 5
(L1) and **0** (SVA), and its cns-consistent set is empty for both. So [B2]'s domain at chr1 is 5 loci
for L1 and **nothing at all for SVA**: this leg proves the wiring, the vocabulary and the file-level
containment, at N=5, on one family. The population exists only at WGS.

**RESULT at `365f352`, job `49703321`: PASS.** Wall 6:03 + 5:55 for the two arms.

| check | result |
|---|---|
| 0. the arms are the code they claim | `gate_cns_disc_cluster` BEFORE=0 AFTER=1; `ledger.record(` 0 and 2; `[TD-LEDGER]` **2 in both arms** |
| 1. every file both arms produce | **99 identical, 1 ledger reconciled** (L1's, identical after removing its 5 B2 rows), 2 characterised, **0 differing, 0 missing** |
| 2. new or removed files | **0** — as expected, both arms already write the ledgers |
| 3. the vocabulary | **18 reasons**, of which exactly **4** are B2's; no B2 reason under another phase and no other reason under B2's |
| 4. the probe's reconciliation | **OK for both families** — the ledger's B2 rows are exactly `picked sources` minus `cluster-consistent`, from different files |
| 5. call counts | **1 / 14 / 2 / 17 in both arms** |

**AND THE FIVE ROWS ARE ALL `no_cns_disc_reads` — the ABSENCE reason, not a threshold.** At chr1 not one
of L1's 5 picked sources had a disc read survive the third bwa pass to the consensus. That is the
absence/threshold split earning its keep on the first run that could test it: **the set-difference
derivation would have reported all 5 as threshold rejections**, because subtracting two site files
cannot see the difference. It is also why chr1 cannot gate [B2]'s thresholds at all — at this scale the
gate never reaches them.

### 4.2 `gate_cns_cluster_ledger_wgs.sbatch` — one arm, and a prediction made before the run

One arm against frozen BASELINE-0, for the reason §4.2 of the previous handoff gives: the correctness
claim is "no call moves" and chr1 establishes that over the whole tree; what chr1 cannot supply is the
population.

**Its check 5 is a prediction, not a measurement.** The funnel probe derived [B2]'s population from two
of the step's own site files — `picked sources` minus `cluster-consistent` — with no ledger in
existence: **135 (L1) and 49 (SVA)** on job `49634290`. The ledger must now say the same numbers, from
a completely different source. If they disagree, the site-file derivation every prior handoff rested on
is what was wrong. The domain arithmetic is sound — `prep_representative_pos` is lossless, one entry per
`(chrm, pos)`, so [B2]'s domain is exactly the picked sources and its drops are 159−24 and 105−56.

Check 6 is the probe's own **set-wise** reconciliation on this arm (not just counts), and check 7 reruns
the FN attribution with [B2] in the vocabulary — reported, not gated, because either answer is
information: [B2]'s domain is only ~160 loci per family, so it may attribute nothing.

### 4.3 The probe now CHECKS what it used to only measure

`probe_decision_funnel.py` compares its set difference against the ledger's B2 rows and **fails loud**
if they disagree. Two derivations from different files agreeing is what makes "the gate is where the
ledger says it is" a measurement rather than a claim. On a run that predates the instrumentation it
prints **"SKIPPED, not passed"** rather than letting an empty phase read as agreement — verified against
job `49634290`'s output.

## 5. The trap this change walked into, and the one it did not

**PHASE ORDER IS SEMANTIC, NOT LEXICAL — and this was live, not hypothetical.**
`"B2_cns_cluster" < "B_source"` is True in ASCII (`'2'` is 0x32, `'_'` is 0x5F), and
`ledger_attribution.attribute` found "the gate that determined this false negative" with `min(rows)` on
the phase **string**. That was correct only while `A < B < C < D < E` held by accident. With B2 present
it would have attributed misses to B2 that gate B rejected first — a wrong attribution that looks
entirely plausible in the output table. `transduction_evidence.PHASE_ORDER` / `phase_rank` make it
explicit; every consumer goes through them; `test_every_reason_has_a_phase` fails if a new phase has no
rank.

**Verified answer-preserving on real output**: the tool reproduces job `49634290`'s numbers exactly
(306 `no_ledger_row` / 8 `low_clip_support` / 8 `no_polya_clip` = 322) after the change.

**Six mutations were run against [B2]'s tests and every one that changes behaviour bites** — the table
is in `tests/test_transduction_evidence.py`'s docstring. The one the differential does *not* catch is
`for_cns_cluster`'s dead-arm default, and that is correct rather than a hole: the caller always passes
`b_transduct` positionally, so the default is unreachable from production.

## 6. What NOT to redo

Plan §10 is the authority and gained two corrections plus one new entry. Everything in the previous
handoff's §5 still stands except where noted.

- **The decision is SIX gate groups behind THREE batched passes**, not five behind two. `[B2]` sits
  behind the disc-to-consensus realignment of gate [B]'s survivors. Do not re-derive this.
- **`[B2]`'s body must stay in `x_clip_disc_filter.py`** — `cns_filter` shares it. Only the predicate
  belongs in `transduction_evidence`, and the ledger must stay **caller-supplied**.
- **Never sort ledger phase names.** Use `phase_rank`. `"B2_cns_cluster" < "B_source"` lexically.
- **`no_cns_disc_reads` is an absence, not a zero count**, and must never be summed with the three
  thresholds beside it.
- **The SVA `decoy_cross_family` zero is a dead arm, not a measurement** (§3). Do not cite it as
  evidence for or against §5's SVA/Alu-competition hazard — that hazard is untested.
- **SVA reference-copy transductions are unreachable today** (§3), and fixing that is call-moving. Do
  not fix it as a tidy-up; it needs an arm and BASELINE-0 will move.
- **The two preserved defects, the two coverage windows, the ledger reasons as public contract, the
  `.slct` two-independent-screens rule, `_extract_transduct_info`'s deadness, and
  `check_seprt_disc_consistency`'s dead `b_transduct`** — all unchanged from the previous handoff's §5.
- **Two of §5's four "interrogate" thresholds remain measured near-inert**: the two-sided polyA/polyT
  veto rejects 1 locus (L1) and 2 (SVA) genome-wide, the large-indel check 1 each.

## 7. NEXT

1. **The SVA decoy fix (§3) — now the biggest known sensitivity item inside the transduction step.**
   Make the decoy test structural (`chrom~start~end~subfamily~rc` has integer coords in fields 2 and 3;
   a consensus record has none) rather than a substring match. **Its own arm**: BASELINE-0 will move,
   and `probe_td_source_classes.py` is the before/after instrument. Two things to settle in the same
   arm: how much larger gate B's SVA domain becomes (210 loci today), and whether the SVA library
   should also get the decoy consensus concatenated in, which is what would make §5's SVA/Alu hazard
   testable at all.
2. **Plan §5's semantic half, still reprioritised onto gate A** (97% of the rejecting). Expect
   *precision*: 95% of BASELINE-0's misses never reach this decision. The `decoy_cross_family` caveat
   the previous handoff attached to this item is now **answered** — see §3, and do not spend the item
   on it.
3. **`NEARBY_CLIP=50`'s distribution is still unmeasured, and it cannot be measured offline.** §5 wants
   it replaced with a measured clip-to-insertion distance distribution; the quantity is
   `abs(anchor_pos - clip_pos)` inside gate [B2]'s evidence loop, and the disc-to-consensus SAM it
   comes from is deleted at the end of the step. So this needs an instrumented run that emits the
   distribution as a diagnostic — gate [B2] is now the natural place to hang it.
4. **The per-locus decision loop's remaining box** — the batched evidence stage (the `--td-manifest`
   driver caller, which has no missing input left). The decision is a library, so the loop's third box
   is a call rather than a rewrite.
5. **The staged-barrier handoff's items 2, 3, 4, 5**, unchanged: reconsider flipping
   `--td-gather-driver` on; the last ~18s (**measure before building** — the +16.3s it is inferred from
   is inside ±18s of per-family variance); `--emit-per-type`; and Stage 3's own gate spec.

## 8. Where things live

| | |
|---|---|
| **the agenda** | `docs/handoff/DOWNSTREAM_OVERHAUL_PLAN.md` (§4 Stage 3's two blocks, §5, §8, §10) |
| the decision | `xtea/transduction_evidence.py` — `TransductionParams`, `SiteEvidence`, `gate_*`, `PHASE_ORDER`/`phase_rank`, `DecisionLedger` |
| gate [B2]'s body | `xtea/x_clip_disc_filter.py:check_seprt_disc_consistency` — shared with `cns_filter`, ledger passed in |
| its unit + seam gate | `tests/test_transduction_evidence.py` (38 tests; the [B2] differential drives the shipping method) |
| its end-to-end gates | `docs/handoff/runs/gate_cns_cluster_ledger_{chr1,wgs}.sbatch`, `gate_golden_downstream.sbatch` |
| the FN attribution | `benchmark/ledger_attribution.py` (`--chrom` for a single-contig run) |
| the funnel reconciliation | `docs/handoff/runs/probe_decision_funnel.py` (now checks B2 rather than inferring it) |
| the SVA finding | `docs/handoff/runs/probe_td_source_classes.py` (libraries + a finished run, seconds) |
| the step | `xtea/x_transduction.py` — `call_candidate_transduction_v3`, `parse_transduction_disc_algnmt` (the decoy test is at `:1160`), `filter_transduction`, `add_new_rcd_with_cutoff` |
| worktrees for this session | `/n/scratch/users/c/cos689/xtea2_wt/td_b2_365f352` (AFTER; `xtea_scan` **copied** from the BEFORE worktree, legitimately — `src/` is byte-identical across the range), BEFORE is `…/td_evidence_48ec4c7` |
| everything else | see `HANDOFF_2026-08-04_transduction_evidence.md` §7 |
