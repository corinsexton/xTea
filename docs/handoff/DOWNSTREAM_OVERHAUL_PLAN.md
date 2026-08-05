# Downstream overhaul — unify the three repeat types, rewrite transduction

> *Provenance: drafted 2026-07-29 as a plan-mode plan file
> (`~/.claude/plans/we-need-to-plan-precious-sketch.md`). **This committed copy is canonical** — edit
> it here, not there, so the reasoning stays version-controlled with the code it governs.*

> **This document is the entry point.** It is written to survive a context clear: every claim that was
> verified is marked as such, every claim that was *not* is flagged, and §10 lists what must not be
> re-derived. Read §1-4 before writing code; read §10 before disbelieving anything here.
>
> **It also outranks the `docs/handoff/` chain on *agenda*.** The handoffs through
> `HANDOFF_2026-07-29_transduction_inputs.md` predate this plan and end on a transduction
> micro-optimization agenda that this plan supersedes (their *measurements* all still stand). If a
> handoff and this plan disagree about what to do next, this plan wins.

**Progress — 2026-07-29/30.** Stage 0 *acquisition* for the germline benchmark is **done**: the HG002
alignment and truth sets are staged and verified at
`/n/data1/hms/dbmi/park/corinne/claude_sandbox/test_data/HG002/` (see its `PROVENANCE.md`, and §7.1 for
the one deviation — the BAM is hs38DH-aligned, not `no_alt`).

**Progress — 2026-08-01. STAGE 0 IS COMPLETE AND ITS GATE PASSES.** Stage 1 is unblocked, and its
scope conditional is resolved: **the read gather is 73-87% of transduction, so the C++ per-site
gather subcommand is the first move.** The one item deliberately left open is `benchmark_truth_set.py`
(Stage 3 prep, design settled in §7.1), and the somatic pair for BASELINE-0 (§6.3 flags its existence
as unverified; the plan puts the somatic gate in Stage 4).

| Stage 0 item | status |
|---|---|
| phase timers (measure 1) | **DONE** `6ce36c7` — 13 phases + total, unconditional, `[TD-PHASE]` to stderr |
| phase breakdown (measure 2) | **DONE** — job `49437990`. **Gather 73-86%**, bwa 12-25%, SAM parse <2%. Stage 1 conditional resolves to the C++ gather; see §4 Stage 0 |
| type-set-dependence fix | **DONE** `de96727` — the BED-plan helper, + `tests/test_mate_in_rep_type_partition.py` |
| `_tprt`'s meaning | **RESOLVED** — §7.1; canonical-mechanism filter, `_all` is the denominator |
| `col5` → `mate_in_rep_clip_cnt` rename | **DONE** — `mate_in_rep_bed_plan`, CLI `--mate-in-rep-families` (+ deprecated alias) |
| BASELINE-0 (both type sets) | **DONE, and the Stage 0 gate PASSES** — jobs `49441918-22`. 3-type L1 188 / Alu 1223 / SVA 82 (merged 1493); L1-only 200; SVA-only 85. Frozen with sha256s under each run's `BASELINE-0/`. See §4 Stage 0 gate result |
| chr1 slices (HG002, SMHT005-3G) | **DONE** — `test_data/HG002/chr1/HG002.chr1.bam` (3.4 GB, 61.2M reads); `test_data/SMHT005-3G/chr1/` (5 CRAMs, 43 GB, `bam_list.txt` written) |
| how transductions are scored | **design settled** — §7.1: derive truth status from the ALT residual (feasibility verified), 5-way error taxonomy, two denominators. `benchmark_truth_set.py` still to implement |

**Progress — 2026-08-01, later. STAGE 1's LEAD ITEM IS BUILT AND ITS GATE PASSES.**
`src/csrc/td_gather.hpp`, a native per-site transduction gather that replaces both Python gathers
and is **byte-identical** to them (clip FASTQ, `all_disc_pos`, disc mate FASTA) on real chr1 data
and in a hermetic multi-contig fixture. The chr1 gate (job `49453032`) is **PASS** on all three
legs: legacy-vs-native agreement 9/9 artifacts, native-vs-native determinism 9/9, and
`test_downstream_from_cns.py` against the NA12878 chr1 golden with the native gather forced on.
**At WGS scale the gather is 15.25× (L1) / 17.82× (SVA)** — 812.0s → 53.3s and 1011.3s → 56.8s,
all three artifacts byte-identical and deterministic (job `49453857`), which puts
`downstream_parallel` at roughly 1.8×. Off by default; `--td-gather-scan` turns it on. See §4
Stage 1 for what was built, the numbers, both gate results, and the one scope item deliberately
deferred to Stage 2 (the manifest).

**Progress — 2026-08-05. THE DECOY DEFAULT IS FLIPPED ON, ITS THREE GATES PASS, AND STAGE 3's GATE NOW
HAS A NO-REGRESSION MODE. Two of §5's open items are closed by measurement, and this plan's own
reasoning is corrected once more — this time by a name collision it never noticed.**

| item | status |
|---|---|
| **`--td-decoy-structural` ON BY DEFAULT** | **DONE** `3d1437a` (+ `f2e0160`), the first change to `xtea/` since `55055ba`. Measured first, flipped second, both scales re-run third — the order §10.1 of the handoff makes an obligation. **Three jobs, all PASS**: no-regression `49791196` (**19 s**), chr1 containment `49792011` (18:19), WGS payoff `49792012` (38:32). WGS reports **SVA reference-copy picked sources 0 → 241** of 346, L1 and Alu byte-identical to BASELINE-0, `L1 = 188` / `Alu = 1223` exact |
| the flip's **cost** | **SVA position precision 93.9% → 90.6%** — three new calls, **all off-truth**, TP 77 and FN 13 unchanged. REPORTED, never gated. This is not explained away by §7.1's "the truth set cannot score transductions as transductions"; that is a reason it is not a *verdict*, not a reason it is free |
| the cost, **ADJUDICATED** `b2a7f41` | **CLOSED, and on a statement about instruments rather than about the caller.** Of the three gained calls: **1 VERIFIED, 2 UNDECIDABLE** — not the "two, both VERIFIED" the handoff had, because the pool's second VERIFIED entry sits at a position BASELINE-0 already called. Three explanations killed: not a windowing near-miss (nearest truth record of *any* family is 250 kb / 1.7 Mb / 2.5 Mb); **chrY:10120435 cannot be on-truth at all** — the SVA truth set has **0 chrY records** of 90, while chrY *is* in scope (Alu has 2), so it is a denominator artifact; and all three are **inter-chromosomal**, so none is the source copy's own reads mis-assigned locally. The two UNDECIDABLEs are SVA_D source windows at **2.7% and 3.8%** divergence, inside the verifier's 5% cutoff. `docs/handoff/runs/probe_sva_flip_offtruth.py` |
| **the flip's UNCOUNTED win** | `chr17:45593489` went `TD_SRC=not_transduction` → `chr17:46252374-46253788~1~orphan`: an **existing** call gained a transduction assertion, **VERIFIED on 9 disc reads**, and it moves **no precision number**. It is also why the verification pool grew by **4** while the call set grew by **3** — a discrepancy carried through three handoffs without being reconciled. Count TD_SRC conversions, not just new positions |
| Stage 3's gate — **no-regression mode** | **DONE** `f10daa7`. `verify_transduction.py --baseline <verdicts.tsv>` gates on **absolute VERIFIED counts per stratum** and refuses, rather than badly performs, two comparisons: a baseline at a different operating point (exit 3), and strata read from columns instead of re-derived. Harness `docs/handoff/runs/td_gate_no_regression.sbatch` refuses a third — see the new §10 entry on arm scope |
| the gate's **primary stratum** | **DONE** `f10daa7` — `family × call class × shape × sidedness`, source class kept beneath with `EMPTY DENOMINATOR` markers. Proven **reporting-only**: the per-entry verdicts are byte-identical to the frozen baseline once the two new columns are stripped. Two-sided **85.2%** vs one-sided **56.2%** on the 3-arm HG002 pool, reproducing 83.9% / 55.0% on master's NA12878 |
| **the gate is cheap** | 19 s on **8 cores, `-p short`**, for three full adjudications against a 43 GB WGS BAM. It was treated as expensive enough to justify freezing a baseline; the baseline is worth freezing for *reproducibility*, not for cost. It can be a per-change gate |
| §5's `ileftmost_cns` 500-vs-5000 | **CLOSED BY MEASUREMENT** `c04e798` — and the premise was wrong: **the live 10× difference is ZERO.** §5 |
| §5's four "interrogate" thresholds | **ALL FOUR now measured.** Nothing on that list is an open lever |
| **§5's SEMANTIC HALF ON GATE [A] — MEASURED** | **DONE** `4b28d53`, job `49836568`, **57 s**, and it needed no pipeline run: the clip SAM survives, so gate [A]'s whole input is on disk. Check 0 reproduces the run TWICE (the ledger's `A_clip` split *and* its survivor file as a set). **Its polyA THRESHOLD is dead** — `low_polya_clip` is unreachable at the shipping cutoff and first fires at 5, so the decision is a **membership** test; **membership, not clip depth, is what binds the gate** (L1 6,056 dropped against 2,802 kept, 594 of the dropped at `n_clip >= 10`); and **gate [A] has two thresholds and one dial**, since `n_polya_cutoff` is derived from `n_clip_cutoff`. §5 |
| **`NEARBY_CLIP`'s "measured inert" is SCOPED to gate [B2]** | **CORRECTED** — the same constant filters the clip parse at `x_transduction.py:865`, which is gate [A]'s domain, and unbounded it is worth **3.7–3.9× gate [A]'s survivors** and 4× its clip records. Not a recommendation to move it; a statement that it is unmeasured where it matters most. §5, §10 |
| the two **"orphan"s** | **DONE** `18aed56` — no longer a measurement but a **mechanism**: one writer creates the SITE and emits a *bare* source, and the other axis reads the *tail*, so the two are **disjoint by construction**. 37/37 on xTea2, **43/43 on master** — inherited. Different words now, in code (`SUBTYPE_ORPHAN_SITE` vs `CALL_CLASS_ORPHAN`) as well as prose; three wrong-population citations fixed including one in `verify_transduction`'s own printed header; **no scored number moved** (re-verified). §10 |

**AND A THIRD CLAIM OF THIS PLAN'S IS CORRECTED, in the same family as the other two but from a new
cause — two fields whose names both say "orphan" and which agree on nothing.** §1's "BASELINE-0's L1
orphans are 0/37 on-truth" is about `SUBTYPE == orphan_or_sibling_transduction`. The Stage 3 gate's
"orphan-class" pool entries are `td_src_call_class == "orphan"`. Cross-tabulated over every
transduction-asserting call in job `49792012`: **all 37** of the SUBTYPE calls are `regular` by
`TD_SRC`, and SVA has **zero** SUBTYPE-orphan calls at all. **The two do not agree on a single call.**
So the handoff §10.1's stated reason for not resting the flip on 2 orphan-class VERIFIED entries cited
evidence about a different population. The *decision* is unaffected — declining was conservative
either way, and the flip rests on structural correctness plus `0 → 241` — but any sentence of the form
"orphan calls are X" now has to name its field. §10.

**Progress — 2026-08-04. STAGE 3's STRUCTURAL HALF IS COMPLETE AND GATED AT BOTH SCALES; ITS TRUTH GATE
NOW HAS A SPEC; AND TWO OF THIS PLAN'S OWN CLAIMS WERE MEASURED WRONG.**

| item | status |
|---|---|
| the decision as a record + named gates (`xtea/transduction_evidence.py`) | **DONE** `48ec4c7` — three legs PASS (chr1 98/98 byte-identical `49634289`, golden 3/3 `49634291`, WGS 4/4 BASELINE-0 identical `49634290`) |
| gate **[B2]** — `check_seprt_disc_consistency` instrumented, the last big un-attributed population | **DONE** `ff7cde1` — three legs PASS (`49703320` golden, `49703321` chr1 99/99, `49703322` WGS) |
| Stage 3's **transduction** gate | **BUILT, BASELINE FROZEN, RULE AGREED, AND AS OF 2026-08-05 IMPLEMENTED — see the block above for the restratification and the `--baseline` no-regression mode.** As of 2026-08-04: **BUILT, BASELINE FROZEN, AND ITS PASS/FAIL RULE IS NOW AGREED (Corinne, 2026-08-05).** `docs/handoff/runs/probe_td_pool.py` (the pool) + `benchmark/verify_transduction.py` (the verifier), `2d38feb`/`9d8e1de`. **Three-arm HG002 pool: 115 entries, 59.1–81.0% over 324 sweep points**, 82 of 115 asserted by all three arms — frozen with sha256s at `/n/scratch/users/c/cos689/xtea2_testdata/td_gate_baseline_hg002/`. Master on HG002 is job `49763841` (1:10:15). **The rule: primary stratum is shape × sidedness, source class kept with empty cells; the margin is absolute VERIFIED counts per stratum and PROVISIONAL; and the decoy fix is justified on structural correctness + step-level recovery, with this gate only showing no regression.** Also: **the source-class stratum was reading disc-pair midpoints as element intervals — see §10** — so the spec's per-family "reference-copy source" counts and its §7 step 3 validation target do not mean what they say. State, the decisions and their reasons: **`HANDOFF_2026-08-04_td_gate_built.md`** (§10, §11) |
| the SVA reference-source defect | **FLIPPED ON BY DEFAULT 2026-08-05 — see the block above; the "off by default" below is superseded.** Originally **FIXED behind `--td-decoy-structural`, off by default** `55055ba` — one shared classification (`transduction_evidence.decoy_tags`), 11 tests of which 4 drive the shipping seam, **four gate jobs all PASS**. **241 reference-copy SVA sources from 0**; L1 and Alu byte-identical to BASELINE-0 with it on; cost one off-truth call; sensitivity unchanged. The *default* stays the frozen behaviour until the Stage 3 transduction gate can score the 241 — and that gate needs this flag on. §4 Stage 3's last block, §10. **AND THE GATE NOW SAYS IT CANNOT SCORE THE 241, because 241 is the wrong level.** The gate's pool is `TD_SRC` in the call set; there the arm contributes **4** reference-copy SVA entries, 3 orphan-class and 1 regular, of which the verifier returns 2 VERIFIED (both orphan) and 2 UNDECIDABLE — the undecidable pair because their source windows are SVA_D copies at 2.7% and 3.8% divergence, and SVA's youth makes that the rule rather than the exception. **The default decision is now MADE (Corinne, 2026-08-05): flip it on, justified on structural correctness plus the step-level `0 → 241`, with the gate only required to show no precision regression** — measure no-regression first, flip second, re-run both gate scales third. `HANDOFF_2026-08-04_td_gate_built.md` §10.1 and §8 NEXT item 3 |
| the decoy **packaging** half of the same defect | **SCOPED, NOT BUILT** — the 45-record decoy set is reconstructible byte-identically from `rep_lib/consensus/` on every build, and `consensus/transduction_decoy.fa` has **zero references in xTea2 or master** while `construct_novel_flanks`' comment claims it is concatenated in. This is what makes §5's SVA/Alu hazard testable at all. §5 |

**The decision is SIX gate groups behind THREE batched passes, not five behind two.** `[B2]` sits behind
the disc-to-consensus realignment of gate [B]'s survivors, a barrier this plan's diagram omitted.

**AND THE SENSITIVITY CEILING IS NOT IN THIS DECISION, confirmed twice.** 306 of BASELINE-0's 322 false
negatives (95%) never reach the transduction step; all 16 attributed misses are gate [A] with an excess
over chance of **~8 L1 calls**; and instrumenting `[B2]` — the second-largest population in the step —
attributed **zero**. So Stage 3's threshold work on transduction is a **precision** story, which is why
its gate spec measures precision and says so out loud.

**TWO CLAIMS IN THIS PLAN WERE WRONG AND ARE CORRECTED IN §10. Both were wrong in KIND, not degree, and
both came from reading a set difference as a decision:**

1. **`check_seprt_disc_consistency` is NOT "the step's largest threshold-driven rejection after gate A".**
   Its population is right (135 of L1's 159 picked sources, 49 of SVA's 105) but **146 of those 184 are
   an ABSENCE of evidence and only 38 are threshold rejections**, so it ranks *below* gate C. Consequence:
   **`NEARBY_CLIP=50` is a small lever** — it can reach 38 rejections plus 80 survivors, not 184.
2. **`decoy_cross_family` firing 144× for L1 and 0 for SVA is not evidence against §5's SVA/Alu hazard.**
   The SVA arm of that gate **cannot fire**, and behind the same line **every reference SVA source flank
   is misread as the consensus decoy and skipped** — measured: **0 of SVA's 105 picked sources are
   reference copies, all 105 are polymorphic**, on every genome build. §5's hazard is *untested*.
   **FIXED LATER THE SAME DAY, behind `--td-decoy-structural` and off by default** (`55055ba`) — the
   *code* half. The *packaging* half is what leaves §5's hazard untested, and it turns out to need no
   library repackaging: the 45-record decoy set is reconstructible byte-identically from
   `rep_lib/consensus/` on every build, and `consensus/transduction_decoy.fa` has zero references
   anywhere in xTea2 or master.

~~**AND A KNOB THIS PLAN HAS NEVER NAMED IS NOW THE LARGEST ONE BETWEEN GATES [B] AND [C]:**
`MIN_DISC_MAPPED_RATIO = 0.7` (`global_values.py:122`) produces all **146** `no_cns_disc_reads` loci —
four times everything `[B2]`'s own thresholds decide.~~ **MEASURED THE SAME DAY AND IT IS WRONG IN THE
SAME WAY THE CLAIM IT CORRECTED WAS** — a population was attributed to a threshold that does not decide
it. `MIN_DISC_MAPPED_RATIO`'s **ceiling is 34 of those 146 loci**, and at 0.5 it recovers **10**;
**112 of the 146 have EVERY disc-to-consensus record come back unmapped from bwa**, which no value of
the cutoff can reach. See §10, and note the enabling fact: **the disc-to-consensus SAM is not deleted
at the end of the step** — every gate passes `--keep-intermediates` — so this was offline and took
seconds (`docs/handoff/runs/probe_disc_cns_mapped_ratio.py`, reconciled exactly against the ledger).

**Progress — 2026-08-03. THE PER-LOCUS LOOP'S FIRST DEPENDENCY IS DISCHARGED: SITE ENUMERATION.**
"Unified site creation" was a named Stage 3 dependency with a size (52 of 1,515 records) and no
mechanism. It has one now, and it is smaller than it looked: **all 52 are already rows in the
raw-disc file the disc adapter writes before `cns_filter` runs**, and the selection that turns that
file into transduction's candidate list needs no BAM, no bwa and no cns evidence. That rule is now
`xtea/site_universe.py`, gated **offline, 4/4 byte-identical** against a real WGS run's own
`.slct`/`.slct2` files (`docs/handoff/runs/gate_site_universe.py`). Consequences: the per-locus
loop's iteration domain is **41,604 union loci** (vs 75,853 per-family visits, 1.82×); the 52 are
*created*, not moved, so no matching window recovers them; **the 7-of-22 ceiling is unchanged**; and
**§3's diagram order is wrong** — family cannot be assigned before evidence is gathered. See §4
Stage 2c's last block.

**Progress — 2026-08-02. STAGE 2a IS DONE AND GATED; THE STAGE 3 BENCHMARK EXISTS AND HAS ALREADY
CORRECTED THIS PLAN THREE TIMES.**

| item | status |
|---|---|
| **Stage 2a** — retire `IS_CALL_SVA` for a family parameter bundle | **DONE** `de529a4`, `xtea/family_params.py`. chr1 agreement gate job `49481391` **PASS**, 9/9 byte-identical on BOTH gather paths against the Stage 1 gate's arms. Found and preserved a dead branch — see §4 Stage 2 |
| **Stage 2b** — the td_gather manifest | **DEFERRED into Stage 2c, with reasons** — §4 Stage 2 |
| the Stage 2c gate's baseline | **DONE** `1c56b8d` — `locus_table` now records every family collision it resolves (`merged/locus_family_conflicts.tsv`) instead of discarding it |
| **Stage 3 prep** — `benchmark/benchmark_truth_set.py` | **DONE** `7f0bb95`. Scores BASELINE-0 today. Its findings are below and two of them are load-bearing |
| **Stage 1's default flip** | **DONE** — gate job `49480961` **PASS** (13/13, 4/4, 6/6 vs BASELINE-0), and `--td-gather-scan` is now the default. `downstream_parallel` **1144.6s → 457.4s (2.50×)** |

**BASELINE-0 against the HG002 truth set, measured** (position window ±200 bp, one-to-one nearest;
**not comparable to the paper's numbers** — hs38DH, §7.1):

| | truth | calls | TP | FN | FP | sensitivity | precision |
|---|---|---|---|---|---|---|---|
| L1 | 197 | 188 | 140 | 57 | 48 | 71.1% | 74.5% (**92.7%** excluding orphan sites — see below) |
| Alu | 1355 | 1223 | 1103 | 252 | 120 | 81.4% | 90.2% |
| SVA | 90 | 82 | 77 | 13 | 5 | 85.6% | 93.9% |
| all | 1642 | 1493 | 1320 | 322 | 173 | **80.4%** | **88.4%** |

**Two findings that change what Stage 3 can gate on. Read both before designing its gate.**

1. **The transduction denominator is 0 (L1) and 2 (SVA), not the ~51-60 §7.1 projected. The
   ALT-residual route cannot gate Stage 3's transduction work on HG002, and this is a property of
   the data, not a mis-set threshold.** The attrition is fully visible (job `49482568`): 197 L1
   truth records → 189 align to the consensus → **20** carry a 3′ residual ≥100 bp → 10 map back at
   MAPQ ≥30 → 7 also align over ≥80% → 4 do not land inside a reference L1 → **0** sit within 10 kb
   of a full-length (≥5000 bp) reference L1.

   Two independent reasons, both measured, and the second is the one that settles it:

   - `hg38_L1.fa.out` holds 978,710 L1 intervals of median length **285 bp**, of which only 10,437
     (1.1%) reach 5000 bp. "Adjacent to a full-length reference L1" is asking for a rare
     coincidence, and the nearest L1 to each rejected residual is a fragment of 207-3233 bp.
   - **Most long 3′ residuals are divergent element body, not transduced sequence.** 11 of the 20
     L1 residuals ≥100 bp map *inside* a reference L1, and **every one is an older subfamily** —
     L1PA6/7/8 at 4-7% divergence, L1MB/MD/ME at 16-25%. `consensus/LINE1.fa` has 4 L1HS-centric
     records, so an older L1's body partly fails to align and surfaces as "residual". The single
     record that used to pass the whole conjunction sat inside an **L1ME3G at 25.2% divergence** —
     it was divergent body next to a full-length L1 by coincidence. Adding the containment test
     took L1 from 1 positive to **0**.

   Loosening the thresholds does not fix this; it admits divergent body.

   **A SECOND, INDEPENDENT METHOD AGREES** (job `49482779`,
   `docs/handoff/runs/probe_alt_vs_genome.py`). Align the truth ALTs to the **genome** instead of
   to the consensus, which cannot have a subfamily-coverage problem because the reference contains
   every subfamily, and read the segmentation. The genome explains the ALTs almost completely —
   median coverage **0.99**, 161/196 L1 and 77/90 SVA above 0.95 — and **195 of 197 L1 ALTs contain
   no unique non-L1 segment of 100 bp or more at all.** One L1 candidate survives and is not
   adjacent to a full-length L1; SVA yields 10 candidates of which 3 are. So: **L1 zero by both
   methods, SVA two or three.**

   **The ALT route is exhausted on HG002. Do not build Stage 3's transduction gate on it and do not
   try a third variant.** What remains is the paper's supplement, if it enumerates transductions, or
   a different sample. A plausible explanation for the absence, offered as hypothesis not
   measurement: the benchmark was curated by verifying "a polyA tail and target-site
   duplication/deletion" per insertion, which selects for clean canonical TPRT products — plausibly
   the ones that do *not* carry a transduced 3' segment.

2. **All 37 of BASELINE-0's L1 orphan-transduction calls are off-truth, and they are 37 of L1's 48
   false positives.** An orphan transduction *site* is where transduced sequence landed **without**
   the element body — it is not an MEI insertion, and this truth set enumerates MEI insertions, so
   an orphan call cannot match whatever the caller did. The split is total rather than noisy, which
   is what makes it a category error and not a quality signal: the 125 L1 calls with no `TD_SRC` are
   118/125 on-truth, the 26 sibling calls 22/26, the 37 orphan calls **0/37**. Quoting L1 precision
   as 74.5% without saying this is misleading; so is quietly reporting 92.7%. The harness prints
   both, as a named stratum.

**Three corrections to §7.1 itself, found by implementing it:**

- **§7.1's residual figures ("median 18 bp, ≥100 bp: 60, ≥200 bp: 51, max 4620") are TOTAL soft
  clip, both ends** — reproduced exactly at 17 / 60 / 51 / 4619. §7.1's own criterion asks for the
  **3′** residual, which is 7 / 22 / 15 before the supplementary-alignment fix and **5 / 20 / 12**
  after. The leading clip is the larger half (40 records ≥100 bp, max 1745) and is 5′ inversion or
  genomic flank in the ALT — the very confounders §7.1 cites as its reason for requiring 3′.
  **Anything sized off the 60 is wrong.**
- **Supplementary alignments are consensus sequence, not residual.** bwa splits a long ALT rather
  than clipping it whenever the far end also hits the consensus, which a 5′ inversion does; 38 of
  189 aligned L1 truth ALTs and 33 of 90 SVA ones carry one. The residual is the ALT minus the union
  of *all* alignment spans.
- **"3′" is the element's 3′, not the ALT's right end.** The ALT is written reference-forward, so a
  minus-strand insertion's 3′ end is its LEFT end. The plausible shortcut — take the primary's
  trailing soft clip, since bwa revcomps SEQ — returns the 5′ end on every minus-strand insertion.

---

## 1. Context

xTea2's downstream phase is three near-identical per-type pipelines (L1, Alu, SVA), each a chain of
`x_TEA_main.py` subprocesses. It is the last unrestructured part of the tool: upstream detection is
already a single C++ `xtea_scan` pass, and the type-agnostic tail (gene/genotype/gVCF) is already
unified behind `locus_table.py`. The middle is still type-major, and it is where nearly all time goes.

**Why now — three measured facts:**

1. **Transduction is ~93% of downstream wall time** (WGS, `-n 8`: L1 1481s, SVA 1772s), and it is the
   least clean code in the tree — six phases, ~15 hardcoded cutoffs each able to veto a call, three bwa
   passes, ~4.4 GB of SAM intermediates, ~108K per-site temp files per type.
2. **The three pipelines process the same loci.** Raw-disc positions are 53,953 / 53,960 / 53,962 for
   L1/Alu/SVA with **53,945 shared by all three** (only 8 L1 positions absent from Alu's), because the
   raw-disc write gate keys on family-*independent* columns (`disc_adapter.py:227-229`). Family is
   decided one step later by the disc cutoff — not by cns, and not by this list.
3. **The redundancy is hidden by parallelism.** L1 and SVA run as concurrent subprocesses at 12 cores
   each, so duplicate gathers cost CPU, not wall time. Restructuring pays only if the shared work runs
   *once at full width* — which requires leaving the per-type-subprocess shape entirely.

**Relationship to `docs/phase_c_unified_downstream.md`** (approved 2026-07-18): that design already
specifies the right architecture — locus as primary key, fold realignment into the scan, unified tail,
type-specific only where biology demands. This plan **supersedes it only in sequencing**. Phase C put
`cns_filter` (move B) first and deferred transduction as "the gnarliest code to fold into a unified
loop; map it carefully before touching" — written before the C++ cns port landed and before transduction
was profiled. Transduction is now the whole cost. **The ordering inverts; the architecture stands.**

**Intended outcome:** one unified, locus-major downstream; transduction rebuilt around an explicit
evidence model with tunable parameters instead of buried constants; judged on sensitivity/precision
against orthogonal truth sets rather than on agreement with master.

## 2. Decisions locked in

| | |
|---|---|
| Implementation | **Python rewrite; C++ only for the hot loop**, as a new `xtea_scan` subcommand |
| Calling logic | **Full redesign permitted.** Agreement with master is explicitly *not* a requirement |
| Scope | **Entire downstream**: cns_filter, transduction, sibling/orphan, postF ×2, gene/genotype/gVCF, **plus somatic/case-control and mosaic** |
| Output | **One sorted VCF**; per-type VCFs behind `--emit-per-type` as a debug/benchmark artifact |
| Validation | **Germline: HG002** (`parklab/xTea_paper` truth sets). **Mosaic: SMHT005 descending colon (3G)** vs the in-house long-read sMEI truth set. **NA12878 retired** |
| Cutover | **Feature flag** `--downstream=legacy|unified`; both runnable for A/B; delete legacy once the benchmark validates |

## 3. Target architecture

```
BAM/CRAM
   │
   ▼  xtea_scan — ONE pass (existing, extended)
      clip+disc detection · realign to all 3 consensi · mate_in_rep count
      NEW: unified per-site read gather for the downstream
   │
   ▼  locus table (one row per candidate locus, per-family support side by side)
      Locus { chrom, pos, clip_L/R, mate_in_rep, disc_L/R,
              support: {L1:…, Alu:…, SVA:…}, assigned_family, … }
   │
   ▼  family assignment — ONE decision per locus (not three parallel runs)
   │
   ▼  per-locus loop, family-branched only where biology demands:
      cns evidence → transduction (L1/SVA) → sibling/orphan → postF
   │
   ▼  type-agnostic tail ONCE: gene → genotype → VCF
   │
   ▼  one sorted VCF  (+ per-type VCFs if --emit-per-type)
```

> **CORRECTION 2026-08-03 — the ORDER in that diagram is wrong, the inversion below is not.** It puts
> family assignment *above* the per-locus loop. For **52 of 1,515** final records there is nothing to
> assign on at that point: their only evidence is the transduction evaluation's own output (§4 Stage
> 2c, measured twice). The order the data supports is **enumerate loci → gather per-family evidence →
> assign family → filter**, i.e. assignment lives in the loop's tail. The locus set itself *can* be
> enumerated up front — that is `xtea/site_universe.py`, and it is 41,604 loci — so only the
> assignment box moves, not the "one pass" shape. Designing the loop off the picture as drawn rebuilds
> the barrier §4 Stage 2c already measured as unable to close (7 of 22, ceiling).
>
> **AND THE FIRST BOX IS BUILT, 2026-08-03 — `xtea/downstream_unified.py`, `--locus-universe`.** The
> enumeration now runs in the driver for every family, at the barrier point (which is forced: the rule
> excludes against `cns_filter`'s own output — see §10). 41,604 loci, and the driver's site list is
> byte-identical to what each family's own step derives, gated at WGS scale offline
> (`gate_locus_universe.py`, PASS) and end to end on chr1.
>
> **ONE MORE CORRECTION TO THE PICTURE, and it is about the middle box, not the order.** "per-locus
> loop" cannot mean per-locus end to end. `call_candidate_transduction_v3` builds one flank reference
> and runs **four batched bwa passes over the whole `.slct` list** before a single per-site record
> exists; one bwa per locus is a different program. The shape the code supports is
>
> ```
> enumerate loci (DONE)  →  batched per-family evidence production  →  per-locus decision  →  assign  →  filter
> ```
>
> with only the third box a loop. That is not a retreat from "one pass": the batched stage is what the
> `--td-manifest` gather already does for all families out of one BAM decode.

The key inversion: **family becomes an output field of one pass, not the input that forks the
pipeline.** `xtea/locus_table.py` already implements the cross-family merge (±50 bp window, priority
L1→SVA→Alu, `two_side_tprt_both` override) and stays — but moves from "reconcile three pipelines
afterwards" to "assign family once, up front."

## 4. Sequencing

Each stage is independently committable and leaves the tree working. **Gates change character at
Stage 3** — that boundary is the plan's most important structural feature.

### Stage 0 — Measure and freeze (prerequisite; blocks everything)

Nothing is rewritten until we know where the time goes and what the calls currently are.

**Measure:**
- **Instrument transduction.** Add `time.time()` phase timers around the four blocks in
  `x_transduction.py:call_candidate_transduction_v3` (flank-library build, read gather, the three bwa
  passes, SAM parsing), gated on the existing `XTEA_TD_PROFILE` env var — which today instruments
  **only** `parse_td_sibling_by_site` (~305 calls), *not* the ~50K-site gather.
- **Get the breakdown.** Two runs of `docs/handoff/runs/transduction_attribution.sbatch` (~25 min each,
  staged cns inputs) to confirm reproducibility and split the 1481s between gather / bwa / parse.

**Freeze:**
- **BASELINE-0 on HG002**, for **`L1,Alu,SVA` and `L1` alone** (and a somatic pair if one exists) —
  single-type and somatic regressions are otherwise invisible for the whole rewrite.
- Cut an **HG002 chr1 slice** and an **SMHT005-3G chr1 slice** for fast iteration.

**Fix, as standalone commits:**
- **The `col5` → `mate_in_rep_clip_cnt` rename** (§6.1) — mechanical.
- **The type-set-dependence bug** (§6.1) — clean signature: **three-type output must not move at all**,
  while **SVA-only output changes**, which is the bug being fixed. Doing this in Stage 0 means every
  later baseline is taken against correct behaviour.

**Acquire (longest lead time in the plan — start first):** — **germline inputs DONE 2026-07-29**,
staged at `/n/data1/hms/dbmi/park/corinne/claude_sandbox/test_data/HG002/` with a `PROVENANCE.md`
that is the authority on what was chosen and why. Summary:
- ~~HG002 truth VCFs and `x_cmp.py` from `parklab/xTea_paper`~~ — **DONE**, pinned to commit
  `2338ece`. They also already existed at `corinne/dev_xTea/truth_set_hg002/` (byte-identical), which
  additionally holds an `all_tes.bed` the paper repo does not ship.
- ~~Download the paper's `HG002_bwa_mem_GRCh38.sorted.cram`~~ — **NOT POSSIBLE, see §7.1.** The
  paper's CRAM was built in-house under a directory that no longer exists and is not hosted anywhere;
  GIAB has no bwa-mem GRCh38 Illumina alignment for HG002 at all. **Substituted** with the NovaSeq
  PCR-free 30x bwa-mem BAM from `brain-genomics-public` (43.7 GB, md5-verified) — which is aligned to
  **hs38DH, not `no_alt`**, an accepted and documented deviation.
- ~~Resolve what `_tprt` means~~ — **DONE 2026-08-01, §7.1.** It is a canonical-TPRT-mechanism filter,
  a strict subset of `_all`; the sensitivity denominator is `_all` (197 L1 / 1355 Alu / 90 SVA).
  **How transductions are scored is still open**, and Stage 3 still blocks on it: the truth sets carry
  no source-element annotation at all (verified), so transduction status has to come from the ALT
  sequence or the supplement.
- Mosaic inputs need **no** acquisition — the five SMHT005-3G CRAMs are local, on our reference, and
  header-compatible (§7.2). They need `--ignore-rg` and a chr1 slice.
- **Still to do in Stage 0:** the phase timers, BASELINE-0 itself (both type sets), the chr1 slices,
  and the two standalone fixes above. Only *acquisition* is complete.

**Gate:** BASELINE-0 frozen per type set; transduction phase breakdown reproducible within ~5%;
three-type calls byte-identical across the rename+fix; `--types SVA` and `--types L1,Alu,SVA` agree on
the SVA call set.
**Rollback:** instrumentation is additive and inert by default; the rename is mechanical; revert either
commit independently.

> ### GATE RESULT 2026-08-01 — PASS. `docs/handoff/runs/stage0_gate_check.sh` reproduces it.
>
> Five HG002 arms, jobs `49441918-22`, post-fix at `9a43b80` and pre-fix at `6ce36c7`, same BAM, same
> `XTEA_BASIC_INFO_SEED=42`, all with `--keep-intermediates`.
>
> | run | LINE1 | ALU | SVA | merged |
> |---|---|---|---|---|
> | post-fix `L1,Alu,SVA` | 188 | 1223 | 82 | **1493** |
> | post-fix `L1` | 200 | — | — | 200 |
> | post-fix `SVA` | — | — | 85 | 85 |
> | pre-fix `L1,Alu,SVA` | 188 | 1223 | 82 | 1493 |
> | pre-fix `SVA` | — | — | 85 | 85 |
>
> - **A — three-type byte-identical across the rename+fix: PASS**, 7 of 7 artifacts (3 per-family VCFs,
>   3 per-type call lists, the merged VCF). The rename touched the C++, the Python, the CLI and the
>   spec, and moved nothing.
> - **B — SVA-only: PASS, output-neutral.** The mechanism moved (1,434,938 positions, monotone) but no
>   call did. Full detail in §6.1; do not restate it as "the fix changed nothing".
> - **C — the partition property: PASS.** `--types SVA` and `--types L1,Alu,SVA` produce a
>   **byte-identical SVA per-type call list**. This is the property §6.2 requires and the one the bug
>   violated.
>
> **`--types` partitions cleanly for L1 too** (checked while verifying C): the per-type L1 call list is
> byte-identical between `--types L1` and `--types L1,Alu,SVA` — 200 in both. The merged-VCF difference
> (200 → 188) is **entirely `locus_table` family reassignment**: all 12 missing positions match an ALU
> call within the ±50 bp merge window (deltas 0 to −41). None is dropped.
>
> **Method note worth keeping:** matching those 12 by *exact position* accounts for only 4 and reports
> 8 as "genuinely dropped" — a false alarm. Any comparison across the family merge must use
> `locus_table`'s ±50 bp window, because the merge can shift the representative position.
>
> **BASELINE-0 artifacts:** `<run>/BASELINE-0/` holds the VCFs plus `PROVENANCE.txt` (commit, job id,
> BAM md5, reference, seed, wall time, sha256 per VCF). Do not read a VCF's `##reference=` header as
> provenance — `x_gvcf.py` hardcodes it to `GRCh38_full_analysis_set_plus_decoy_hla.fa` regardless of
> the reference actually used.

> **Stage 1's scope is deliberately conditional on this measurement.** If the gather dominates, the C++
> subcommand is the right first move. If bwa or SAM parsing dominates, the first move is collapsing the
> alignment passes and streaming the SAMs. **Do not pre-commit the C++ boundary before the numbers
> exist** — the claimed ~2× from a shared gather is bounded below by ~390s of 1481s and could be
> anywhere from 26% to 80% of the step.
>
> ### RESOLVED 2026-08-01 — the gather dominates. Take the C++ gather branch.
>
> Job `49437990`, `docs/handoff/runs/transduction_phase_timing.sbatch`, live tree at `de96727`, the
> staged NA12878 cns inputs, L1+SVA concurrent at `-n 8`, two repeats.
>
> | | gather | bwa (3 passes) | SAM parse | TOTAL |
> |---|---|---|---|---|
> | **L1** | **857.7s / 86.2%** | 115.4s / 11.6% | 18.3s / **1.8%** | 995.6s |
> | **SVA** | **898.2s / 72.9%** | 310.1s / 25.2% | 18.1s / **1.5%** | 1231.9s |
>
> Split further: `clip_gather` is the single biggest phase in both (577.6s L1 / 622.8s SVA),
> `disc_gather` next (280.1s / 275.5s). `flank_build`, `flank_bwa_index`, `disc_cns_reselect`,
> `disc_cns_bwa_realign`, `disc_cns_sam_parse`, `filter_transduction` and `genotype_features` are each
> **under 6s** — together under 1% of the step.
>
> **Gate met.** Repeat-to-repeat spread 2.4% (L1) and 0.6% (SVA), inside the ~5% requirement; the two
> repeats' `candidate_disc_filtered_cns2.txt` are **byte-identical**, so the step is deterministic. The
> 13 phases account for **99.6%** of the function, and the function is **95-96%** of the whole
> `--transduction` step (whole-step wall 17:28.8 L1 / 21:21.4 SVA) — so there is no material
> unmeasured remainder hiding a different answer.
>
> **Three consequences, two of which correct this plan:**
>
> 1. **The C++ per-site gather subcommand is the right Stage 1 move**, and by a wide margin — the
>    gather is 73-86%, far above even the optimistic end of the 26-80% range this note worried about.
> 2. **SAM parsing is not a cost. It is 1.5-1.8%.** §5's "streaming the SAMs" and the framing of the
>    ~2.4 GB clip SAM parse as a target are **refuted**: parsing all ~4.4 GB of SAM intermediates costs
>    ~18s per family. Do not spend Stage 1 effort there. (Writing/reading those files may still be worth
>    deleting for I/O and tidiness reasons — but not for time.)
> 3. **bwa is secondary but asymmetric**: 11.6% for L1, 25.2% for SVA — SVA's realign passes cost 2.7×
>    L1's, almost all of it in `clip_bwa_realign` (160.1s vs 85.4s) and `disc_bwa_realign` (144.0s vs
>    29.2s). §5's "target two alignment passes, not three" is worth *more* for SVA than for L1, but it
>    is a Stage 3 question and cannot be the Stage 1 lead.
>
> **Amdahl ceiling, so Stage 1 is not oversold:** removing the gather *entirely* caps the step at
> **7.2× (L1)** and **3.7× (SVA)**. The separate ~1.9× cross-family sharing win multiplies against the
> gather portion only, and §6.2 already records that a single-type run gets none of it.
>
> **Do not compare these totals against the plan's "L1 1481s / SVA 1772s" as a speedup.** Different
> scope (that is the whole step from a WGS run's `step_timings.tsv`, this is one function on staged cns
> inputs) and different code (the live tree carries `fb1d766`'s BAM cache; the figures in §1 predate
> it). No A/B was run, and per §1 of the 07-30 handoff none should be.

### Stage 1 — Mechanics only (agreement-gated)

Restructure *how* work happens without changing *what* is decided. **This is the last stage where a
byte-diff is a valid gate**, so extract maximum value from it.

- New `xtea_scan` subcommand for the unified per-site gather, reusing two patterns already proven in
  `src/csrc/cns_emit.hpp`: `FetchRegion{chrom,beg,end,own_from}` coalescing with the exact ownership
  rule (`:502-523` — a read belongs to a region iff its start ≥ the previous region's end), and
  **candidate-major emission** (`:798-843`).
- Eliminate the per-site temp-file round-trip and the per-chromosome `Pool` churn
  (`x_clip_disc_filter.py:351-358` creates a fresh Pool per chromosome, ~25 fork storms; each task
  writes two tiny files the parent re-reads and unlinks).
- Manifest-driven, mirroring `cns_manifest.txt` (`run_xtea.py:995-996`, parsed `xtea_scan.cpp:532`).

**Gate:** `tests/test_downstream_from_cns.py` (VCF cols 1-8) plus a chr1 diff vs BASELINE-0 —
**byte-identical calls required** — plus a two-run determinism check.
**Rollback:** ~~legacy remains the default; the flag selects the new gather~~ — **the native
gather became the DEFAULT 2026-08-02** once the WGS flip gate passed; `--no-td-gather-scan` is the
rollback. See the gate result below.

> ### DEFAULT-FLIP GATE 2026-08-02 — PASS, and `--td-gather-scan` IS NOW THE DEFAULT
>
> Job `49480961`, `docs/handoff/runs/stage1_flip_wgs.sbatch`, HG002 **WGS**, `L1,Alu,SVA`, two arms
> at `58d2203` from the frozen `stage1_flip_wt` worktree, 24 cores each, sequential.
> This is the run the Stage 1 handoff left open: *"one WGS run whose calls match a BASELINE-0 arm —
> the chr1 gate does not substitute, being 17 calls."*
>
> | leg | result |
> |---|---|
> | arms ran different paths | LEGACY 2 `clip_gather` / 0 `native_gather`; NATIVE the reverse |
> | AGREEMENT LEGACY vs NATIVE | **PASS** — 13/13 (per-type call lists, sibling lists, postF_hc lists, per-type VCFs, merged VCF) |
> | arm LEGACY vs frozen BASELINE-0 VCFs | **PASS** — 4/4. *No drift across the five commits since `9a43b80`* |
> | arm NATIVE vs frozen BASELINE-0 VCFs | **PASS** — 4/4. **The flip criterion** |
> | per-type call lists vs BASELINE-0 | **PASS** — 6/6 for each arm |
> | counts | all three: L1 188 / Alu 1223 / SVA 82, merged 1493 — identical to BASELINE-0 |
>
> **And the numbers the flip is for**, both arms at n=24 on one node — these supersede the
> projections below, which were extrapolated from the standalone oracle at different thread counts:
>
> | | legacy | native | speedup |
> |---|---|---|---|
> | L1 gather | 664.6s (85.5% of the function) | **50.4s** (29.5%) | **13.2×** |
> | SVA gather | 741.1s (74.8%) | **52.9s** (17.3%) | **14.0×** |
> | L1 `call_candidate_transduction_v3` | 777.2s | **170.9s** | 4.5× |
> | SVA `call_candidate_transduction_v3` | 991.1s | **304.8s** | 3.3× |
> | **`downstream_parallel`** | **1144.6s** | **457.4s** | **2.50×** |
> | whole run | 2888s (48:08) | **2166s (36:06)** | 1.33× |
>
> `downstream_parallel` at **2.50×** beats the ~1.8× this plan projected. The whole-run figure is
> the honest headline and it is limited by detection: `scan_pass` is ~1500s of both runs and this
> plan does not touch it.
>
> **What flipping the default also does, recorded so it is not a surprise:** it makes the multi-BAM
> divergence reachable. The Python gather truncates one shared `all_disc_pos` per input BAM, so with
> N>1 BAMs only the last one's positions survive; `td_gather.hpp` gathers each BAM against its own
> anchors. That is a fix, but the mosaic path now runs it with no baseline of its own — a Stage 4
> concern (§6.4, §7.2), not a blocker here. `--no-td-gather-scan` restores the old behaviour, and
> the three gate scripts under `docs/handoff/runs/` now pass it explicitly to build their legacy arm.

> ### GATE RESULT 2026-08-01 — PASS. `docs/handoff/runs/stage1_gate_chr1.sbatch` reproduces it.
>
> Job `49453032`, HG002 chr1, `L1,Alu,SVA`, `--keep-intermediates`, `XTEA_BASIC_INFO_SEED=42`,
> three arms at commit `49154dc` from the frozen `stage1_wt` worktree: **A** legacy gather, **B**
> native, **C** native again.
>
> - **AGREEMENT (A vs B): PASS, 9 of 9 artifacts byte-identical** — 3 per-type
>   `candidate_disc_filtered_cns2.txt`, 2 sibling lists (Alu produces none), 3 per-type VCFs, the
>   merged VCF.
> - **DETERMINISM (B vs C): PASS**, the same 9 of 9.
> - **`tests/test_downstream_from_cns.py`: PASS** with `XTEA_SCAN_TD_GATHER` forced on — i.e. the
>   native gather also reproduces the NA12878 chr1 *master* golden's downstream, a second sample and
>   the only external baseline in the gate.
> - All three arms: L1 1 / Alu 14 / SVA 2, merged 17. Whole-run 419s / 344s / 356s.
>
> **How the chr1 A/B substitutes for "vs BASELINE-0".** BASELINE-0 is a WGS artifact and cannot be
> diffed against a chr1 run; an A/B at one commit is the same statement about the port with the
> sample and commit confounds removed. BASELINE-0 remains Stage 2's comparator, where calls are
> allowed to move.
>
> **Read this gate's weakest link honestly: 17 calls.** That is normal for a chr1 slice — the
> NA12878 chr1 master golden is 32 (L1 0 / Alu 27 / SVA 5), and §7.1 already notes chr1 cripples
> L1 — but it means 17 identical calls is a *wiring* proof, not a coverage proof. The load-bearing
> evidence for this change is that the gather's own output FILES are byte-identical: 3733 chr1
> sites (89,716 clip records / 60,246 disc anchors / 41,859 mates) and 49,042 L1 / 53,569 SVA sites
> over 24 contigs at WGS scale (`stage1_gather_oracle_wgs.sbatch`). Identical bytes into bwa is
> what makes identical calls a consequence rather than a coincidence.

> **Ordering determinism is load-bearing.** `cns_emit.hpp:798-817` documents that bwa's tie-break among
> equal-scoring placements depends on read index in the input, so reordering the FASTQ silently changes
> which repeat copy multi-mapping clips land on. Candidate-major emission order must be preserved
> exactly, or Stage 1 cannot be agreement-gated at all.

> ### BUILT 2026-08-01 — `src/csrc/td_gather.hpp`, byte-identical to the Python gather
>
> One `xtea_scan --td-sites` invocation replaces **both** Python gathers —
> `collect_clipped_reads` (per-chromosome `Pool`, one forked task and two temp files per site) and
> `collect_disc_reads` → `extract_mate_reads_by_name` (whole-genome scan in 50 Mb bins, re-parsing
> the entire disc-pos file once per bin). Reuses the two patterns this section names: the
> `FetchRegion` merge with its exact ownership rule (`cns_emit.hpp:502-523`) and candidate-major
> emission (`:798-843`). Doing both halves in one process is safe because the disc half never
> depended on anything computed between them — `collect_disc_reads` ignores its own site-list
> argument and reads the file the clip half wrote (`x_clip_disc_filter.py:392`).
>
> **THREE ORDERING RULES, and they are not the same rule.** Getting any of them wrong passes a
> record-level comparison and fails the byte-level one:
>
> 1. The clip FASTQ is candidate-major in the **site file's** order — contigs by first appearance,
>    positions in file order within a contig. `sf_raw_slct` is **not** position-sorted (the cns
>    candidates are written first, then the raw-disc lines), so copying `emit_cns_multi`'s
>    position-sorted fan would have been wrong here. That is why this is a separate collector
>    rather than a parameter on `emit_cns_multi`.
> 2. `n_cnt` (read-name field 5) keeps master's **gaps**: master builds the right-clip name and
>    increments the counter *before* that side's length/polyA/quality gates, so a failing right
>    clip burns an index. Nothing reads field 5 — reproducing it is what makes `cmp` a valid
>    oracle, and index-burning right clips are therefore held as no-emit events.
> 3. The disc FASTA is **BAM-header** contig order (decoy names skipped), coordinate order within
>    a contig, and per mate read its insertion sites in disc-pos order.
>
> **Measured attribution** (chr1 L1 site list, 3733 sites, `NA12878.chr1.bam`, **2 threads** so it
> is apples-to-apples with the legacy `Pool` at `n_jobs=2`):
>
> | | legacy | native |
> |---|---|---|
> | clip half | 24.0s | collect 19.0s + fan 0.1s |
> | disc half | 106.5s | mate-gather 11.1s |
> | **total** | **130.5s** | **30.2s** |
>
> **The two halves are not equally improved, and the ratio INVERTS with scale** (equal thread
> counts on both sides of each row):
>
> | | clip half (legacy → collect+fan) | disc half (legacy → mate-gather) |
> |---|---|---|
> | chr1, 3733 sites, 2 threads | 24.0s → 19.1s = **1.26×** | 106.5s → 11.1s = **9.6×** |
> | WGS L1, 49,042 sites, 12 threads | 615.0s → 23.8s = **25.8×** | 197.0s → 29.1s = **6.8×** |
>
> The chr1 row alone reads as "the clip half is index-fetch bound, so its win comes from threads
> rather than from C++" — and **the WGS row refutes that reading**: same threads, same windows,
> 25.8×. The clip half's cost is per-site fork-and-temp-file overhead, invisible at 3733 sites and
> dominant at 49,042. The disc half moves the other way: its old cost was a fixed genome scan
> (~106s at 207 sites and at 3733, 197s at 49K) while the native mate gather does work
> proportional to the number of mates.
>
> **End to end, in the pipeline** (HG002 chr1, gate job `49453032`, arms A vs B, `-n 6` per type
> under `--parallel-downstream`):
>
> | | legacy gather | native gather | `..._transduction_v3` | the whole step |
> |---|---|---|---|---|
> | L1 (2815 sites, extnd 956) | 75.1s (clip 38.5 + disc 36.7) | **21.8s (3.44×)** | 83.6s → 30.5s (2.74×) | 100.7s → **46.9s (2.15×)** |
> | SVA | 82.1s (clip 45.6 + disc 36.5) | **22.4s (3.67×)** | 101.6s → 41.8s (2.43×) | 115.4s → **54.8s (2.11×)** |
>
> Three levels, because quoting the wrong one oversells this: the gather 3.4-3.7×, the instrumented
> function 2.4-2.7×, the **step** 2.1× — the step also pays automatic-parameter generation, the
> annotation tree and the re-selection, none of which moved. **All chr1 figures, and they
> understate the port**: the WGS block below measures the gather at 15-18×.
>
> Whole chr1 run 6:59 → 5:44, which is the honest end-to-end figure and is small for a reason
> worth remembering: on chr1 the **detection** `scan_pass` is 203.9s of the run, half of it, and
> this plan does not touch it.
>
> ### MEASURED AT WGS SCALE 2026-08-01 — 15-18×, and it CORRECTS the chr1-based projection
>
> Job `49453857`, `docs/handoff/runs/stage1_gather_oracle_wgs.sbatch`, `NA12878.wgs.bam`, the run's
> own `.slct` site lists, extnd 753, **12 threads for both arms**, peak RSS 0.94 GB:
>
> | | sites | legacy | native | speedup |
> |---|---|---|---|---|
> | **L1** | 49,042 over 24 contigs | **812.0s** (clip 615.0 + disc 197.0) | **53.3s** (collect 22.0 + fan 1.8 + mate 29.1) | **15.25×** |
> | **SVA** | 53,569 over 24 contigs | **1011.3s** (clip 817.0 + disc 194.3) | **56.8s** (collect 23.2 + fan 2.1 + mate 31.0) | **17.82×** |
>
> **All three artifacts byte-identical in both families** — clip FASTQ 177,935,469 / 180,890,118
> bytes, `all_disc_pos` 90,338,623 / 94,115,888, disc FASTA 235,340,814 / 245,224,351 — over
> 1.46M / 1.49M clip records, 1.04M / 1.08M disc anchors and 592K / 629K gathered mates. **And
> deterministic**: two native runs per family, byte-identical on all three.
>
> **This supersedes the "1.5-1.7×" projection that used to be here, which was built on the chr1
> 3.4× ratio and understated the win by an order of magnitude on the half that matters.** Why the
> chr1 ratio was wrong: the legacy clip gather's cost is per-site overhead — one fork and two temp
> files each — so it scales with the site count, and chr1 has 13× fewer sites. 6.4 ms/site at 3733
> sites, 12.5 ms/site at 49,042. The native collect is 22.0s for the same windows.
>
> **Recomputed step-level effect** (mixing thread counts, so ±: the 857.7s/898.2s gather figures
> are from job `49437990` at `-n 8`, these native figures are at 12 threads): L1's step
> 1480.7s → **~676s (2.2×)**, SVA's 1772.3s → **~931s (1.9×)**, and `downstream_parallel`
> 1884s → **~1040s (1.8×)** — the SVA chain becomes the critical path. Against the Amdahl ceiling
> of 7.2× on the *function*, the native gather realises 5.2× of it (995.6s → ~191s).
>
> **Suggestive, not established: the legacy clip gather does not scale with threads.** 577s at
> `-n 8` (job `49437990`) versus 615s at 12 threads here, same input root. Consistent with
> fork/temp-file/I-O binding rather than CPU, and with the native collect scaling instead
> (19.0s at 2 threads on chr1's 3733 sites, 22.0s at 12 on 49,042). The two runs' site lists were
> regenerated between them, so do not quote it as a controlled comparison.
>
> **Deferred to Stage 2, deliberately: the manifest.** A `cns_manifest.txt`-style fan is what lets
> ONE collection serve all three types, and today the three types' transduction steps are separate
> concurrent subprocesses with separate site lists — there is nothing to fan across. The collector
> already takes an arbitrary site list and nothing in it is type-specific, so unifying is a caller
> change in Stage 2, not a rewrite here. §1's point stands: the cross-family sharing win only
> materialises when the shared work runs once at full width.
>
> **Three divergences, all documented in the header, all unreachable from this gate** — each differs
> only on input the legacy path cannot process at all: with N>1 input BAMs each is now gathered
> against its own anchors (legacy truncates one shared `all_disc_pos` per BAM, so only the last
> one's positions survive and every sample's mates are gathered against them); an exact duplicate
> site line is collapsed rather than crashing; a read with no base qualities skips the clip quality
> gate rather than raising `TypeError`.

### Stage 2 — Unify the family dimension (predicted movement)

- One pass over loci with per-family support side by side; assign family once via `locus_table.py`.
- Replace the three per-type subprocess chains with a single per-locus loop; family becomes a branch
  (`if family in (L1, SVA): transduction`), not a fork.
- Retire `global_values.IS_CALL_SVA` global state for an explicit per-family parameter bundle carrying
  `flank` (None for Alu), `leftmost_cns` (500 / **5000 for SVA**), `annotation_extend` (0 / **200 for
  SVA**).
- Emit one VCF; per-type VCFs behind `--emit-per-type`.

**Gate:** calls may move, but **only in a predicted direction** — borderline multi-family loci may flip
family. Every moved call must be explainable as a family-assignment flip; report the per-family delta vs
BASELINE-0 and hand-inspect the flips. Anything else is a bug.

> ### 2a DONE 2026-08-02 — `de529a4`, and it surfaced a branch that has never run
>
> `xtea/family_params.py`: a `FamilyParams` bundle carrying `has_flank`, `annotation_extend`,
> `leftmost_cns` and the two SVA cns-path relaxations, passed to `XClipDiscFilter` /
> `XTransduction` / `XOrphanTransduction` instead of flipped process-wide.
> `global_values.IS_CALL_SVA` and `turn_on_sva()` are **deleted**, not merely unused — a leftover
> attribute would let a stale `global_values.IS_CALL_SVA = True` assignment (one existed, in
> `docs/handoff/runs/td_gather_oracle.py`) keep parsing while silently doing nothing.
>
> **Gate: PASS.** Job `49481391`, `docs/handoff/runs/refactor_gate_chr1.sbatch` — HG002 chr1, both
> gather paths, 9/9 artifacts byte-identical to the Stage 1 gate's arms A and B in each case. The
> validity condition that script cannot check itself was checked by hand:
> `git diff 49154dc..HEAD -- xtea/ src/` is the gather port plus this refactor and nothing else.
>
> **THE FINDING, and why the obvious cleanup is a Stage 3 change.** `--sva` is passed to the
> **cns_filter step only** — `build_type_steps` matches the golden `run_xTEA_pipeline.sh`, which
> puts it on -C/-D/-N and not on transduction/sibling/postF. So `IS_CALL_SVA` was **False inside the
> transduction subprocess even for SVA**, and SVA's `leftmost_cns = 5000` has never executed, in
> xTea2 or in master. Resolving the bundle from `--rtype` instead — which looks like tidying — wakes
> it, and that is not a family-assignment flip:
>
> > `pick_full_length_insertion` keeps a candidate whose cns cluster start is
> > `0 < start < leftmost_cns`. The longest record in `consensus/SVA.fa` is **3183 bp**, so a 5000
> > threshold admits *every* candidate with a positive cluster start. The polymorphic-flank library
> > stops being a filtered subset and becomes the whole candidate list, and every downstream bwa
> > realignment goes with it.
>
> Preserved dormant, with `tests/test_family_params.py` pinning the resolution rule by reading the
> argv `build_type_steps` actually emits — so waking it has to be a decision, made in Stage 3
> against the truth set, rather than a side effect.

> ### 2b DEFERRED into 2c 2026-08-02 — the manifest, and the arithmetic for not building it yet
>
> `HANDOFF_2026-08-01_stage1_gather.md` §4 lists the manifest as the next item. It is being folded
> into 2c rather than built standalone, for two reasons that are worth writing down because the
> handoff's own framing already implies the first:
>
> - **The wall-clock case is ~zero.** L1's and SVA's transduction steps are separate concurrent
>   subprocesses, so their two native gathers (53.3s and 56.8s at WGS scale) **overlap**. A shared
>   collection would still be ~57s of wall time; what it saves is CPU. §1 already says the sharing
>   win "only materialises when the shared work runs once at full width", and full width is what 2c
>   creates.
> - **The caller half is throwaway work.** Hoisting the gather means the driver must compute the
>   step's `iextnd` (from `automatic_gnrt_parameters`) and run `re_slct_with_clip_raw_disc_sites`
>   itself, then hand each type pre-gathered artifacts — re-implementing step logic that 2c deletes.
>
> The C++ side is genuinely reusable and is what 2c should build: generalise `td_gather.hpp` from one
> `TdSites` to a vector of them, build the fetch plan over the **union** of the entries' windows, and
> emit per entry. The union is nearly free here — 53,945 of ~54K sites are shared by all three
> families — so one decode serves every family.
>
> ### 2b's C++ SIDE BUILT 2026-08-02 — `run_td_gather_multi`, `xtea_scan --td-manifest`
>
> `td_gather.hpp` now takes a `std::vector<TdGatherEntry>`; `run_td_gather` (the single-family form
> `XClipDiscFilter.native_td_gather` still calls) is a one-entry vector into it and **nothing else**,
> so "the multi-entry version cannot change the one-entry output" is structural rather than tested.
> The CLI mirrors `--cns-emit-manifest`: `--td-manifest F`, one line per family,
> `tag sites extnd clip_fastq disc_fa [disc_pos]`, `-` for an output the entry does not want.
>
> **The three things that make sharing invisible**, each with a test:
>
> - collection runs over the **union** site index at the **largest** entry's `extnd`, so `ev` is a
>   superset of any solo run's — which is safe only because the header's existing invariant holds:
>   *holding is permissive, emission is exact.* The fan's two `ref_end <= wbeg` tests are the exact
>   half, and `test_entries_may_have_different_extnd` is what would catch losing either.
> - order survives the wider plan: regions own reads starting at or after the previous region's end,
>   so `ev` is ascending however the regions merge, the `stable_sort` stays a no-op, and equal-start
>   reads keep the BAM order **bwa's tie-break rides on**.
> - the mate gather merges its region plan over the union of the entries' wanted mate positions;
>   different boundaries, same ownership rule, so each entry's FASTA is the same byte sequence. Each
>   stripe now writes one part file **per entry**, and the mate's sequence is decoded once however
>   many entries claim it.
>
> **The manifest is stricter than the cns one on purpose**: a short line is a hard error, not a
> silent skip. A dropped field there means a family quietly does not run, which on a WGS gate reads
> as "that family found nothing". Two entries sharing an output path are rejected up front too.
>
> **GATES, both PASS.**
>
> - **WGS scale, job `49516920`** (`docs/handoff/runs/td_manifest_wgs_gate.sbatch`): L1 and SVA's
>   real `sf_raw_slct` lists — 34,197 and 39,452 sites at `extnd=954` — run solo, then both from one
>   manifest. **6 of 6 products byte-identical**: clip FASTQ 189,654,477 / 194,726,108 bytes, disc
>   FASTA 242,849,015 / 258,527,335, disc-pos 92,938,486 / 98,907,461. No part files left behind.
>   The unit tests are a 40-site chr1 slice; the union plan, the per-entry fan and the per-entry
>   mate-gather striping only meet thousands of merged regions and ~1 GB FASTAs here.
> - **chr1 end-to-end, job `49516538`** (`refactor_gate_chr1.sbatch`): the production single-entry
>   path now runs through `run_td_gather_multi`, so this is the leg that matters most. **9/9
>   identical on the legacy arm and 9/9 on the native arm** vs `stage1_gate_chr1_49453032`, call
>   counts unchanged (L1=1 Alu=14 SVA=2 merged=17).
>
> **The CPU number:** solo L1 50s + solo SVA 53s = 103s of gather, versus **56s** for both together.
>
> **No production caller yet, deliberately** — hoisting the gather to the driver is the throwaway
> half described above. `--td-manifest` is what the per-locus loop will call. Anyone measuring this
> before then will find nothing: the two gathers are already concurrent processes, so the saving is
> CPU, not wall clock, until the downstream is one job at full width.
>
> **A gate-writing note that cost a submission** (`49516776`, FAIL): `local t=$1 d="$OUT/solo_$t"`
> expands before `local` assigns, so `$t` was a leaked global from an earlier loop; one arm
> overwrote the other's output and the gate blamed the code. `refactor_gate_chr1.sbatch` already
> documents this trap from job `49452897`. The counter that kept it from being a silent pass was
> `products compared: N (expected 6)` — copy that into any gate comparing a variable number of
> things.

> ### 2c's GATE BASELINE EXISTS 2026-08-02 — `1c56b8d`
>
> "Calls may move, but only as family-assignment flips" cannot be checked against a population
> nobody counted: a 2c run that drops a call for an unrelated reason looks exactly like one that
> reassigned it. `locus_table.LocusTable` now records every collision it resolves — winner, loser,
> and whether load-order **priority** or a **`two_side_tprt_both` override** decided it — and
> `run_xtea` writes `merged/locus_family_conflicts.tsv` and logs the summary. Diagnostic only: the
> merged hc_pf list and the type map are byte-identical whether or not the sidecar is written, and
> a test asserts that directly.
>
> **MEASURED on BASELINE-0 (HG002 WGS) 2026-08-02 — 22 conflicts, and this IS Stage 2c's gate:**
>
> | winner | loser | rule | n |
> |---|---|---|---|
> | Alu | L1 | `tprt_override` | **12** |
> | SVA | Alu | `priority` | 6 |
> | L1 | SVA | `priority` | 2 |
> | L1 | Alu | `priority` | 1 |
> | Alu | SVA | `tprt_override` | 1 |
>
> It reconciles exactly: 200 L1 + 1230 Alu + 85 SVA per-type calls − 22 = **1493** merged, the
> BASELINE-0 count, at 188 / 1223 / 82. The 12 `Alu over L1` rows are precisely the 12 the Stage 0
> gate found by hand when comparing `--types L1` (200) against `--types L1,Alu,SVA` (188) — so that
> anecdote is now the whole population, per rule, and the winner and loser positions are within
> 41 bp of each other (median 0; 12 of the 22 are at the identical position).
>
> **So Stage 2c's gate has a number in advance rather than after the fact: at most ~22 loci are
> currently contended, and a 2c run that moves substantially more than that — or that moves a locus
> not on this list — is not doing family reassignment, it is doing something else.**

> ### 2c.1/2c.2/2c.3 BUILT 2026-08-02 — and the premise this section was written on is wrong
>
> **`scan_clip_realign.tsv` is not per-family evidence.** The decomposition above says "that table
> already carries `chrom pos family nL nR mate_in_rep` — three rows per position, i.e. per-family
> support side by side… a read of that plus the disc side, not a new computation." The three rows
> exist and each comes from its own bwa (`realign.hpp` runs one per family), but **the numbers they
> carry are the same number**. Measured on BASELINE-0 WGS:
>
> | | |
> |---|---|
> | positions with clip support in all three families | 8,621,063 of 8,622,517 |
> | position-pairs whose `nL/nR/mate_in_rep` differ across families | **1,147 of 17,244,556** |
> | positions with ≥5 clip support | 155,615, of which **504** differ across families |
> | resulting clip candidate lists | 577,293 / 577,312 / 577,347 — the same positions |
>
> So there is nothing to decide family on at the clip stage. **Family-discriminating evidence first
> appears on the disc side**, where the mate is tested against each family's own annotation:
> `candidate_list_from_disc.txt` is 2,244 L1 / 9,468 Alu / 1,183 SVA, union **12,223**, with only
> **672** rows (5.2%) claimed by more than one family. That file — not the realign TSV — is
> `cns_filter`'s actual input and the real starting point for a locus evidence table.
>
> **How far up front the decision can go is fixed by one fact:** the rule's `two_side_tprt_both`
> override — which decides 13 of the 22 conflicts — is written **by `cns_filter`**
> (`x_clip_disc_filter.py:1845`). Nothing earlier has it. So `cns_filter`'s output is the earliest
> file the existing rule can read at all, and deciding there is a *relocation*; deciding at the
> candidate list would be a different rule, adopted without evidence to support it.
>
> **Do not expect throughput from Stage 2c.** Family contention is small at every level: 672 of
> 12,895 disc candidates, 16 of 2,796 cns rows, 22 of 1,515 postF_hc calls. Partitioning the work by
> family saves 0.6–5% of the downstream, and after Stage 1 the downstream is only ~21% of the run
> (`downstream_parallel` 457s of 2,166s; `scan_pass` ~1,500s). **The reason to do 2c is that family
> stops being decided three chains later by whichever list a call happened to survive on** — not
> speed. Anyone who picks this up expecting a speedup will conclude, correctly, that it isn't worth
> doing, and will be answering the wrong question.
>
> **What was built:**
>
> - **2c.2 (the lift).** `locus_table.load_rows` / `LocusTable.assign_rows` — the rule now takes
>   pre-parsed claims, so one copy of it serves both ends of the pipeline. `assign_type` is defined
>   as `load_rows + assign_rows`, and a test asserts the two agree, because "moved, not rewritten" is
>   the entire safety argument.
> - **2c.3 (the relocation).** `locus_table.assign_and_filter` resolves family over the three
>   `candidate_disc_filtered_cns.txt` and rewrites each to its own winners (originals kept at
>   `.pre_family_assign`). `run_xtea.downstream_waves` splits the per-family chain into
>   `cns_filter` → **barrier** → `transduction…postF_hc`, using `build_type_steps`' existing
>   `start_from`/`stop_after`, so the argv cannot drift from the single-wave chain.
>   `--family-assign-early`, **default OFF**, no-op for a single `--types` entry.
>
> **GATE RUN — job `49502115` (after `49498881` died on the companion-file bug below). The result
> is ZERO moved calls, and the reason matters more than the number.**
>
> | | OFF | ON |
> |---|---|---|
> | merged VCF | 1493 | **1493** |
> | per family (L1 / Alu / SVA) | 188 / 1223 / 82 | **188 / 1223 / 82** |
> | call set `(chrom, pos, family)` | — | **byte-identical to OFF** |
> | vs frozen BASELINE-0 VCFs | **4/4 identical** | — |
> | conflicts resolved at the barrier | — | 16 |
> | conflicts resolved at the merge | 22 | **15** |
>
> **The barrier does not REPLACE the merge — both run.** So it can change output only where it
> removes a row that would have *won* at the merge, and it never does. Of its 16 decisions, **7 are
> also seen by the merge and all 7 produce the identical winner, loser and rule** (verified locus by
> locus); the other 9 the merge never sees. 22 − 7 = 15, exactly the ON arm's merge count.
>
> **The 15 the merge still resolves are not reachable from cns at all.** The losing family has no
> cns row there because transduction and sibling create its record *after* the barrier has run — 12
> of the 15 are the `Alu over L1 tprt_override` population. **So Stage 2c cannot end with one
> decision point until site creation is also unified.** That is a Stage 3 dependency this section
> did not previously know about: the barrier is a second decision point, not a relocated one, and
> the merge's resolution cannot simply be deleted.
>
> **A previous version of this block predicted 4 moved loci, all toward truth. That prediction was
> wrong and `probe_early_assign.py` is fixed.** Two methodological errors, both worth recognising
> because they are easy to repeat:
>
> - it compared the two resolutions as **alternatives** when both actually run;
> - it matched each late locus to the **nearest early locus of any family**, so a late L1 call at a
>   position where L1 had no cns row got paired with a nearby Alu or SVA row and reported as a
>   family flip. The barrier left L1 at **300 of 300 rows** — it never took anything from L1.
>
> The probe now compares the two **conflict populations** instead, and reproduces the real run
> exactly: 7 seen by both with 0 disagreements, 9 barrier-only, 15 merge-only.
>
> **What the gate is worth even at zero movement:** the relocation is now *proven* answer-preserving
> on real WGS data rather than argued. That is what would license removing the merge's resolution
> later — the two decision points agree everywhere they overlap.
>
> **The bug the first run found**, and the reason a chr1 gate would not have: a step can open an
> input it was never handed. See §10's entry; `SVA:transduction` died with a bare
> `KeyError: 31329020` because the barrier filtered `candidate_disc_filtered_cns.txt` and left the
> matching row in `…​.high_confident`, which `x_TEA_main.py:551` derives from its own `-i`.
>
> ### SITE CREATION MEASURED 2026-08-02 — the Stage 3 dependency now has a number
>
> §2.4 named unified site creation as a Stage 3 dependency but never sized it.
> `docs/handoff/runs/probe_site_creation.py` does, offline, from the gate's kept intermediates
> (no BAM, no run). Where each family's final postF_hc record **first exists**, ±50 bp one-to-one:
>
> | | postF_hc | cns | cns2 (transduction) | sibling |
> |---|---|---|---|---|
> | L1 | 200 | 150 | **50** | 0 |
> | Alu | 1230 | 1230 | 0 | 0 |
> | SVA | 85 | 83 | **2** | 0 |
>
> **52 of 1,515 per-family final records do not exist when the barrier runs** — all L1 or SVA, all
> created by the **transduction** step, **none by sibling**. (Sibling's own output file contributes
> nothing directly; the postF_hc chain reads `candidate_disc_filtered_cns2.txt.high_confident`.)
>
> And on the conflict population: **exactly 7 of the 22 have both sides present at cns.** The other
> 15 have a winner or a loser transduction had not yet invented — 12 loser-invented (the
> `Alu over L1 tprt_override` block), 2 winner-invented, 1 both.
>
> **That 7 is the gate's own number, derived from the opposite direction.** The gate counted which
> conflicts the barrier *actually* resolved (7 both, 0 disagreements); this counts which ones it
> *could have* and gets the same 7. So **the ceiling on any decision point placed before
> transduction is 7 of 22**, and it is a site-creation limit, not a rule or a window one. Widening
> `MERGE_WINDOW`, moving the barrier, or changing the priority order cannot raise it.
>
> **DECIDED 2026-08-02 — `--family-assign-early` stays OFF, and this is not a "revisit later".**
> The flag is proven answer-preserving and buys nothing on its own, because it is a *second*
> decision point rather than a relocated one. The only case for flipping it is that it is a
> prerequisite for deleting the merge's resolution — so **flip it in the same change that makes it
> load-bearing (unified site creation), not before.** Turning it on earlier costs a barrier and
> keeps the merge. Recorded at the flag itself in `run_xtea.py` so it is not re-litigated from the
> help text alone.
>
> **Still open in 2c**, deliberately not built: the single per-locus loop that replaces the three
> subprocess chains, and `--emit-per-type`. The chains are still three, and there are still two
> decision points. **The per-locus loop is the item to know about**: §8 lists it under *Rewrite*
> alongside `downstream_unified.py` / `transduction_evidence.py`, which are Stage 3 modules, and
> §2.4's finding adds unified site creation to its dependencies. It is not a Stage 2-sized change
> that happens to be left over; it is the front edge of Stage 3, and sequencing it as such is the
> honest reading of what this stage measured.

> ### SITE ENUMERATION IS SOLVED 2026-08-03 — `xtea/site_universe.py`, and the 52 come from a file
> ### that already exists at the barrier
>
> "Unified site creation" was a named dependency with a size (52 of 1,515) and no mechanism. Two
> offline measurements on the 2c gate's OFF arm give it one, and the first kills the obvious
> alternative reading of the probe:
>
> - **The 52 are created, not moved.** Every postF_hc position is present at **distance 0** in
>   `candidate_disc_filtered_cns2.txt` — the post-cns chain does not re-anchor coordinates at all —
>   and the nearest cns row to each of the 52 is **≥ 3,293 bp away, median 2.4 Mb**. So the probe's
>   ±50 window was not hiding movement; there is none to hide. (Worth stating because the probe's own
>   docstring justifies that window by "the chain shifts a call's coordinate", which is true of
>   `postF` re-anchoring but not of this hop.)
> - **All 52 are already rows in `candidate_list_from_disc.txt.clip_sites_raw_disc.txt`** — the
>   raw-disc file `disc_adapter.run_disc_from_scan` writes **before `cns_filter` runs** — and all 52
>   survive into `<raw_disc>.slct`, the selection `re_slct_with_clip_raw_disc_sites` derives from it.
>
> So site creation is **a selection over a file that already exists at the barrier**, needing no BAM,
> no bwa and no cns evidence — only raw_disc, the cns candidate list (proximity exclusion) and the
> family's rmsk annotation. The rule is now `xtea/site_universe.py`
> (`write_slct_file` / `select_new_sites`), with `XTransduction.re_slct_with_clip_raw_disc_sites`
> reduced to a delegation, so one copy serves the step and any driver-side caller.
>
> **GATE: PASS, 4/4 byte-identical, offline** — `docs/handoff/runs/gate_site_universe.py`, no BAM and
> no sbatch, ~2 min. It reruns the rule with the arguments the 2c run itself used (parsed out of that
> run's `*.steps.log`, not retyped) and byte-compares against the files that run produced, for **both**
> production callers, which differ in every argument that could hide a bug:
>
> | | rows | bytes |
> |---|---|---|
> | L1 `.slct` (transduction caller: 300 candidates + 33,897 new) | 34,197 | **identical**, 793,811 |
> | L1 `.slct2` (sibling caller: `b_with_candidate=False`, extend = mean IS) | 33,927 | **identical**, 790,127 |
> | SVA `.slct` (292 + 39,160) | 39,452 | **identical**, 915,815 |
> | SVA `.slct2` | 39,158 | **identical**, 911,479 |
>
> The cutoff is the trap the gate removes: `b_automatic` is True by default, so the step ignores
> `--nd 5` and uses `Parameters.get_par_by_cov(31.816)[1]` = **4**. The gate derives it from the run's
> own cached basic-info file for that reason. Unit side: `tests/test_site_universe.py` (+11), whose
> centre is a randomized differential test against a **verbatim copy of the deleted loop** — which is
> how the half-open proximity window below became a verified property rather than a reading.
>
> **THE LOCUS UNIVERSE HAS A NUMBER: 41,604.** Union over the three families of (cns candidates ∪
> new raw-disc sites), exact position: 41,604 loci — 7,910 claimed by one family, 33,139 by two, 555
> by all three — where the three per-family chains between them visit 75,853. **1.82×.** That is the
> per-locus loop's iteration domain, and nobody had it; the 53,945-of-54K figure earlier in this
> section is about the td gather's fetch *windows*, a different thing.
>
> **WHAT THIS DOES NOT CHANGE — read before hoping.** Knowing a site exists is not knowing which
> family claims it. `.slct` is a **34-40K-row superset per family**, and the per-family evidence at
> those rows only comes into being when the transduction evaluation runs. A decision point placed
> before that evaluation still cannot resolve the 15 remaining conflicts, so **the 7-of-22 ceiling is
> unchanged and `--family-assign-early` still stays OFF.** What moves is the *iteration domain*: the
> per-locus loop can enumerate its loci before gathering any evidence, which is what "one pass over
> loci" in §3 requires.
>
> **AND IT CORRECTS §3's ORDERING.** The diagram puts *family assignment* before the per-locus loop.
> For 52 of 1,515 records there is nothing to assign on at that point — their only evidence is the
> transduction evaluation's own output. The order the data supports is **enumerate loci → gather
> per-family evidence → assign family → filter**, with assignment inside the loop's tail rather than
> ahead of it. Designing the loop off §3's picture would build the barrier that §2.4 already measured
> as unable to close.
>
> **The named caller this unblocks:** 2b.4 above says hoisting the td gather to the driver would mean
> "re-implementing step logic the per-locus loop deletes", and names
> `re_slct_with_clip_raw_disc_sites` as that logic. It is no longer step logic. ~~What the driver still
> lacks for a `--td-manifest` caller is the step's `iextnd` (from `automatic_gnrt_parameters`).~~
> **DISCHARGED 2026-08-03 — the driver already had it.** The step's `iextnd` is
> `max(400, int(mean_is + 3*std_var))` (`x_TEA_main.py:473`, `:497-499`), and `run_xtea`'s `cns_extnd`
> is that second term already; the 400 floor never binds at WGS insert sizes.
> `downstream_unified.transduction_extnd()` is the one copy. **So the `--td-manifest` driver caller has
> no missing input left**: site list, extnd, output paths.

> ### THE ENUMERATION NOW RUNS IN THE DRIVER 2026-08-03 — `xtea/downstream_unified.py`, `--locus-universe`
>
> The block above ends with a rule that *could* be called for every family up front and no caller
> doing it. There is one now, and building it moved three things.
>
> **1. WHERE it can run was not free, and the docstrings had it wrong.** The rule excludes raw-disc
> rows near a cns candidate, and that candidate list is `cns_filter`'s **output** — so the universe is
> computable before **transduction**, not before `cns_filter`. Both `site_universe.py` and
> `x_transduction.py` said "available before cns_filter runs"; corrected in the tree, and in §10.
> The consequence is convenient rather than costly: the earliest legal point is exactly where the
> Stage 2c barrier already sits, so `--locus-universe` reuses `downstream_waves`' split. It does
> **not** imply `--family-assign-early`, and it does not inherit that flag's single-family exemption.
> **Where both are on, the universe runs first** — the assignment barrier rewrites each family's cns
> list to its own winners, so the other order would enumerate a post-assignment universe and change
> what the number means without saying so.
>
> **2. What it is worth today: a byte comparison instead of an argument.** The flag is diagnostic on
> purpose. Each family's transduction step still derives its own `.slct` through the same one copy of
> the rule, so the driver writes a *replica* to `merged/locus_universe/` and no call can move — and
> after any run, `cmp <type>/…raw_disc.txt.slct merged/locus_universe/<type>.slct` must be identical.
>
> **GATES, both PASS.**
>
> - **WGS scale, offline, ~5s** — `docs/handoff/runs/gate_locus_universe.py`. It differs from
>   `gate_site_universe.py` in the one way that matters: that gate hands the rule **the arguments a
>   run's own step log records**, deliberately, so it cannot drift into testing the driver. This one
>   builds every argument **the driver's way** — `REP_TYPES` layout, coverage-derived cutoff, the
>   step's own 450 boundary extend, 225 min copy — and **fails on any mismatch with the run's**, then
>   byte-compares. Both halves are needed: agreeing bytes with a disagreeing annotation path would
>   mean the rule is insensitive to it *on this data*, not that the driver is right. Result: 4/4
>   arguments match for all three families; L1 34,197 rows / 793,811 B and SVA 39,452 / 915,815 B
>   byte-identical; Alu correctly has no step `.slct` to compare.
> - **chr1 end to end** — `gate_locus_universe_chr1.sbatch`, job `49547574`: ON vs OFF, every VCF /
>   cns2 / sibling list identical, plus the replica comparison with a **flanked-family** counter (2,
>   not 3). The wiring is what the offline gate cannot see, and the wave split is where a real chain
>   fails — job `49498881` died at `SVA:transduction` for exactly that reason.
>
> **The trap that gate exists to catch:** transduction is handed `annotation1`
> (`{g}_L1.fa.out`), `cns_filter` is handed `annotation`
> (`{g}_L1_larger_500_with_all_L1HS.out`). Both exist, both look like an L1 annotation, and the wrong
> one changes the rule's third rejection.
>
> **3. A number that reframes the 7-of-22 ceiling: 38,815 of the 41,604 loci are reachable only by
> creation** — no family has a cns row there. Quote it as a fact about the iteration **domain**, never
> as a call count (the evaluation rejects almost all of them; only 52 become records). What it says is
> that the created population *dominates* the universe rather than being a tail, so the ceiling is
> structural rather than a small-numbers accident.
>
> **AND A CORRECTION TO §3's MIDDLE BOX, found by reading the step it replaces.**
> `call_candidate_transduction_v3` builds one flank reference and runs **four batched bwa passes over
> the whole `.slct` list** before any per-site record exists, so "per-locus loop" cannot mean per-locus
> end to end. The shape is *enumerate → batched per-family evidence production → per-locus decision*,
> only the third box a loop. Recorded at the diagram. The loop's record type also already exists:
> `collect_ins_info`'s 11 fields plus `_collect_gntp_features`' 10, and a **created** site's record is
> the same shape as any other's — which answers the open question "what does per-family evidence look
> like at a locus only `.slct` knows about" without a probe.
>
> **Still open, and now the whole of it:** the batched evidence stage (the `--td-manifest` caller,
> which has no missing input left) and the per-locus decision loop.

### Stage 3 — The evidence-model rewrite (truth-set gated)

Agreement-gating stops being coherent here, because the decisions themselves change. See §5.
**The truth-set harness must be working before this stage starts.**

> ### THE DECISION IS NOW A RECORD AND NAMED GATES 2026-08-03 — `xtea/transduction_evidence.py`
>
> §5's first sentence is built: the transduction step's veto cascade is one evidence record feeding
> named gates, with every threshold a parameter and every drop a reason. **Answer-preserving** —
> the thresholds did not move, only their location, so this is the *structural* half of Stage 3 and
> it is still agreement-gated. Changing a number is the next half and that one is truth-set gated.
>
> **`decide()` IS FIVE GATE GROUPS, NOT ONE FUNCTION, AND THAT IS THE DATA FLOW.** §5 draws one
> call because it assumes a locus's evidence is all available at once. It is not, and no refactor
> makes it so: **two of the five groups sit behind batched BAM fetches over sets the earlier groups
> define.**
>
> ```
> clip SAM parse -> [A] clip     ----.
> disc SAM parse -> [B] source ---> [C] candidate
>                                     |  (ReadDepth over the SURVIVING candidates only)
>                                     v [D] sanity
>                                     |  (the genotyper over the SURVIVING sites only)
>                                     v [E] genotype
> ```
>
> `filter_transduction`'s `calc_coverage_of_two_regions` and `_collect_gntp_features`' genotyper each
> open the BAM for the set handed to them. Flattening this would compute coverage and genotype
> features for the whole `.slct` list — 34-40K loci per family instead of a few thousand. So **"one
> decision function" means one place per threshold, not one call**, and anyone designing the
> per-locus loop off §5's picture should read this first.
>
> **WHAT IT BUYS TODAY: the reason ledger.** `<out>.decision_ledger.tsv` per flanked family, one row
> per dropped locus with `phase`, `reason` and the numbers that decided it. §5 says `reason` "is what
> makes the benchmark interpretable, because every missed truth call gets attributed to the specific
> gate that dropped it" — **BASELINE-0's 322 false negatives have no such attribution today**, and
> the information existed only as `print()` calls in the step's stdout. Diagnostic: nothing reads the
> file, and the 14 reason strings are the module's only public contract (rename one and a join
> against an older run's ledger silently misses).
>
> **The 14 reasons at this commit — 18 after gate [B2], see the next block** — by phase, and four of
> them are §5's own "interrogate" items:
> `A_clip`: `low_clip_support`, `no_polya_clip`, `low_polya_clip`.
> `B_source`: `decoy_cross_family` (§5's three-way decoy tally), `low_disc_support`,
> `no_dominant_source` (the 0.45 dominance and 0.5 uniqueness ratios, both kept).
> `C_candidate`: `no_clip_evidence`, `no_disc_polya`, `no_disc_cluster`.
> `D_sanity`: `same_family_repeat`, `coverage_island`, `two_sided_polya`, `two_sided_polyt` (§5's
> "blunt tandem-dup filter — make conditional on source strength").
> `E_genotype`: `large_indel_excess`.
> **Added 2026-08-04, `B2_cns_cluster`**: `no_cns_disc_reads` (an ABSENCE, not a threshold),
> `no_cns_disc_cluster`, `low_two_sided_cns_disc`, `low_one_sided_cns_disc` (unreachable today).
>
> **TWO DEFECTS ARE PRESERVED ON PURPOSE**, marked `DEFECT-PRESERVED` in the module and pinned by a
> test that fails if either is "fixed":
>
> 1. **The cluster test sits inside the per-read loop** (`dominant_side_cluster`). It runs once per
>    read, and only the last read's evaluation survives for the side it picked while the other side
>    keeps whatever an earlier iteration left there — so a locus can carry `b_l_trsdct` from read 3
>    next to boundaries from read 9, and "both sides true" is reachable **only** through that stale
>    flag. `collect_ins_info` reads exactly that pair to decide which side's cns interval to blank,
>    so hoisting the test out of the loop MOVES CALLS. It is a Stage 3 change with a truth-set
>    measurement behind it, not a cleanup.
> 2. **`_pick_trsdct_rescued`'s second return value was both discarded and wrong** — its
>    accumulation sat inside `if ins_chrm not in m_n_source_picked`, so it recorded only the first
>    locus per chromosome. Not rebuilt; the ledger's `no_dominant_source` rows are that population,
>    correctly.
>
> **A DOCUMENTATION DEFECT FIXED IN PASSING, worth knowing because both windows are gated on:**
> `filter_transduction` passes `focal_win = LOCAL_COV_WIN` (**200**) and `focal_win2 =
> COV_ISD_CHK_WIN` (**900**), and its comments described them the other way round. So `i_lcov1/i_rcov1`
> are the 200 bp window, `i_lcov2/i_rcov2` the 900 bp one, all four are tested against `i_max_cov`,
> and **cov2 is what `collect_ins_info` reports as the site's local coverage**. Only the comment
> changed — swapping the constants would move calls.
>
> **GATED PASS AT BOTH SCALES on one commit (`48ec4c7`), three legs**: chr1 before/after **98/98
> byte-identical** (job `49634289`), the NA12878 golden downstream gate **3/3 families** (job
> `49634291`), and WGS **4/4 frozen BASELINE-0 VCFs identical** with 188/1223/82/1493 calls and
> `downstream_parallel` 463.3s against the 462.7 / 511.75s of the two previous no-flag WGS runs (job
> `49634290`). **The refactor's cost is not detectable at WGS.**
>
> **WHERE THIS DECISION ACTUALLY SPENDS ITS REJECTIONS (WGS, job `49634290`) — this is new, and it
> reprioritises §5's "interrogate" list:**
>
> | phase | L1 | SVA |
> |---|---|---|
> | A_clip | **20,332** (`low_clip_support` 14,276, `no_polya_clip` 6,056) | **27,590** (16,286, 11,304) |
> | B_source | 676 (`low_disc_support` 424, `decoy_cross_family` 144, `no_dominant_source` 108) | 105 (92, 0, 13) |
> | C_candidate | 13 | 41 |
> | D_sanity | **1** | **2** |
> | E_genotype | **1** | **1** |
>
> **Gate A does 97% of the rejecting; gates D and E reject one or two loci each in a whole genome.**
> So the two-sided polyA/polyT veto and the large-indel sanity check — two of §5's four interrogate
> items — cannot be worth meaningful sensitivity at this population, and the leverage is gate A's
> clip-support and polyA-signature cutoffs. Also: `decoy_cross_family` fires **144 times for L1 and 0
> for SVA**, the opposite of what §5's SVA-is-SINE-VNTR-Alu argument predicts. Look at that before
> touching the decoy tally.
>
> **AND THE SENSITIVITY CEILING IS NOT HERE. `benchmark/ledger_attribution.py` on that run: 306 of
> BASELINE-0's 322 false negatives (95%) have NO ledger row — they never reached this decision.** All
> 16 attributed misses are in gate A; none in B/C/D/E. With the tool's own base-rate control (a drop
> lands within 200 bp of a *true positive* 12.1% of the time for L1), the L1 enrichment is 2.17× and
> the **excess over chance is ~8 calls**; SVA's 1.18× is nothing. **Quote the excess, never the raw
> 16.** The missing 322 are upstream, in `cns_filter` and what reaches it — so Stage 3's threshold work
> on transduction is a *precision* story, not a sensitivity one, and that is a change of expectation
> worth having before the numbers move.
>
> **GATES.** See §4 Stage 3's gate block above for the truth-set gate that the *next* half needs;
> this half's are agreement gates, and the chr1 one is the strongest in the repo: **every file both
> arms produce must be byte-identical**, ~98 of them, not just the VCFs. That threshold is available
> because a dry run over an earlier gate's two arms (job `49599837`) established that 98 files under
> the sample directory are already identical between two arms of the *same* commit, and the only two
> that are not are characterised — `cns_manifest.txt` (embeds each arm's own output paths) and
> `scan_disc_pos.tsv` (a **permutation** only; the multithreaded scan writes rows in completion
> order, so the gate compares it sorted).

> ### AND IT IS SIX GATE GROUPS, NOT FIVE — GATE [B2] 2026-08-04
>
> `check_seprt_disc_consistency`'s decision is now
> `transduction_evidence.gate_cns_disc_cluster` (commit `ff7cde1`), which closes the largest hole
> in the attribution: it drops **135 of L1's 159 picked sources and 49 of SVA's 105** at WGS, it is
> where §5's `NEARBY_CLIP=50` lives, and until now those loci were indistinguishable in the ledger
> from loci that never reached the transduction decision at all.
>
> **[B2] SITS BEHIND A THIRD BATCHED PASS THE DIAGRAM ABOVE OMITTED.** `re_select_disc_reads`
> writes only the disc reads of [B]'s survivors and they are realigned to the CONSENSUS before
> [B2] can see them — the same barrier argument as the two BAM fetches, one stage earlier. So the
> shape is **six ordered gate groups with three batched passes between them**:
>
> ```
> clip SAM parse -> [A] clip ---------------------------------.
> disc SAM parse -> [B] source                                |
>                        |  (B's survivors' disc reads         |
>                        |   realigned to the CONSENSUS)       |
>                        v [B2] cns cluster -----> [C] candidate
>                                                      |  (ReadDepth over survivors)
>                                                      v [D] sanity
>                                                      |  (the genotyper over survivors)
>                                                      v [E] genotype
> ```
>
> **[B2] IS THE ONE GATE WHOSE PREDICATE DID NOT MOVE, and the reason is a constraint, not a
> compromise.** `check_seprt_disc_consistency` lives in `x_clip_disc_filter.py` and `cns_filter`
> calls it twice (`:1483`, `:1655`), so its body stays where both callers reach it and only its
> decision moved. **The ledger is caller-supplied**: the transduction caller passes one, cns_filter's
> two pass nothing, a None ledger records nothing. That is what keeps the instrumentation off a path
> the refactor deliberately did not touch.
>
> **FOUR REASONS, AND THE FIRST IS NOT A THRESHOLD.** `no_cns_disc_reads` is an ABSENCE — no disc
> read of that locus survived the third bwa pass — and §4.7's rule holds: an absence and a threshold
> are never added together. `no_cns_disc_cluster` (neither side clusters) and
> `low_two_sided_cns_disc` (both cluster, total support below the cutoff) are the live thresholds.
> `low_one_sided_cns_disc` is the `b_transduct=True` arm, **dead in all three callers** and
> implemented anyway, so turning it on later gets its attribution for free.
>
> **PHASE ORDER IS NOW SEMANTIC, NOT LEXICAL, AND THIS WAS A LIVE TRAP.**
> `"B2_cns_cluster" < "B_source"` in ASCII (`'2'` is 0x32, `'_'` is 0x5F), and
> `ledger_attribution.attribute` took `min(rows)` on the phase STRING to find the gate that
> determined a false negative. That was correct only while A < B < C < D < E held by accident; with
> B2 present it would have attributed misses to B2 that gate B rejected first — a wrong attribution
> that looks entirely plausible in the output table. `transduction_evidence.PHASE_ORDER` /
> `phase_rank` make it explicit and every consumer goes through them. Verified answer-preserving:
> the tool reproduces job `49634290`'s numbers exactly (306 / 8 / 8 = 322) after the change.
>
> **ITS GATES, and the chr1 one is stronger here than its parent was.** Both arms write ledgers now,
> so the ledger cannot be excluded as a new file — it is the one artifact this change may alter, and
> the check is exact: **AFTER's ledger with every B2 row deleted must be byte-identical to
> BEFORE's**, which asserts in one comparison that rows were added and that nothing else moved, not
> a row, not a reason, not an order. chr1's B2 domain is 5 loci for L1 and **empty for SVA**, so that
> leg proves wiring at N=5 on one family; the WGS leg carries a **prediction made before the run** —
> the site-file derivation said 135/49, and the ledger must now say the same, from a different
> source. `docs/handoff/runs/gate_cns_cluster_ledger_{chr1,wgs}.sbatch`, plus
> `gate_golden_downstream.sbatch`, which is the NA12878 golden downstream gate finally committed as
> a script rather than retyped per session.

> ### THE FIRST THRESHOLD-CHANGING ITEM IS IN: THE SVA DECOY FIX 2026-08-04 — `55055ba`
>
> Not a threshold, strictly: a **classification**, and the first change in Stage 3 that moves calls
> rather than relocating a decision. §10's `decoy_cross_family` entry is the finding; this is what was
> done with it.
>
> **`transduction_evidence.decoy_tags(name, b_structural)` + `is_flank_record`** replace a substring
> test that was written inline **twice** (`x_transduction.py:1160` primary hit, `:1285` XA tags), and a
> third copy — exact equality — survives in the dead `parse_clip_algnmt_rescued` (`:1051`), which is
> left alone because it is one of §8's zero-caller functions. **A flank id is
> `chrom~start~end~subfamily~rc` with integer coordinates in fields 2 and 3** (polymorphic flanks
> included: `gnrt_flank_regions_of_polymerphic_insertions` writes `start == end`) and a consensus id
> (`LINE1`, `ALU_AluYb8`, `SVA_D`) has none of that. It fails **closed**: enough tildes but unparseable
> coordinates is still classified by the tags.
>
> ~~**IT SHIPS OFF, AND THE DEFAULT IS THE FROZEN BEHAVIOUR RATHER THAN THE CORRECT ONE.**~~
> **FLIPPED ON 2026-08-05 (`3d1437a`); the paragraph below describes the pre-flip state and is kept
> because its reasoning is what the flip had to satisfy.** It did: `b_structural=False` is still the
> pre-fix expression verbatim and still pinned, the "gate both arms then flip on the numbers" shape was
> followed exactly, and the three jobs are in the 2026-08-05 progress block.
>
> `run_xtea --td-decoy-structural` → `XTEA_TD_DECOY_STRUCTURAL`, resolved once per step because the
> step is its own `x_TEA_main` subprocess. `b_structural=False` is the pre-fix expression **verbatim**,
> pinned in `tests/test_transduction_evidence.py` against a copy of it — anything shipping off has to
> be provably inert or its default is not a rollback. This is the `--td-gather-scan` shape: gate both
> arms, then flip on the numbers. Flipping it means **re-freezing BASELINE-0**, and every gate under
> `docs/handoff/runs/` that compares against it has to restate its expectations in the same commit.
>
> **THREE THINGS THE FLIP CHANGED THAT ARE NOT "the default value", and they are the parts to know:**
>   1. **`run_xtea` is now the SINGLE WRITER of the policy** and sets the env var on every path, `"1"`
>      or `"0"`. Pre-flip it wrote only on opt-in so an inherited `"0"` could not fake an opt-in; with
>      the fix on by default the risk reverses — an inherited value could fake an opt-**out** and
>      silently restore the defect.
>   2. **The env var's own rule is UNCHANGED: absent still means OFF** (`x_transduction.py:78`, and the
>      pre-existing test is untouched, which is evidence rather than assertion). The default is a
>      **driver** policy with a blast radius of one entry point, because **a dozen frozen harnesses
>      under `docs/handoff/runs/` invoke `x_TEA_main` directly** with baselines measured under the
>      substring test. If absent meant ON, one flip would silently move all of them. **Do not "fix"
>      this asymmetry without re-measuring those baselines.**
>   3. **BASELINE-0 was NOT re-frozen; `--no-td-decoy-structural` is what reproduces it**, and the chr1
>      gate proves that byte-for-byte over the whole output tree (100 identical / 2 tolerated / **0
>      differing**). So the paragraph's "flipping means re-freezing BASELINE-0" was avoided rather than
>      satisfied: comparisons against BASELINE-0 become before/after, with the opt-out as the "before".
>      **The `AFTER_OFF` arm of the chr1 gate now has to pass `--no-td-decoy-structural`** — left alone
>      it would have been `AFTER_ON` under another name and the gate would have failed on a difference
>      it caused itself.
>
> **L1 IS UNTOUCHED BY CONSTRUCTION, AND THAT IS LOAD-BEARING.** `TD_DECOY_LINE` is `"LINE1"`, not
> `"L1"`, so no L1 flank id collides (the probe measures 0 in both L1 libraries) and Alu never runs
> this step. So the WGS gate's check 1 is **two-sided**: L1 and Alu byte-identical to BASELINE-0 *with
> the fix on*, SVA not — a moved L1 file falsifies the premise, an unmoved SVA file means the flag
> reached nothing. L1's 21,158 ledger rows are gated too, which is finer than its VCF.
>
> **WHAT IT DOES NOT DO:** the hg38 SVA flank library still carries no decoy records, so SVA's arm of
> `gate_decoy_cross_family` still has nothing to count and §5's SVA/Alu hazard is **still untested**.
> That half is now scoped exactly — see §5, and note that `consensus/transduction_decoy.fa` has zero
> references anywhere.
>
> **It also unblocks the third arm of the Stage 3 transduction gate's candidate pool** (below), which
> is "xTea2 with the SVA decoy test made structural" and which has no other source of reference-copy
> SVA sources.
>
> **Its gates:** `docs/handoff/runs/gate_decoy_structural_{chr1,wgs}.sbatch`. chr1 runs **three arms** —
> BEFORE, AFTER_OFF (byte-identical over the whole tree: the seam claim the unit differential cannot
> make) and AFTER_ON (wiring, not population: chr1's gate B has 5 L1 survivors and **0 SVA**). WGS runs
> one arm with the flag ON against frozen BASELINE-0, whose check 5 is the payoff — the probe's SVA
> reference-copy picked sources, **0 of 105** today on every build and on master.
>
> **AND ITS WGS EFFECT WAS PREDICTED BEFORE THE RUN, from the pre-fix run's own SAM.** The flag
> changes only how a hit's reference name is classified — after the alignment, and without touching
> the bwa reference or the reads — so job `49703322`'s `temp_transduction_disc.sam` is the same file
> the fixed run produces, and replaying the shipping method over it with the flag both ways is the
> fix's exact effect on gate [B] from a different source (`docs/handoff/runs/probe_decoy_fix_delta.py`;
> the OFF arm reproduces the run's picked sources and full gate-[B] ledger split on both families,
> which is what licenses reading the ON arm):
>
> **SVA picked sources 105 → 346, of which reference copies 0 → 241; orphan loci 26 → 87; gate [B]'s
> SVA domain 210 → 1,432 loci. L1 identical to the last row, `decoy_cross_family` included.** What it
> cannot say is what the run is for: everything after gate [B] is untouched by a replay, so 346 is an
> input to the cascade, not a call count.
>
> **WGS RESULT at `11b0e61`: PASS on both arms** (`49724055` 20 cores, `49724056` 24 cores — and the two
> are byte-identical to each other, which retires the core-count caveat their scheduling introduced).
> **L1 and Alu are byte-identical to BASELINE-0 WITH THE FIX ON**, L1's ledger identical at 21,158 rows
> and by reason; SVA gains **3 calls** (82 → 85, three positions, none lost) and its ledger goes 27,787
> → 29,001. **The payoff: SVA reference-copy picked sources 241 of 346, from 0 of 105.**
> `downstream_parallel` 494.9s at 24 cores, inside the no-flag band. And **every number of the offline
> prediction landed exactly** — 346 picked, 241 reference, 1,086 gate-[B] drops, 1,432 domain, 87 orphans.
>
> **THE COST IS ONE CALL, NOT THREE, AND THE BENEFIT IS NOT MEASURABLE ON HG002.** All three new calls
> are `one_side_and_half_transduction` with a `TD_SRC` naming a reference SVA copy — the fix working —
> and **two of the three are orphan SITES** (`~orphan`, `SVLEN=0`), which are not MEI insertions and
> which this truth set therefore cannot contain, exactly as §4 Stage 3 measured for all 37 of L1's.
> `benchmark_truth_set.py`'s orphan stratum misses them because it keys on
> `subtype == "orphan_or_sibling_transduction"` and these carry a different `SUBTYPE`; fix that key
> before quoting the 93.9% → 90.6%. SVA *sensitivity* is unchanged (85.6%), all-families sensitivity
> unchanged (80.4%), and one false negative moves from `no_ledger_row` to `low_disc_support` — a locus
> that previously had no ledger row at all because its hits were discarded before any gate saw them.
>
> **SO THE DEFAULT STAYS OFF FOR A REASON THAT IS NOT CAUTION.** What the 241 newly reachable sources
> are *worth* is exactly what §7.1's exhausted ALT route cannot say (denominator 2-3), and the
> instrument that can is the Stage 3 transduction gate below — which needs this flag **on** to build its
> third pool arm. Build `benchmark/verify_transduction.py`, score the 241, then flip; and when flipping,
> re-freeze BASELINE-0 and restate every gate's expectations in the same commit.
>
> **chr1 RESULT at `11b0e61`, job `49721312`: PASS.** Walls 7:51 / 5:59 / 5:57. AFTER_OFF is
> **100 identical + 2 characterised, 0 differing, 0 missing** against BEFORE. With the fix ON, exactly
> **one** artifact in the whole tree moves and it is **SVA's decision ledger** — L1 0 files, Alu 0
> files, no call anywhere (1 / 14 / 2 / 17 in all three arms). And the ledger delta is the fix arriving:
> SVA's rows go **2,281 → 2,315**, all of it new gate-**B** decisions on loci whose hits used to be
> discarded at the decoy branch before any gate saw them — **`low_disc_support` 2 → 35, plus one
> `no_dominant_source`**. So at chr1 the fix reaches **34 loci and every one is rejected for thin
> evidence**, which is what a 5-loci-of-L1, 0-loci-of-SVA gate-B domain predicts. The population is
> WGS's. NA12878 golden downstream is **PASS** at the same commit (job `49721733`), exactly the four
> documented L1 carry-forward notes and the Alu genotype flips.

**Gate:** per-family sensitivity/PPV against the HG002 truth set, no regression beyond a margin
calibrated on BASELINE-0; BASELINE-0 delta reported for attribution but **not** as pass/fail.

> ### THE TRANSDUCTION HALF OF THAT GATE — SPEC, NOT YET AGREED 2026-08-04
>
> **Full spec: `docs/handoff/HANDOFF_2026-08-04_stage3_td_gate_spec.md`.** The line above works for the
> call set as a whole (BASELINE-0 scores 80.4% / 88.4%) and **does not work for transductions**, because
> §7.1's truth-side denominator is **0 for L1 and 2-3 for SVA** and that route is exhausted. Corinne's
> direction (2026-08-04) supplies the missing input: **truth from a run of master xTea on HG002.**
>
> **Master enumerates the hypotheses; the reads return the verdict; the metric is precision per stratum.**
> That distinction is load-bearing rather than pedantic — scoring against master's *calls* would
> contradict §2's locked decision (*"agreement with master is explicitly not a requirement"*), and it
> would **reject the SVA decoy fix**, because master shares that bug (`master:xtea/x_transduction.py:1127`,
> same substring test, same `hg38_FL_SVA_flanks_3k.fa`) and its own SVA output carries **zero
> reference-copy sources** — measured 2026-08-04 on master's NA12878 WGS run (job `48825977`): **92
> picked SVA sources, 0 reference copies**, against L1's 135 with 84. Using master only as a candidate
> pool avoids both.
>
> **And the third arm now exists**: the SVA decoy test is structural behind `--td-decoy-structural`
> (`55055ba`), so that arm is a flag rather than a branch to write.
>
> **Why precision is the right target and not a consolation** — measured, not assumed: 306 of 322 FNs
> never reach this step, all 16 attributed misses are gate [A] with an excess over chance of ~8 calls, and
> gate [B2] attributes **zero** (job `49703322`). Sensitivity lives upstream in `cns_filter`.
>
> **AND THE THIRD ARM IS NOW A MEASURED QUANTITY, NOT A HOPE: 241 reference-copy SVA sources** (WGS
> `49724055`/`49724056`, `--td-decoy-structural`), against **0 of 105** in BASELINE-0 (same sample) and **0 of 92** in master's own
> NA12878 output (different sample, same signature). That is the pool contribution this gate could not have without the fix, and it is ~2.6× the
> 92 polymorphic-only SVA sources master's NA12878 run supplies.
>
> **The pool is THREE arms** — master-on-HG002 ∪ BASELINE-0 ∪ **xTea2 with the SVA decoy test made
> structural**. The third is not optional: without it the pool has no reference-copy SVA sources at all.
> Measured on master's existing NA12878 WGS run (job `48825977`): **115 transduction-asserting calls**
> (86 L1, 29 SVA, 0 Alu) against the ALT route's 0-3. ~~HG002's pool size is **unmeasured**.~~
> **MEASURED 2026-08-04: HG002's three-arm pool is 115 entries** (master on HG002 contributes 104
> transduction-asserting calls — 72 L1 + 32 SVA + 0 Alu), frozen with sha256s. Note **86 of master's
> NA12878 L1 assertions are the BARE shape**, so the "86 L1" above is disc spans, not source intervals
> — see §10.
>
> **`docs/handoff/runs/master_serial_wgs.sbatch` already works** — repoint `BAM`/`SAMPLE` from NA12878,
> and run **master's own pipeline end to end**, because §10 records that master's transduction figures
> come in two non-interchangeable flavours and only the own-pipeline one samples master's call set.
>
> **The deliverable is `benchmark/verify_transduction.py`**, not the run: per pooled `(locus, source)`
> pair, decide from the BAM whether reads at the locus carry sequence from `[Y_end, Y_end + flklen]` that
> does **not** align to the family consensus over ≥80%, rejecting source windows that sit inside a
> same-family repeat — **which is exactly the trap §7.1 hit**, where 11 of 20 L1 residuals were divergent
> element body rather than transduced sequence. Three verdicts; `UNDECIDABLE` never folded into either.
>
> ~~**Open, and listed in the spec's §8**~~ — **ALL THREE NOW ANSWERED, two of them fully:**
>   * **HG002's pool size: 115** (three arms, ±200 bp one-to-one). Also measured: **82 of 115 asserted
>     by all three arms**, 8 by two, 25 by one; master contributes 21 unique, the decoy arm 4,
>     **BASELINE-0 zero**.
>   * **The margin: absolute VERIFIED counts per stratum, and PROVISIONAL** (Corinne, 2026-08-05) — not
>     a percentage tolerance, because the strata that matter hold 1-4 entries where one reclassification
>     swings a percentage by 100 points. Provisional until an arm exists that is **neither master nor
>     xTea2**: three-arm agreement is a floor on independence, not a measure of it, because two of the
>     three arms differ by one flag. Implemented as `verify_transduction.py --baseline`.
>   * **The sibling `POS` question: LEANING, NOT SETTLED.** The verifier searches both sides and reports
>     which supported, so the question is now answerable from output rather than argument. On the 3-arm
>     HG002 pool: **D (downstream, `[Y_end, Y_end+W]`) 518 supporting reads vs U (its upstream mirror)
>     298** — about 64/36. Leaning downstream, but **not the lopsided split that would settle it**, and
>     by the verifier's own stated rule a near-even split means the window is catching something
>     orientation-independent and the *claim* should be loosened rather than the parameter tuned. Treat
>     as open, and note it is per-read rather than per-entry, which is the weaker unit.
>
> **And one thing the spec's §7 step 3 asked for DOES NOT EXIST**: its designated validation target,
> "the 33 L1 reference-copy sources master reports on NA12878", are the BARE shape — disc-pair
> midpoints. The verifier was validated instead on the stratum that does exist, **26 of 31 two-sided
> entries confirmed (83.9%)**. §10.

### Stage 4 — Somatic, mosaic, tail, legacy deletion

- Fold case-control into the new parameter model (§6.3). **No somatic regression gate exists today** —
  building one is part of this stage and is a precondition for touching the tumor thresholds.
- Land mosaic as a parameter profile plus one scan field (§6.4); gate on the SMHT005-3G truth set (§7.2).
- Collapse the tail invocations; delete the legacy per-type path.

## 5. The evidence model (Stage 3)

Replace the six-phase veto cascade with **one explicit per-site evidence record feeding one decision
function**, so every threshold is a named parameter rather than a buried constant.

```
SiteEvidence:
  clip:      per-side counts, polyA/polyT counts and ratios
  disc:      per-side dominant source, unique vs non-unique counts, anchor positions
  source:    family, source id, cross-family conflict flags
  cluster:   per-side consensus cluster extents
  sanity:    local coverage, large-indel ratio, repeat-context divergence

decide(evidence, params) -> (call | None, reason)
```

`reason` is not decoration — it is what makes the benchmark interpretable, because every missed truth
call gets attributed to the specific gate that dropped it.

**Load-bearing, keep:** polyA/polyT signature (the transcriptional hallmark), disc uniqueness
(MAPQ ≥ 25), source dominance ratio (0.45), unique/non-unique ratio (0.5), repeat-divergence guard (5%),
coverage and large-indel sanity.

**Interrogate — likely implementation artifacts:** `NEARBY_CLIP=50` (arbitrary proxy — replace with a
measured clip-to-insertion distance distribution); `ileftmost_cns` 500-vs-5000 (unexplained 10×
asymmetry — replace with an explicit full-length criterion); the two-sided polyA/polyT veto (blunt
tandem-dup filter — make conditional on source strength); the three-way decoy tally.

> **THREE OF THOSE FOUR ARE NOW MEASURED NEAR-INERT — reprioritise before touching any of them, and
> do not spend an item on the first three.** The two-sided polyA/polyT veto rejects **1 locus (L1) and
> 2 (SVA) in a whole genome** and the large-indel check 1 each (§4 Stage 3's WGS table), so neither can
> be worth meaningful sensitivity.
>
> **AND `NEARBY_CLIP=50` IS THE THIRD, measured 2026-08-04 by sweeping the shipping gate itself**
> (`docs/handoff/runs/probe_nearby_clip.py`, offline, seconds — see §10 on the SAM that was thought
> deleted). Swept **0 → ∞** through `check_seprt_disc_consistency` over job `49703322`'s own inputs,
> gate [B2]'s verdicts change for **exactly one locus in the whole genome**, and only at ∞: L1 keeps
> 24 survivors at every value, SVA keeps 56 at every value up to 500 and 57 above the largest observed
> distance. The probe reproduces the run's own ledger B2 split exactly at 50 (L1 113 + 22, SVA 33 + 16),
> which is what makes the sweep a measurement.
>
> **This list asked for "a measured clip-to-insertion distance distribution", and here it is — it says
> the proxy is FINE.** Over the disc records clipped on exactly one side (64 of L1's 296, 178 of SVA's
> 505) the distance `|anchor_pos - ins_pos|` is sharply **bimodal**: L1 median 3 with p75 **411**, SVA
> median 129 with p75 358 and p90 787. 50 sits in the empty gap between the modes, so moving it has
> almost nothing to move. (`prep_representative_pos` is the identity — `clip_pos` IS the insertion
> position, `x_transduction.py:710` — which is why this needs only the SAM's read names.)
>
> ~~**So the only one of the four still worth interrogating is `ileftmost_cns` 500-vs-5000**~~ —
> **MEASURED 2026-08-05, AND THE PREMISE ABOVE IS WRONG: THE LIVE 10× ASYMMETRY IS ZERO.** This list
> called it "unexplained", but the code explains it and says it is *dormant*.
> `family_params.py:18` — "`leftmost_cns` was therefore **always 500** — SVA's 5000 branch";
> `x_TEA_main.py:311` resolves the bundle from `from_sva_flag(options.sva)`, and
> `run_xtea.build_type_steps` passes `--sva` to **cns_filter alone** (matching the golden
> `run_xTEA_pipeline.sh`), so the transduction step resolves `GENERIC` = **500**. **Both families run
> 500 today.** The question was never "why do they differ" but "should SVA's dormant 5000 be woken",
> and waking it is a documented one-word change.
>
> Swept through the **shipping** `pick_full_length_insertion`, identical across three independent
> finished runs (`docs/handoff/runs/probe_leftmost_cns.py`):
>
> | `ileftmost_cns` | L1 (of 300) | SVA (of 292) |
> |---|---|---|
> | 500 | **58 · 19.3%** ← live | **162 · 55.5%** ← live |
> | 2000 | 69 · 23.0% | 289 · 99.0% |
> | 5000 | 167 · 55.7% | **292 · 100.0%** ← dormant canonical |
> | 10000 | 300 · 100.0% | 292 · 100.0% |
>
> `family_params`' claim — a 5000 threshold "admits essentially the whole candidate list" — is
> **CONFIRMED, and is now a measurement rather than an inference from a sequence length**: waking it
> takes SVA 162 → 292, **+130 candidates (+80%)** into the polymorphic-flank library and every
> realignment built from it. The consensus-length mechanism checks out independently — `SVA.fa`'s
> longest record is **3183 bp** (< 5000, so SVA saturates), `LINE1.fa`'s is **6342 bp** (> 5000, so L1
> does not, reaching 100% only at 10000).
>
> **NEW, and this is the part worth an item: at the value actually in use the two families are nowhere
> near equivalent — L1 19.3% against SVA 55.5%.** "Both run 500" is not "both are filtered the same";
> SVA's cns clusters sit much closer to the 5′ end. **Explain that before deciding anything about
> waking 5000.** Whichever way it goes, it needs a truth-set-scored two-arm gate, and §7.1's obstacle
> applies: the derived transduction denominator is 0 (L1) and 2 (SVA), so only SVA *position*
> sensitivity/precision is scorable. `family_params`' own words: "not to flip blind".
>
> **So all four of this list's items are now measured and none is an open lever.** The leverage overall
> is **gate A**, which does 97% of the rejecting. Expect *precision*, not sensitivity: 95% of
> BASELINE-0's misses never reach this decision.

> ### GATE A IS NOW MEASURED TOO — 2026-08-05, job `49836568`, and it took 57 SECONDS
>
> The semantic half was filed against gate [A] through four handoffs as though it needed a
> pipeline run. It did not: **`tmp/transduction/temp_transduction_clip.sam` survives a finished
> run**, so gate [A]'s entire input — the clip realignments *and* the polyA map, which the same
> parse builds — is on disk. `docs/handoff/runs/probe_clip_gate.py` drives the shipping parse
> (`x_transduction.py:281`) and the shipping gate (`:283`), and reproduces the run twice over: the
> ledger's `A_clip` reason split **and** the run's own `.unique_trsdct_half_clip` file as a *set*
> of loci.
>
> **1. `low_polya_clip` IS UNREACHABLE, so the polyA THRESHOLD decides nothing.** The polyA record
> is created and incremented in the same block (`x_transduction.py:929-975`), so every entry that
> exists has `sum(rcd[:6]) >= 1` — measured: **0 zero-sum entries** of 4,604 (L1) and 5,261 (SVA).
> At the coverage-derived cutoff of 3 the derived `n_polya_cutoff` is 0.5, so the comparison needs
> a sum of zero. The reason first fires at cutoff **5**. This list keeps the polyA signature as
> "load-bearing" and it *is* — but as a **membership test**, and the named parameter is dead.
>
> **2. Membership, not clip depth, is what binds gate [A].** Of loci with `n_clip >= 3`: L1 **6,056
> dropped for polyA-map absence against 2,802 kept** (68%), SVA **11,304 against 3,233** (78%). And
> 594 / 810 of the dropped carry `n_clip >= 10`, one of them 969 — strong clip support discarded for
> the absence of a single polyA-ish read. The kept population's `n_clip` is only modestly higher
> (median 5 vs 4), so this is not clip depth under another name.
>
> **3. Gate [A] has TWO thresholds and ONE dial — a §5-shaped defect in the parameterisation.**
> `n_polya_cutoff = n_clip_cutoff/2 - 1` is derived inside `TransductionParams`, so the sweep is
> single-variable only up to 4; at 5 the polyA threshold wakes and rejects 286 (L1) / 364 (SVA)
> loci the shipping configuration keeps. **"Raise the clip cutoff" and "raise the polyA cutoff"
> are not separable changes today**, which is exactly the "every threshold a named parameter"
> ambition failing at the one gate that matters most.
>
> **4. The decoy-consensus asymmetry reproduces from the CLIP side** — independent corroboration
> for the packaging item below. The clip path does **not** skip decoy hits (the skip is commented
> out, `x_transduction.py:872-875`), so L1 counts **3,992 of 99,245** qualified clip records
> (4.02%) as support from `>LINE1`/`>SVA`/`>ALU` consensus records, on 529 of its 2,802 survivors;
> SVA counts **0 of 128,262**, because hg38's SVA flank library has no decoy records concatenated
> in. Same family asymmetry the disc side's `decoy_cross_family` shows (144 L1 / 0 SVA), reached by
> a different route. **No survivor rests on decoy support alone**, so it is not manufacturing loci
> — it is a family-asymmetric evidence channel, and the concatenation would make it symmetric.
>
> **AND IT CORRECTS `NEARBY_CLIP`'s "measured inert" verdict by scoping it — see §10.** That
> measurement was inside gate [B2]; the same constant filters this parse at `x_transduction.py:865`
> and is worth **3.7–3.9× gate [A]'s survivors** unbounded. Not a recommendation to move it.

**Unify sibling / orphan / regular transduction into one call class with a typed source.** These are
epistemological categories, not biologically distinct events — the evidence shape is identical (clip
polyA signature + one-sided disc source, or two-sided anchors on a common region). Model as
`source ∈ {reference_element, polymorphic_element, inferred_region}`, one decision function for all
three, distinction surviving only as an output label. This collapses `x_orphan_transduction.py` into the
main path.

**Alignment passes: target two, not three.** Today: clip→flank library, disc→flank library,
disc→consensus. ~~The flank library already has consensus concatenated in
(`x_transduction.py:186`)~~ — **that is a COMMENT, not behaviour, and it is false for SVA on hg38
(measured 2026-08-04).** `construct_novel_flanks` concatenates the family flank library with the
polymorphic flanks it just built and nothing else; whether any decoy consensus is present is purely a
property of how that build's library was packaged (both L1 libraries: 45 records; hg38 SVA: none). And
`rep_lib/consensus/transduction_decoy.fa` — 3 records, exactly the `>LINE1`/`>SVA`/`>ALU` the libraries
carry — has **zero references in xTea2 or in master**. So the third pass's marginal contribution is
position refinement plus a cross-family decoy tally that is dead on the SVA/hg38 path; measure whether
it changes any call before keeping it, and do not size that measurement off the comment.

> **On merging the per-family flank libraries into one reference:** a design pass argued this is safe
> because the families are evolutionarily distinct. **That reasoning is wrong for SVA** — SVA is
> literally SINE-VNTR-**Alu**, carrying Alu-derived sequence, so SVA/Alu mapping competition is a real
> hazard. The existing three-way decoy tally (`x_transduction.py:1103-1124`), which explicitly counts
> wrong-family hits, is evidence the original authors hit exactly this. Treat library merging as an
> **experiment with a measured outcome**, not a refactor, and run it only after the truth-set gate exists.
>
> > **AND THE DECOY TALLY DOES NOT CURRENTLY TEST THIS FOR SVA — read §10 before using it as
> > evidence either way.** The decoy test is a *substring* match on the flank record's name
> > (`TD_DECOY_SVA = "SVA"`), the SVA flank library has no decoy records concatenated in at all, and
> > every one of its 5,454 reference flank records is named `…~SVA_D~0L` and so is misread as the SVA
> > consensus and skipped. Measured: **0 of SVA's 105 picked sources at WGS are reference copies; all
> > 105 are polymorphic.** So the hazard above is untested, and the decoy tally's SVA arm is dead
> > code rather than a measurement of zero competition.
> >
> > **HALF OF THAT IS FIXED (2026-08-04) AND THE OTHER HALF NOW HAS AN EXACT RECIPE.** The substring
> > test is now `transduction_evidence.decoy_tags(name, b_structural)`, and `--td-decoy-structural`
> > makes reference SVA sources reachable — call-moving, so it ships off; see §10. That is the *code*
> > half. The *packaging* half — an hg38 SVA flank library with no decoy records, so SVA's arm of
> > `gate_decoy_cross_family` has nothing to count — is what leaves this hazard untested, and it is
> > now fully determined rather than open:
> >
> > **The decoy set is 45 records and every one of them is reconstructible from `rep_lib/consensus/`
> > on any build, byte-identically** (measured, offline): `transduction_decoy.fa`'s 3 (`>LINE1`,
> > `>SVA`, `>ALU`) + 35 `>ALU_Alu*` whose sequences match `consensus/ALU.fa` + 6 `>SVA_[A-F]` from
> > `consensus/SVA.fa` + `>HERV-K` from `consensus/HERV.fa`. Both L1 libraries carry all 45 and the
> > hg19 SVA library carries 44 of them (no HERV-K); hg38's SVA and chm13's SVA and L1 carry none.
> > So making the hazard testable needs **no library repackaging and no build-specific data** — it is
> > a concatenation in `construct_novel_flanks`, which is where this file's own comment
> > (`x_transduction.py:207`) already claims it happens. It is still CALL-MOVING and still needs its
> > own arm.

## 6. Cross-cutting design requirements

### 6.1 Kill the name `col5`, and fix the bug behind it

`col5` names a **column position**, not a concept, and leaks that non-name into the C++, the Python, the
CLI (`--col5-families`), the logs, and the spec. What it actually is:

> **the number of clipped reads at a site whose mate lands inside any selected repeat family's
> annotation** — a repeat-anchored support signal. Master called it `mate_within_rep_cnt` and computed it
> Alu-only (non-deterministically — whichever type won a shared-`pub_clip` race). xTea2 deliberately
> changed it to a union over families, each family's boundary extend baked into its BED (Alu/L1 = 0,
> **SVA = 200**), emitted identically for every family at a site. **Consumed only by the SVA cns path**
> (`IS_CALL_SVA` branches in `x_clip_disc_filter.py`).

The concept *is* documented: `docs/clipdisc_contract_spec.md` §4 **DC1** is thorough, including the
measured justification (at golden SVA sites 111639406 and 168451690 the Alu-only count is 0 while 8 and
6 reads have mates in SVA). Only the *name* fails to connect.

**Rename to `mate_in_rep_clip_cnt`**, CLI `--mate-in-rep-families` (keep `--col5-families` as a
deprecated alias), across `src/csrc/xtea_scan.cpp` (`clip_matein`, the `--l1-bed/--alu-bed/--sva-bed`
help), `src/csrc/realign.hpp:316,363-369`, `xtea/scan_adapter.py`, `xtea/run_xtea.py:701-715,829-848`,
and the spec. Document it in the evidence record, not as "the 5th column of a file."

> **DONE 2026-08-01.** Two notes for anyone reading the old vocabulary:
> - **The C++ identifiers were already right** — `clip_matein`, `mate_in_rep()`, `matein_merged`,
>   `rep_idx`. Only comments, the `--help` label, and one local in `realign.hpp` actually said `col5`,
>   so the C++ side was a comment/help change, not an identifier rename.
> - **`col5` deliberately survives in three places**, all of them explaining the retired name (two
>   comments plus the deprecated flag's help/warning). `--col5-families` still works — several committed
>   `docs/handoff/runs/*.sbatch` pass it, and those are historical run records that should not be
>   rewritten. Passing both flags with different values is a hard error; passing the alias warns.

**The bug found while checking this.** `run_xtea.py:837` iterates **requested** types and skips those
outside `--col5-families`, so the BED set is `requested_types ∩ col5_families`:

| invocation | BEDs generated | `mate_in_rep` value |
|---|---|---|
| `--types L1,Alu,SVA` (default families `alu,sva`) | ALU + SVA | mate-in-(Alu ∪ SVA) |
| `--types SVA` | SVA only | mate-in-**SVA only** |
| `--types L1` | none (`rep_any` false) | **0 at every site** |

This value gates the SVA cns consistency check, and union ≥ subset always, so **`--types SVA` can emit
fewer SVA calls than the same sample run as `--types L1,Alu,SVA`.** A live correctness bug, not a
redesign artifact: asking for less must not change what you get for what you did ask for.

**Fix:** the contributing family set is a property of the *signal definition*, not of which families were
requested. Build those BEDs from `REP_TYPES` for every family in `--mate-in-rep-families` regardless of
`--types` (annotation paths resolve independently of the requested set — small change). Add the
`--types X` vs `--types L1,Alu,SVA` comparison as a permanent regression test.

> ### MEASURED 2026-08-01 — the fix is correct, and OUTPUT-NEUTRAL on HG002
>
> This section predicted the fix's signature as "three-type output does not move at all; **SVA-only
> output changes**". The second half **did not happen on HG002**, and that is recorded here rather
> than quietly dropped.
>
> The two runs demonstrably computed different things (jobs `49441922` pre-fix, `49441920` post-fix):
>
> | | BEDs built | positions with nonzero `mate_in_rep` | sum | SVA calls |
> |---|---|---|---|---|
> | pre-fix `--types SVA` | SVA only (5,827 iv) | 61,548 | 83,445 | **85** |
> | post-fix `--types SVA` | ALU (1,238,490 iv) + SVA | 1,493,594 | 1,774,091 | **85** |
>
> The value changed at **1,434,938 of 8,622,501** positions — a 21× increase in total signal —
> including **38 of the 85 called SVA sites**. And the SVA call set is **byte-identical**.
>
> **Both facts matter.** The fix is right: `post >= pre` at *every one* of the 8.6M positions (0
> violations), which is the union-≥-subset property the argument rests on, now measured rather than
> asserted. And it removes a real type-set dependence. But on this sample the SVA cns consistency
> check was **not the binding constraint** at any called site, so nothing flipped.
>
> **Consequences:**
> - Do not describe this fix as "changed nothing" — the signal it feeds moved by 21×.
> - Do not describe it as "recovered N SVA calls" either. On HG002 it recovered none.
> - Whether it *ever* moves a call is sample-dependent and still unmeasured elsewhere. The direction
>   is at least pinned: monotone upward, so it can only ever **add** SVA calls, never remove them.
> - `stage0_gate_check.sh` gate B was rewritten to test the **mechanism** (does `mate_in_rep` change,
>   monotonically) rather than call movement. The genuine failure mode is the mechanism *not* moving —
>   that would mean the pre/post distinction never took effect, which would also make gate A's
>   "identical" result meaningless.

### 6.2 Running a single type

The efficiency case rests partly on cross-family sharing, so single-type behaviour is an explicit design
requirement, not an afterthought:

| source of gain | `L1,Alu,SVA` | `L1` only |
|---|---|---|
| shared read gather across families (~1.9× fewer gathers) | **yes** | no — nothing to share |
| one family-assignment pass instead of three | **yes** | degenerate (L1 or nothing) |
| kill per-site temp files + per-chromosome `Pool` churn | yes | **yes** — per-type overhead |
| collapse/stream the bwa passes and SAM intermediates | yes | **yes** |
| unified tail invoked once | yes | already once |

**L1-only gets meaningfully faster, but not by the ~1.9× headline** — never quote that number for a
single-type run.

Requirements, each needing a test:
- The unified pass must **degrade cleanly to N=1 family**: no path may require all three families'
  evidence to emit a call; family assignment with one candidate reduces to "L1 or nothing", not an empty
  competition.
- **`--types` must partition cleanly** — a family's call set must not depend on which other families were
  requested. §6.1 is the one known violation.
- **Alu-only** must skip transduction entirely rather than run it empty (`REP_TYPES["Alu"]["flank"] is
  None` — a non-autonomous SINE with no internal Pol II promoter cannot read through; correct, not a gap).

### 6.3 Somatic / case-control

Fits the locus-major model *better* than germline does, but has two wrinkles and no safety net.

**Today:** two steps per type after `postF_all` (`run_xtea.py:386-406`). `case_control` re-examines the
**control** BAM list at the tumor's called loci, filtering the raw `postF_all` list; `case_control_hc`
(`--somatic_hc`) intersects that with the type's own `hc_pf` list. Both take **no `--rtype`** and no
`--tumor/--purity` — they are **already family-agnostic**. `type_stop_after` becomes `case_control_hc`,
and the shared tail keys off that instead of `postF_hc`.

**Where tumor-specific behaviour actually lives** — upstream, and it must survive the parameter-bundle
refactor: the raw-disc retention threshold (`disc_adapter.py:218-220`, `cutoff/2` tumor vs
`int(cutoff*3/4)+1` germline — **this feeds transduction**, so getting it wrong silently moves
transduction/orphan calls), and `--purity`/`--tumor` through cns_filter/transduction/sibling, with
`--tumor` alone through both postF steps.

**Consequences:**
1. **The unified gather must cover the control BAM(s) too**, or the somatic stage re-opens BAMs per locus
   and reintroduces the pattern being deleted. Reuse `cns_emit.hpp`'s per-file `ftag` machinery
   (`:536-537`, stripped on emit at `:855`); do not invent a second mechanism.
2. **The `hc_pf` intersection must be re-expressed over the unified call set** — the per-type `hc_pf`
   intermediate disappears, so "intersect with this type's HC list" becomes "filter the unified HC set"
   with family as a field. A simplification, not a port.
3. Somatic *gains* from unification: one control-BAM pass over the merged locus set instead of three.

**Gap:** no somatic regression gate exists. Build one in Stage 4. **Not verified: whether a
case-control test pair exists on scratch** — if none does, source or simulate one before Stage 4.

### 6.4 Mosaic mode

**Already planned — do not re-plan.** `docs/handoff/MOSAIC_PLAN.md` (164 lines, committed `46cba5b`,
2026-07-20) is the prepared design, backed by two code-reading research passes. Read it first.
`xtea/x_mosaic_calling.py` exists (9 KB, `MosaicCaller(XPostFilter)`) and `--postFmosaic` is wired
(`x_TEA_main.py:114-116, 726-729`). Upstream has an `xtea_mosaic` branch. What is *missing* is the clip
step: master's `--mosaic` clip path was deleted in trim step 4, and `x_TEA_main.py:708-713` already
records the intended replacement — a `mosaic_adapter.py` reading a polyA tally from `xtea-scan`.

| what mosaic needs | how it lands here |
|---|---|
| per-clipped-read polyA tally, accumulated per position | **one extra cheap per-read fact in the scan pass we are already extending** |
| `PEAK_WINDOW` 100 → 30 bp | a parameter, not a code path — **already shared with `--somatic`**, in scope |
| AF cutoffs, hardcoded `(0.35, 0.02)` per type (`x_mosaic_calling.py:173-179`) | exactly what the parameter bundle is for; AF becomes an evidence field |
| polyA-dominant-both-sides rejection, same-family repeat rejection, blacklist dedup | in-memory per-locus filters, **zero BAM re-reads** |

So mosaic is **a parameter profile plus one scan field, not a fourth pipeline** — true *only because* the
rewrite turns scattered constants into parameters, which makes now the cheapest moment and later
disproportionately expensive.

Stated priority is **speed at 300×** ("prohibitively slow on 300x data… memory can scale up… I can give
many cores"). `MOSAIC_PLAN.md`'s cost model: the bottleneck is not the polyA check (trivial) but doing
*any* per-read Python work at 10× depth — the exact problem `xtea-scan` solved. This rewrite does that
folding anyway, so mosaic gets the win for free.

`--joint` (cross-sample singleton filter, `x_joint_calling.py`) stays **out of scope** — a post-hoc set
operation on finished per-sample lists, not a statistical joint caller; addable later.

## 7. Data and validation

### 7.1 Germline: HG002

**NA12878 is retired.** Cost, stated so it is a choice: the frozen NA12878 comparators (master 50:18,
port 54:39, transduction attribution, chr1 goldens) become historical context, not gates. Since the
rewrite abandons master agreement anyway that is mostly notional, but two things must be rebuilt —
**BASELINE-0 on HG002** (Stage 1's byte-diff becomes "HG002 before vs after", still a real gate) and a
**fast iteration loop** (an HG002 chr1 slice; without it every gate is a ~50 min job and iteration dies).
**The HG002 chr1 slice EXISTS as of 2026-08-01** (job `49440267`):
`test_data/HG002/chr1/HG002.chr1.bam` + `.bai`, 3.4 GB, 61,165,483 reads, `quickcheck`-clean.
No flag filtering was applied and the 3,366-contig `@SQ` is unpruned, both so the slice keeps
predicting the full run rather than quietly becoming an easier problem.

**Truth sets** — `parklab/xTea_paper/benchmark/HG002/`, "coordinate based on hg38" per the repo README,
matching our `GRCh38_no_alt`: `L1_all.vcf`, `L1_tprt.vcf`, `Alu_all.vcf`, `Alu_tprt.vcf`,
`SVA_tprt_all.vcf`. Verified characteristics:

- **Assembly-derived, phased.** ALT carries the full inserted sequence; INFO is
  `QNAME=NA24385#S10-H1-000000F;QSTART=…;QSTRAND=±` (NA24385 = HG002) with a trailing haplotype column
  (`H1`, or `H1,H2`).
- **Headerless** — no `##` or `#CHROM` lines. The parser must not assume a well-formed VCF header.
- Record counts, **measured 2026-08-01** on the staged files (this corrects an earlier "`L1_tprt.vcf`
  has ~26 records" in this section, which was simply wrong):

  | | `_all` | `_tprt` | tprt positions absent from `_all` | `_all` positions absent from `_tprt` |
  |---|---|---|---|---|
  | L1 | 197 | 184 | **0** | 13 |
  | Alu | 1355 | 1257 | **0** | 98 |
  | SVA | — | 90 (`SVA_tprt_all.vcf`) | | |

  197 + 1355 + 90 = **1642**, which matches the benchmark size the paper states, so the `_all` sets
  (plus SVA's single combined file) are the paper's full benchmark.

**RESOLVED 2026-08-01 — `_tprt` is a canonical-mechanism filter, not "transduction-positive."**
It is a strict *subset* of `_all`: every `_tprt` position appears in `_all` and none is new, with
byte-identical INFO structure (`QNAME`/`QSTART`/`QSTRAND`) and an insertion-size distribution
indistinguishable from `_all`'s (L1 mean 2604 vs 2647 bp). A transduction-positive set would be a
small minority carrying *extra* non-repeat sequence, so it would skew longer and would not nest this
cleanly. The paper's HG002 benchmark was manually curated by verifying "a polyA tail and target-site
duplication/deletion" per insertion — the molecular hallmarks of TPRT — and xTea's own
`two_side_tprt_both` (`locus_table.py`, `x_post_filter.py:756-786`) uses TPRT the same way: both
breakpoints resolved and consistent. Reading `_tprt` as "transduction-positive" is refuted.
Source: Chu *et al.* 2021, *Nat Commun* 12:3836 (PMC8219666).

**Consequence for scoring:** the sensitivity **denominator is the `_all` set** (197 L1 / 1355 Alu /
90 SVA). `_tprt` is usable as a secondary, higher-confidence stratum — reporting sensitivity on both
is nearly free and separates "we miss hard calls" from "we miss well-formed calls" — but it is not the
headline denominator.

**Also: these truth sets carry no explicit transduction/source-element annotation** — **verified
2026-08-01**, not assumed: no `transduct*`, `TD=`, `SOURCE=`, `ORPHAN` or `SIBLING` token appears
anywhere in any of the five VCFs, and `_tprt` carries no field `_all` lacks. So transduction status
**must** be derived from the ALT sequence (the non-repeat co-mobilized portion) or taken from the
supplement; there is no annotation to fall back on. Transductions are the point of the rewrite, so
**settle how they are scored before Stage 3** — a call with the right position but the wrong source
element is a different error from a missed call, and position-only matching cannot see the difference.

#### How transductions get scored — design settled 2026-08-01, implementation still to do

**What xTea already emits.** The VCF carries `TD_SRC` (`x_gvcf.py:148,385`, from record field 23) and
`SUBTYPE`. Measured on the NA12878 WGS VCF, `TD_SRC` takes exactly three shapes:

| shape | example | n |
|---|---|---|
| `not_transduction` | — | 1438 |
| source **region** | `chr8:51817625-51818695`; start may be unknown → `chrX:-1-45277389` | ~60 |
| sibling **point** | `chrY:9591830~sibling` | ~46 |

`SUBTYPE` separately distinguishes `orphan_or_sibling_transduction` (46),
`one_side_and_half_transduction` (3) and the `two_side_tprt*` / `one_half_side*` non-transduction
classes. **Parsing hazard:** the region delimiter is `-` and an unknown start is literally `-1`, so
`chrX:-1-45277389` must not be split naively on `-`.

**Deriving truth transduction status — feasibility VERIFIED, and then MEASURED, and the two do not
agree.** Aligning all 197 `L1_all.vcf` ALT sequences against `consensus/LINE1.fa` (24,702 bp over 4
subfamily entries) with `bwa mem` and measuring the primary alignment's soft-clip residual:

| non-L1 residual | count of 197 | 3′ ONLY, supplementary hits excluded |
|---|---|---|
| median | 18 bp (polyA + TSD — no transduction) | **5 bp** |
| ≥ 100 bp | **60** | **20** |
| ≥ 200 bp | **51** | **12** |
| max | 4620 bp | 2959 bp |

> **The left column is TOTAL soft clip, both ends, and it is not what this section's own criterion
> asks for.** Reproduced exactly by `benchmark/benchmark_truth_set.py` (17 / 60 / 51 / 4619). The
> right column is the 3′ residual the criterion actually specifies, with supplementary alignments
> counted as explained consensus sequence rather than as residual. The leading clip is the LARGER of
> the two halves — 40 records ≥100 bp, max 1745 bp — and it is 5′ inversion or genomic flank in the
> ALT, exactly the confounders named two paragraphs below as the reason to require 3′. **Nothing
> sized off the 60 is right.** SVA's figures are median 11, ≥100 bp: 14, ≥200 bp: 11, max 929.

So ~10% — not ~26-30% — of truth L1 insertions carry substantial non-consensus sequence at the 3′
end, and those residuals are long enough (100-2959 bp) to attempt a mapping back to the reference
and yield a **source locus derived independently of xTea**. That is the truth-side transduction
annotation the VCFs lack, and §4 Stage 2's progress block records how far it actually gets.

**Residual is a candidate, not a transduction.** A soft-clip residual can also be a 5′ inversion (a
known L1 structure), an assembly artefact, or plain genomic flank included in the ALT. To call a truth
record transduction-positive, require all three: the residual is at the **3′** end, it **maps
uniquely** to the reference, and it sits **adjacent to a full-length reference L1** — that
conjunction is the 3′-transduction signature. Anything failing it is `TRANSDUCTION_INDETERMINATE` and
leaves the transduction denominator, rather than silently counting as a negative.

**Error taxonomy — the whole point is that these are different errors:**

| outcome | meaning |
|---|---|
| `TP_POS_AND_SOURCE` | position within tolerance **and** source locus within tolerance |
| `TP_POS_WRONG_SOURCE` | position matches, source disagrees — a *source* error, not a miss |
| `TP_POS_SOURCE_MISSED` | position matches, we said `not_transduction`, truth is transduction-positive |
| `FP_TRANSDUCTION` | we assert a transduction where truth has none |
| `FN` | truth transduction with no call at all |

**Two denominators, never conflated:** position sensitivity over all truth records (197 / 1355 / 90
per above), transduction sensitivity over only the truth-derived transduction-positive subset
(~51-60 for L1, `TRANSDUCTION_INDETERMINATE` excluded). Quoting one number against the wrong
denominator is the specific failure this taxonomy exists to prevent.

> ### IMPLEMENTED 2026-08-02 — `benchmark/benchmark_truth_set.py`, and the conjunction yields 1
>
> `derive-td` (truth VCFs + bwa + consensus + reference + rmsk → a cached
> `truth_transductions.tsv`) and `score` (calls + truth + that table → the report). Both the
> `TD_SRC` parser and the taxonomy above are implemented; `tests/test_benchmark_truth_set.py` pins
> the parsers, the residual and the matcher, and it caught two real bugs before the first real run.
>
> **The conjunction this section specifies yields a transduction denominator of 0 for L1 and 2 for
> SVA** (job `49482568`): 197 → 189 align to the consensus → **20** with a 3′ residual ≥100 bp → 10
> at MAPQ ≥30 → 7 also ≥80% aligned → 4 not inside a reference L1 → **0** within 10 kb of a
> full-length (≥5000 bp) reference L1. See §4 Stage 2's progress block for why loosening the
> thresholds is the wrong response: 11 of the 20 residuals are divergent element body (every one an
> older subfamily, L1PA6-8 or L1M*), which the L1HS-centric 4-record consensus cannot cover, so a
> looser criterion admits noise rather than recovering signal. **This route cannot gate Stage 3's
> transduction work on HG002.**
>
> **Two design points the taxonomy above does not have, both forced by real output:**
>
> - **`UNADJUDICATED_SIBLING`.** A sibling transduction's source element is by definition ABSENT
>   from the reference, so the adjacency test at the heart of this derivation cannot be satisfied
>   even when the call is right. Folding sibling claims into `FP_TRANSDUCTION` charges the caller
>   for the derivation's blind spot — it turned 17 of BASELINE-0's 82 SVA calls into false
>   positives. They are reported, not scored.
> - **The orphan stratum, which is 37 of L1's 48 position-level false positives.** An orphan
>   transduction *site* is where transduced sequence landed WITHOUT the element body; it is not an
>   MEI insertion, and this truth set enumerates MEI insertions, so an orphan call cannot match
>   whatever the caller did. The split is total, not noisy — 125 L1 calls with no `TD_SRC` are
>   118/125 on-truth, 26 sibling calls 22/26, **37 orphan calls 0/37** — which is what makes it a
>   category error rather than a quality signal. L1 precision is 74.5% with them and 92.7% without;
>   the report prints both.
>
> Alu carries no transductions at all (non-autonomous, §6.2), and `derive-td` **refuses**
> `--family Alu` rather than emitting an all-negative table that would read as a measurement.

> ### SETTLED 2026-08-02 — HG002 cannot gate transduction, and the search stops here
>
> `HANDOFF_2026-08-02` §4 NEXT item 1 asks: before Stage 3, decide how transductions get scored at
> all, and check whether the paper's supplement enumerates them. **It does not, and the truth files
> themselves answer it without needing the supplement.** A third independent line of evidence, after
> the ALT-vs-consensus derivation (job `49482568`) and the ALT-vs-genome cross-check (`49482779`):
>
> **The truth VCFs' own INFO vocabulary is a GIAB SV benchmark's, plus assembly provenance, and
> contains no transduction field.** The complete key set across all five files is
> `BREAKSIMLENGTH ClusterIDs ClusterMax* Dist* ExactMatchIDs HG003_GT HG004_GT HG2/3/4count
> *calls *exactcalls MendelianError MultiTech* NumClusterSVs NumExactMatchSVs NumTechs* REFWIDENED
> REPTYPE SVLEN SVTYPE TRall TRgt100 TRgt10k segdup sizecat` + `QNAME QSTART QSTRAND`
> (NA24385 haplotype-contig provenance) + `xTEA_specific` (SVA only, **5 records**). These are
> insertion-site records with sequence. There is no source-locus column, and there is no place one
> was dropped from — the schema never had one.
>
> **So the decision, and Stage 3 should not revisit it:**
>
> 1. **Stage 3's gate is position-level, per family, plus the named strata** (orphan / sibling /
>    no-`TD_SRC`), exactly as `benchmark/benchmark_truth_set.py score` already reports it.
>    Transduction *calls* are reported and stratified; they are **not** scored as right or wrong.
> 2. **Transduction correctness needs a different sample**, and acquiring one is Stage 4 work, not a
>    Stage 3 blocker. Nothing in Stage 3's evidence-model rewrite is unable to proceed without it.
> 3. **Do not derive transduction truth from the ALT sequence again.** Two methods have been run,
>    they agree (L1 zero by both; SVA two or three), and this INFO-key check is a third. The most
>    likely reason is a property of the curation, not of the data: the benchmark was built by
>    verifying a polyA tail and a target-site duplication per insertion, which selects for clean
>    canonical TPRT products — plausibly the ones that do not carry a transduced 3′ segment.
>
> One consequence worth carrying into any SVA number quoted from this set: `SVA_tprt_all.vcf`
> contains **5 `xTEA_specific` records**, i.e. 5 of its 90 truth records came from xTea itself.
> That is 5.6% of the SVA denominator, and it is mildly circular. It does not change any conclusion
> at the precision this benchmark supports, but an SVA sensitivity quoted to the tenth of a percent
> is quoting through it.

> ### ACQUIRED 2026-07-29 — this subsection's fetch instruction could not be followed as written
>
> Everything below about *why* to prefer fetching over realigning still stands. What changed is that
> the specific file it names does not exist to be fetched:
>
> - **The paper's CRAM is gone and unobtainable.** `run_tools/` carries the absolute path
>   `/n/data1/hms/dbmi/park/simon_chu/projects/XTEA/xTEA_paper/HG002/`; `simon_chu/projects/`
>   survives but `XTEA/` does not, and the file is not hosted publicly under that or any other name.
> - **GIAB has no bwa-mem GRCh38 Illumina alignment for HG002.** Per
>   `data_indexes/AshkenazimTrio/`, the only GRCh38 Illumina alignment indexes are
>   `Illumina300X_wgs_novoalign` and `Illumina_2x250bps_novoalign`. A different aligner is a *larger*
>   confounder than the one accepted below.
> - **The promising local file is a dangling symlink.**
>   `vinay/datasets/for_acb/HG002_30x_2_hg38_bwamem.sorted.bam` points into
>   `park/DATA/GIAB/AshkenazimTrio/Illumina_short/`, which no longer exists. Every other local HG002
>   BAM is the `1mb_bin` subset. The "no usable full local HG002" conclusion holds, but for this
>   reason rather than the 20 MB-subset reason given below.
>
> **What was used instead:** `HG002.novaseq.pcr-free.30x.dedup.grch38.bam` (43.7 GB, md5
> `7378e27ba9a623ccb5eb1a4353aaebd9`) from
> `gs://brain-genomics-public/research/sequencing/grch38/bam/novaseq/wgs_pcr_free/30x/` — bwa-mem
> 0.7.17 + MarkDuplicates 4.1.2.0, `SM:HG002`, coordinate-sorted. **Aligned to hs38DH (3,366 `@SQ`),
> not `no_alt`**, so it carries exactly the alt-aware MAPQ confounder the bullets below warn about.
> Accepted because the faithful alternative is a FASTQ download plus a bwa-mem run for comparability
> this plan does not gate on (§2, §7.3), and because the pipeline is already proven on a 3,366-contig
> BAM — `NA12878.wgs.bam` is bwa-mem 0.7.15 against GRCh38 **with** alt/decoy, also 3,366 `@SQ`.
>
> **Consequence to honor:** our absolute sensitivity/FDR may run below the paper's for alignment
> reasons alone. **Never quote them side by side without saying so.** If a published-number comparison
> becomes the point, rebuild on `no_alt` — the FASTQs are in the same bucket under
> `research/sequencing/fastq/novaseq/wgs_pcr_free/30x/`, and the paper's exact reference is still local
> and already bwa-indexed at
> `/n/data1/hms/dbmi/park/simon_chu/projects/data/reference/GRCh38/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna`.
>
> Full reasoning, header dumps and dead ends: `test_data/HG002/PROVENANCE.md`.

**Fetch, do not realign.** The paper's `run_tools/` commands show the other callers running on
**`HG002_bwa_mem_GRCh38.sorted.cram`** against a `GRCh38/GCA_000001405.15…` reference. Get that file and
pass its reference:
- **Realignment is unnecessary** — GRCh38 primary-chromosome coordinates are identical whether or not
  alt/decoy contigs are present, so hg38 annotations apply directly. `validate_headers` enforces
  identical `@SQ` only **across multiple input files**, irrelevant for a single CRAM, and the pipeline
  already handles contig-rich references (NA12878.wgs.bam has 3,366).
- **It removes a confounder** — xTea filters on MAPQ (clip ≥ 12, disc ≥ 20) and alt-aware alignment
  lowers MAPQ in exactly the repeat-homologous regions that matter. Aligning ourselves would make our
  numbers incomparable to the paper's for reasons unrelated to our code.

So this is a **download**: CRAM + `.crai` + reference, then verify with `samtools quickcheck` and a
spot-check that annotation contig names resolve. Fall back to FASTQ + `bwa` only if unobtainable. GIAB
hosts HG002 Illumina WGS under
`ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/data/AshkenazimTrio/HG002_NA24385_son/`
(`NIST_Illumina_2x250bps/` confirmed present; **the exact path of the bwa-mem GRCh38 CRAM is NOT
verified** — locate it before budgeting).

Local HG002 that is **not** usable — both checked 2026-07-29, do not re-chase:
`…/ACB_miniproject_data/bams_for_test/HG002_30x_2_hg38_bwamem.sorted.1mb_bin.w_rg.bam` is **20 MB**, a
1 Mb-bin test subset; and `vinay/datasets/for_acb/HG002_30x_2_hg38_bwamem.sorted.bam`, which looks
full-genome and correctly aligned from its name, is a **dangling symlink** into a deleted
`park/DATA/GIAB/` tree. Long reads *are* local
(`/n/data1/hms/dbmi/park/yuwei/alignment/` — PacBio, HiFi, ONT-UL): not the benchmark input, but useful
for eyeballing individual transduction calls, since the inserted flank is directly visible in a long read.

### 7.2 Mosaic: SMHT005 descending colon (3G)

HG002's germline truth set cannot validate mosaic calls, but a **real orthogonal somatic-MEI truth set
exists in-house**, which supersedes the master-parity gate entirely (so the awkwardness of gating mosaic
on master agreement while germline explicitly is not, disappears):

**Truth set:** `/n/data1/hms/dbmi/park-smaht_dac/DATA/Globus_DataShare/MEI/p25/filtered_calls/v1.0.0/SMHT005-SMaHT-sMEI-v1.0.0.vcf`

- 55 PASS records, **all `SVTYPE=INS:ME:LINE1`**, all `TYPE=Single-tissue` (i.e. all somatic/mosaic).
- **15 carry `TISSUE=3G`** (Colon, Desc) — the descending-colon benchmark set.
- **Long-read derived** (`LR_CALLER=PALMER,longcallD`) — genuinely independent of short-read xTea, not
  another short-read caller's opinion. The strongest kind of gate available here.
- Scoring-ready INFO: `MEINFO=L1HS,strand,start,end` (gives element extent, so "right site, wrong 5′
  truncation" is checkable), `TSD`/`TSD_LEN`, `POLYA`, `PHASING`, `HETSNP`, `BASEINFO`.
- **24 per-tissue genotype columns** (`GT:PB_DP:PB_AD:ONT_DP:ONT_AD:ILL_DP:ILL_AD`) — the 23 non-carrier
  tissues are a **built-in specificity panel**: a call emitted in a non-carrier tissue is a false
  positive, with per-tissue depth to interpret it.
- Insertion sizes 112–1400 bp; Illumina depth in 3G ~350–480×, matching mosaic's high-depth target.

**Stated assumption, recorded because it is material.** For all 15 3G calls the long-read support is
**1–2 alt reads out of 9–52** (VAF ≈ 2–8%) and **`ILL_AD` alt is 0 at every one**. That AD field almost
certainly cannot represent a 112–1400 bp insertion from ~150 bp reads — a spanning-read AD is
structurally 0 at this size — so it is *not* evidence Illumina lacks signal, and xTea's evidence is
clip + discordant pairs, a different measurement. **We proceed assuming the short-read signal is
recoverable.** If early runs disagree, the honest fallback is to report mosaic **specificity** (via the
23-tissue panel) and treat sensitivity as bounded by what the data supports, rather than quoting a
sensitivity that was never achievable.

**Limitation: L1-only** — no Alu, no SVA, so mosaic Alu/SVA cannot be scored from it. `MEI/` also holds a
`Raw/` subdirectory and other donors; check those before claiming any mosaic result beyond L1.

**Input data — local, on our reference, no conversion.** Five CRAMs:
`…/DATA/GCC_Broad/SMHT005/SMHT005-3G/illuminaNovaseq_bulkWgs/seq_data/` (`SMAFI5HKGSER` 180.5 GB,
`SMAFI6WNA1XM` 142.6 GB) and `…/DATA/GCC_UW/SMHT005/SMHT005-3G/illuminaNovaseq_bulkWgs/seq_data/`
(`SMAFI4852JH2` 94.9 GB, `SMAFI4D3DFD4` 78.3 GB, `SMAFIH9MAGML` 64.3 GB). Use **all** of them, pooled via
the existing multi-BAM `bam_list`. SMaHT is aligned to
`/n/data1/hms/dbmi/park-smaht_dac/ref/GRCh38_no_alt/hg38_no_alt.fa` — **the same reference xtea2 runs on.**

**Pass `--ignore-rg`; header compatibility is otherwise VERIFIED** (Broad `SMAFI5HKGSER` vs UW
`SMAFI4852JH2`):

| check | result |
|---|---|
| `@RG` `SM` | **differ** — `SMATSBL826FN` vs `SMATSKLUWAXC` → `--ignore-rg` required |
| `@SQ` contig count | 195 in both |
| `@SQ` `SN`+`LN`+`M5` | **byte-identical** (`e677cdd9…`) once `UR:` is stripped |
| only difference | `UR:` — a per-run CWL staging path |

`validate_headers` (`xtea_scan.cpp:408-451`) compares count/names/lengths and **passes**. All five CRAMs
came from the same pipeline (`sentieon_bwamem_202308.01_GRCh38` in every filename) and matching `M5`
confirms the same reference *content*. `--ignore-rg` is the only flag needed; its docstring already names
this case ("as with GCC-split SMaHT data").

**Watch-item:** `cns_emit.hpp`'s per-file `ftag` path — what keeps same-named reads from different CRAMs
cross-pairing — has never run at five-CRAM scale.

**Data volume drives the mosaic speed requirement: ~561 GB across five CRAMs.** For scale, germline WGS
takes ~50 min on a single 39 GB BAM, and combined depth here is ~350–480×. So **the scan pass alone is a
multi-hour job with CRAM decode on the critical path.** Size mosaic timing runs accordingly, and treat a
chr1 mosaic slice as mandatory for iteration.

> **The chr1 mosaic slice EXISTS as of 2026-08-01** (job `49440268`,
> `docs/handoff/runs/chr1_slices.sbatch smht005`):
> `/n/data1/hms/dbmi/park/corinne/claude_sandbox/test_data/SMHT005-3G/chr1/` — all five CRAMs sliced,
> **43 GB total** (7.7% of the 561 GB source, i.e. chr1's genome fraction), each `.crai`-indexed and
> `samtools quickcheck`-clean, with `bam_list.txt` already written in the order the pipeline wants.
> Kept as **CRAM, not BAM**, deliberately: the mosaic input *is* CRAM, and the per-file `ftag` machinery
> flagged in the watch-item below has never run at five-file scale — converting would quietly change the
> thing under test. The `@SQ` dictionary is left unpruned (195 contigs).
> **`--ignore-rg` confirmed still required on the slices:** the two Broad files carry `SM:SMATSBL826FN`
> and the three UW files `SM:SMATSKLUWAXC`, exactly as the full CRAMs do — slicing does not touch `@RG`.

The same tissue also has PacBio HiFi (`pacbioHifi_bulkWgs`) and ONT (`ontStd_r10_bulkWgs`) — the right
way to adjudicate individual disputed calls.

### 7.3 The harness, and `x_cmp.py`

Nothing about `x_cmp.py` is algorithmically special: a window match (`--extnd`, default 200; the README
example uses 100), TP/FP/FN, `F1 = 2tp/(2tp+fp+fn)`, `sensitivity = tp/(tp+fn)`, `FDR = fp/(fp+tp)`, plus
annotated dumps. We could write better in an afternoon.

What it has is **comparability**: its quirks *are* the definition the published numbers were computed
under. Its `cmp_dict()` is benchmark-major (iterate truth, ask whether any call falls in the window),
which does not obviously enforce one-to-one matching — and whether one call may satisfy two truth sites,
and whether matching is greedy or nearest, moves sensitivity/FDR by percentage points. It also carries a
high-confidence-region notion, and HC restriction changes denominators.

**Write our own harness; keep `x_cmp.py` as a calibration reference.** `benchmark/benchmark_truth_set.py`
does what theirs cannot: one-to-one greedy matching, per-family scoring, transduction-aware scoring
(source element, not just position), stratification by insertion size and local coverage, tolerance
sweeps at 50/100/200 bp, and a provenance record (truth-set commit, pipeline commit, command line) per
run. Then run `x_cmp.py` once on the same inputs and reconcile: agreement within noise means ours is
trustworthy *and* on the published axis; divergence is itself the finding and must be explained before
any metric is quoted.

`x_cmp.py` lives at the root of `parklab/xTea_paper`, self-contained (`sys`, `os`, `optparse`,
`subprocess`) — a one-file recovery. **It was deleted from this repo in an earlier trim**, which is why
it must be fetched rather than found.

Calibrate targets on what the current code achieves, not on asserted thresholds. GIAB genome
stratifications are local at `/n/data1/hms/dbmi/park/SOFTWARE/REFERENCE/GIAB/genome_stratifications/`
(symlinked from `…/park-smaht_dac/ref/GIAB/`) if HC-region restriction is wanted.

### 7.4 Test-suite changes

`tests/test_downstream_from_cns.py` and `tests/test_e2e_concordance.py` encode agreement with the current
call set — keep strict through Stage 2, then convert to format-validity + truth-set metrics. Keep
`test_locus_table.py`, `test_gvcf_multi.py`, `test_transduction_core_split.py`. Add: unit tests for
`decide()` and the evidence record (the redesign's point is that this becomes unit-testable), a
VCF-format validation test, the `--types X` vs `--types L1,Alu,SVA` partition test (§6.1/6.2), and the
missing somatic gate (§6.3).

## 8. Module layout

**Rewrite:** `xtea/x_transduction.py` (evidence model, library form, no subprocess entry);
`xtea/x_orphan_transduction.py` (fold into the above); `xtea/run_xtea.py` (per-locus loop replaces
`build_type_steps`' three chains); `xtea/x_clip_disc_filter.py` (drop the gather + Pool machinery at
`:80-400`; keep `call_MEIs_consensus` as a library call).

**New:** ~~the locus universe~~ — **DONE 2026-08-03: `xtea/site_universe.py`**, the site-creation rule
lifted out of `XTransduction.re_slct_with_clip_raw_disc_sites` (which is now a delegation).
Enumerates every locus a family's chain can create a record at, with no BAM and no cns evidence, so
the per-locus loop can know its iteration domain (41,604 union loci) before gathering anything.
Tests: `tests/test_site_universe.py`. Gate: `docs/handoff/runs/gate_site_universe.py` (offline, 4/4
byte-identical). See §4 Stage 2c for what it does and does not unlock;
`xtea/downstream_unified.py` — **STARTED 2026-08-03**: the locus universe half is built
(`enumerate_family_sites`, `LocusUniverse`, `transduction_extnd`), driven by `run_xtea`'s
`--locus-universe` at the `cns_filter → transduction` barrier, diagnostic until the loop
consumes it. Tests: `tests/test_downstream_unified.py`. Gates:
`docs/handoff/runs/gate_locus_universe.py` (offline, WGS, argument derivation + byte compare)
and `gate_locus_universe_chr1.sbatch` (wiring: ON vs OFF calls identical). Still to come in
this module: the locus **evidence** table and the per-locus decision loop — and note the
batched-bwa constraint in §3 before shaping them;
`xtea/transduction_evidence.py` — **BUILT 2026-08-03, extended 2026-08-04** (`SiteEvidence`,
`TransductionParams`, `DecisionLedger`, `PHASE_ORDER`/`phase_rank`, and the gates; `decide()` turned
out to be **six gate groups**, see §4 Stage 3's two blocks for why that is the data flow rather than a
compromise). The transduction step's whole decision cascade now lives here as named thresholds and
named reasons, answer-preserving; six methods were deleted from `x_transduction.py`. **[B2] is the
exception**: `check_seprt_disc_consistency` is shared with `cns_filter`, so its body stays in
`x_clip_disc_filter.py` and only `gate_cns_disc_cluster` moved here, recording into a
**caller-supplied** ledger. Tests: `tests/test_transduction_evidence.py` (differentials against
verbatim copies of the deleted bodies — and for [B2], against the SHIPPING method, because there the
method is the seam). Gates: `docs/handoff/runs/gate_transduction_evidence_{chr1,wgs}.sbatch` and
`gate_cns_cluster_ledger_{chr1,wgs}.sbatch` (the latter's chr1 leg compares the ledgers exactly:
AFTER minus its B2 rows must be byte-identical to BEFORE);
`xtea/mosaic_adapter.py` (§6.4); `benchmark/benchmark_truth_set.py`; ~~a new `xtea_scan`
subcommand~~ — **DONE 2026-08-01: `src/csrc/td_gather.hpp`**, driven by `xtea_scan --td-sites`,
called from `XClipDiscFilter.native_td_gather`, selected by `--td-gather-scan` →
`XTEA_SCAN_TD_GATHER`. Tests: `tests/test_td_gather.py`. Gate: `docs/handoff/runs/stage1_gate_chr1.sbatch`.

**Reuse, do not reinvent:** `xtea/interval_index.py` (`GenomeIntervalIndex` — sorted-array +
prefix-max-end, 105× faster to build than the AVL tree it replaced); `xtea/x_alignments.py:44-75`
(`get_cached_bam`/`get_cached_bam_meta`, pid-keyed, fork-safe); `xtea/locus_table.py`;
`src/csrc/cns_emit.hpp` patterns; `xtea/scan_adapter.py` / `disc_adapter.py`.

**Delete — five dead functions only** (verified zero callers): `call_transduction_with_unknown_source`,
`call_high_confident_transduction`, `update_existing_candidate_sites` (all three empty stubs),
`collect_realign_reads`, `parse_clip_algnmt_rescued`.

**Do NOT expect to delete modules.** An earlier revision of this plan listed seven (`x_cmp.py`,
`union_find_set.py`, `x_insert_size.py`, `x_population_AF.py`, `x_ref_TE_deletion.py`,
`x_somatic_panel_normal.py`, `x_sv.py`). **All seven are already absent** — removed in an earlier trim; a
stale `xtea/__pycache__/x_cmp.cpython-310.pyc` is the only trace. The original "verified zero references"
claim checked references but not existence.

## 9. Environment and conventions

- **The interactive session has 2 cores — never time anything in-session; always `sbatch`.**
- Turnarounds: chr1 ≈ 5 min, germline WGS ≈ 45-60 min, transduction-only ≈ 25 min, **mosaic full-depth
  multi-hour**.
- Slurm: `-A park` with `park|short|medium|long`, or `-A park_contrib` with `park|priopark`. Never mix.
  Never throttle arrays with `%N`.
- **Never put comma-bearing values in `sbatch --export`** — it splits on every comma and the run looks
  successful. Type lists go as positional args.
- R needs `module load gcc/14.2.0 R/4.4.2` first. Figures use Nimbus Sans.

## 10. Verified facts and corrections — do not re-derive

Checked against the code/data this session. Several came from design passes that had them **wrong**:

- `MIN_POLYMORPHIC_SOURCE_DIST` is **1000**, not 100 (`global_values.py:166`).
- **`TD_SRC` HAS FIVE SHAPES, AND TWO OF THEM USED TO PARSE TO `None`** (fixed `2d38feb`, verified
  2026-08-04 against master's NA12878 and both HG002 arms). The grammar is
  `{chrom}:{start}-{end}~{rc}` followed by either `~orphan` or `~{cluster_l}-{cluster_r}`, plus the
  two tail-free shapes `{chrom}:{start}-{end}` and `{chrom}:{pos}~sibling`. **The source region is the
  first `~`-field after the colon; everything after it is provenance.** `{source}~orphan` and
  `{source}~{l}-{r}` both returned `None` — indistinguishable from `not_transduction` — dropping 7 of
  master's 115 assertions, 3 of BASELINE-0's and 5 of the decoy arm's, and reaching the scorer as
  `TP_POS_SOURCE_MISSED`. It moves exactly one reported BASELINE-0 number: L1's `FP_TRANSDUCTION`,
  0 → 1. **Never `rsplit("-", 1)` the whole field** — the tail is another interval.
- **`TD_SRC` HAS TWO WRITERS AND ITS *BARE* SHAPE IS TWO DISCORDANT-PAIR MIDPOINTS, NOT A SOURCE
  INTERVAL** (verified 2026-08-04 from both codebases at identical lines). `x_transduction.py:1252`
  writes `"{chrm_fl_L1}:{source_start}-{source_end}~{sourc_rc}"` from a **flank record** and
  `call_candidates` always appends a tail, so a real source region is never bare.
  `x_orphan_transduction.py:371` writes a bare `"{src_chrm}:{i_start}-{i_end}"` whose coordinates its
  own comments at `:352-353` call the **"mid-position of left discordant pair"** and **"mid-position of
  right discordant pair"**. **So on the bare shape `start < end` says NOTHING about copy status** —
  `benchmark_truth_set.td_src_source_class` returns `disc_span` for it, which is not a subtype of
  reference_copy / polymorphic / unknown_start. Measured corroboration: bare spans are 1–2814 bp
  (median 28) against a ~6 kb full-length L1; **none of master's 33 bare intervals lie inside or within
  5 kb of any source region** in `unique_trsdct_disc_only`; the tailed shape's first interval reaches
  6030 bp. **Consequence: only 7 of master's 115 NA12878 assertions and 3 of BASELINE-0's 90 have a
  flank-derived source interval**, so any per-family "reference-copy source" count taken from a VCF is
  reading disc spans. **The SVA decoy measurement is NOT affected** — `0 of 105 → 241 of 346` is taken
  from `unique_trsdct_disc_only` column 3, which is the flank-derived id.
- **WITHIN THE BARE SHAPE, `start == -1` MEANS THE SECOND SIDE WAS EXCUSED BY REPEATS — NOT MERELY
  MISSING — AND IT VERIFIES ~29 POINTS WORSE.** Measured twice: master's NA12878 pool two-sided
  (`start ≥ 0`) **83.9%** against one-sided **55.0%**; the 3-arm HG002 pool **85.2%** against
  **56.2%**. Two samples, same direction, same size.
  **PINNED TO A LINE 2026-08-05** (`f3dd6d0`), replacing the earlier "strongly indicated" reading.
  `-1` is `mid_pos`'s initialiser in `_form_unique_cluster` (`x_orphan_transduction.py:1064`), returned
  unchanged when a side has too few disc pairs (`:1072`) or no dominant source chromosome (`:1082`) —
  **but a failed side alone does not produce a surviving record.** `b_consistent` (`:596-619`) requires
  one side to have formed, and for the other it calls `_check_form_repeat_or_same_cluster`, which
  returns `True, -1` **only at `:1120`, when that side's disc mates are predominantly INSIDE REPEATS**
  (`:1126` returns a real position; `:1128` does not survive). So a one-sided entry is one whose second
  side is **repeat-driven**, which is exactly what the verifier's repeat test exists to distrust — the
  ~29-point gap has a mechanism rather than being a correlation. Byte-identical in master (`:1093`,
  `:1148`, `:1155`). **And `END` cannot be `-1`**: the writer takes `i_end = max(...)` (`:359-363`) and
  a formed cluster's position is the median of a non-empty list
  (`_position_of_same_chrm_form_cluster`, `:1086-1096`), so `parse_td_src`'s `100--1` branch is
  unreachable for this writer rather than merely unobserved.
- **THERE ARE TWO DIFFERENT "ORPHAN"s AND THEY ARE DISJOINT BY CONSTRUCTION** — measured first,
  then **pinned to one line of each of two writers 2026-08-05 (`18aed56`), which makes it a
  prediction rather than a property of one run.** The mechanism: `SUBTYPE ==
  orphan_or_sibling_transduction` is set by exactly **one** writer, `x_orphan_transduction.py:378`,
  on a site the caller **CREATED** — and that writer emits a **bare** source with **no tail**
  (`:371`). `td_src_call_class` reads the **tail**, and "no tail" maps to `regular`. So every
  created orphan-or-sibling SITE is `regular` by `TD_SRC` whatever its biology, and every
  `~orphan`-tailed source came from the other writer (`transduction_evidence.py:847`, on a gate [A]
  survivor whose disc anchors bracket the insertion point — `select_orphans`, which is explicitly
  not a gate). Confirmed in **two codebases**: 37 of 37 on xTea2's WGS decoy arm (`49792012`) and
  **43 of 43** on master's HG002 (`49763841`), so the collision is **inherited**. In both runs every
  `~orphan`-tailed call carries a different `SUBTYPE` and the `tailed` shape.
  **THE WORDS, now different in code as well as in prose:** a **"created orphan-or-sibling SITE"**
  (`benchmark_truth_set.SUBTYPE_ORPHAN_SITE`, which exists so that a bare string comparison cannot
  reintroduce the confusion) versus an **"orphan-TAILED source"** (`CALL_CLASS_ORPHAN`). Also fixed:
  three places cited the wrong population's evidence — `td_src_call_class`'s docstring, its test's
  docstring, and **`verify_transduction`'s PRIMARY STRATUM header, which printed the 0/37 above a
  table strata'd by the tail**. Re-scored to prove it moved no number: TP 140, FN 57, FP 48, 74.5%,
  92.7% over 151, all identical. **And a test pins the mechanism, so §11's bare-shape writer fix
  will FAIL it the moment that writer gains a tail** — at which point the two populations stop
  being disjoint and every claim keyed on either field has to be re-read.
  The original measurement, kept because it is the evidence: `SUBTYPE ==
  orphan_or_sibling_transduction` is what `benchmark_truth_set`'s precision carve-out keys on and what
  §1's **"BASELINE-0's L1 orphans are 0/37 on-truth"** measured. `td_src_call_class == "orphan"` is the
  `~orphan` tail in `TD_SRC` and what the Stage 3 gate strata by. **All 37** of L1's SUBTYPE-orphan
  calls are `regular` by `TD_SRC`; SVA has **zero** SUBTYPE-orphan calls, so the carve-out never fires
  for SVA at all. **The two fields do not agree on a single call in that run.** Consequences: any claim
  of the form "orphan calls are X" must name its field; the 0/37 is *not* evidence about TD_SRC
  orphans (which is what the handoff's §10.1 cited it as, for a decision that is nonetheless
  unaffected); and SVA's post-flip precision drop is **not** explainable as the category L1 excludes.
  The carve-out is deliberately **left keyed on SUBTYPE** — that is the population 0/37 was measured
  on, so re-keying it would move a reported precision number. Changing it is a decision, not a cleanup.
- **A CALL-COUNT DELTA UNDERSTATES A FLAG THAT CHANGES `TD_SRC` ON EXISTING CALLS** (measured
  2026-08-05, `b2a7f41`). The decoy flip's call set went 82 → 85 SVA, a strict superset. But its
  *verification pool* went 90 → 94, and the fourth entry is `chr17:45593489` — a position **both**
  arms call, whose `TD_SRC` went from `not_transduction` to `chr17:46252374-46253788~1~orphan`. So
  the flip gained a transduction assertion on an already-called site, **VERIFIED on 9 disc reads**,
  at zero precision cost. Three handoffs carried the 4-vs-3 discrepancy without reconciling it.
  **When a change can rewrite a field on a call it did not create, diff the FIELD and not only the
  position set** — and note the corollary for the truth-set score: a TD_SRC conversion is invisible
  to a position metric in both directions.
- **THE SVA TRUTH SET HAS ZERO chrY RECORDS, SO A chrY SVA CALL IS OFF-TRUTH BY CONSTRUCTION**
  (0 of 90; L1 also 0 of 197). **This is not a contig exclusion** — the Alu truth set carries 2 chrY
  records, so chrY is in scope and the absence is family-specific. Whether it reflects a real
  absence or SVA-specific incompleteness on chrY is not decidable from the truth set, but either way
  no SVA call there can match, so it is a denominator artifact rather than a caller error. This is
  one of the three calls behind the flip's 93.9% → 90.6%. Relevant beyond this item: with **90**
  truth records and a **5**-FP budget, one SVA call is worth ~1.2 precision points, so any SVA
  precision delta smaller than a few points is within the reach of a single unmatchable position.
- **GATE [A]'s polyA THRESHOLD IS UNREACHABLE, SO THE DECISION IS A MEMBERSHIP TEST** (measured
  2026-08-05 at WGS, job `49836568`). `m_polyA`'s per-locus record is created and one of its six
  counters incremented in the same block (`x_transduction.py:929-975`), so **every entry that exists
  has `sum(rcd[:6]) >= 1`** — measured 0 zero-sum entries of 4,604 (L1) and 5,261 (SVA). At HG002's
  coverage-derived `n_clip_cutoff = 3` the derived `n_polya_cutoff` is **0.5**, so `< 0.5` needs a
  sum of zero and `low_polya_clip` cannot fire; it first appears at cutoff **5**. So `n_polya_cutoff`
  is a named parameter that decides **nothing** at the shipping operating point, and gate [A]'s two
  live decisions are a clip-count cutoff and **polyA-map membership**. Membership is the binding one:
  of loci with `n_clip >= 3`, L1 drops **6,056** for absence against 2,802 kept and SVA **11,304**
  against 3,233, with 594 / 810 of the dropped carrying `n_clip >= 10` (max 969). **And the two
  thresholds share one dial** — `n_polya_cutoff = n_clip_cutoff/2 - 1` — so above 4 a clip-support
  sweep is not single-variable. `R_NO_POLYA_CLIP` writes **no detail**, which is why the ledger
  cannot answer any of this and the probe can.
- **A FINISHED RUN'S CLIP SAM SURVIVES, so gate [A] is measurable offline** — `temp_transduction_clip.sam`
  (2.46 GB L1, 1.6 GB SVA at WGS) is under `tmp/transduction/`, and `parse_transduction_clip_algnmt`
  builds **both** of gate [A]'s inputs from it, the clip maps and the polyA map. This is the same
  discovery as the cns SAM's (which was believed deleted, see the `MIN_DISC_MAPPED_RATIO` entry): the
  step's inputs outlive it, and a gate filed as "needs an instrumented run" for four handoffs took
  **57 seconds** on 2 cores. **Check the intermediates before scoping a run.** `probe_clip_gate.py`
  also shows the strongest available check 0: a probe that reproduces the run's own *output file* as
  a set, not merely its counts.
- **A SINGLE-ARM POOL IS NOT COMPARABLE TO THE 3-ARM FROZEN BASELINE, and the difference is arm scope,
  not caller behaviour.** The frozen baseline is the union over master + BASELINE-0 + decoy, 115
  entries, of which **master contributes 21 uniquely**. A single-arm pool holds ~90, so comparing one
  against it "loses" verified entries in nearly every stratum for reasons that have nothing to do with
  the change under test. Use the frozen baseline for the gate's **precision range**; use two single-arm
  pools against each other for a **one-flag delta**. This is the same error as the two above — a number
  taken from the wrong population — and `td_gate_no_regression.sbatch` refuses it in a comment.
- **A DROPPED HYPOTHESIS IS INVISIBLE TO A PER-STRATUM VERIFIED COUNT**, because it is absent from both
  sides of it. "No stratum lost a verified entry" and "no hypothesis was removed" are different claims
  needing different checks — hence check 2a. Measured for the flip: the ON pool is a **strict superset**
  of the OFF pool, all 90 entries verbatim plus exactly 4, so on that pair a VERIFIED loss could only
  have come from nondeterminism, which makes the no-regression check *also* a determinism check.
- **A PROBE THAT REIMPLEMENTS A SHIPPING PREDICATE WILL BE WRONG, and here is the case in numbers.**
  `probe_leftmost_cns.py`'s first version carried its own `_reach_5_prime_end` *and* its own candidate
  parser while its docstring claimed to import them. `XIntemediateSites.load_in_candidate_list_str_version`
  stores `fields[2:]` — **dropping chrom and pos** — so the shipping `rcd[17..20]` are the line's
  columns 20-23, the four cluster fields. Reading the whole line put `rcd[17]` two columns early and
  fed plain **counts** into the cluster test (`0 < 4 < 500` passes), so L1 at 500 read **50.3% instead
  of 19.3%** — plausible, 2.6× wrong, and on its way into a handoff. Import the predicate; a probe with
  its own copy of a rule can agree with a scorer that disagrees with it.
- **LINE-NUMBER CITATIONS IN THIS REPO DRIFT, and this is the second time.** The two `TD_SRC` writers
  were cited as `x_transduction.py:1196` and `x_orphan_transduction.py:365` in ten places across three
  files plus this plan; they are at **`:1252`** and **`:371`**, and the midpoint comments at
  **`:352-353`** not `:345-346`. Corrected `ad109e4`. When citing a line, grep the statement rather
  than trusting an offset — these particular citations *are* the evidence for the bare-shape argument,
  so one landing 56 lines away is one a reader abandons.
- **THE CALL CLASS AND THE SOURCE CLASS ARE INDEPENDENT AXES; "orphan" does NOT imply
  "polymorphic".** Master's 4 SVA orphans are polymorphic (`start == end`); its 3 L1 orphans are
  reference copies spanning ~6 kb. The call class is in the tail (`td_src_call_class`), the source
  class in `start` vs `end` (`td_src_source_class`).
- **A POLYMORPHIC SOURCE IS VERIFIABLE FROM THE REFERENCE.** What it lacks there is the element
  **body**; the transduced sequence comes from downstream of its 3′ end, which is ordinary genomic
  sequence — literally the window `gnrt_flank_regions_of_polymerphic_insertions` cuts its flank from,
  which is *why* it writes `start == end`. Master's 4 SVA orphan-polymorphic entries verify 4/4. Do
  not re-derive "polymorphic means unverifiable"; it was written down once and it was wrong.
- **THE STEP'S PICKED SOURCES AND THE CALL SET'S `TD_SRC` ARE DIFFERENT UNIVERSES, NOT A FILTERING
  CHAIN.** `x_transduction.py:238` keys `…unique_trsdct_disc_only` on the insertion positions of the
  step's `--input2`, which the golden's chain sets to `candidate_list_from_disc.txt.clip_sites_raw_disc
  .txt` — the **raw disc-site universe, 36,687 rows for L1** — while the VCF is a ~200-call file at
  clip-refined positions. Measured 2026-08-04: **104 of its 111 L1 loci are exact matches in that raw
  file**, and only **4 of the VCF's 59 non-sibling transduction calls** land within ±200 bp of any
  picked row. So the counts (L1 159 vs 39; SVA BASELINE-0 105 vs 1; SVA `--td-decoy-structural` 346 vs
  5) are not a survival rate and must not be quoted as one — an earlier version of this entry called
  them "attrition" and that was wrong. **Quote the level you mean.** A picked-source count sizes the
  *step*, which is where Stage 3's thresholds act; `TD_SRC` sizes what the call set *asserts*, which is
  what a precision gate scores. Reading one as the other is how the decoy fix's contribution to the
  gate's pool was recorded as 241 when it is 4.
- **SOURCE-INTERVAL LENGTH DOES NOT PREDICT THE VERDICT — checked, negative.** Mean 214 bp for
  REFUTED against 232 for VERIFIED on master's L1 reference-copy stratum, and VERIFIED includes
  intervals of 2, 5 and 6 bp. Most L1 reference-copy source intervals are tiny (median ~25 bp), so
  `start < end` versus `start == end` is a much thinner distinction than "a full-length copy versus a
  point". Do not re-propose interval length as a filter.
- **THE STEP'S CUTOFFS ARE NOT `--nfclip/--nfdisc`. They are derived from COVERAGE, and at HG002 WGS
  they are 3 and 4, not 4 and 5** (verified 2026-08-04 by a replay that reproduced 94 of 105 picked
  SVA sources until it used the right number). `x_TEA_main` sets `b_automatic = True` unless `--user`,
  which nothing in this repo passes, so both cutoffs come from
  `x_parameter.Parameters.get_par_by_cov(coverage)`: HG002 WGS is **31.816×** →
  `m_cov_par[30] = (3, 4, 1)` → **`n_clip_cutoff` 3, `n_disc_cutoff` 4, `n_half_disc` 2.0, decoy
  cutoff 4.0**, and `n_polya_cutoff` 0.5 rather than 1.0. `run_xtea.py`'s locus-universe docstring
  records the disc half of this ("wrong by one at HG002 (4, not 5)"); nothing else did, and
  `tests/test_transduction_evidence.py`'s header called 4/5 "the production ones" until this was
  measured. **Any threshold question sized off `--nd 5` is sized off a run that never happened**, and
  a different sample's coverage moves every one of these numbers.
- The polyA signature is **5 consecutive** bases — `AAAAA`/`TTTTT`/`AATAA`/`TTATT` (`x_polyA.py:55`) —
  not 3.
- `docs/handoff/master_golden_seeded/` **does not exist**; do not cite it as a golden fixture.
- **All seven "dead modules" are already absent** (§8). Checking references is not checking existence.
- ~~The claim that transduction's bottleneck is *not* the bwa passes is **unverified**~~ — **MEASURED
  2026-08-01 and CONFIRMED** (job `49437990`): the bottleneck is the read **gather** (73-86%), not
  bwa (12-25%). And **SAM parsing is 1.5-1.8%, not a cost at all** — the "2.4 GB SAM parse" framing
  is refuted. See §4 Stage 0.
- **Master's transduction figures come in two non-interchangeable flavours:** master's own pipeline did
  L1 448s / SVA 587s on its *own* smaller inputs (36,687 L1 raw sites); master's *code* on the *port's*
  inputs took 1126s / 1379s. The whole-run comparator was master's 50:18.
- **Do not compare raw-disc files across families by whole line** — the `nL`/`nR` columns are
  family-specific, so line-wise diffs report ~66K spurious differences. Compare `cut -f1,2`.
- Filtering `raw_disc` to the cns candidate list is **not** a safe optimization: on WGS L1, 64 of 322
  output sites (20%) and **all 59** orphan/sibling calls are absent from the 300-site candidate list.
  Transduction's wide pre-cns input is load-bearing.
- The `.slct` site set is **not** a subset of `raw_disc` (206 L1 / 151 SVA positions are new), because
  the re-selection can introduce positions.
- SVA structurally contains Alu-derived sequence — cross-family flank competition is real (§5).

Added 2026-07-29 during germline acquisition:

- **The paper's `HG002_bwa_mem_GRCh38.sorted.cram` cannot be obtained** — built in-house under a
  deleted directory, not hosted publicly, and **GIAB hosts no bwa-mem GRCh38 Illumina alignment for
  HG002 at all** (only novoalign). Do not go looking again; see §7.1.
- **`NA12878.wgs.bam` is aligned to GRCh38 *with* alt/decoy** — bwa-mem 0.7.15, NYGC CCDG, **3,366
  `@SQ`**. The plan previously cited its 3,366 contigs only as evidence the pipeline tolerates
  contig-rich references; it is also the reason an hs38DH HG002 is the internally consistent choice.
- **The HG002 truth sets were already local** at `corinne/dev_xTea/truth_set_hg002/`, byte-identical
  to `parklab/xTea_paper@2338ece` for all five VCFs. That directory also holds **`all_tes.bed`**
  (39 KB), which the paper repo does not ship — identify it before writing
  `benchmark/benchmark_truth_set.py`, it may already be what the harness needs.
- The truth VCFs are confirmed **headerless** (no `##`, no `#CHROM`) on the staged files — the §7.1
  warning is verified, not assumed.

Added 2026-08-01:

- **`_tprt` is a canonical-TPRT-mechanism filter and a strict subset of `_all`** — measured, not
  inferred: 0 `_tprt` positions are absent from `_all`, in either family. **The sensitivity denominator
  is `_all`.** Do not re-litigate this from the filename (§7.1).
- **This section previously stated `L1_tprt.vcf` has ~26 records. It has 184.** The correct counts are
  L1 197/184, Alu 1355/1257, SVA 90 (`_all`/`_tprt`). Anything sized off the old figure is wrong.
- **The truth sets carry no transduction or source-element annotation of any kind** — grepped for
  `transduct*`, `TD=`, `SOURCE=`, `ORPHAN`, `SIBLING` across all five VCFs: nothing. Do not go looking
  for a field to score transductions against; there isn't one.
- **`all_tes.bed` in `corinne/dev_xTea/truth_set_hg002/` is not a transduction annotation.** It is
  exactly **1642** rows of bare `chrom start end` with `end = start + 1` — i.e. the union of the three
  `_all`-flavour sets' positions (197 + 1355 + 90), carrying **no family column and no annotation of any
  kind**. Useful as a harness convenience only; it adds nothing the VCFs lack, and it cannot tell you
  which family a position belongs to. The earlier suggestion that it "may already be what the harness
  needs" is **half right**: right that it is a combined all-family position list, wrong if read as
  carrying annotation.
- **col5 was type-set-dependent and is now fixed** (`de96727`, §6.1). The value ceased to depend on
  `--types`; `col5_bed_plan()` does not take the requested types as an argument, by design.
- **`Error happen at merge clip and disc feature step: <contig> not exist` IS NOT AN ERROR.**
  `x_intermediate_sites.py:514` prints it when a contig carries clip candidates but no *discordant*
  candidates — the normal situation on alt/unplaced contigs. **Pre-existing and not HG002-specific:**
  the NA12878 WGS run emits **4,395** of them. HG002 emits **87,039** across **1,625** contigs (1,448
  `chrUn`, 135 `_alt`, 40 `_random`, plus `chrM`/`chrEBV`) because its BAM is alt-rich and the sample
  is male. **No primary chromosome is affected and the run completes normally** (the L1-only arm
  finished in 22:57 with 200 L1 calls). Do not chase it, and do not read it as the hs38DH deviation
  biting. It is worth fixing *eventually* only because 87K lines saying "Error" is where a genuine
  error goes to hide — that log was 8 MB.
- **`f=$(ls MISSING 2>/dev/null | head -1)` is a landmine under `set -euo pipefail`**, and
  `docs/handoff/runs/port_wgs.sbatch` still contains it. `2>/dev/null` hides ls's *message*, not its
  exit code 2; `pipefail` carries that through the pipe into the assignment and `set -e` kills the
  script. It has never fired there only because `merged/` always existed. It killed
  `baseline0_hg002.sbatch` job `49440169` *after* a fully successful 22-minute pipeline run. Use a
  glob helper. Related: `grep -c` exits 1 on a zero count while still printing `0`, so
  `grep -c ... || echo 0` prints `0` twice.
- **`run_xtea.py` deletes everything but the merged VCF on success** (`cleanup_successful_run`).
  Any comparison needing per-type, pre-merge call lists — which includes the Stage 0 partition gate —
  **must pass `--keep-intermediates`**, or the files are gone before you look.

Added 2026-08-01 while building the Stage 1 gather:

- **The legacy disc gather's cost is nearly independent of how many sites it was given.**
  `extract_mate_reads_by_name` splits the whole genome into 50 Mb bins and fetches each bin in
  full, re-parsing the entire disc-pos file once per bin: **106s for 207 chr1 sites and 106s for
  3733** (2 threads). Do not reason about it as per-site work, and do not benchmark it on a small
  site list expecting a small number.
- **Do not size the gather port from a chr1 measurement.** The legacy clip gather's cost is
  per-site overhead (one fork, two temp files), so it scales with the site count: **1.26×** for the
  clip half on chr1's 3733 sites, but the whole gather is **15.25× (L1) / 17.82× (SVA)** on the
  real WGS lists at 49,042 / 53,569 sites (job `49453857`, both arms at 12 threads). The first
  version of this plan's Stage 1 note projected 1.5-1.7× off the chr1 ratio and was wrong by an
  order of magnitude on that half.
- **`[TD-PHASE]` lands in `<sample>_<TYPE>.steps.log`, not `<sample>_xtea2.log`, whenever
  `--parallel-downstream` is on.** The earlier note (which says `_xtea2.log`) is right only for the
  serial path. A harness that greps only `_xtea2.log` reports an empty timing table and looks like
  the timers are missing.
- **An exact duplicate `(chrom,pos)` line in the transduction site list CRASHES the legacy gather** —
  `FileNotFoundError` out of `collect_clipped_reads_of_given_list`'s merge loop
  (`x_clip_disc_filter.py:373`), because the second Pool task overwrites the first's temp file and
  the merge unlinks it after the first read. Observed, not theorised. `re_slct_with_clip_raw_disc_sites`
  is what normally prevents it (raw-disc lines within `TWO_SITES_CLOSE_DIST` of a candidate are
  dropped), so a hand-built site list is where it bites.
- **The transduction gather is broken for multi-BAM input** and has been all along:
  `collect_clipped_reads_of_given_list` writes ONE fixed `all_disc_pos.disc_pos`, truncated per
  input BAM, so with N inputs only the last BAM's disc positions survive and every sample's mates
  are then gathered against them. No baseline exists (all germline runs are single-BAM; the
  multi-CRAM mosaic path is Stage 4). `td_gather.hpp` gathers each BAM against its own anchors.
- **`local arm=$1 armout="$OUT/$arm"` does not work under `set -u`.** The shell expands the whole
  command line before `local` assigns anything, so `$arm` is still unset and the script dies. It
  killed the first Stage 1 gate submission (job `49452897`) before arm A. Split the assignments.

Added 2026-08-02 while doing Stage 2a and the Stage 3 benchmark:

- **`--sva` reaches the cns_filter step and nothing else**, matching the golden
  `run_xTEA_pipeline.sh` (-C/-D/-N only). So the transduction and sibling steps have always run
  with the non-SVA parameters even for SVA, and **SVA's `leftmost_cns = 5000` has never executed**,
  in xTea2 or in master. Do not "fix" this as tidying: `consensus/SVA.fa`'s longest record is
  **3183 bp**, so 5000 admits every candidate with a positive cluster start and the
  polymorphic-flank library stops being a filtered subset. `tests/test_family_params.py` pins the
  resolution rule; §4 Stage 2 has the reasoning.
- **`consensus/*.fa` are multi-copy libraries, not single consensus sequences** — SVA 88 records
  (longest 3183 bp), LINE1 4 (6342), ALU 47 (490). A threshold compared against a "consensus
  length" has to be compared against the longest RECORD.
- **§7.1's ALT-residual figures were total soft clip, not the 3′ residual** — see §7.1's own
  correction. 3′-only is 20 records ≥100 bp for L1, not 60.
- **bwa splits a long query rather than soft-clipping it when the far end also hits the reference.**
  38 of 189 aligned L1 truth ALTs and 33 of 90 SVA ones carry a supplementary alignment to the
  consensus (5′ inversion is the usual cause). Any "residual = the primary's soft clip" calculation
  counts that second consensus hit as non-consensus sequence.
- **The ALT of a minus-strand insertion has its 3′ end on the LEFT.** VCF ALTs are reference-forward.
  Taking the primary SAM record's trailing soft clip — which looks right, since bwa revcomps SEQ —
  returns the element's 5′ end on every minus-strand insertion.
- **`hg38_L1.fa.out` is fragment-dominated**: 978,710 intervals, median length **285 bp**, only
  10,437 (1.1%) at or above 5000 bp. An "adjacent to a full-length reference L1" test against it is
  far stricter than it reads.
- **Most long 3′ ALT residuals are divergent L1 body, not transduced sequence.** 11 of the 20 L1
  truth residuals ≥100 bp map INSIDE a reference L1, and every one is an older subfamily — L1PA6/7/8
  at 4-7% divergence, L1MB/MD/ME at 16-25%. `consensus/LINE1.fa` is 4 L1HS-centric records, so an
  older element's body partly fails to align and surfaces as "residual". Any ALT-residual method on
  this consensus inherits the confound; loosening its thresholds admits body, not transductions.
- **BASELINE-0's L1 orphan-transduction calls are 0/37 on-truth and are 37 of L1's 48 position-level
  FPs.** Orphan transduction sites are not MEI insertions and the HG002 truth set does not
  enumerate them. Never quote L1 precision on this truth set without saying which stratum:
  74.5% overall, 92.7% excluding orphan sites.

Added 2026-08-02 (Stage 2c):

- **`scan_clip_realign.tsv`'s `nL`/`nR` are NOT family-discriminating**, despite being produced by
  one bwa per family (`realign.hpp`). Of 17,244,556 cross-family position-pairs, **1,147 differ**;
  of the 155,615 positions with ≥5 clip support, **504** differ. The three clip candidate lists are
  the same positions (577,293 / 577,312 / 577,347). §4 Stage 2c was originally written assuming
  this file carries "per-family support side by side" — it does not, and any locus evidence table
  built on it is three copies of one column.
- **Family-discriminating evidence starts at `candidate_list_from_disc.txt`**, where the mate is
  tested against each family's own annotation: 2,244 L1 / 9,468 Alu / 1,183 SVA, **union 12,223,
  only 672 rows claimed by more than one family**. This — not the realign TSV — is `cns_filter`'s
  input and the real basis for a locus evidence table.
- **`two_side_tprt_both` is written by `cns_filter`** (`x_clip_disc_filter.py:1845`), so
  `candidate_disc_filtered_cns.txt` is the **earliest** file the family-resolution rule can read.
  "Move the decision up front" cannot mean the candidate list without also meaning "adopt a
  different rule".
- **Family contention is small at every level and Stage 2c is not a speedup**: 672 of 12,895 disc
  candidates, 16 of 2,796 cns rows, 22 of 1,515 postF_hc calls. Partitioning by family saves
  0.6–5% of a downstream that is itself ~21% of the run after Stage 1. Do the stage for the
  architecture, not the clock.
- **The 12 `Alu over L1 tprt_override` conflicts are not contended at cns time at all.** The losing
  L1 record is created later, by transduction/sibling — 35 of the 1,493 called loci have no
  cns-level counterpart. So moving the decision earlier reproduces those 12 answers without the
  override, rather than losing them, which is the opposite of the obvious prediction.
- **The HG002 truth VCFs' INFO vocabulary is a GIAB SV benchmark's** (`ClusterIDs`, `ExactMatchIDs`,
  `NumTechs`, `TRall`, `segdup`, `sizecat`, `HG003_GT`/`HG004_GT`, …) plus `QNAME`/`QSTART`/`QSTRAND`
  assembly provenance. There is no transduction field and no schema slot one could have been dropped
  from. Third independent confirmation; §7.1 has the resulting decision. Do not go to the paper's
  supplement for this.
- **`SVA_tprt_all.vcf` contains 5 `xTEA_specific` records** — 5.6% of the 90-record SVA truth
  denominator came from xTea itself. Harmless at this benchmark's precision, but an SVA sensitivity
  quoted to a tenth of a percent is quoting through it.
- **A step can open an input it was never handed.** `x_TEA_main.py:551` builds
  `sf_candidate_list + HIGH_CONFIDENT_SUFFIX` from its own `-i`, so the transduction step reads
  `candidate_disc_filtered_cns.txt.high_confident` without that path appearing anywhere in the argv
  `build_type_steps` emits. Anything that rewrites a step's `-i` file must rewrite that companion
  too, or `update_high_confident_callset` (`x_transduction.py:632`) — which looks every companion
  row up in the transduction output with **no membership test** — dies on a bare `KeyError: <pos>`.
  This is how job `49498881` failed on `SVA:transduction`: the barrier dropped SVA's `chr6:31329020`
  (lost to Alu by `tprt_override`) from the candidate list and left it in the companion.
  `locus_table.CNS_COMPANION_SUFFIXES` is where that list lives. The other two companions,
  `.before_calling_transduction` and `.before_filtering`, are written and consumed **inside**
  cns_filter (genotype features) and must not be touched.
- **The companion is joined on `(fields[0], fields[1])`, not on `call_position`.**
  `update_high_confident_callset` reads `int(fields[1])` directly, while the locus clustering uses
  `min(lclip, rclip)` from `fields[2]`/`fields[3]`. Filtering the companion on the clustering key
  reintroduces the same mismatch it was meant to remove.
- **The family barrier and the postF_hc merge both run — the barrier does not replace the merge.**
  Measured (job `49502115`): the barrier resolves 16 conflicts, the merge then resolves 15, and the
  output is byte-identical to the un-barriered arm. A change at the barrier can only move a call if
  it drops a row that would have WON at the merge. Do not reason about the barrier as though the
  merge were gone.
- **12 of the 22 family conflicts are unreachable before transduction.** The losing family has no
  cns row at those loci; transduction and sibling create its record afterwards. So no decision point
  placed before transduction can resolve them, and **Stage 2c cannot collapse to one decision until
  site creation is unified too** — a Stage 3 dependency.
- **Do not compare two family assignments by nearest-locus matching across families.** A call at a
  position where its family has no row at the earlier stage will pair with a *different* family's
  nearby row and read as a flip. This produced a confident, wrong prediction of 4 moved calls
  (`probe_early_assign.py`, since fixed). Compare the two **conflict populations** on the winner
  position instead.
- **`run_xtea.log()` writes to stdout, not to `<sample>_xtea2.log`.** A gate grepping the latter for
  a driver-level line finds nothing and reports the feature as not having run — job `49502115`
  claimed "arm ON did not take the barrier" about a run whose barrier had resolved 16 conflicts.
  Driver lines are in whatever the caller redirected stdout to (`driver.log` in the gate scripts);
  `<sample>_xtea2.log` carries per-STEP output.
- **The per-type VCFs are at `merged/<bam_basename>_{LINE1,ALU,SVA}.vcf`**, named after the BAM and
  not the sample, since move A.2 made gVCF fan out per record from one shared tail. A glob built on
  `<sample>/<TYPE>/<sample>_*.vcf` matches nothing — and a comparison of two empty sets passes.
- **The td gather's "holding is permissive, emission is exact" invariant is now load-bearing, not
  just a nicety.** `run_td_gather_multi` collects over the union of the entries' site windows at
  the *largest* entry's `extnd`, so a narrower entry is always fanned out of a collection holding
  reads it must not emit. Anything that tightens the hold test in `collect_stripe`, or loosens the
  fan's two `ref_end <= wbeg` tests, silently corrupts multi-entry output while leaving the
  single-entry path correct — which is exactly the failure a one-entry test cannot see.
- **A shared BAM decode cannot reorder the gather, and the reason is the ownership rule, not the
  sort.** A fetch region owns reads starting at or after the previous region's end, so `ev` comes
  out ascending for *any* merge boundaries; the `stable_sort` that follows is a no-op and equal-start
  reads keep BAM order. That is what lets the fetch plan widen to a union without moving the FASTQ
  bytes bwa's tie-break depends on. Do not "simplify" it into a sort by anything but `ref_start`.
- **Sharing the td gather buys CPU, not wall clock, until the downstream is one job.** L1's and SVA's
  transduction steps are already separate concurrent subprocesses, so their two native gathers
  (53.3s and 56.8s at WGS scale) overlap; one shared gather is still ~57s of wall time. Measuring
  `--td-manifest` against two solo runs today and concluding it "doesn't help" is measuring the
  wrong thing.
- **52 of 1,515 per-family final records are created AFTER cns_filter**, all by the transduction
  step (50 L1, 2 SVA), none by sibling — measured, `probe_site_creation.py`. That is why exactly 7
  of the 22 family conflicts have both sides present at cns, and why **7 of 22 is the hard ceiling
  for any decision point placed before transduction.** It is a site-creation limit, not a rule or
  window one: widening `MERGE_WINDOW`, moving the barrier earlier, or reordering the priority
  cannot raise it.
- **Those 52 records are *created*, not moved, and they are already rows in the raw-disc file.**
  Every postF_hc position sits at distance **0** in `candidate_disc_filtered_cns2.txt` (the post-cns
  chain does not re-anchor at all), and the nearest cns row to each of the 52 is ≥ 3,293 bp away,
  median 2.4 Mb — so "created" is literal, and no wider matching window recovers them. All 52 are
  present in `candidate_list_from_disc.txt.clip_sites_raw_disc.txt`, written **before `cns_filter`
  runs**, and in the `.slct` selection derived from it — which needs no BAM, no bwa and no cns
  evidence. The rule is `xtea/site_universe.py`; `gate_site_universe.py` reproduces all four real
  `.slct`/`.slct2` files byte-for-byte offline. **Do not re-derive any of this by running the
  pipeline.**
- **Site enumeration being solved does NOT raise the 7-of-22 ceiling.** `.slct` is a 34-40K-row
  superset per family; the per-family evidence at those rows only exists after the transduction
  evaluation. What it unlocks is the per-locus loop's iteration domain — **41,604 union loci** vs
  75,853 per-family visits, 1.82× — not an earlier family decision.
- **The transduction/sibling steps ignore `--nd`.** `b_automatic` is True unless `--user_specific`, so
  the disc cutoff is `Parameters.get_par_by_cov(coverage)[1]` — **4** at HG002's 31.8×, not the `--nd 5`
  on the command line. Any offline reproduction of a step that reads `n_disc_cutoff` must derive it
  from `tmp_basic/basic_cov_is_rlth_info.txt`, not from the argv.
- **`XIntemediateSites.are_sites_close` implements a HALF-OPEN window**, because it walks
  `range(pos1-iwin, pos1+iwin)`. In the site-creation rule that makes the exclusion window on the
  raw-disc site's offset from a cns candidate `(-100, +100]`: exactly 100 bp to the *left* of a
  candidate is kept, exactly 100 bp to the *right* is suppressed. Symmetrising it changes which sites
  transduction re-selects.
- **`local a=$1 b="$x/$a"` is a bug in every gate script in this repo**, twice now (jobs `49452897`,
  `49516776`). Bash expands the whole `local` line before assigning any of it, so `$a` there is
  whatever the *caller* had. Use separate `local` statements, and do not name a loop variable the
  same as a helper's first parameter. The second occurrence made a passing change report FAIL.

Added 2026-08-03 (the driver-side locus universe):

- **The locus universe is computable before TRANSDUCTION, not before `cns_filter`.** The
  site-creation rule excludes raw-disc rows near a cns candidate, and that candidate list *is*
  `cns_filter`'s output. raw_disc and the rmsk annotation do precede `cns_filter`; the cns list does
  not, and it is the input that fixes where a driver-side pre-pass can run — the same point the
  Stage 2c barrier occupies. `site_universe.py`'s and `x_transduction.py`'s "available before
  cns_filter runs" said otherwise and are corrected in the tree. Placing the pre-pass on that
  sentence would have it excluding against a file that does not exist yet.
- **Transduction is handed `annotation1`, not `annotation`.** `{g}_L1.fa.out`, not
  `{g}_L1_larger_500_with_all_L1HS.out` (which is `cns_filter`'s). Both exist, both are plausible L1
  annotations, and the wrong one changes the rule's third rejection — so any driver-side caller must
  take `REP_TYPES[t]["annotation1"]`. `gate_locus_universe.py` fails on a mismatch rather than
  trusting the bytes to notice.
- **The transduction step's `iextnd` is `max(400, int(mean_is + 3*std_var))`** (`x_TEA_main.py:473`,
  `:497-499`), which is `run_xtea`'s `cns_extnd` plus a 400 floor that never binds at WGS insert sizes
  (HG002: 954). Both read the same cached `tmp_basic/basic_cov_is_rlth_info.txt`, so they cannot
  disagree. `downstream_unified.transduction_extnd()` is the one copy. **This closes the last gap 2b.4
  named for a `--td-manifest` driver caller** — the driver now has the site list *and* the extnd.
- **38,815 of the 41,604 universe loci are reachable only by creation** — no family has a cns row
  there (`gate_locus_universe.py`). Quote it as a fact about the iteration DOMAIN, never as a call
  count: the transduction evaluation rejects almost all of them and only 52 become records. What it
  says is that the created population *dominates* the universe rather than being a tail, which is why
  the 7-of-22 ceiling is structural rather than a small-numbers accident.
- **THE COST OF A cns_filter→transduction BARRIER IS NOT THE WORK DONE AT IT — it is that the barrier
  waits for Alu.** Measured at WGS, job `49550925`, and the critical path decomposes to the second:

  | | OFF | GATHER |
  |---|---|---|
  | `cns_filter` L1 / SVA / **Alu** | 57.3 / 62.5 / **124.1s** | 55.0 / 61.1 / **115.5s** |
  | SVA `transduction` (the long pole) | 348.1s | 292.0s |
  | `downstream_parallel` | **462.7s** | **505.3s** |

  OFF's critical path is `SVA cns 62.5 + td 348.1 + sib 49.9 + postF 2.1 = 462.6` (measured 462.7).
  GATHER's is `wait-for-Alu-cns 115.5 + universe 2.7 + gather 42.3 + SVA td 292.0 + sib 49.5 +
  postF 3.1 = 505.1` (measured 505.3). So the +42.6s breaks down as **+52.9 barrier, +2.7 universe,
  +42.3 gather, −56.0 transduction**.

  **Alu's `cns_filter` is ~2× L1's and SVA's**, and in the one-wave chain L1 and SVA start transduction
  as soon as *their own* `cns_filter` finishes, overlapping Alu's remaining ~60s. Any barrier there
  destroys that overlap. **This retro-explains the Stage 2c gate's unattributed +40s** (job `49502115`,
  OFF 2126s → ON 2166s): `--family-assign-early` uses the same barrier and paid the same cost.
- **So the gather hoist itself is a small WIN at WGS and the barrier is the loss.** −56.0s off the long
  pole for +42.3s serialized = **−13.7s**, plus **~66s of CPU returned** (51.49 + 57.22 → 42.3).
  On chr1 it was ~+0 (22s of overlapped gather for 19.6s serialized) with ~25s of CPU returned.
  Total run wall is NOT the place to read this: it moved 2263s → 2228s only because `scan_pass` varied
  1567.3 → 1500.8s, which the change cannot touch.
- **AND THE +52.9s IS TAKEN BACK — `run_xtea.StagedBarrier`, the default since 2026-08-03.** The join
  was never required, because the three things that run at this point have three different dependency
  sets and only one needs every family. Each family's `.slct` is built from **its own** cns list and
  raw_disc only (`enumerate_family_sites`), so enumeration needs **nobody** and runs the moment that
  family's `cns_filter` returns. The gather needs the **flanked** families (L1, SVA) and is needed by
  only them — Alu takes x_TEA_main's copy-through branch and never calls
  `call_candidate_transduction_v3` at all, so it never reads the gathered FASTQ and must not wait for
  it. Only the *union table* and the *family-assign* barrier need all three, and the union table still
  blocks nobody because nothing in the run reads it. **`--family-assign-early` therefore still joins,
  correctly** — its rule compares every family's cns rows against every other's.
  - **Same work, same order, fewer waits.** The ordering that mattered is preserved for free: every
    family enumerates before it registers and the assignment cannot be elected until all have
    registered, so all three enumerations still precede the assignment.
  - **The staged interleaving is the one production already runs.** With the flags off, the one-wave
    chain already has L1 and SVA in transduction while Alu is in `cns_filter`. Staging reproduces
    that, which is why it is answer-preserving rather than merely gated as such.
  - `--global-barrier` restores the join, and exists so the numbers above stay reproducible as a
    one-flag A/B rather than a checkout.
  - Two things staging introduces and the code handles explicitly: the driver gather now runs while a
    flankless family may still be in `cns_filter` (brief oversubscription, the same trade as the
    transduction core split), and `XTEA_TD_GATHER_PRECOMPUTED` is set **before any worker starts**
    rather than on the gather thread, because mutating `os.environ` concurrently with `fork`/`exec` is
    the one race staging would otherwise add. Setting it early is free: the step's test is
    `env AND both files present`, so a family the gather skipped still falls back to gathering itself.
  - **A staged barrier has to handle a death that a join got for free.** Under a join nothing waits on
    a family that died; here it does, so every failure path calls `abort`, which latches the failure
    and releases the waiters. Without it a non-zero `cns_filter` hangs the run inside the executor's
    shutdown instead of failing it.
  - **GATED PASS at both scales on one commit (`d565445`):** chr1 6/6 (job `49599837`) and WGS 6/6
    (job `49599838`) — `OFF == BASELINE-0` 4/4, GLOBAL and STAGED each 9/9 against OFF, the universe
    artifacts 3/3 identical between the barrier arms, 188/1223/82/1493 calls in all three arms.
  - **What it bought, measured (job `49599838`):** the flankless family's wait at the barrier went
    **45.8s → 0.1s** against a 45.49s gather, and `downstream_parallel` went OFF 511.75 / GLOBAL
    557.98 / STAGED 516.59 — **41.4s of the join's 46.2s recovered, ~90%**.
  - **THAT NODE WAS ~10% SLOWER than job `49550925`'s** (OFF 511.75s vs 462.7s; `cns_filter` 65.5 /
    155.1 / 70.6 vs 57.3 / 124.1 / 62.5). **Only within-job arm deltas are comparable between these
    two jobs.** Do not read one job's absolutes against the other's.
  - **The residual is +4.8s against OFF, not the −13.7s the hoist alone predicts**, and the likely
    cause is the gather contending with Alu's still-running `cns_filter` for cores the join handed it
    exclusively. **Likely, not measured**: the same run shows GLOBAL's L1/SVA `cns_filter` at
    83.4/88.9s against OFF's 65.5/70.6s — +18s in an arm where nothing touches `cns_filter` — so
    per-family variance is ±~18s and the +16.3s this is inferred from sits inside it. The only direct
    measurement of the wait is the 45.8 → 0.1s.
- **`--locus-universe` and `--family-assign-early` share the wave split but nothing else.** Either
  flag alone opens the `cns_filter → transduction` barrier (`downstream_waves`), and the universe does
  NOT inherit the assignment barrier's single-family exemption. **Order matters where both are on:
  the universe runs first**, because the assignment barrier *rewrites* each family's cns list to its
  own winners — enumerating after it would silently produce a post-assignment universe.
- **`call_candidate_transduction_v3` is set-oriented by construction**, so "per-locus loop" cannot
  mean per-locus end to end: it builds one flank reference and runs four batched passes over the
  whole `.slct` list (bwa index, clip realign, disc realign to flank, disc realign to consensus)
  before any per-site record exists. The shape the code supports is **enumerate → batched per-family
  evidence production → per-locus decision**, with only the third box a loop. Designing one bwa per
  locus is a different program, not a refactor.
- **The per-locus evidence record already exists in today's code, and a created site's is the same
  shape as any other's.** `collect_ins_info` assembles `(n_l_disc, lcns_start, lcns_end, n_r_disc,
  rcns_start, rcns_end, n_cns_clip, n_polyA_clip, i_lcov2, i_rcov2, s_src)` and
  `_collect_gntp_features` the 10-field genotype tuple; `add_new_rcd_with_cutoff` writes a created
  record from exactly those, which is why "what does evidence look like at a locus only `.slct` knows
  about" has the same answer as everywhere else. That is the loop's record type — do not invent one.
- **Every `driver.log` line carries seconds since driver start** (`[xtea2 +  123.4s] …`), added with
  the staged barrier. `step_timings.tsv` records *durations*, so two rows that overlap and two that
  abut are indistinguishable in it — and under `--parallel-downstream` overlap is the whole question.
  Nothing greps the `[xtea2]` prefix (checked across `docs/handoff/runs/` and `tests/`); every reader
  matches on the message, so gate scripts are unaffected unless they want the stamp.
- **`LocusUniverse`'s input map's iteration order IS the output file's column order.** `rows()` emits
  each locus's `families`/`claims` pairs in `self.families` order, i.e. the order
  `build_locus_universe`'s map was built in. **Job `49596361` failed on exactly this**: chr1, staged
  vs global, an identical 3,272 loci with identical claims, and `SVA,L1` where the global arm wrote
  `L1,SVA` — because the staged barrier populated the map in `cns_filter`-completion order. Every
  call was identical, so **only a comparison of the universe artifacts themselves could see it**;
  that is why the gate compares them separately from the calls. The canonical order is the driver's
  family order. It is not sorted inside `LocusUniverse` on purpose — sorting would silently rewrite
  every table already gated.
- **A concurrency change's failure mode is usually a byte, not a call.** The staged barrier moved no
  call in any arm at either scale, and its one real defect was invisible to every call-level check.
  When staging or reordering work, diff the *diagnostic* artifacts too, not only the deliverables.
Added 2026-08-04 (measured on the refactor's own WGS gate, job `49634290`):

- **THE TRANSDUCTION STEP IS TWO INDEPENDENT SCREENS OVER `.slct`, NOT A FUNNEL.** The clip gather
  writes `all_disc_pos` for the **whole** `.slct` list, and `collect_disc_reads_of_given_list` reads its
  site-list argument into a variable it never uses — so the disc evidence domain is **not** a subset of
  the clip-selected set. Subtracting one from the other (what "funnel" invites) yields a meaningless
  number. Measured: L1 2,802 clip-selected and 835 with disc evidence but only **157 with both**; SVA
  3,233 and 210 with **32 both**. That intersection is the real domain of everything after gate B, and
  **34,197 L1 sites become 4 records**. The 678 L1 / 178 SVA loci with disc evidence but no clip
  selection are why the reason `no_clip_evidence` exists. `docs/handoff/runs/probe_decision_funnel.py`
  recomputes all of this by set difference and asserts five domain assumptions.
- **Gate A does 97% of the rejecting; gates D and E reject one or two loci in a whole genome.** L1
  20,332 of 21,023 ledger rows are A_clip, SVA 27,590 of 27,739; `D_sanity` is 1 (L1) and 2 (SVA),
  `E_genotype` 1 each. **So the two-sided polyA/polyT veto and the large-indel check cannot be worth
  meaningful sensitivity** — two of §5's four "interrogate" items are near-inert at this population,
  and the leverage is gate A's clip-support and polyA-signature cutoffs.
- **`decoy_cross_family` fires 144 times for L1 and ZERO for SVA — and the reason is NOT biology.
  ANSWERED 2026-08-04: the SVA arm of that gate cannot fire, and a second, larger defect sits behind
  the same line of code.** The decoy test in `parse_transduction_disc_algnmt`
  (`x_transduction.py:1160`) is a **substring** test on the flank record's name, with
  `TD_DECOY_SVA = "SVA"`. Two independent consequences, both measured:

  1. **The hg38 SVA flank library contains no decoy records at all** — and this half is a LIBRARY
     PACKAGING difference that **varies by build**, so state the build. `hg38_FL_SVA_flanks_3k.fa` is
     5,454 flank records and nothing else, while **`hg19_FL_SVA_flanks_3k.fa` carries all 44**
     (`>ALU`, 35 `>ALU_Alu*`, `>LINE1`, `>SVA`, 6 `>SVA_[A-F]`) and chm13's SVA carries none. On the
     builds with none, `gate_decoy_cross_family`'s SVA arm reads `n_decoy_alu`/`n_decoy_l1` and no
     code path can increment either, so **BASELINE-0's zero is a dead arm and §5's SVA/Alu-competition
     hazard is UNTESTED, not refuted** — the previous version of this entry read the zero as evidence
     against §5 and that was wrong. **And chm13's L1 library carries no decoy records either**, so on
     chm13 even L1's cross-family tally is dead. Only hg38/hg19 L1 have the full 45.
  2. **Every reference SVA source flank is misclassified as the SVA consensus decoy and skipped — and
     THIS half is a code defect present on every build.** The SVA library's records are named
     `chr13~51077120~51078923~SVA_D~0L`, so `"SVA" in hit_flank_id` is True for all of them; line 1160's
     branch counts the hit into `m_SVA` and `continue`s before the source is extracted. Collision rate:
     **hg19 99.1% (5,114 of 5,158), hg38 100% (5,454), chm13 100% (6,072); L1 collides on none.** So
     **SVA transduction can only ever be called from POLYMORPHIC flanks** (`~polymerphic~`, which
     `gnrt_flank_regions_of_polymerphic_insertions` writes with `start == end`). Measured on job
     `49634290`: of L1's 159 picked sources 114 are reference copies and 45 polymorphic; of SVA's 105,
     **0 are reference copies and 105 are polymorphic.** A 5,454-record library built to be the
     reference SVA source set, every hit to it discarded.

  **THE CORRECT FORM EXISTS IN THE SAME FILE, BUT IN A DEAD FUNCTION** — so it is evidence the exact
  form was intended, not a live precedent. `parse_clip_algnmt_rescued` (`:1051`) tests `hit_flank_id ==
  TD_DECOY_LINE` — exact equality, not `in` — over the same four constants (it omits HERV), and it is
  one of §8's five zero-caller functions (verified: 0 non-def references). `check_algnmt_tag_decoy` (`:1285`) uses the
  substring form on XA-tag names, but its effect there is **benign**: an L1 flank id matches none of the
  four tags, and on an SVA run the mis-tallied hits land in `n_decoy_sva`, which SVA's own arm does not
  read. So the fix is `:1160` plus `:1285` for consistency, and the shape is already in the file.

  **It arrived from master byte-identical** (`46520ab`: same four constants, same substring test), so
  xTea2 did not break it — and that is now a MEASUREMENT on master's own output rather than a reading
  of master's source. The probe's run half over master's NA12878 WGS run
  (`master_golden_wgs_seeded/NA12878`, job `48825977`): **SVA 92 picked sources, 0 reference copies, 92
  polymorphic**, against **L1 135 with 84 reference copies**. A different sample and a different
  pipeline, same signature.

  **THE CODE HALF IS FIXED, 2026-08-04, BEHIND A SWITCH THAT SHIPS OFF** (`55055ba`). The
  classification is `transduction_evidence.decoy_tags(name, b_structural)` plus `is_flank_record`, one
  function asked by both the primary-hit test (`x_transduction.py:1160`) and `check_algnmt_tag_decoy`'s
  XA-tag copy of it (`:1285`), so one run cannot hold two definitions of "decoy". `b_structural=False`
  is the pre-fix expression verbatim — pinned against a copy of it in
  `tests/test_transduction_evidence.py`, because anything shipping off has to be provably inert or its
  default is not a rollback. `run_xtea --td-decoy-structural` sets `XTEA_TD_DECOY_STRUCTURAL`, which
  crosses the x_TEA_main subprocess boundary the same way `XTEA_SCAN_TD_GATHER` does.

  **The default is the FROZEN behaviour, not the correct one, and that is deliberate**: turning it on
  moves SVA calls, so BASELINE-0 and every gate frozen against it stay valid and the flip is its own
  measured decision. How much it unlocks was **not derivable from existing artifacts** — gate B's SVA
  domain was 210 loci, all through polymorphic flanks — which is why it needed a run: see §4 Stage 3's
  gate results for what the two gates measured.

  **Two things the fix does NOT change, and one of them bounds what it can ever prove.** L1 is
  untouched by construction (0 collisions in both L1 libraries), which is what makes an L1-identical
  WGS arm a check on the fix rather than a coincidence. And the **packaging** half is still open: an
  hg38 SVA flank library with no decoy records means SVA's arm of `gate_decoy_cross_family` still has
  nothing to count, so §5's SVA/Alu hazard remains untested even with the fix on. That half now has an
  exact recipe — the 45-record decoy set is reconstructible byte-identically from `rep_lib/consensus/`
  on any build, and `consensus/transduction_decoy.fa` (3 records) has zero references anywhere in
  xTea2 or master. See §5.

  `docs/handoff/runs/probe_td_source_classes.py` reproduces both halves in seconds from the flank
  libraries plus a finished run, and is the before/after instrument for that arm. It now **imports**
  both shipping predicates instead of writing its own copy, so a collision count it prints is a
  statement about the step; it reproduces every number in the table above against them.
- **`NEARBY_CLIP=50` lives in `check_seprt_disc_consistency` — and this entry's claim that it is "the
  step's largest threshold-driven rejection after gate A" was WRONG IN KIND. CORRECTED 2026-08-04 by
  instrumenting it (gate `[B2]`, `ff7cde1`, WGS job `49703322`).** The 135 of L1's 159 picked sources
  and 49 of SVA's 105 are the right *population*, but they are **not** threshold rejections:

  | | L1 | SVA | total |
  |---|---|---|---|
  | `no_cns_disc_reads` — an ABSENCE of evidence | 113 | 33 | **146** |
  | `no_cns_disc_cluster` — a threshold | 22 | 16 | **38** |

  So the threshold count is **38 genome-wide, not 184**, and the corrected ranking after gate A is
  `B_source` 781, `C_candidate` 54, **[B2]'s thresholds 38**, `D_sanity` 3, `E_genotype` 2 — [B2] sits
  **below gate C**. `NEARBY_CLIP` cannot touch the absences (they are decided before [B2] sees the
  locus), so its reachable population is **38 rejections + 80 survivors**. It is a small lever.
  **The obstacle this entry named was real and shaped the solution**: the function is shared with
  `cns_filter`, so its body stays in `x_clip_disc_filter.py` and only its predicate moved, into
  `transduction_evidence.gate_cns_disc_cluster`, with a **caller-supplied** ledger.

  **AND `NEARBY_CLIP` IS NOT MERELY A SMALL LEVER, IT IS MEASURED INERT — swept 0 → ∞ through the
  shipping gate, it moves ONE locus in the whole genome, and only at ∞.** That "38 rejections + 80
  survivors" is a bound on what the constant could in principle reach; the realised effect at every
  value from 0 to 500 is **zero** (L1 24 survivors, SVA 56), and 57 for SVA only above the largest
  observed distance. `docs/handoff/runs/probe_nearby_clip.py` drives
  `check_seprt_disc_consistency` itself at each value over job `49703322`'s own inputs and
  **reproduces that run's ledger B2 split exactly at 50** (113 + 22 / 33 + 16), which is what makes it
  a measurement rather than a model. The reason is in the distribution §5 asked for: over the disc
  records clipped on exactly one side, `|anchor_pos - ins_pos|` is **bimodal** — L1 median 3 with p75
  411, SVA median 129 with p75 358 — and 50 sits in the gap between the modes. Two facts that make
  this computable at all: `prep_representative_pos` is the **identity** (`x_transduction.py:710`, so
  `clip_pos` IS `ins_pos`, despite the name), and both quantities are in the cns SAM's read names.
  ~~**So this item is closed. Do not open an instrumented run for it.**~~

  > **SCOPE CORRECTED 2026-08-05, AND THE CORRECTION IS LARGE: EVERY SENTENCE ABOVE IS ABOUT GATE
  > `[B2]`, AND THE SAME CONSTANT IS THE BIGGEST LEVER MEASURED SO FAR AT GATE `[A]`.**
  > `NEARBY_CLIP` has a **second** use the entry above never names: `x_transduction.py:865` skips a
  > clip realignment whose `|ref_mpos - ori_insertion_pos|` exceeds it, **before** the unmapped and
  > clip-qualified tests — which is gate [A]'s *domain*, and gate [A] does 97% of the step's
  > rejecting. Measured through the shipping parse and the shipping gate at WGS
  > (`docs/handoff/runs/probe_clip_gate.py`, job `49836568`, both halves of check 0 reproduced):
  >
  > | at gate [A], WGS | L1 at 50 | L1 unbounded | SVA at 50 | SVA unbounded |
  > |---|---|---|---|---|
  > | qualified clip records | 99,245 | **415,429** (4.2×) | 128,262 | **519,097** (4.0×) |
  > | domain (loci) | 23,134 | 31,789 | 30,823 | 38,212 |
  > | survivors at the shipping cutoff | 2,802 | **10,339** (3.7×) | 3,233 | **12,548** (3.9×) |
  >
  > **This is NOT a recommendation to raise it**, and two reasons why not are in the measurement:
  > part of the growth is the same reads re-assigned to nearby loci rather than new evidence, and
  > gate [A]'s survivors still face [B], [B2], [C], [D] and [E] — where 306 of BASELINE-0's 322
  > false negatives already fail to appear at all. What it says is that **"measured inert" is a
  > claim about [B2] and gate [A] is unmeasured**, so the bimodal distance distribution quoted
  > above — which was taken over the *disc-to-consensus* records — is not the distribution this
  > second use consults. Any move needs a truth-set-scored two-arm run.
- ~~**`MIN_DISC_MAPPED_RATIO = 0.7` IS THE LARGEST UNEXAMINED KNOB BETWEEN GATES [B] AND [C]** …
  is what produces the **146** `no_cns_disc_reads` loci … **Interrogate this before `NEARBY_CLIP`.**~~
  **MEASURED 2026-08-04 AND THE CLAIM IS WRONG IN KIND — it is a SMALLER lever than `NEARBY_CLIP`, not
  a larger one.** This entry made the same mistake as the one directly above it: it attributed a
  population to the threshold nearest to it. `parse_disc_algnmt_consensus` has **two** filters and they
  are not comparable — an unmapped skip and a mapped-fraction cutoff — and the population belongs
  almost entirely to the first:

  | | L1 | SVA | total |
  |---|---|---|---|
  | `no_cns_disc_reads` loci | 113 | 33 | **146** |
  | …of which **every** record came back unmapped from bwa — **unreachable at any cutoff** | 89 | 23 | **112** |
  | …of which ≥1 record mapped and failed the ratio — **the cutoff's ceiling** | 24 | 10 | **34** |
  | loci actually recovered at 0.5 (from the shipping 0.7) | +4 | +6 | **+10** |
  | loci lost at 0.9 | −4 | −8 | −12 |

  So its reachable population is **34**, and the realistic movement is **10 loci genome-wide**. It is
  not the largest knob in this region; the largest thing here is that **90.2% (L1) and 84.9% (SVA) of
  disc-to-consensus alignment records come back unmapped**, which is not a threshold at all. (It is
  still the *bigger* of the two knobs in this region: `NEARBY_CLIP` was swept the same way and moves
  **one** locus — see the entry above.)

  **THE ENABLING CORRECTION, and it unblocks more than this entry: the disc-to-consensus SAM is NOT
  deleted at the end of the step.** The handoff chain said it was, which is why both this knob and
  `NEARBY_CLIP` were filed as needing an instrumented run. It is deleted by
  `run_xtea.cleanup_successful_run` — at the end of the RUN, and only without `--keep-intermediates`,
  which **every gate under `docs/handoff/runs/` passes**. So
  `<fam>/tmp/transduction/temp_transduction_disc_cns.sam` (1.4 MB per family at WGS) is sitting in every
  gate's output tree, and so are the 393 MB disc-to-flank SAM and the 2.46 GB clip-to-flank SAM.
  `docs/handoff/runs/probe_disc_cns_mapped_ratio.py` reads the first, imports the shipping predicate
  rather than copying it, and **reconciles exactly against the ledger** (113 and 33 by two derivations,
  a SAM and a TSV). `NEARBY_CLIP` is measurable the same way, from the same file.

  **And one more thing that file says, which is a perf observation and not a sensitivity one:**
  **56% of what the third bwa pass aligns was already unmapped against the flank library** (2,784 of
  L1's 4,969 focal reads), and those reads contribute **3 of the 456 records that survive to gate
  [B2]**. So the pass and its parse handle ~2.2× the reads that can matter — but filtering them would
  move 3 records, so it is call-moving at the margin rather than a free win.
- **Gate [B2] attributes ZERO false negatives** (WGS job `49703322`): the FN table is unchanged at 306
  `no_ledger_row` / 8 `low_clip_support` / 8 `no_polya_clip` = 322. So [B2] is not a sensitivity gate
  either, and §4.6's conclusion — the ceiling is upstream of this decision — survives its instrumentation.
- **`check_seprt_disc_consistency`'s `b_transduct` parameter is DEAD IN ALL THREE CALLERS** and always
  has been. The two `cns_filter` call sites omit it (default False) and the transduction one passes
  `False` explicitly; `git log -S` finds no commit that ever passed True, so it arrived dead from
  master. Consequence for instrumenting that function: only one behaviour to reason about. Turning the
  arm on is a **call-moving** change — a per-side minimum of `ndisc_cutoff/2` for one-sided clusters,
  and since the transduction caller passes the already-halved cutoff that is a quarter of `--nd`. Its
  effect is not measurable from existing artifacts, which record only chrom/pos for survivors. **The
  arm is now IMPLEMENTED in `gate_cns_disc_cluster` and switched off** (`b_one_sided_cns_min`, False
  in every caller), with the reason `low_one_sided_cns_disc` in the vocabulary and unreachable, so
  turning it on is a one-line params change that arrives with its own attribution.
- **THE SENSITIVITY CEILING IS NOT IN THE TRANSDUCTION DECISION. 306 of BASELINE-0's 322 false
  negatives (95%) never reach it** (`benchmark/ledger_attribution.py`). All 16 attributed misses are in
  gate A; none in B/C/D/E. **And quote the EXCESS, not the raw 16**: a dropped locus lands within the
  200 bp match window of a *true positive* 12.1% of the time for L1, so ~6.9 of the 15 L1 attributions
  are chance — enrichment 2.17×, excess **~8 calls**. SVA's 1.18× (1 observed, 0.8 expected) is
  nothing. So Stage 3's transduction threshold work is a **precision** story, not a sensitivity one.

Added 2026-08-03 (the decision refactor — `transduction_evidence.py`):

- **The transduction decision cannot be one function, and the reason is THREE BATCHED PASSES, not
  code shape.** `filter_transduction`'s `calc_coverage_of_two_regions` runs over the loci that
  survived the candidate gate and `_collect_gntp_features`' genotyper over the sites that survived the
  sanity gate, so the decision is necessarily *gate → fetch → gate → fetch → gate*. Flattening it
  means computing coverage and genotype features for the whole 34-40K-row `.slct` list per family.
  **And there is a third barrier one stage earlier, added to this entry 2026-08-04**:
  `re_select_disc_reads` writes only the disc reads of gate [B]'s survivors and they are realigned to
  the CONSENSUS before gate [B2] can see them. So it is **six gate groups, three batched passes** —
  this entry said two fetches and five groups. §4 Stage 3's blocks have the diagram. **Do not
  re-derive this as a design question.**
- **A LEDGER PHASE NAME'S LEXICAL ORDER IS NOT ITS CASCADE ORDER, and this was a live trap rather
  than a hypothetical.** `"B2_cns_cluster" < "B_source"` is True in ASCII (`'2'` 0x32 before `'_'`
  0x5F), and `ledger_attribution.attribute` found "the gate that determined this false negative" with
  `min(rows)` on the phase STRING — correct only while A < B < C < D < E held by accident. Use
  `transduction_evidence.PHASE_ORDER` / `phase_rank`; never `sorted()` on a phase. Any new phase name
  must be added to `PHASE_ORDER`, and `test_every_reason_has_a_phase` fails if it is not.
- **`filter_transduction`'s two coverage windows are 200 (cov1) and 900 (cov2), and its comments said
  the opposite.** `focal_win = LOCAL_COV_WIN = 200`, `focal_win2 = COV_ISD_CHK_WIN = 900`; all four
  numbers are tested against `i_max_cov` and **cov2 is the one `collect_ins_info` reports as local
  coverage**. The comment is fixed; the constants are not touched, because swapping them moves calls.
- **Two defects in the decision are preserved deliberately and pinned by a test**: the cluster test
  inside `call_out_candidate_transduct`'s per-read loop (a locus can carry one side's flag from read 3
  next to the other side's boundaries from read 9, and "both sides true" is reachable only that way —
  `collect_ins_info` reads that pair), and `_pick_trsdct_rescued`'s discarded second map, whose
  accumulation sat inside an `if chrm not in map` and recorded only the first locus per chromosome.
  Neither is a cleanup: the first moves calls, the second no longer exists.
- **`_extract_transduct_info` is a SIXTH dead function** (zero callers, verified 2026-08-03) on top of
  the five §8 lists. Left in place — deleting it belongs to the §8 sweep, not to a decision refactor.
- **`scan_disc_pos.tsv` differs between two runs of the SAME commit as a PERMUTATION only** — the
  multithreaded scan writes rows in completion order. `cns_manifest.txt` differs because it holds each
  run's own absolute output paths. **Everything else — 98 files under the sample directory — is
  byte-identical between two arms of one commit** (measured on job `49599837`). That is what makes a
  full-tree byte comparison a usable gate for a pure refactor, and it is the strongest agreement gate
  available here; use it (`gate_transduction_evidence_chr1.sbatch`) rather than comparing only VCFs.
- **A concurrency check must not be order-dependent, or it passes for reasons the code does not
  control.** The tempting staged-barrier assertion — "did some family cross the barrier before the
  slowest `cns_filter` landed" — holds or fails depending on which family happens to be slowest and
  on whether the gather outruns the gap behind it. The structural one is the flankless family's own
  gap, from its enumeration line to its wave-2 announcement: under staging it contains only the
  union-table write, under the join it must contain the entire gather, **on any ordering**. That is
  `gate_staged_barrier_chr1.sbatch` check 5, and the difference is why it is check 5 and not the
  first thing tried.

Added 2026-08-04 (gate [B2] — the cns-cluster ledger):

- **chr1 CANNOT gate [B2]'s thresholds, and the reason is stronger than "small N".** chr1's gate B
  survivors are 5 (L1) and **0** (SVA), so [B2]'s domain is 5 loci on one family — and measured on job
  `49703321`, **all five drops are `no_cns_disc_reads`**, the ABSENCE reason. Not one of the five had a
  disc read survive the third bwa pass to the consensus, so at chr1 the gate never reaches a threshold
  at all. chr1 gates [B2]'s wiring, vocabulary and file-level containment; WGS gates its thresholds.
- **THE ABSENCE/THRESHOLD SPLIT IS NOT PEDANTRY — the set-difference derivation gets it wrong.**
  `picked sources` minus `cluster-consistent` is the right *population* but says nothing about why, and
  on chr1 it would have reported 5 threshold rejections where the true answer is 5 absences. Any
  statement of the form "gate [B2] rejected N loci" must come from the ledger's reasons, not from
  subtracting two site files.
- **A ledger is not excludable from a gate once BOTH arms write one.** The parent refactor could treat
  it as a new file; a change that only *adds rows* cannot. The usable check is exact rather than
  tolerant: **AFTER's ledger with the new phase's rows deleted must be byte-identical to BEFORE's**,
  which asserts in one comparison that rows were added and that nothing else moved — including row
  ORDER, since a new phase's rows land between two existing phases'. `gate_cns_cluster_ledger_chr1.sbatch`
  check 1. Measured there: 99 files identical, 1 ledger reconciled, 0 differing.
- **`prep_representative_pos` is lossless** — one entry per `(chrm, pos)`, no filtering — so gate
  [B2]'s domain is exactly the picked-source set and its drop count is `picked − cluster-consistent`
  with no slack. That is what makes the site-file number a real prediction for the ledger to hit
  rather than an approximation.
- **`grep -q X && { …; exit 1; }` is a landmine under `set -e`** and it is the *negative* preflight
  form that bites: when the grep correctly finds nothing, the AND-list's status is nonzero and the
  script dies on the case it was written to accept. Write it as `if grep -q …; then … fi`. Third shell
  footgun recorded in these gates, after `local a=$1 b="$x/$a"` and `f=$(ls MISSING | head -1)`.
