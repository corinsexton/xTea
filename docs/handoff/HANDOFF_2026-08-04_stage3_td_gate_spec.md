# xTea2 handoff — 2026-08-04, Stage 3's transduction gate: the spec

**This is the entry point for ONE decision: how Stage 3's transduction work gets gated.** It is not a
state handoff — for state read `HANDOFF_2026-08-04_sva_decoy_fix.md`, and for agenda read
`docs/handoff/DOWNSTREAM_OVERHAUL_PLAN.md`.

> **STATUS: SPEC, NOT AGREED — but §2, §3 and §7 steps 1-3 ARE NOW BUILT AND MEASURED.** Corinne
> directed (2026-08-04) that Stage 3's transduction gate should assume **truth from a run of master
> xTea on HG002**, which is the input the plan was missing — until then §7.1's conclusion was "the ALT
> route is exhausted; what remains is the paper's supplement or a different sample." This document is
> the design that direction implies. The pass/fail thresholds in §5-6 are still the part most likely
> to need Corinne's judgement rather than mine.
>
> **Read `HANDOFF_2026-08-04_td_gate_built.md` for what was built and the four numbers in this
> document that measurement changed.** The corrections are also inline below, each marked
> **[CORRECTED]**, so this document cannot be read straight without meeting them. The largest: §2's
> "241 reference-copy SVA sources" is measured at the *step's* level and the pool contribution is
> **4** — see §2's table and §5.

**The one-sentence version: master enumerates the hypotheses, the reads return the verdict, and the
metric is precision per stratum — not sensitivity, and not agreement with master.**

---

## 1. Why this gate has to exist, and why it can only measure precision

Stage 3 changes decisions, so agreement-gating stops being coherent (plan §4 Stage 3). Its stated gate
is "per-family sensitivity/PPV against the HG002 truth set" — which **works for the call set as a whole**
(BASELINE-0 scores 80.4% / 88.4%) and **does not work for transductions**, because the truth-side
denominator is **0 for L1 and 2-3 for SVA**, established by two independent methods and declared
exhausted (plan §7.1, jobs `49482568` / `49482779`).

**Precision is the right target anyway, and this is measured rather than a consolation.** The reason
ledger says the sensitivity ceiling is not in this decision at all:

- **306 of BASELINE-0's 322 false negatives (95%) never reach the transduction step.**
- All 16 attributed misses are in gate **[A]**; with the base-rate control the excess over chance is
  **~8 L1 calls** and SVA's is ~0.
- **Gate [B2] attributes zero** (WGS job `49703322`) — instrumenting the step's second-largest
  population did not move the FN table.

So sensitivity for MEI calling lives upstream in `cns_filter`; what the *transduction decision* controls
is whether the transductions it asserts are real. That is precision, and a precision gate is matched to
the question.

## 2. The pool — three arms, and the third one is not optional

The pool is the union of `(insertion locus X, source Y)` pairs from the `TD_SRC` INFO field.

| arm | why it is in the pool |
|---|---|
| **master on HG002** | the new hypotheses; the only arm not descended from xTea2 |
| **BASELINE-0** (xTea2, frozen) | the reference point the plan already gates deltas against |
| **xTea2 with the SVA decoy test made structural** | **without it the pool contains ZERO reference-copy SVA sources**, so the gate cannot score the fix that matters most — see §5. **THIS ARM NOW EXISTS AND IS MEASURED**: `run_xtea --td-decoy-structural` (`55055ba`, off by default), and on HG002 WGS it contributes **241 reference-copy SVA sources** against **0 of 105 in BASELINE-0 (same sample)** and **0 of 92 in master's own NA12878 output (different sample, same signature)** — jobs `49724055`/`49724056`, both PASS. Its run tree is on disk, so this arm needs no new compute. **[CORRECTED] The 241 is the STEP's level, not the pool's.** It is counted in `unique_trsdct_disc_only`, the transduction step's own hypotheses before downstream filtering; this section defines the pool on `TD_SRC` in the call set, and there the same arm contributes **4** reference-copy SVA sources against BASELINE-0's 0 — of which **3 are orphan-class and 1 is regular**. Counted on both: L1 159 picked → 39 non-sibling `TD_SRC`, SVA 346 → 5, BASELINE-0's SVA 105 → 1. **These are NOT nested populations, so the ratio is not a filtering rate** — `x_transduction.py:238` keys the picked-source file on the step's `--input2`, the RAW disc-site universe (**36,687 rows for L1** against a 287-row call list; 104 of its 111 L1 loci are exact matches there, and only 4 of the VCF's 59 non-sibling transduction calls fall within ±200 bp of a picked row). Which makes the point stronger: 241 cannot be the pool contribution under *any* filtering story. 241 sizes the step, where Stage 3's *thresholds* act; 4 sizes what the call set asserts, which is what a precision gate scores (§1). `probe_td_pool.py` prints the reconciliation on every arm |

**Merge loci at ±200 bp, one-to-one**, the same window and matching discipline `benchmark_truth_set.py`
uses (`DEFAULT_POS_WINDOW`), so one call cannot satisfy two pooled entries. Plan §7.3 records that
greedy-vs-nearest and one-to-one-vs-not "moves sensitivity/FDR by percentage points"; the harness's
choice is the one to keep.

> **[CORRECTED] "One call cannot satisfy two pooled entries" also has to hold WITHIN an arm, and the
> harness's matcher cannot give you that.** `match_one_to_one` enforces one-to-one only *between* two
> sets, so a single arm's entries are never matched against each other and a caller that emits the same
> locus twice puts two entries in the pool. Master's NA12878 arm has one such pair —
> `chr9:87033714` and `chr9:87033895`, **181 bp apart with the identical source `chr3:-1-32916262`** —
> which the verifier then refuted twice, as though it were two failed hypotheses. `probe_td_pool.py`
> now collapses intra-arm duplicates (same locus within the window **and** the same source) and reports
> the count. **Master's pool is 114 entries; the 115 below is a CALL count, which is a different thing
> and stays 115.**

### 2.1 `TD_SRC`'s shapes, measured — not guessed

> **[CORRECTED] There are FIVE shapes, not four, and the parser handled THREE of them.** Both
> `{source}~orphan` and `{source}~{cluster_l}-{cluster_r}` returned `None` from
> `benchmark_truth_set.parse_td_src` — indistinguishable from `not_transduction` — so 7 of master's
> 115 NA12878 assertions, 3 of BASELINE-0's and 5 of the decoy arm's were invisible, and reached the
> scorer as `TP_POS_SOURCE_MISSED` for calls that named their source perfectly well. The grammar, from
> both writers: `transduction_evidence.call_candidates` writes `{trsdct_source}~{l}-{r}` or
> `{trsdct_source}~orphan`, and `trsdct_source` is itself `chrom:start-end~rc` — so the source is the
> **first** `~`-field after the colon and everything after it is provenance. Fixed by stripping the
> tail (`2d38feb`). On BASELINE-0 it moves exactly one reported number: L1's `FP_TRANSDUCTION`, 0 → 1.
>
> **[CORRECTED — THE LARGEST ONE] THE TABLE BELOW'S `reference-copy source` COLUMN IS READING THE
> WRONG FIELD, and §7 step 3's validation target does not exist.** `TD_SRC` has TWO writers.
> `x_transduction.py:1196` (the transduction step) writes
> `"{chrm_fl_L1}:{source_start}-{source_end}~{sourc_rc}"` from a **flank record** — a real source
> element — and `call_candidates` always appends a tail, so that shape is never bare.
> `x_orphan_transduction.py:365` (the sibling / orphan module) writes a **bare**
> `"{src_chrm}:{i_start}-{i_end}"` whose coordinates its own comments at `:345-346` call the
> **"mid-position of left discordant pair"** and **"mid-position of right discordant pair"**.
> Byte-identical in master at the same lines. So on the bare shape `start < end` is a disc-evidence
> span and says **nothing** about copy status.
>
> Measured, three ways: bare spans are 1–2814 bp (median 28) where a full-length L1 is ~6 kb; **none
> of master's 33 bare intervals lie inside or within 5 kb of any source region** in
> `unique_trsdct_disc_only`; and in the tailed shape the first interval reaches 6030 bp, which is what
> the field looks like when it really carries a source region.
>
> `td_src_source_class` now returns a fourth value, **`disc_span`**, not a subtype of the others. So
> master's L1 row is **56 regular `disc_span` + 3 orphan reference-copy + 27 sibling**, and only
> **7 of its 115 assertions** have a flank-derived source interval at all (BASELINE-0: **3 of 90**).
> **§4's source-class stratum is therefore defined for 7 of 115 entries** — a far harder limit than
> §6.2's. §7 step 3's "the 33 L1 reference-copy sources … have real intervals" is false; what the
> verifier actually validated on is the two-sided bare stratum, 26 of 31 = 83.9%.
>
> **The SVA decoy finding is NOT affected**: `0 of 105 → 241 of 346` is measured on
> `unique_trsdct_disc_only` column 3, which IS `x_transduction.py:1196`'s flank-derived id.
>
> **[CORRECTED] "orphan" does not imply "polymorphic", and the table below therefore undercounts L1's
> reference copies.** Master's 4 SVA orphans are polymorphic (`start == end`); its 3 L1 orphans are
> reference copies spanning ~6 kb (`chr4:79966908-79972933~0~orphan`). Call class and source class are
> **independent axes** — which is why §4's two strata are two columns. (An intermediate version of
> this note said L1's reference copies are 36 rather than 33; with the block above, the answer is
> **3** — the tailed orphans — and the other 33 were never reference copies.)
>
> **[CORRECTED] A polymorphic source IS verifiable.** What it lacks in the reference is the element
> BODY; the transduced sequence comes from downstream of its 3′ end, which is ordinary genomic
> sequence — literally the window `gnrt_flank_regions_of_polymerphic_insertions` cuts its flank from,
> which is *why* it writes `start == end`. Master's 4 SVA orphan-polymorphic entries verify 4/4.

`x_gvcf.py:148` declares the field and `:385` writes it; **master emits it identically** (`master:xtea/x_gvcf.py:147`, `:321`) and its real output confirms it. Counted on master's existing
NA12878 WGS run (job `48825977`, the run that already completed at 50:18):

| family | calls | assert a transduction | sibling | reference-copy source | polymorphic | unknown-start |
|---|---|---|---|---|---|---|
| L1 | 203 | **86** | 27 | **33** | 0 | 23 |
| Alu | 1262 | **0** | — | — | — | — |
| SVA | 73 | **29** | 25 | **0** | 4 (orphan) | 0 |

**115 transduction-asserting calls per WGS sample**, against the ALT route's 0-3. That is the whole case
for this direction.

Four shapes, all already handled by `benchmark_truth_set.py:207`'s parser:

- `not_transduction` — excluded.
- `chrom:START-END` — a source region. **`START == END` means a polymorphic source**, `START < END` a
  reference copy (`gnrt_flank_regions_of_polymerphic_insertions` writes polymorphic flanks with both
  coordinates equal to `ins_pos`).
- `chrom:-1-END` — unknown start. **Benign for this gate**: transduced sequence comes from *downstream of
  the element's 3′ end*, so `END` is the coordinate the verifier needs and the missing start costs
  nothing. 23 of master's 86 L1 assertions are this shape.
- `chrom:POS~sibling` — a position, not a region. **Use xTea's own convention: `[POS-100, POS+100]`**,
  which is literally what `x_orphan_transduction.py:1190-1191` does when it consumes the field. Do not
  invent a window.
- `chrom:POS-POS~0~orphan` also appears (master's SVA). Orphan, polymorphic.

## 3. The verifier — this is the actual deliverable

`benchmark/verify_transduction.py`. For each pooled `(X, Y)`, decide from the **BAM alone** whether the
insertion at X carries sequence from downstream of Y. No caller in the loop; that is what makes the
result truth rather than agreement.

1. **The source's transduced window** is `[Y_end, Y_end + W]` with **`W = flklen = 3000`** — the step's
   own flank length (`run_xtea.py`'s golden defaults `--nfclip 4 --nfdisc 5 --flklen 3000`), so the gate
   and the caller mean the same thing by "downstream". For a sibling, `Y_end = POS + 100` per §2.1.
2. **The reads at X**: clipped reads whose clipped part does not map at X, plus discordant mates, inside
   the scorer's ±200 bp window of X.
3. **Align those parts twice** — to (a) the source window's reference sequence and (b) the family
   consensus (`consensus/*.fa`). A read supports the transduction only if it aligns to **(a)** and
   **not** to **(b)** over ≥80%. This is the "explained by consensus is not residual" discipline §7.1
   already had to learn on the truth side, applied to calls.
4. **Reject a source window that is itself a same-family repeat.** Use `annotation1` (`{g}_L1.fa.out`,
   the index transduction is actually handed — plan §10) and drop windows sitting inside a
   low-divergence copy, because there (a) and (b) are the same sequence. **This is the exact trap §7.1
   hit**: 11 of 20 L1 residuals ≥100 bp mapped *inside* a reference L1 and every one was an older
   subfamily, so "residual" was divergent body. A verifier without this test will confirm transductions
   that are element body.
5. **Three verdicts, and the third is not a rounding error.** `VERIFIED` (≥N supporting reads),
   `REFUTED` (0), `UNDECIDABLE` (source window in a same-family repeat, or local coverage below a
   floor). **`UNDECIDABLE` is never folded into either** — the same rule the ledger enforces between
   `no_ledger_row` and "unattributed", and between `no_cns_disc_reads` and a threshold. Getting that
   wrong is how a set difference produced 184 "threshold rejections" where there were 38.

   > **[CORRECTED] There are FOUR states.** This item names `VERIFIED (≥N)` and `REFUTED (0)` and
   > leaves `0 < support < N` unnamed. `SUB_THRESHOLD` is its own state in the implementation:
   > folding it into `REFUTED` lets the choice of N manufacture false positives, folding it into
   > `VERIFIED` lets N do nothing. It is also the clearest reason §3.1 is right — every
   > `SUB_THRESHOLD` entry becomes `VERIFIED` at a lower N. Neither it nor `UNDECIDABLE` is ever in
   > the precision denominator.

### 3.1 The verifier's own thresholds are a hazard — report a sweep, not a point

> **[CORRECTED] The parameter that actually decides the answer is not in this section's list.** Over
> master's NA12878 pool (114 entries, 324 sweep points), `N`, the consensus fraction, the source
> fraction and the coverage floor each move precision by **≤ 2 points** across their whole range. The
> **repeat-divergence cutoff** of §3.4 moves it **10 points** (67.1% → 77.4%) and moves `UNDECIDABLE`
> from 24 to 52. It is swept too, and it is the one to argue about.
>
> On the HG002 decoy arm it is not influential but *decisive*: the fix's only **regular**-class SVA
> reference-copy entry is `UNDECIDABLE` because its source is a **2.7%-divergent SVA_D**, and a second
> is 3.8%. SVA is a young family, so its reference copies sit at low divergence and this test
> disqualifies most of them **by construction**. See §5.

`N`, the MAPQ floor, the 80% consensus fraction and the coverage floor are **as arbitrary as the
constants this gate exists to interrogate**. A single operating point would just relocate the problem.
So the verifier reports its verdict **as a function of them** over a small sweep, and the gate quotes the
range. For calibration the codebase's own values are the natural starting points:
`MINIMUM_DISC_MAPQ = 20`, `F_MIN_TRSDCT_DISC_MAP_RATION = 0.65`, `MIN_DISC_MAPPED_RATIO = 0.7`.

## 4. Scoring — stratified, because the classes are not the same claim

Per family × call class × source class:

- **call class**: regular (region source) / **sibling** / **orphan**. Not cosmetic — BASELINE-0's L1
  orphan calls are **0/37 on-truth** against GIAB, because an orphan transduction site is where
  transduced sequence landed *without* the element body and is therefore not an MEI insertion at all
  (plan §1). Scoring them together with regular calls hides a category error.
- **source class**: reference copy / polymorphic / unknown-start.
- **metric**: `precision = VERIFIED / (VERIFIED + REFUTED)`, with `UNDECIDABLE` reported beside it and
  **never in the denominator**.

## 5. Pass/fail — and the one rule that keeps the gate from rejecting the fix it exists to enable

- Gate on **no regression in `VERIFIED` count and no regression in precision**, per family, beyond a
  margin calibrated on the pool's own three-arm agreement.
- **The SVA decoy fix is scored as GROWTH in `VERIFIED`, not as a precision delta.** Recovered
  reference-copy SVA sources are new true positives. A gate built the naive way — master's calls as
  truth — would score every one of them as a **false positive and reject the fix**, because master
  shares the bug (`master:xtea/x_transduction.py:1127` is the identical substring test and it reads the
  identical `hg38_FL_SVA_flanks_3k.fa` via `gnrt_pipeline_local_v38.py:803`). Master's own SVA output
  proves it: every source region it reports is `start == end`, plus 25 siblings, and **zero
  reference-copy sources**. L1 is the control at 33.

  > **[CORRECTED] THIS RULE IS RIGHT AND IT HAS ALMOST NOTHING TO SCORE — the decision is Corinne's.**
  > At the pool's own level the fix contributes **4** reference-copy SVA entries, not 241 (§2), and the
  > verifier adjudicates them: **2 VERIFIED (both orphan-class), 2 UNDECIDABLE** — the two undecidable
  > ones because their source windows sit in SVA_D copies at 2.7% and 3.8% divergence. The single
  > **regular**-class entry is one of the undecidable pair. So:
  >
  > - the GROWTH claim is supportable at **+2 VERIFIED, both orphan** — and §4 says orphan sites are
  >   not MEI insertions at all, so that is the weakest class to rest a fix on;
  > - **§6.2's warning that SVA's regular class has N ≈ 0 survives the fix rather than being retired
  >   by it**, which is the opposite of what this section assumed;
  > - and the blocker is §3.4's repeat test interacting with SVA's youth, not the fix.
  >
  > Three ways out, none of them mine to pick: loosen the divergence cutoff for SVA specifically and
  > say so; accept orphan-class evidence as sufficient for the growth claim; or score the fix at the
  > *step's* level (the 241) with a separate, non-precision argument. **Decide before the pool is
  > frozen** — §7 step 5.
- **BASELINE-0 delta reported for attribution, not pass/fail** — the plan's existing Stage 3 rule,
  unchanged.

## 6. What this gate does NOT claim — put this in its docstring

1. **Not sensitivity.** The pool is bounded by what some caller found; it cannot say what all three
   missed, and therefore cannot state how many transductions HG002 has.
2. **SVA gets almost no *regular* transduction truth.** Master's SVA is 25 sibling + 4 orphan + **0**
   reference-copy, so SVA's regular-class precision has N≈0 until the decoy-fixed arm contributes. A
   per-family table that does not say this invites the reader to trust an empty cell.
   **[CORRECTED] It stays N≈0 AFTER the decoy-fixed arm contributes.** The arm adds exactly one
   regular-class SVA reference-copy entry and the verifier returns `UNDECIDABLE` on it (§3.1, §5). The
   empty cell is still empty, and `verify_transduction.py` prints `EMPTY DENOMINATOR -- this cell says
   nothing` under it rather than a `0.0%` a reader could mistake for a rate.
6. **[ADDED] Not a verdict where the source window is a young same-family copy.** §3.4 is a real
   epistemic limit, not a tuning knob: for SVA it disqualifies most reference-copy sources because SVA
   is young. The verifier reports those as `UNDECIDABLE` with the divergence and the subfamily name, so
   the limit is visible per entry.
3. **Alu contributes nothing**, by construction — no transduction flank, so the step never runs. Biology,
   not a gap.
4. **NOT master agreement.** Master enumerates; reads decide. This is what keeps the design inside plan
   §2's locked decision — *"agreement with master is explicitly not a requirement"* — which scoring
   against master's calls directly would contradict.
5. **Not comparable to the paper's numbers** — hs38DH vs `no_alt` (plan §7.1).

## 7. Order of operations, and the prerequisite that will bite if ignored

1. ~~**Repoint and run master on HG002.**~~ **IN FLIGHT / DONE 2026-08-04 — job `49763841`, and it is
   NOT a repoint of `master_serial_wgs.sbatch`.** That script opens `rm -rf "$OUT"` over the NA12878
   tree this whole document cites, so a repointed invocation that forgot `MASTER_OUT` would delete the
   evidence. `docs/handoff/runs/master_wgs_hg002.sbatch` is a separate file with a refusal guard on
   `OUT`. It is still **master's OWN pipeline end to end** — the golden's per-type
   `run_xTEA_pipeline.sh` files (themselves the output of master's `gnrt_pipeline_local_v38.py`)
   transformed in place by `sed`, so exactly one thing differs: the BAM. §10's two non-interchangeable
   flavours are respected. **The golden's `--lc/--rc/--nfclip/--nfdisc` needed no repointing**:
   `x_TEA_main` sets `b_automatic = True` unless `--user`, which overrides them from the BAM's own
   coverage (`e38d0f5`), so master self-adjusts to HG002. Out:
   `/n/scratch/users/c/cos689/xtea2_testdata/master_wgs_hg002/HG002`. L1 took 30:05.
2. **Measure HG002's pool size.** The 115 above is **NA12878**; HG002's is unknown. **The tool now
   exists** — `docs/handoff/runs/probe_td_pool.py`, which reproduces §2.1's NA12878 table exactly as a
   ledger check and BASELINE-0's as a second. The two HG002 arms already on disk measure **90**
   (BASELINE-0) and **94** (decoy-structural), merged pool **94**; master's arm is step 1.
3. ~~**Build the verifier**~~ **DONE — but THIS STEP'S PREMISE IS FALSE.** `benchmark/verify_transduction.py`
   (`9d8e1de`). "The 33 L1 reference-copy sources … have real intervals" is wrong: they are the bare
   shape, i.e. disc-pair midpoints (§2.1's largest `[CORRECTED]` block). The verifier is validated on
   the stratum that does exist — the two-sided bare entries, over master's own NA12878 WGS BAM: **26
   VERIFIED / 5 REFUTED / 1 SUB_THRESHOLD / 1 UNDECIDABLE — 83.9%**. The whole pool -- 114 entries, see §2.1's intra-arm note -- runs at
   **65.4% – 79.0%** over 324 sweep points. Per stratum: L1 regular reference-copy 83.9%, L1 regular
   unknown-start **55.0%**, L1 sibling 83.3%, SVA orphan-polymorphic 100% (4/4), SVA sibling **41.2%**.
   An unknown start costs 29 points and SVA's siblings are the weakest class in the table — which is
   §4's argument, measured.
4. ~~**The SVA decoy fix must exist as a branch BEFORE the pool is frozen.**~~ **DONE 2026-08-04** —
   it is a flag, not a branch (`--td-decoy-structural`, `55055ba`), and its WGS arm has already run:
   `/n/scratch/users/c/cos689/xtea2_testdata/gate_decoy_structural_wgs_49724056/HG002.wgs`, with 241
   reference-copy SVA sources in `TD_SRC`. So this prerequisite costs zero compute now; do not re-run
   it. See `HANDOFF_2026-08-04_sva_decoy_fix.md` §5.2.
5. **Freeze the pool, score it — that is the gate's baseline.** Then Stage 3's threshold work is
   gateable.

## 8. Open questions I could not settle from the code

- ~~**HG002's pool size**~~ **ANSWERED for two of three arms**: BASELINE-0 **90** transduction-asserting
  calls, the decoy-structural arm **94**, merged pool **94** with 4 unique to the decoy arm and 0 unique
  to BASELINE-0. Master's arm is §7 step 1. The earlier note here — "346 picked sources of which 241 are
  reference copies" — is the STEP's level, not the pool's; see §2's `[CORRECTED]` row.
- **Whether a sibling's `POS` is the source element's 3′ end or its insertion point.** §2.1's `±100`
  convention is what xTea's own consumer does, which makes it the right *default*, but it is a window
  around a point rather than an oriented 3′ end — so for siblings the verifier searches both directions
  or accepts a looser claim. Worth resolving before N is calibrated, because 25 of master's 29 SVA
  assertions are siblings.
  **PARTLY ANSWERED, FROM DATA.** The verifier searches both `[end, end+W]` and `[start-W, start]` and
  reports which side supported. Over master's NA12878 pool: **D 676 / U 295**. Downstream dominates at
  70/30, so the oriented reading is the better default — but `--one-sided` would discard **30%** of the
  support, which is too much to drop silently. The generalisation matters too: this is not only a
  sibling question, because `TD_SRC` never carries the source element's strand for *any* class.
- **What margin makes §5's "no regression" a real test** rather than a rubber stamp. Needs the pool's
  three-arm agreement measured first — which is step 5, so the margin is an output of the baseline, not
  an input to it. **Two arms are measured and they agree on 90 of 94** — but they are both xTea2 and one
  is the other plus a flag, so that agreement is a floor on what a real margin can be, not the margin.
- **[ADDED] Whether §3.4's repeat test should differ by family.** It is the only swept parameter that
  moves the answer (§3.1) and for SVA it is decisive rather than influential (§5). A family-specific
  cutoff would be defensible on biology — SVA is young — and indefensible as a knob tuned until the fix
  passes. This needs a reason stated before a value is chosen.
- **[ADDED] Why L1's *regular unknown-start* stratum scores 55.0% against reference-copy's 83.9%.**
  Both are L1, both regular, and the only difference is that the caller could not name the source's
  start. 29 points is too large to leave unexplained: either the window placement degrades when only
  `end` is known, or the `-1` shape marks weaker calls. The verdicts TSV has the 22 entries.

## 9. Where things live

| | |
|---|---|
| **the agenda** | `docs/handoff/DOWNSTREAM_OVERHAUL_PLAN.md` (§4 Stage 3's gate block carries this spec's load-bearing half) |
| **what was built against this spec** | `docs/handoff/HANDOFF_2026-08-04_td_gate_built.md` — read it for the four numbers here that measurement changed |
| the state handoff | `docs/handoff/HANDOFF_2026-08-04_sva_decoy_fix.md` |
| **the pool builder** | `docs/handoff/runs/probe_td_pool.py` — §2.1's table per arm, the N-arm merge, the picked-vs-pool reconciliation, two ledgers |
| **the verifier** | `benchmark/verify_transduction.py` — built, and §7 step 3's validation PASSES at 83.9% |
| **master on HG002** | `docs/handoff/runs/master_wgs_hg002.sbatch`, job `49763841` → `/n/scratch/users/c/cos689/xtea2_testdata/master_wgs_hg002/HG002`. **Not** a repoint of `master_serial_wgs.sbatch`, and §7 step 1 says why |
| master's NA12878 harness | `docs/handoff/runs/master_serial_wgs.sbatch` (job `48825977` PASS, 50:18). **It opens `rm -rf "$OUT"` over the NA12878 tree** — never invoke it without `MASTER_OUT` set |
| master's NA12878 output | `/n/scratch/users/c/cos689/xtea2_testdata/master_serial_out/NA12878/` |
| the truth-set harness | `benchmark/benchmark_truth_set.py` (`TD_SRC` parser at `:206` — now five shapes, plus `td_src_call_class` / `td_src_source_class`, which are §4's two strata as functions) |
| the SVA finding | `docs/handoff/runs/probe_td_source_classes.py`, and plan §10's two entries |
