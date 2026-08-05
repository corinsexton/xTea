# xTea2 handoff — 2026-08-04, Stage 3's transduction gate: the spec

**This is the entry point for ONE decision: how Stage 3's transduction work gets gated.** It is not a
state handoff — for state read `HANDOFF_2026-08-04_sva_decoy_fix.md`, and for agenda read
`docs/handoff/DOWNSTREAM_OVERHAUL_PLAN.md`.

> **STATUS: SPEC, NOT AGREED.** Corinne directed (2026-08-04) that Stage 3's transduction gate should
> assume **truth from a run of master xTea on HG002**, which is the input the plan was missing — until
> then §7.1's conclusion was "the ALT route is exhausted; what remains is the paper's supplement or a
> different sample." This document is the design that direction implies. **Nothing here is built**, and
> the pass/fail thresholds in §6 are the part most likely to need Corinne's judgement rather than mine.

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
| **xTea2 with the SVA decoy test made structural** | **without it the pool contains ZERO reference-copy SVA sources**, so the gate cannot score the fix that matters most — see §5. **THIS ARM NOW EXISTS AND IS MEASURED**: `run_xtea --td-decoy-structural` (`55055ba`, off by default), and on HG002 WGS it contributes **241 reference-copy SVA sources** against **0 of 105 in BASELINE-0 (same sample)** and **0 of 92 in master's own NA12878 output (different sample, same signature)** — jobs `49724055`/`49724056`, both PASS. Its run tree is on disk, so this arm needs no new compute |

**Merge loci at ±200 bp, one-to-one**, the same window and matching discipline `benchmark_truth_set.py`
uses (`DEFAULT_POS_WINDOW`), so one call cannot satisfy two pooled entries. Plan §7.3 records that
greedy-vs-nearest and one-to-one-vs-not "moves sensitivity/FDR by percentage points"; the harness's
choice is the one to keep.

### 2.1 `TD_SRC`'s shapes, measured — not guessed

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

### 3.1 The verifier's own thresholds are a hazard — report a sweep, not a point

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
- **BASELINE-0 delta reported for attribution, not pass/fail** — the plan's existing Stage 3 rule,
  unchanged.

## 6. What this gate does NOT claim — put this in its docstring

1. **Not sensitivity.** The pool is bounded by what some caller found; it cannot say what all three
   missed, and therefore cannot state how many transductions HG002 has.
2. **SVA gets almost no *regular* transduction truth.** Master's SVA is 25 sibling + 4 orphan + **0**
   reference-copy, so SVA's regular-class precision has N≈0 until the decoy-fixed arm contributes. A
   per-family table that does not say this invites the reader to trust an empty cell.
3. **Alu contributes nothing**, by construction — no transduction flank, so the step never runs. Biology,
   not a gap.
4. **NOT master agreement.** Master enumerates; reads decide. This is what keeps the design inside plan
   §2's locked decision — *"agreement with master is explicitly not a requirement"* — which scoring
   against master's calls directly would contradict.
5. **Not comparable to the paper's numbers** — hs38DH vs `no_alt` (plan §7.1).

## 7. Order of operations, and the prerequisite that will bite if ignored

1. **Repoint and run master on HG002.** `docs/handoff/runs/master_serial_wgs.sbatch` already works —
   change `BAM` and `SAMPLE` from NA12878 to HG002. **Master's OWN pipeline end to end**, because §10
   records that "master's transduction figures come in two non-interchangeable flavours": master's own
   pipeline on its own smaller inputs (36,687 L1 raw sites) versus master's *code* on the port's inputs.
   Only the former samples master's call set. ~50 min (job `48825977` was 50:18).
2. **Measure HG002's pool size.** The 115 above is **NA12878**; HG002's is unknown.
3. **Build the verifier**, and validate it where the answer is checkable: the 33 L1 reference-copy
   sources master reports on NA12878 have real intervals, so a verifier that cannot confirm a healthy
   fraction of those is wrong before it is applied to anything.
4. ~~**The SVA decoy fix must exist as a branch BEFORE the pool is frozen.**~~ **DONE 2026-08-04** —
   it is a flag, not a branch (`--td-decoy-structural`, `55055ba`), and its WGS arm has already run:
   `/n/scratch/users/c/cos689/xtea2_testdata/gate_decoy_structural_wgs_49724056/HG002.wgs`, with 241
   reference-copy SVA sources in `TD_SRC`. So this prerequisite costs zero compute now; do not re-run
   it. See `HANDOFF_2026-08-04_sva_decoy_fix.md` §5.2.
5. **Freeze the pool, score it — that is the gate's baseline.** Then Stage 3's threshold work is
   gateable.

## 8. Open questions I could not settle from the code

- **HG002's pool size** — only NA12878's 115 is measured. **Partially answered for the third arm**:
  its HG002 SVA contribution is 346 picked sources of which 241 are reference copies, and 87
  orphan-labelled loci (job `49724056`).
- **Whether a sibling's `POS` is the source element's 3′ end or its insertion point.** §2.1's `±100`
  convention is what xTea's own consumer does, which makes it the right *default*, but it is a window
  around a point rather than an oriented 3′ end — so for siblings the verifier searches both directions
  or accepts a looser claim. Worth resolving before N is calibrated, because 25 of master's 29 SVA
  assertions are siblings.
- **What margin makes §5's "no regression" a real test** rather than a rubber stamp. Needs the pool's
  three-arm agreement measured first — which is step 5, so the margin is an output of the baseline, not
  an input to it.

## 9. Where things live

| | |
|---|---|
| **the agenda** | `docs/handoff/DOWNSTREAM_OVERHAUL_PLAN.md` (§4 Stage 3's gate block carries this spec's load-bearing half) |
| this session's state | `docs/handoff/HANDOFF_2026-08-04_cns_cluster_ledger.md` |
| master's proven WGS harness | `docs/handoff/runs/master_serial_wgs.sbatch` (NA12878 today; job `48825977` PASS, 50:18) |
| master's existing output | `/n/scratch/users/c/cos689/xtea2_testdata/master_serial_out/NA12878/` — **no HG002 master run exists on disk** |
| the truth-set harness | `benchmark/benchmark_truth_set.py` (`TD_SRC` parser at `:207`, `DEFAULT_POS_WINDOW`) |
| the SVA finding | `docs/handoff/runs/probe_td_source_classes.py`, and plan §10's two entries |
| the verifier | `benchmark/verify_transduction.py` — **to build** |
