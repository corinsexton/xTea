# xTea2 handoff — 2026-08-04, the Stage 3 transduction gate is BUILT, and its source-class stratum was reading the wrong field

**This is the current entry point for state.** It is NOT the entry point for agenda — the plan is
(`docs/handoff/DOWNSTREAM_OVERHAUL_PLAN.md`), and for the gate's *design* the entry point is
`HANDOFF_2026-08-04_stage3_td_gate_spec.md`, whose corrections are all inline there marked
`[CORRECTED]` as well as summarised here.

> **CORINNE DECIDED ALL FOUR OPEN QUESTIONS ON 2026-08-05 — see §10.** The gate's pass/fail rule is no
> longer open, which is what the spec's "SPEC, NOT AGREED" banner was waiting on. What remains is
> implementation, and §8 NEXT is rewritten around it.

> Supersedes `HANDOFF_2026-08-04_sva_decoy_fix.md`. Its §0 agenda item **1** — "decide the default"
> for `--td-decoy-structural` — is **still open and is now better informed rather than answered**: the
> gate that was supposed to justify flipping it can adjudicate **2** of the fix's 4 pool entries, both
> orphan-class, and the one entry that matters most is `UNDECIDABLE` for a structural reason. §4. Its
> items 2, 3 and 4 stand untouched.

**Three things, and the third is the one to argue about.**

1. **The gate's two programs exist and the verifier works** — `probe_td_pool.py` builds the pool,
   `verify_transduction.py` adjudicates it from the BAM with no caller in the loop, and it confirms
   **26 of 31** two-sided entries where read support should exist. §2, §3.
2. **`TD_SRC`'s BARE shape is not a source interval — it is two discordant-pair midpoints**, written by
   the sibling module, and reading it as an element interval invented a "reference-copy" stratum that
   does not exist. It relabels **56 of master's 86 L1 assertions** and it is the biggest thing in this
   handoff. **§6.**
3. **The number admitting the decoy-fixed arm to the pool was measured at the wrong level.** 241 →
   **4**, of which 2 are adjudicable and both orphan-class. Not a reason to abandon the fix; a reason
   the *gate* cannot be what justifies it without a decision from Corinne. §4.
4. **`TD_SRC` had five shapes and the truth-set parser handled three**, silently. §5.

---

## 0. State at handoff

| | |
|---|---|
| branch | `xTea2` |
| commits this session | `2d38feb` the parser fix + `probe_td_pool.py` + `master_wgs_hg002.sbatch`; `9d8e1de` `verify_transduction.py`; `3bb4311`/`2260737`/`0303dda` three corrections to my own work (see §4, §6, and §1's intra-arm note); `9361e9b` the frozen baseline. **Nothing under `xtea/` was touched** — the caller is byte-identical to `55055ba`, so BASELINE-0 and every frozen gate result still stand |
| tests | **218 passed / 31 skipped** in the fast subset (was 213/31); five new in `tests/test_benchmark_truth_set.py` — three pinning the parser, two pinning `td_src_shape` and the refusal to call a disc span a reference copy. The skips are the C++ gates with no `XTEA2_SCAN_BIN` |
| in flight | **nothing.** Job `49763841` COMPLETED at **1:10:15** (L1 30:05 / Alu 17:28 / SVA 22:38, MaxRSS 39.9 GB). Out: `/n/scratch/users/c/cos689/xtea2_testdata/master_wgs_hg002/HG002` |
| **the gate's baseline** | **spec §7 step 5 is DONE** — `/n/scratch/users/c/cos689/xtea2_testdata/td_gate_baseline_hg002/` (pool TSV, verdicts TSV, both reports, `SHA256SUMS`). §2.1 |
| **the pass/fail rule** | **AGREED, Corinne 2026-08-05 — §10.** All four open questions decided; nothing about this gate now waits on a judgement call |
| what is NOT done | §8's items 1-3, which are pure implementation of §10's decisions: restratify, build the no-regression gate, then flip `--td-decoy-structural` on by default. **No compute is blocked** |

## 1. What the two programs are, and the one discipline they share

Neither of them parses `TD_SRC` or matches loci with its own code. `probe_td_pool.py` imports
`parse_td_src`, `td_src_call_class`, `td_src_source_class`, `match_one_to_one` and
`DEFAULT_POS_WINDOW` from `benchmark/benchmark_truth_set.py`; `verify_transduction.py` reads the pool
TSV the first one writes. That is the same rule `probe_td_source_classes.py` follows for
`transduction_evidence`, and it exists because **a probe with its own copy of a rule can agree with a
scorer that disagrees with it.** A pool built here and a score computed there cannot disagree about
which loci are the same locus.

- **`docs/handoff/runs/probe_td_pool.py`** — per arm, spec §2.1's table stratified by §4's two axes;
  across arms, the pool merged at ±200 bp **one-to-one** through the harness's own matcher; the
  picked-source reconciliation of §4 below; and a pool TSV. Two ledgers, both **PASS**: master's
  NA12878 (115 assertions in 6 strata) and BASELINE-0 (90 in 6).
  - Matching is per **(family, source)**, not per locus. Two arms calling one locus with *different*
    sources are asserting two independently verifiable things, and folding them would silently pick a
    winner.
  - **It also collapses each arm's OWN duplicates, which the first version did not.** Spec §2 requires
    that "one call cannot satisfy two pooled entries", and `match_one_to_one` enforces that only
    *between* two sets — an arm's entries are never matched against each other. Master's NA12878 arm
    had one such pair: `chr9:87033714` and `chr9:87033895`, **181 bp apart with the identical source**
    `chr3:-1-32916262`, both of which the verifier then refuted separately. One hypothesis, counted and
    adjudicated twice. **So master's pool is 114 entries, not 115** — the 115 is the *call* count, and
    the per-arm table still reports it as such, which is why both ledgers still PASS.
- **`benchmark/verify_transduction.py`** — spec §3. Its five decisions and their citations are in its
  docstring; the two worth repeating here are in §3 below.

## 2. The validation — and the spec's designated validation target does not exist as described

Spec §7 step 3: "the L1 reference-copy sources master reports on NA12878 **have real intervals**, so a
verifier that cannot confirm a healthy fraction of those is wrong before it is applied to anything."
**That premise is false, and §6 below is why.** Those 33 entries are the BARE `TD_SRC` shape, whose
coordinates are discordant-pair midpoints, not a source element's interval. There is no
"reference-copy source" stratum of any size in master's NA12878 output to validate against — it has
**7** entries with a flank-derived source interval, total, across both families.

**The verifier is nonetheless validated in substance, on the stratum that actually exists.** Run over
master's own NA12878 WGS BAM, 114-entry pool, at N=3 / f_cns=0.80 / f_src=0.65 / cov≥5 / rep_div<5:

| stratum | VERIFIED | REFUTED | SUB_THR | UNDEC | precision |
|---|---|---|---|---|---|
| **L1 regular disc_span** | **37** | **14** | 2 | 2 | **72.5%** |
| — of it, two-sided (`start ≥ 0`) | 26 | 5 | 1 | 1 | **83.9%** |
| — of it, one-sided (`start = -1`) | 11 | 9 | 1 | 1 | **55.0%** |
| L1 sibling | 5 | 1 | 0 | 21 | 83.3% |
| SVA orphan polymorphic | 4 | 0 | 0 | 0 | 100.0% |
| SVA sibling | 7 | 10 | 3 | 5 | **41.2%** |
| L1 orphan reference_copy | 0 | 0 | 0 | 3 | *empty denominator* |

Whole pool (114 entries): **65.4% – 79.0%** across 324 sweep points.

**26 of 31 two-sided entries confirmed is the "healthy fraction" the step was asking for**, so the
verifier does what it was built to do — it finds real read support linking X to the coordinates the
caller named. What it cannot do is confirm *reference-copy sources*, because the pool barely contains
any.

**And this resolves §8's open question about the 29-point gap.** It was framed as "unknown-start
versus named-interval", i.e. about whether the caller could name a source start. It is not:
**both are the same bare shape, and `start = -1` means one of the two discordant-pair midpoints was
never established — one-sided evidence.** One-sided verifies at 55.0%, two-sided at 83.9%. That is a
usable finding about the caller: one-sided disc support for a transduction is ~30 points less
reliable. (The reading of `-1` follows from the writer taking `min(src_pos_start, src_pos_end)`; it is
not yet pinned to a line that assigns `-1`, so treat the *mechanism* as strongly indicated and the
*split* as measured.)

**The strata are not one population, which is §4 of the spec measured rather than argued.** SVA's
siblings — **25 of master's 29 SVA assertions** — are the weakest class in the table at 41.2%. A
per-family number would have averaged that away.

## 2.1 The gate's baseline — spec §7 step 5, the three-arm HG002 pool

`master` (job 49763841) + `BASELINE-0` + `--td-decoy-structural` (job 49724056), merged at ±200 bp
one-to-one, **115 entries**. Frozen with sha256s at
`/n/scratch/users/c/cos689/xtea2_testdata/td_gate_baseline_hg002/`.

**Precision 59.1% – 81.0% across 324 sweep points.** Per stratum at the midpoint (N=3, f_cns=0.80,
f_src=0.65, cov≥5, rep_div<5):

| stratum | VERIFIED | REFUTED | SUB_THR | UNDEC | precision |
|---|---|---|---|---|---|
| L1 regular disc_span | 32 | 11 | 3 | 4 | 74.4% |
| SVA sibling | 9 | 7 | 5 | 5 | 56.2% |
| L1 sibling | 2 | 2 | 2 | 20 | 50.0% |
| SVA orphan polymorphic | 3 | 0 | 0 | 0 | 100.0% |
| **SVA orphan reference_copy** | **2** | **0** | 0 | 1 | **100.0%** ← the decoy fix's contribution |
| SVA regular polymorphic | 2 | 0 | 0 | 0 | 100.0% |
| L1 regular polymorphic | 0 | 1 | 0 | 0 | 0.0% (n=1) |
| L1 orphan reference_copy / L1 regular reference_copy / **SVA regular reference_copy** | 0 | 0 | 0 | 1 / 2 / **1** | *empty* |

**Arm agreement, which is what §8's margin question needs**: **82 of 115 entries are asserted by all
three arms**, 8 by two, 25 by one. Master contributes 21 unique, the decoy arm 4, **BASELINE-0 zero**.

**§7 step 2 is answered**: HG002's master pool is **104 transduction-asserting calls** (72 L1 + 32 SVA
+ 0 Alu), against NA12878's 115.

**And the decoy signature is now confirmed on the gate's own sample, third independent time.** Master
on HG002: **SVA 72 picked sources, 0 reference copies**; L1 134 with 79. Same signature as master's
NA12878 (0 of 92) and BASELINE-0 (0 of 105). So §5's rule — never score the fix against master's calls
— rests on a measurement of HG002 itself, not a transfer from another sample.

## 2.2 Timing and memory, master vs xTea2, same sample and allocation

Not what this run was for, but it is the first like-for-like measurement on HG002 and it reverses the
plan's NA12878 figure.

| | xTea2 BASELINE-0 (job 49441918) | master (job 49763841) |
|---|---|---|
| **wall, 24 cores, one job, 3 types** | **45:31** | **1:10:15** |
| **peak RSS** | **~18.7 GB** | **39.9 GB** |
| front end | one `xtea-scan` pass, all 3 types: 24:50 | clip+disc **per type**: 13:38 + 11:24 + ~9:00 |
| downstream | 3 types **concurrent**: 16:56 | serial: L1 16:27, Alu ~8:00, SVA ~13:30 |
| per type | — | L1 30:05, Alu 17:28, SVA 22:38 |

**xTea2 is 1.54× faster on 0.47× the memory** — and the C++ `td_gather` (15.25× L1 / 17.82× SVA, job
49453857) was **off** in BASELINE-0, so that is the figure *before* its largest known win.

Two things not to over-read. The plan records master 50:18 serial against xTea2 1:18:23 on **NA12878**
— master ahead; this is a different sample *and* a much later xTea2 (the scan pass and parallel
downstream landed between), so the reversal is real but not a like-for-like contradiction of that
figure. And master's 24-core scaling is poor on its own terms: L1 went 37:15 at 8 cores (NA12878
golden) to 30:05 at 24, ~20% for 3× the cores.

## 3. Two things the spec did not have, both found by running it

## 4. The 241, and why it does not size what it was cited for

Spec §2 admits the decoy-structural arm to the pool on the strength of **"241 reference-copy SVA
sources"**. That number is right, and it is measured in `unique_trsdct_disc_only`. The same section
defines the pool on `TD_SRC` in the call set. Measured on both:

| | picked (locus,source) pairs | non-sibling `TD_SRC` in the VCF |
|---|---|---|
| L1 (both arms, identical) | 159 | 39 |
| SVA, BASELINE-0 | 105 | 1 |
| SVA, `--td-decoy-structural` | 346 | 5 |

**These two columns are not nested, so the ratio between them is not a filtering rate — and I had it
as one.** `x_transduction.py:238` keys the picked-source file on the insertion positions of the step's
`--input2`, which the golden's chain sets to `candidate_list_from_disc.txt.clip_sites_raw_disc.txt`:
the **raw disc-site universe, 36,687 rows for L1** against a 287-row call list. Measured — **104 of its
111 L1 loci are exact matches in that raw file**, and only **4 of the VCF's 59 non-sibling transduction
calls** land within ±200 bp of any picked row, because raw disc anchors and clip-refined call positions
are different coordinates. So these are differently-scoped populations in different coordinate systems,
and "159 → 39, 4.1×" was the wrong frame even though both counts are right.

**Which makes the conclusion stronger, not weaker.** 241 cannot be the arm's pool contribution under
*any* filtering story, because the two are not the same universe.

So at the pool's own level the fix contributes **4** reference-copy SVA entries against BASELINE-0's
0 — **3 orphan-class, 1 regular.** Named, because four entries is few enough to name:

```
SVA chr17:45593489  orphan   <- chr17:46252374-46253788    VERIFIED    (9 disc)
SVA chr20:45799810  orphan   <- chr8:100891476-100893344   VERIFIED    (15 clip + 3 disc)
SVA chrY:10120435   orphan   <- chr20:30846336-30848045    UNDECIDABLE (SVA_D @ 3.8%)
SVA chr12:39235233  regular  <- chr2:197618052-197619758   UNDECIDABLE (SVA_D @ 2.7%)
```

**Neither number is wrong. They answer different questions**, and this is the same class of error plan
§10 already records twice: reading a set difference as a decision, here reading an internal tally as a
pool contribution. 241 sizes the STEP, which is the right figure for Stage 3's *threshold* work
because that is where the thresholds act. 4 sizes what the CALL SET asserts, which is what a
*precision* gate scores — spec §1 defines the question as "whether the transductions it asserts are
real".

**The consequence for spec §5, which is Corinne's decision and not mine.** §5's rule — score the fix
as GROWTH in `VERIFIED`, never as a precision delta against master — is right and it now has almost
nothing to score:

- the GROWTH claim is supportable at **+2 VERIFIED, both orphan-class** — and §4 of the spec says
  orphan sites are not MEI insertions at all (BASELINE-0's L1 orphans are 0/37 on-truth), so that is
  the weakest class to rest a fix on;
- **§6.2's warning that SVA's regular class has N ≈ 0 survives the fix rather than being retired by
  it**, which is the opposite of what §5 assumed;
- and the blocker is **structural, not tuning**: SVA is a young family, so its reference copies sit at
  low divergence, and §3.4's repeat test disqualifies them *by construction* — where the source window
  and the family consensus are the same sequence, no read separates transduction from element body.
  That is the exact trap §7.1 hit with 11 of 20 L1 residuals.

Three ways out, none of them mine to pick: **(a)** loosen the divergence cutoff for SVA specifically
and state the biological reason before choosing the value, **(b)** accept orphan-class evidence as
sufficient for the growth claim, or **(c)** score the fix at the step's level (the 241) with a
separate, explicitly non-precision argument. **Decide before the pool is frozen** — spec §7 step 5.

## 5. `TD_SRC` has five shapes; the parser handled three

Both `{source}~orphan` and `{source}~{cluster_l}-{cluster_r}` returned `None` from
`benchmark_truth_set.parse_td_src`, which is **indistinguishable from `not_transduction`**. That
dropped 7 of master's 115 NA12878 assertions, 3 of BASELINE-0's and 5 of the decoy arm's, and reached
the scorer as `TP_POS_SOURCE_MISSED` for calls that had named their source perfectly well.

The grammar, from both writers rather than from one: `transduction_evidence.call_candidates` writes
`"{trsdct_source}~{trsdct_l}-{trsdct_r}"` for a regular call and `"{trsdct_source}~orphan"` for an
orphan, and `trsdct_source` is itself `chrom:start-end~rc`. Master's
`x_orphan_transduction.py:1203` consumes the same grammar and says so in its own comment. So the
source is the **first** `~`-field after the colon and everything after it is provenance — fixed by
stripping the tail, not by adding two special cases. The tail is also where the **call class** lives,
which is now `td_src_call_class`.

**What it moves on BASELINE-0: exactly one reported number.** L1's `FP_TRANSDUCTION` 0 → 1, the
previously-invisible `chr6:13190783` whose source is `chr17:66641176-66641176~0~66641176-66641893`.
Every position sensitivity and precision figure is byte-identical (L1 140/57/48, Alu 1103/252/120,
SVA 77/13/5). Verified by running the scorer with the fix and with it stashed.

## 6. `TD_SRC`'s BARE shape is two discordant-pair midpoints, not a source interval

**This is the session's most consequential finding and it invalidates a stratum, including one I had
already reported.** There are two writers putting different things in the same field:

| writer | shape | what the interval IS |
|---|---|---|
| `x_transduction.py:1196` — the transduction step | always **tailed** (`~rc~{l}-{r}` or `~orphan`) | `s_source = "{chrm_fl_L1}:{source_start}-{source_end}~{sourc_rc}"`, from a **flank record** — a real source element |
| `x_orphan_transduction.py:365` — the sibling / orphan module | **bare**, no tail | `s_src = "{src_chrm}:{i_start}-{i_end}"` where `:345-346`'s own comments read **"mid-position of left discordant pair"** and **"mid-position of right discordant pair"** |

Byte-identical in master at the same line numbers. So a bare shape's `end - start` is how far apart two
disc-pair midpoints are, and **`start < end` says nothing about whether the source is a reference
copy.**

Three measurements made this visible before the code was read, any one of which should have been
enough:

- bare spans are 1–2814 bp (median **28** master, **52** BASELINE-0) where a full-length L1 is ~6 kb;
- **none of master's 33 bare intervals lie inside, or even within 5 kb of, any source region** in
  `unique_trsdct_disc_only`;
- in xTea2's *tailed* shape the first interval reaches **6030 bp** — element scale — showing what the
  field looks like when it really does carry a source region.

**What it costs.** `td_src_source_class` now returns a fourth value, `disc_span`, which is **not** a
subtype of the other three. Re-measured:

| arm | entries with a flank-derived source interval | of how many assertions |
|---|---|---|
| master NA12878 | **7** (3 L1 orphan reference-copy + 4 SVA orphan polymorphic) | 115 |
| BASELINE-0 HG002 | **3** (1 L1 reference-copy, 1 L1 polymorphic, 1 SVA polymorphic) | 90 |

So §4's second stratum — the axis the whole SVA decoy story is told in — **is defined for 7 of 115
entries and 3 of 90.** That is a far harder limit than §6.2's "SVA's regular class has N ≈ 0"; it is
"the source-class axis has N ≈ 0 almost everywhere". Spec §2.1's per-family table (`reference-copy
source: L1 33, SVA 0`) is reading disc spans in that column, and its §7 step 3 validation target —
"the 33 L1 reference-copy sources … have real intervals" — does not exist.

**What it does NOT touch, and this matters:** the SVA decoy finding is unaffected. `0 of 105 → 241 of
346 reference copies` is measured by `probe_td_source_classes.py` on `unique_trsdct_disc_only` column
3, which is `x_transduction.py:1196`'s **flank-derived** source id. That measurement is on exactly the
right field, and the defect it found is in the function that writes it.

### 6.1 Two smaller corrections to my own work, from the same pass

1. **"orphan" does not imply "polymorphic".** Master's 4 SVA orphans are polymorphic (`start == end`);
   its 3 L1 orphans are reference copies spanning ~6 kb. Both are the *tailed* shape, so both labels
   are sound. Call class and source class are independent axes.
2. **A polymorphic source IS verifiable, and `probe_td_pool` said otherwise in its first version.**
   What it lacks in the reference is the element **body**; the transduced sequence comes from
   downstream of its 3′ end, which is ordinary genomic sequence — literally the window
   `gnrt_flank_regions_of_polymerphic_insertions` cuts its flank from, which is *why* it writes
   `start == end`. Master's 4 SVA orphan-polymorphic entries verify **4/4**.

## 7. Checked and NEGATIVE — recorded so it is not re-proposed

**Source-interval length does not separate VERIFIED from REFUTED.** The five refuted L1
reference-copy entries have intervals of 1, 16, 20, 66 and 965 bp, which looked like a pattern —
degenerate intervals producing spurious transductions. It is not one: mean length is **214 bp for
REFUTED against 232 for VERIFIED**, and the verified set includes intervals of 2, 5 and 6 bp. Most L1
reference-copy source intervals are tiny (median ~25 bp), so `start < end` versus `start == end` is a
thinner distinction than "a full-length copy versus a point" — worth knowing when reading the source
class, but it does not predict the verdict.

## 8. NEXT — all four decisions are made, so this is implementation

**In this order, because 1 and 2 are preconditions for 3.**

1. **Restratify the verifier on shape × sidedness** (§10.2). Primary table
   `family × call class × shape × sidedness`; source class kept beneath with `EMPTY DENOMINATOR`
   markers. `verify_transduction.py` already has both facts per entry — shape from
   `benchmark_truth_set.td_src_shape(raw)`, sidedness from `src_start == -1` — so this is reporting,
   not new measurement. The pool TSV format does not change.
2. **Implement §5's no-regression gate** (§10.4): a `--baseline <verdicts.tsv>` mode on
   `verify_transduction.py` that compares **absolute VERIFIED counts per stratum** and fails if any
   stratum loses entries. Print `PROVISIONAL` in its own output — the margin is not validated until an
   arm exists that is neither master nor xTea2.
3. **Then flip `--td-decoy-structural` on by default** (§10.1), and **in this order**: measure
   no-regression against the frozen baseline *first*, flip second, re-run the chr1 containment gate and
   the WGS payoff gate third. `docs/handoff/runs/gate_decoy_structural_{chr1,wgs}.sbatch` already
   exist and both PASSED with the flag on — what changes is that the flag becomes the default, so the
   *off* arm becomes the one that needs a reason to exist. Note this is the first change to `xtea/`
   since `55055ba`; BASELINE-0 stays frozen and comparisons against it become before/after.
4. **§11's bare-shape writer fix** — scoped, not built. Read §11 before starting it, particularly the
   trap: the disc midpoints are what the verifier anchors on today and must stay available.
5. **Pin `-1`'s assignment in the bare shape to a line**, closing §2's one-sided reading. The split is
   measured (55.0% vs 83.9%); only the mechanism is inferred.
6. **Whether §3.4's repeat cutoff should differ by family** is now explicitly NOT the route to
   justifying the decoy fix (§10.1). It remains a real question about the verifier — it is the only
   swept parameter that moves the answer, 66.2% → 78.4% across 2–10% — but it must be argued on
   biology, in isolation, and never in the same change as a fix whose pass it would affect.
7. The previous handoff's items 2, 3 and 4, unchanged: plan §5's semantic half on gate **[A]**, the
   decoy-consensus concatenation, and `ileftmost_cns` 500-vs-5000.

## 10. The four decisions — Corinne, 2026-08-05

Recorded here because each one closes something the spec left open, and because the *reason* matters
more than the choice: three of the four are choices to avoid tuning a measurement until it agrees.

### 10.1 The decoy fix is justified on structural correctness, and the gate only has to show no regression

**Decision:** do not rest the flip on the gate's ability to score the recovered sources. Justify it on
(a) the defect itself — a substring test that misclassifies **100% of hg38's 5,454 SVA flank records**
as decoy consensus is wrong regardless of how many calls it moves — and (b) the step-level recovery,
**0 → 241 reference copies**, which is measured on `unique_trsdct_disc_only` column 3, the
flank-derived field, i.e. the right field. The gate's role shrinks to **confirming no precision
regression**, which it can do.

**Why this and not the alternatives.** Accepting the 2 orphan-class VERIFIED entries as the GROWTH
claim would rest a shipping default on the class plan §1 says is not an MEI insertion at all
(BASELINE-0's L1 orphans are 0/37 on-truth). Loosening §3.4's divergence cutoff for SVA is defensible
on biology — SVA is young — but it is *the only swept parameter that moves the answer*, so moving it
here would read as a knob tuned until the fix passed. **That is the failure mode this project has
already documented twice** (plan §10's two "read a set difference as a decision" entries).

**What it obliges:** measure no-regression *before* flipping, then flip, then re-run both gate scales.
Not "flip and check".

### 10.2 Shape × sidedness becomes the primary stratum; source class stays with empty cells

**Decision:** report `family × call class × shape × sidedness` as the main table, and keep
`family × call class × source class` beneath it with explicit `EMPTY DENOMINATOR` markers.

**Why:** the source-class axis is defined for **7 of 115** pool entries (§6), while shape × sidedness
is fully populated *and* discriminating — two-sided 74.4%, one-sided 55.0%, sibling 50–56%. Dropping
source class entirely would delete the axis the SVA decoy story is told in and that a future arm could
populate; keeping it silent would let a reader trust an empty cell. Keeping both costs nothing.

### 10.3 The bare-shape writer is documented now and its fix is scoped, not built

**Decision:** consumers are already shape-aware (`td_src_shape`, `disc_span`); plan §10 records the
finding; the writer stays as it is for now.

**Why:** making `x_orphan_transduction.py:365` emit a flank-derived source region would give the
source-class axis real content — and it is a Stage 3-sized change to a shipping path that will move
calls. It gets the same treatment the decoy fix got: **a flag, two scales, its own arm.** See §11.

### 10.4 The margin is absolute per-stratum VERIFIED counts, and it is provisional

**Decision:** gate on **absolute VERIFIED counts per stratum** — no stratum loses entries — rather
than a percentage tolerance. Mark the margin **provisional** until an arm exists that is not master
and not xTea2.

**Why:** the strata that matter most have 1–4 entries (SVA regular reference_copy has **1**), where a
percentage is noise and a single reclassification swings it by 100 points. And three-arm agreement —
82 of 115 — cannot set a real tolerance, because two of the three arms differ by **one flag**; it is a
floor on independence, not a measure of it. Recorded as provisional so nobody later quotes it as
validated.

## 11. Scoped, not built — the bare-shape writer fix

Per §10.3. What it would be, so the next session does not re-derive the shape of it:

- **The change:** `x_orphan_transduction.py:365` writes `"{src_chrm}:{i_start}-{i_end}"` from two
  discordant-pair midpoints. It would instead resolve the flank record those disc reads hit and emit
  `"{chrom}:{source_start}-{source_end}~{rc}"` — the same form `x_transduction.py:1196` produces — so
  `TD_SRC` means one thing across both writers.
- **Why it is worth doing:** it would move the source-class axis from 7 of 115 entries to most of them,
  which is the difference between §4's second stratum being a real axis and being a footnote.
- **Why it is not done here:** it changes a shipping field that postF, the gVCF and every consumer
  read. Behind a flag, off by default, with a chr1 containment gate and a WGS payoff gate — exactly
  `55055ba`'s pattern.
- **The trap to avoid:** the midpoints are not useless. They are what the verifier currently anchors
  its window on, and that anchoring verifies at 83.9% two-sided. Any change must keep them available
  (a second field, or the existing `~{l}-{r}` tail) rather than replacing them.

## 12. Where things live

| | |
|---|---|
| **the agenda** | `docs/handoff/DOWNSTREAM_OVERHAUL_PLAN.md` |
| **the gate's design** | `docs/handoff/HANDOFF_2026-08-04_stage3_td_gate_spec.md` — corrections inline, marked `[CORRECTED]` |
| the previous state handoff | `docs/handoff/HANDOFF_2026-08-04_sva_decoy_fix.md` |
| the pool builder | `docs/handoff/runs/probe_td_pool.py` |
| the verifier | `benchmark/verify_transduction.py` |
| master on HG002 | `docs/handoff/runs/master_wgs_hg002.sbatch`, job `49763841` PASS 1:10:15 → `/n/scratch/users/c/cos689/xtea2_testdata/master_wgs_hg002/HG002` |
| **the gate's frozen baseline** | `/n/scratch/users/c/cos689/xtea2_testdata/td_gate_baseline_hg002/` — pool TSV, verdicts TSV, both verify reports, `SHA256SUMS` |
| **a trap in the neighbouring harness** | `docs/handoff/runs/master_serial_wgs.sbatch` opens `rm -rf "$OUT"` and its default `OUT` is the NA12878 tree the gate spec cites. Never invoke it without `MASTER_OUT`. `master_wgs_hg002.sbatch` has a refusal guard for exactly this |
| the truth-set harness | `benchmark/benchmark_truth_set.py` — `parse_td_src` (five shapes), `td_src_call_class`, `td_src_source_class` |
