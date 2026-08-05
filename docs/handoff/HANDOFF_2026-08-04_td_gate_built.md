# xTea2 handoff — 2026-08-04, the Stage 3 transduction gate is BUILT, and four of its spec's numbers moved

**This is the current entry point for state.** It is NOT the entry point for agenda — the plan is
(`docs/handoff/DOWNSTREAM_OVERHAUL_PLAN.md`), and for the gate's *design* the entry point is
`HANDOFF_2026-08-04_stage3_td_gate_spec.md`, whose corrections are all inline there marked
`[CORRECTED]` as well as summarised here.

> Supersedes `HANDOFF_2026-08-04_sva_decoy_fix.md`. Its §0 agenda item **1** — "decide the default"
> for `--td-decoy-structural` — is **still open and is now better informed rather than answered**: the
> gate that was supposed to justify flipping it can adjudicate **2** of the fix's 4 pool entries, both
> orphan-class, and the one entry that matters most is `UNDECIDABLE` for a structural reason. §4. Its
> items 2, 3 and 4 stand untouched.

**Three things, and the third is the one to argue about.**

1. **The gate's two programs exist and the spec's own validation passes.** `probe_td_pool.py` builds
   the pool; `verify_transduction.py` adjudicates it from the BAM with no caller in the loop. On the
   stratum spec §7 step 3 names as checkable — master's L1 reference-copy sources on NA12878 — it
   returns **26 VERIFIED / 5 REFUTED / 1 SUB_THRESHOLD / 1 UNDECIDABLE, 83.9%**. §2, §3.
2. **`TD_SRC` had five shapes and the truth-set parser handled three**, silently. §5.
3. **The number admitting the decoy-fixed arm to the pool was measured one level too deep, and at the
   pool's own level the fix has almost nothing for this gate to score.** 241 → **4**, of which 2 are
   adjudicable and both are orphan-class. Not a reason to abandon the fix; a reason the *gate* cannot
   be the thing that justifies it without a decision from Corinne. §4.

---

## 0. State at handoff

| | |
|---|---|
| branch | `xTea2` |
| commits this session | `2d38feb` the parser fix + `probe_td_pool.py` + `master_wgs_hg002.sbatch`; `9d8e1de` `verify_transduction.py`. **Nothing under `xtea/` was touched** — the caller is byte-identical to `55055ba`, so BASELINE-0 and every frozen gate result still stand |
| tests | **216 passed / 31 skipped** in the fast subset (was 213/31); three new in `tests/test_benchmark_truth_set.py`, all pinning the parser. The skips are the C++ gates with no `XTEA2_SCAN_BIN` |
| in flight | **job `49763841`** — master's own pipeline on HG002 WGS, spec §7 step 1. L1 finished in 30:05; Alu and SVA to go. Out: `/n/scratch/users/c/cos689/xtea2_testdata/master_wgs_hg002/HG002` |
| what is NOT done | spec §7 step 5 — freeze the three-arm pool and score it. Two arms are measured; the third is the job above |

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

## 2. The validation, which is the reason to believe anything else

Spec §7 step 3: "the L1 reference-copy sources master reports on NA12878 have real intervals, so a
verifier that cannot confirm a healthy fraction of those is wrong before it is applied to anything."
Run over master's own NA12878 WGS BAM, 114-entry pool, at N=3 / f_cns=0.80 / f_src=0.65 / cov≥5 /
rep_div<5:

| stratum | VERIFIED | REFUTED | SUB_THR | UNDEC | precision |
|---|---|---|---|---|---|
| **L1 regular reference_copy** | **26** | **5** | 1 | 1 | **83.9%** ← the checkable class |
| L1 regular unknown_start | 11 | 9 | 1 | 1 | **55.0%** |
| L1 sibling | 5 | 1 | 0 | 21 | 83.3% |
| L1 orphan reference_copy | 0 | 0 | 0 | 3 | *empty denominator* |
| SVA orphan polymorphic | 4 | 0 | 0 | 0 | 100.0% |
| SVA sibling | 7 | 10 | 3 | 5 | **41.2%** |

Whole pool (114 entries): **65.4% – 79.0%** across 324 sweep points.

**The strata are not one population, which is §4 of the spec measured rather than argued.** An unknown
start costs **29 points** against a named interval — same family, same call class, the only difference
being that the caller could not name the source's start. And SVA's siblings, which are **25 of
master's 29 SVA assertions**, are the weakest class in the table. A per-family number would have
averaged all of that away.

## 3. Two things the spec did not have, both found by running it

### 3.1 The parameter that decides the answer is not in §3.1's list

§3.1 names N, the MAPQ floor, the 80% consensus fraction and the coverage floor. Over the NA12878
pool, **each of those moves precision by ≤ 2 points across its whole range.** The one that moves it is
§3.4's **repeat-divergence cutoff**: 67.1% → 77.4%, and `UNDECIDABLE` from 24 to 52. So the verifier
sweeps it too, and it is the one to argue about.

### 3.2 There are four states, not three

§3.5 names `VERIFIED (≥N)` and `REFUTED (0)` and leaves `0 < support < N` unnamed. `SUB_THRESHOLD` is
its own state: folding it into `REFUTED` lets the choice of N manufacture false positives, folding it
into `VERIFIED` lets N do nothing. Neither it nor `UNDECIDABLE` is ever in the precision denominator —
the same rule the decision ledger enforces between `no_ledger_row` and "unattributed".

### 3.3 §8's orientation question, answered partly and from data

`TD_SRC` never carries the source element's strand — for **any** class, not just siblings — so
"downstream" is unoriented. The verifier searches `[end, end+W]` and `[start-W, start]` and reports
which side supported: **D 676 / U 295**. Downstream dominates at 70/30, so the oriented reading is the
better default, but `--one-sided` would discard 30% of the support. Too much to drop silently.

## 4. The 241, and why it does not size what it was cited for

Spec §2 admits the decoy-structural arm to the pool on the strength of **"241 reference-copy SVA
sources"**. That number is right, and it is measured in `unique_trsdct_disc_only` — the transduction
step's own hypotheses, **before** the downstream filtering that produces the call set. The same
section defines the pool on `TD_SRC` in the call set. Measured on both:

| | picked (locus,source) pairs | non-sibling `TD_SRC` in the VCF | attrition |
|---|---|---|---|
| L1 (both arms, identical) | 159 | 39 | 4.1× |
| SVA, BASELINE-0 | 105 | 1 | 105× |
| SVA, `--td-decoy-structural` | 346 | 5 | 69× |

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

## 6. Two corrections to my own work, both caught before they shipped anywhere

1. **"orphan" does not imply "polymorphic".** Master's 4 SVA orphans are polymorphic (`start == end`);
   its 3 L1 orphans are reference copies spanning ~6 kb. Call class and source class are **independent
   axes**, which is why they are two columns. Spec §2.1's L1 reference-copy count is therefore **36**,
   not 33 — the three tailed orphans were invisible.
2. **A polymorphic source IS verifiable, and `probe_td_pool` said otherwise in its first version.**
   What a polymorphic source lacks in the reference is the element **body**; the transduced sequence
   comes from downstream of its 3′ end, which is ordinary genomic sequence — literally the window
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

## 8. NEXT

1. **Job `49763841` lands → freeze the three-arm pool and score it.** That is spec §7 step 5 and the
   gate's baseline. `probe_td_pool.py --arm master=… --arm BASELINE-0=… --arm decoy=… --pool-tsv …`
   then `verify_transduction.py`. Answers spec §8's HG002 pool size for the last arm.
2. **§4's decision — how the decoy fix gets scored.** Blocks the `--td-decoy-structural` default,
   which is the previous handoff's item 1. Needs Corinne.
3. **Why L1's unknown-start stratum scores 55.0% against reference-copy's 83.9%.** 29 points between
   two strata that differ only in whether the caller could name the source's start is too large to
   leave unexplained. The verdicts TSV has the 22 entries.
4. **Whether §3.4's repeat cutoff should differ by family** — the only swept parameter that moves the
   answer, and decisive for SVA. Defensible on biology; indefensible as a knob tuned until the fix
   passes. Needs a stated reason before a value.
5. The previous handoff's items 2, 3 and 4, unchanged: plan §5's semantic half on gate **[A]**, the
   decoy-consensus concatenation, and `ileftmost_cns` 500-vs-5000.

## 9. Where things live

| | |
|---|---|
| **the agenda** | `docs/handoff/DOWNSTREAM_OVERHAUL_PLAN.md` |
| **the gate's design** | `docs/handoff/HANDOFF_2026-08-04_stage3_td_gate_spec.md` — corrections inline, marked `[CORRECTED]` |
| the previous state handoff | `docs/handoff/HANDOFF_2026-08-04_sva_decoy_fix.md` |
| the pool builder | `docs/handoff/runs/probe_td_pool.py` |
| the verifier | `benchmark/verify_transduction.py` |
| master on HG002 | `docs/handoff/runs/master_wgs_hg002.sbatch`, job `49763841` |
| **a trap in the neighbouring harness** | `docs/handoff/runs/master_serial_wgs.sbatch` opens `rm -rf "$OUT"` and its default `OUT` is the NA12878 tree the gate spec cites. Never invoke it without `MASTER_OUT`. `master_wgs_hg002.sbatch` has a refusal guard for exactly this |
| the truth-set harness | `benchmark/benchmark_truth_set.py` — `parse_td_src` (five shapes), `td_src_call_class`, `td_src_source_class` |
