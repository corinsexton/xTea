# xTea2 handoff — 2026-08-05 (second session), gate [A] is measured and three claims are corrected

**This is the current entry point for state.** The agenda entry point is still
`docs/handoff/DOWNSTREAM_OVERHAUL_PLAN.md`; the gate's design entry point is still
`HANDOFF_2026-08-04_stage3_td_gate_spec.md`.

> Supersedes `HANDOFF_2026-08-05_decoy_default_flipped.md`. **Its §6 NEXT items 1 (first half), 2,
> 6 and 7 are DONE.** Item 1's second half — the decoy-consensus concatenation — is what remains and
> is now the top item. Items 3, 4 and 5 are unchanged: §11's writer fix stays scoped per §10.3, the
> repeat cutoff stays explicitly-not-the-route per §10.1, and the third arm is still blocked on
> nothing but itself.
>
> **§10's four decisions (2026-08-04) remain the standing agreement. Nothing here reopens one.**

**Four things, and three of them correct a claim this project has been repeating.**

1. **Gate [A] is measured — the gate that does 97% of the rejecting — and it needed 57 seconds, not
   a pipeline run.** Its polyA *threshold* is unreachable, so the decision is a membership test. §2.
2. **`NEARBY_CLIP`'s "measured inert, do not open an instrumented run" is a claim about gate [B2]
   only.** At gate [A] the same constant is worth ~4× the clip evidence. §3.
3. **The two "orphan"s are disjoint BY CONSTRUCTION, not by coincidence** — pinned to one line of
   each of two writers, and 43/43 in master, so the collision is inherited. §4.
4. **The flip's off-truth cost is adjudicated, and its cleanest win was never counted.** Of the
   three gained calls it is 1 VERIFIED and 2 UNDECIDABLE, not two VERIFIED — and a fourth,
   already-called site gained a VERIFIED transduction assertion at zero precision cost. §5.

---

## 0. State at handoff

| | |
|---|---|
| branch | `xTea2` |
| commits this session | `265e71d` the plan's 2026-08-05 block; `4b28d53` the gate [A] probe + harness; `18aed56` the two-orphan disambiguation; `3f20bf3` the plan's gate [A] + `NEARBY_CLIP` scope correction; `b2a7f41` the off-truth adjudication; `ec52a81` this handoff; `e8bfc1c` the truth-set facts (§5.3) |
| files touched | **8, in `benchmark/`, `docs/` and `tests/` only** — verified `git diff --name-only 61ffad1..HEAD -- xtea/` is empty |
| tests | **234 passed / 30 skipped** (`python3 -m pytest -m "not slow"`, system 3.9.25). One new: `test_the_two_orphans_are_disjoint_by_construction` |
| **no `xtea/` change this session** | every commit is in `benchmark/`, `docs/` or `tests/`. The last `xtea/` change is still `53ca400`, so **BASELINE-0 and every frozen harness baseline are untouched** and no gate needs re-running |
| in flight | **nothing.** Job `49836568` reported PASS in 57 s |
| jobs | `49836568` gate [A] at WGS, **PASS**, **57 s**, `-c 2 -p short -A park`. Artifacts `/n/scratch/users/c/cos689/xtea2_testdata/probe_clip_gate_wgs_49836568/` |
| what is NOT done | the **decoy-consensus concatenation** (§6 item 1) — the last piece of the previous handoff's item 1, and now the top of the list |

## 1. What made this session cheap, and it is a rule rather than luck

Both of this session's measurements were filed as needing runs. Neither did.

- **Gate [A]** was filed against four handoffs as the "semantic half", implicitly a pipeline change.
  Its entire input is `tmp/transduction/temp_transduction_clip.sam`, which **survives a finished
  run** — and `parse_transduction_clip_algnmt` builds *both* of gate [A]'s inputs from it, the clip
  maps and the polyA map. 57 seconds on 2 cores.
- **The off-truth question** (§5) needed two VCFs, the truth set, and the verdicts file already on
  disk.

This is the **second** time the step's intermediates have outlived a handoff's belief that they were
gone — `probe_nearby_clip` found the disc-to-consensus SAM the same way. It is now a §10 entry:
**check the intermediates before scoping a run.**

## 2. Gate [A], measured — and its named polyA threshold decides nothing

`docs/handoff/runs/probe_clip_gate.py` drives the shipping parse (`x_transduction.py:281`) and the
shipping gate (`transduction_evidence.select_clip_sites`, `:283`). Nothing is reimplemented, for the
reason plan §10 records in numbers: `probe_leftmost_cns`' first version carried a private copy of a
shipping predicate and read 2.6× wrong.

**Check 0 is two checks, and the second is stronger than anything the earlier probes used.** At the
shipping operating point the probe must reproduce (a) the run's ledger `A_clip` reason split and
(b) the run's own `.unique_trsdct_half_clip` output file **as a set of loci**. Both, both families:
L1 `14276`/`6056` and 2,802 loci; SVA `16286`/`11304` and 3,233. So the sweep is a measurement, and
`i_flank_length` is *verified* irrelevant to gate [A] rather than assumed to be.

### 2.1 `low_polya_clip` is unreachable, so this is a membership test

`m_polyA`'s record is created and one of its six counters incremented **in the same block**
(`x_transduction.py:929-975`), so every entry that exists has `sum(rcd[:6]) >= 1`. Measured: **0
zero-sum entries** of 4,604 (L1) and 5,261 (SVA). At HG002's coverage-derived `n_clip_cutoff = 3`
the derived `n_polya_cutoff` is **0.5**, so the comparison needs a sum of zero and the reason cannot
fire. It first appears at cutoff **5**.

Plan §5 keeps the polyA signature on the "load-bearing" list. **It is load-bearing — and it is a
boolean.** The named parameter is dead at this operating point, which is a different thing from the
signature being unimportant, and the two were indistinguishable in the ledger until now.

### 2.2 Membership, not clip depth, is what binds the gate

Of loci with `n_clip >= 3`:

| | dropped for polyA-map ABSENCE | kept | of the dropped, `n_clip >= 10` |
|---|---|---|---|
| L1 | **6,056** (68%) | 2,802 | 594 (9.8%), max **969** |
| SVA | **11,304** (78%) | 3,233 | 810 (7.2%), max 707 |

The kept population's `n_clip` is only modestly higher (median 5 vs 4, p90 16 vs 9), so **this is
not clip depth under another name** — a locus with 969 qualified clip realignments is dropped for the
absence of a single polyA-ish read. Whether that is right is a biology question, and it is the one
gate [A]'s semantic half actually contains. `R_NO_POLYA_CLIP` writes **no detail**, which is why the
ledger cannot say any of this.

### 2.3 Two thresholds, one dial

`n_polya_cutoff = n_clip_cutoff/2 - 1` is derived inside `TransductionParams`, so a clip-support
sweep is single-variable only up to 4; at 5 the polyA threshold wakes and rejects 286 (L1) / 364
(SVA) loci the shipping configuration keeps. **"Raise the clip cutoff" and "raise the polyA cutoff"
are not separable changes today** — which is §5's "every threshold a named parameter" ambition
failing at the gate that matters most. Separating them is a one-line change to a constructor and a
Stage 3 decision, not a cleanup: it would move calls.

### 2.4 The decoy asymmetry reproduces from the clip side

The clip path does **not** skip decoy hits — the skip is commented out at `x_transduction.py:872-875`
— so hits to the `>LINE1`/`>SVA`/`>ALU` consensus records count as clip *support*:

| | decoy-consensus clip records | survivors with any | with **only** that |
|---|---|---|---|
| L1 | **3,992** of 99,245 (4.02%) | 529 of 2,802 | **0** |
| SVA | **0** of 128,262 | 0 | 0 |

SVA's zero is the packaging half of the same defect: hg38's SVA flank library has no decoy records
concatenated in. This is the **same family asymmetry** the disc side's `decoy_cross_family` shows
(144 L1 / 0 SVA), reached by an independent route — so it is corroboration for §6 item 1 rather than
a new argument for it. **No survivor rests on decoy support alone**, so this channel is not
manufacturing loci; it is family-asymmetric evidence.

## 3. `NEARBY_CLIP` — the scope correction, and it is large

Plan §10's entry ended **"So this item is closed. Do not open an instrumented run for it."** That
sentence is struck. Every measurement behind it was taken **inside gate [B2]**, where the constant
moves one locus genome-wide. The same constant has a **second** use the entry never named:
`x_transduction.py:865` drops a clip realignment whose `|ref_mpos - ori_insertion_pos|` exceeds it,
**ahead of** the unmapped and clip-qualified tests — which is gate [A]'s domain.

| at gate [A], WGS | L1 @ 50 | L1 unbounded | SVA @ 50 | SVA unbounded |
|---|---|---|---|---|
| qualified clip records | 99,245 | **415,429** (4.2×) | 128,262 | **519,097** (4.0×) |
| domain (loci) | 23,134 | 31,789 | 30,823 | 38,212 |
| survivors | 2,802 | **10,339** (3.7×) | 3,233 | **12,548** (3.9×) |

**This is not a recommendation to raise it, and two reasons are inside the measurement.** Part of
the growth is the same reads re-assigned to nearby loci rather than new evidence — at unbounded
distance the *kept* population's median `n_clip` goes 5 → 9 while the locus count rises only 37%, so
existing loci are gaining support as well as new loci appearing. And gate [A]'s survivors still face
[B], [B2], [C], [D] and [E], where 306 of BASELINE-0's 322 false negatives already fail to appear at
all. What this says is that **"measured inert" is a claim about [B2], and gate [A] is unmeasured** —
including the bimodal distance distribution that justified 50, which was taken over the
*disc-to-consensus* records and is not the population this second use consults.

## 4. The two "orphan"s: a mechanism, not a measurement

The previous handoff measured that `SUBTYPE == orphan_or_sibling_transduction` and
`td_src_call_class == "orphan"` agree on no call in job `49792012`. They are **disjoint by
construction**, and it is one line of each of two writers:

- `SUBTYPE == orphan_or_sibling_transduction` is set by exactly **one** writer,
  `x_orphan_transduction.py:378`, on a site the caller **CREATED** — and that writer emits a **bare**
  source with **no tail** (`:371`).
- `td_src_call_class` reads the **tail**, and "no tail" maps to `regular`.

So every created orphan-or-sibling SITE is `regular` by `TD_SRC` whatever its biology, and every
`~orphan`-tailed source came from the other writer (`transduction_evidence.py:847`, on a gate [A]
survivor whose disc anchors bracket the insertion point — `select_orphans`, explicitly not a gate).
Confirmed **37 of 37** on xTea2's WGS decoy arm and **43 of 43** on master's HG002 (`49763841`), so
the collision is **inherited, not xTea2's**.

**The words, now different in code as well as in prose:** a **"created orphan-or-sibling SITE"**
(`benchmark_truth_set.SUBTYPE_ORPHAN_SITE`, which exists so a bare string comparison cannot
reintroduce the confusion) versus an **"orphan-TAILED source"** (`CALL_CLASS_ORPHAN`).

**Three places were citing the wrong population's evidence, all fixed**: `td_src_call_class`'s
docstring, its test's docstring, and — the worst, because it is user-facing —
`verify_transduction`'s **PRIMARY STRATUM header**, which printed the SUBTYPE 0/37 above a table
strata'd by the tail.

**No scored number moved, and that is verified rather than asserted**: re-scored L1 on `49792012`'s
own VCF — TP 140, FN 57, FP 48, 74.5%, 92.7% over 151 calls, all identical. The carve-out stays
keyed on `SUBTYPE` per §7.2 of the previous handoff.

> **The new test pins the MECHANISM, so §11's bare-shape writer fix will FAIL it** the moment that
> writer gains a tail — at which point the two populations stop being disjoint and every claim keyed
> on either field has to be re-read. Failing there is the intended alarm, not a broken test.

## 5. The flip's off-truth cost, adjudicated — and one VERIFIED entry is not what it was thought to be

`docs/handoff/runs/probe_sva_flip_offtruth.py`. Three explanations were checkable and three are dead:

- **Not a windowing near-miss.** The scorer matches at ±200 bp; the nearest truth record of **any**
  family is **250 kb** (chr12), **1.7 Mb** (chr20), **2.5 Mb** (chrY). These are positions the truth
  set says nothing about, not shifted versions of records it has.
- **`chrY:10120435` cannot be on-truth at all.** The SVA truth set has **0 chrY records** of 90 (L1
  also 0 of 197) — and this is **not** a contig exclusion, because the Alu truth set carries 2 chrY
  records. So it is a **denominator artifact**, not a caller error. Whether SVA's chrY zero is a real
  absence or SVA-specific incompleteness is not decidable from the truth set; either way no SVA call
  there can match.
- **Not the source copy's own reads.** All three are **inter-chromosomal** (chr12←chr2, chr20←chr8,
  chrY←chr20), which is the FP mode this kind of call usually has.

**AND THE READ VERDICTS ARE 1 VERIFIED / 2 UNDECIDABLE, NOT "two, both VERIFIED".** The previous
handoff's §6 item 7 had it as two. Of the three *gained calls*:

| call | reads | why |
|---|---|---|
| `chr20:45799810` | **VERIFIED** | support 18 (15 clip + 3 disc), source divergence 9.0% |
| `chr12:39235233` | UNDECIDABLE | `repeat_window(SVA_D@2.7%)` |
| `chrY:10120435` | UNDECIDABLE | `repeat_window(SVA_D@3.8%)` |

The pool's second VERIFIED entry is `chr17:45593489`, at a position **BASELINE-0 already called** —
so it was never part of the precision cost. Both UNDECIDABLEs are disqualified by §3.4's 5% repeat
cutoff, exactly as §4 predicted for young SVA. All three are also one-sided in clip (`RCLIP=0`), the
stratum §5 measures ~29 points worse.

### 5.1 The uncounted win, and the 4-vs-3 discrepancy it explains

`chr17:45593489` went `TD_SRC=not_transduction` → `chr17:46252374-46253788~1~orphan`. **An existing
call gained a transduction assertion, VERIFIED on 9 disc reads, and it moves no precision number**
because the position was already there.

That is also the reconciliation of something three handoffs carried without noticing: the
verification pool grew **90 → 94** while the call set grew **82 → 85**. The fourth entry is this
conversion. **When a change can rewrite a field on a call it did not create, diff the field and not
only the position set** — and note the corollary: a `TD_SRC` conversion is invisible to a position
metric in *both* directions.

### 5.2 So the item closes on a statement about instruments

Against a 90-record truth set with a 5-FP budget, one SVA call is ~1.2 precision points. The
93.9% → 90.6% is three calls, of which **one cannot be matched by construction**, one is
**read-VERIFIED**, and one is **undecidable by both instruments**. That is neither "the flip is free"
nor "the flip costs precision", and it should be quoted as the split rather than as the delta.

### 5.3 Two facts about the truth set itself, one of which §5 depends on

**The truth set exists in two directories and they are byte-identical** — sha256 equal for all five
VCFs across `claude_sandbox/test_data/HG002/truth` (the one every harness names in its own output
header) and `dev_xTea/truth_set_hg002`. §5's chrY argument would be much weaker if it were a
property of one copy; it is not.

**`all_tes.bed` — in the second copy only — is NOT a reference-TE annotation**, and the name invites
the assumption. Such a file would be directly useful for the repeat-context questions in §6 item 5.
Measured, it is exactly the **union of the three truth VCFs' positions** (1642 rows, every one at
distance **0** from a truth POS), so it holds nothing the VCFs lack — and it writes the 1-based VCF
`POS` straight into the BED **start** field with `end = start + 1`, so a BED-aware tool reads it one
base off. For a real element annotation use the RepeatMasker files under `RMSK`.

**And one truth position is claimed by two families**: `chr19:3742665` is in both `Alu_all.vcf` and
`SVA_tprt_all.vcf`. The per-family denominators are therefore not disjoint, and the scorer's
`--- all families ---` line — built by summing per-family TP/FN/FP — counts it twice. One record of
1642, so nothing reported moves. **What it is matters more than what it costs**: an insertion the
truth set could not assign to one of Alu or SVA, which is the SVA-is-SINE-VNTR-**Alu** confusability
plan §5 names as a mapping hazard, appearing in the truth set itself. Treat a single Alu/SVA
disagreement as an open question, not as a caller error.

## 6. NEXT

**In this order. Nothing below is blocked on compute, and nothing needs a gate re-run.**

1. **The decoy-consensus concatenation** — the last piece of the previous handoff's item 1, and now
   the top. The recipe is exact and build-independent (plan §5: 45 records, all reconstructible from
   `rep_lib/consensus/` byte-identically, and `consensus/transduction_decoy.fa` has **zero
   references** in xTea2 or master while `construct_novel_flanks`' comment claims it is concatenated
   in). §2.4 above adds a second reason to do it and a second instrument to measure it with: the
   asymmetry is now measured on the **clip** side as well as the disc side, so the concatenation has
   a before/after on both. Still call-moving; still needs its own arm, two scales, `55055ba`'s
   pattern.
2. **Gate [A]'s two semantic decisions, in this order.** They are decisions, not measurements, and
   §2 is the measurement they need:
   - **Separate `n_polya_cutoff` from `n_clip_cutoff`** (§2.3). One dial for two thresholds is the
     §5 defect; separating them moves calls, so it is a flag and an arm.
   - **Decide whether polyA-map membership should be a veto at all** (§2.2), given that it drops 68%
     / 78% of the clip-qualified population including 594 / 810 loci with `n_clip >= 10`. Expect
     *precision*: 306 of BASELINE-0's 322 FNs never reach this gate.
3. **`NEARBY_CLIP` at gate [A] is now an open question** (§3), and it is the largest lever measured
   in the step. It needs a truth-set-scored two-arm run, and the arm is a one-line env change rather
   than a code change. **Do not raise it on §3's numbers alone** — read §3's two reasons why the 4×
   is not 4× more evidence.
4. **§11's bare-shape writer fix** — still scoped, not built, per §10.3. §4 above adds a tripwire:
   the new disjointness test fails when this lands, deliberately.
5. **§3.4's repeat cutoff by family** — still explicitly NOT the route to justifying anything about
   the decoy fix (§10.1). §5 raises its stakes again: **two of the flip's three new calls are
   UNDECIDABLE only because of it**, at 2.7% and 3.8% against a 5% cutoff. That makes it the single
   parameter standing between "adjudicated" and "unadjudicated" for this population — which is
   exactly why moving it here would look like tuning until the answer arrived. Argue it on biology,
   in isolation, in its own change.
6. **A third arm that is neither master nor xTea2** — still the only thing that can take the
   no-regression margin out of PROVISIONAL (§10.4). Nothing on this list unblocks it.

## 7. Where things live

| | |
|---|---|
| **the agenda** | `docs/handoff/DOWNSTREAM_OVERHAUL_PLAN.md` |
| **the gate's design** | `docs/handoff/HANDOFF_2026-08-04_stage3_td_gate_spec.md` (corrections inline, `[CORRECTED]`) |
| the previous state handoff | `docs/handoff/HANDOFF_2026-08-05_decoy_default_flipped.md` |
| **the standing decisions** | `HANDOFF_2026-08-04_td_gate_built.md` §10 — four decisions, unchanged |
| **gate [A]'s instrument** | `docs/handoff/runs/probe_clip_gate.py` + `probe_clip_gate_wgs.sbatch` — job `49836568`, PASS, 57 s |
| **the off-truth instrument** | `docs/handoff/runs/probe_sva_flip_offtruth.py` — no job needed |
| the two orphans | `benchmark/benchmark_truth_set.py` — the block above `CALL_CLASS_*`, `SUBTYPE_ORPHAN_SITE`, `td_src_call_class`; the test is `tests/test_benchmark_truth_set.py::test_the_two_orphans_are_disjoint_by_construction` |
| gate [A]'s code | `xtea/transduction_evidence.py:355` (`gate_clip_evidence`) and `:376` (`select_clip_sites`); its parse is `xtea/x_transduction.py:851` |
| the verifier + the no-regression gate | `benchmark/verify_transduction.py` (`--baseline`) |
| the gate's frozen baseline | `/n/scratch/users/c/cos689/xtea2_testdata/td_gate_baseline_hg002/` (+ `SHA256SUMS`) — **untouched this session** |
| the truth set | `/n/data1/hms/dbmi/park/corinne/claude_sandbox/test_data/HG002/truth` — SVA is **90** records with **0 on chrY**. A second copy at `dev_xTea/truth_set_hg002` is **byte-identical** for all five VCFs, so §5's chrY claim is about the truth set and not a copy. That copy's `all_tes.bed` is **not** a reference-TE annotation — it is the truth union with a 1-based POS in a BED start field. See §10 |
| **a trap in the neighbouring harness** | `docs/handoff/runs/master_serial_wgs.sbatch` opens `rm -rf "$OUT"` and its default `OUT` is the NA12878 tree the gate spec cites. Never invoke it without `MASTER_OUT` |
| running the tests | `python3 -m pytest -m "not slow"` on the system 3.9.25 — **234 passed / 30 skipped** |
