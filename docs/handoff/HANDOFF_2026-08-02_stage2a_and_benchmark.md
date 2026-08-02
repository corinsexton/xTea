# xTea2 handoff — 2026-08-02 (Stage 2a done and gated; the Stage 3 benchmark exists and bites)

**This is the current entry point for state. It is NOT the entry point for agenda — the plan is.**

> **Read `docs/handoff/DOWNSTREAM_OVERHAUL_PLAN.md` first.** It still outranks this chain on agenda.
> Everything factual below is also written into the plan at the section that owns it — the plan is
> the durable copy, this is the index. Supersedes `HANDOFF_2026-08-01_stage1_gather.md`, whose §4
> "NEXT" items 1-3 are what this session addressed (item 1 deferred with reasons, item 2 started,
> item 3 delivered).

**One-line status: `global_values.IS_CALL_SVA` is gone, replaced by an explicit `FamilyParams`
bundle, and the chr1 agreement gate PASSES 9/9 on both gather paths. `benchmark/benchmark_truth_set.py`
exists and scores BASELINE-0 today — and building it corrected §7.1 three times and settled that the
transduction gate as specified has a denominator of ZERO for L1, for a reason no threshold fixes.
Stage 2c's gate now has a number in advance:
22 contended loci. Stage 1's WGS default-flip gate was in flight at write time (job `49480961`).**

## 0. State at handoff

| | |
|---|---|
| branch / HEAD | `xTea2` @ `5c68481` — 6 commits this session, `de529a4`..`5c68481` |
| working tree | clean |
| tests | **137 passed, 1 skipped** with `XTEA2_SCAN_BIN` set, excluding `test_downstream_from_cns.py` (the slow chr1 golden, run inside the gate jobs). The skip is `test_e2e_concordance`, needs `XTEA2_OUTPUT_DIR`. Was 108 passed before this session |
| in flight | **job `49480961`** — the Stage 1 WGS default-flip gate. §1 says what to do with each outcome |
| the agenda | `docs/handoff/DOWNSTREAM_OVERHAUL_PLAN.md` §4 Stage 2c, then Stage 3 |

## 1. The one job still running, and what to do with it

`docs/handoff/runs/stage1_flip_wgs.sbatch`, job `49480961`, output
`/n/scratch/.../stage1_flip_wgs_49480961{,.out}`. It is the gate the previous handoff left open:
*"what remains before flipping `--td-gather-scan` [to default] is one WGS run whose calls match a
BASELINE-0 arm — the chr1 gate does not substitute, being 17 calls."*

**Two arms, not one**, and the reason matters if you re-run it: BASELINE-0 was frozen at `9a43b80`
and HEAD is later, so a lone native arm against it conflates "the gather is faithful" with "nothing
else moved in between". Arm `LEGACY` isolates commit drift; arm `NATIVE` vs `LEGACY` isolates the
gather; arm `NATIVE` vs BASELINE-0 is the flip criterion. Arm LEGACY finished in 48:08.

- **PASS** → `--td-gather-scan` may become the default. Note what flipping *also* does: it makes the
  multi-BAM divergence reachable, so the mosaic path would then run the fixed behaviour with no
  baseline of its own. Stage 4 concern, not a blocker.
- **arm LEGACY differs from BASELINE-0** → the gather is not implicated; something in `937f0e6`..
  `58d2203` moved, and the diff is small enough to read (`git diff 9a43b80..58d2203 -- xtea/ src/`).
- **arm NATIVE differs from arm LEGACY** → the gather is implicated at WGS scale despite being
  byte-identical on its own artifacts at WGS scale (job `49453857`). That would be surprising and
  the artifact oracle is where to look first.

## 2. What was built

### 2a — the family parameter bundle (`de529a4`), gate PASS

`xtea/family_params.py`. `FamilyParams` carries `has_flank`, `annotation_extend`, `leftmost_cns`
and the two SVA cns-path relaxations; it is passed to `XClipDiscFilter` / `XTransduction` /
`XOrphanTransduction` rather than flipped process-wide. `global_values.IS_CALL_SVA` and
`turn_on_sva()` are **deleted**, not left unused — a surviving attribute would let a stale
`global_values.IS_CALL_SVA = True` assignment keep parsing while doing nothing (one existed, in
`docs/handoff/runs/td_gather_oracle.py`; fixed).

**Gate: job `49481391`, `docs/handoff/runs/refactor_gate_chr1.sbatch` — PASS.** HG002 chr1, both
gather paths, 9/9 artifacts byte-identical to the Stage 1 gate's arms A and B respectively. The
validity condition that script cannot check itself was checked by hand:
`git diff 49154dc..HEAD -- xtea/ src/` is the gather port plus this refactor and nothing else.

**The finding this surfaced, which is the reason to do refactors like this at all.** `--sva` is
passed to the **cns_filter step only** — `build_type_steps` matches the golden
`run_xTEA_pipeline.sh`, which puts it on -C/-D/-N and not on transduction/sibling/postF. So
`IS_CALL_SVA` was False inside the transduction subprocess even for SVA, and **SVA's
`leftmost_cns = 5000` has never executed, in xTea2 or in master.** It is preserved dormant. Waking
it is not tidying: `consensus/SVA.fa`'s longest record is 3183 bp, so a 5000 threshold makes
`_reach_5_prime_end` true for every candidate with a positive cluster start and the
polymorphic-flank library stops being a filtered subset. `tests/test_family_params.py` pins the
resolution rule by reading the argv the driver emits, so the change has to be deliberate — and it
belongs in Stage 3, against the truth set.

### 2b — the manifest, DEFERRED into 2c with arithmetic

The previous handoff's NEXT item 1. Not built, on purpose:

- **The wall-clock case is ~zero.** L1's and SVA's transduction steps are concurrent subprocesses,
  so their native gathers (53.3s and 56.8s at WGS scale) **overlap**. Sharing saves CPU, not time.
  §1 of the plan already says the sharing win "only materialises when the shared work runs once at
  full width", and full width is what 2c creates.
- **The caller half is throwaway.** Hoisting the gather means the driver computes the step's
  `iextnd` (from `automatic_gnrt_parameters`) and runs `re_slct_with_clip_raw_disc_sites` itself —
  re-implementing step logic 2c deletes.

The C++ half is reusable and is what 2c should build: generalise `td_gather.hpp` from one `TdSites`
to a vector, plan fetches over the **union** of the entries' windows, emit per entry. The union is
nearly free — 53,945 of ~54K sites are shared by all three families.

### 2c's gate baseline (`1c56b8d`) — 22 contended loci

`locus_table.LocusTable` now records every collision it resolves (winner, loser, and whether
load-order **priority** or a **`two_side_tprt_both` override** decided it); `run_xtea` writes
`merged/locus_family_conflicts.tsv`. Diagnostic only — the merged hc_pf list and the type map are
byte-identical whether or not the sidecar is written, and a test asserts that.

Measured on BASELINE-0:

| winner | loser | rule | n |
|---|---|---|---|
| Alu | L1 | `tprt_override` | **12** |
| SVA | Alu | `priority` | 6 |
| L1 | SVA | `priority` | 2 |
| L1 | Alu | `priority` | 1 |
| Alu | SVA | `tprt_override` | 1 |

It reconciles exactly: 200 L1 + 1230 Alu + 85 SVA per-type calls − 22 = **1493** merged at
188/1223/82. The 12 `Alu over L1` rows are precisely the 12 the Stage 0 gate found by hand comparing
`--types L1` (200) with `--types L1,Alu,SVA` (188). **So Stage 2c's gate has a number before the
work starts: a run that moves substantially more than ~22 loci, or moves one not on this list, is
not doing family reassignment.**

### Stage 3 prep — `benchmark/benchmark_truth_set.py` (`7f0bb95`, `1c56b8d`)

Two subcommands. `derive-td` (truth VCFs + bwa + consensus + reference + rmsk → a cached
`truth_transductions.tsv`) depends only on the truth set and the reference, so it is an artifact you
regenerate when those change, not per scored run — `docs/handoff/runs/derive_truth_td.sbatch` is its
harness, output `/n/scratch/.../truth_td_2222ec6/`. `score` runs in seconds and needs no external
tool.

**BASELINE-0 scored** (±200 bp, one-to-one nearest; **not comparable to the paper's numbers** —
hs38DH, §7.1, and the report says so on every run):

| | truth | calls | TP | FN | FP | sensitivity | precision |
|---|---|---|---|---|---|---|---|
| L1 | 197 | 188 | 140 | 57 | 48 | 71.1% | 74.5% / **92.7%** excluding orphan sites |
| Alu | 1355 | 1223 | 1103 | 252 | 120 | 81.4% | 90.2% |
| SVA | 90 | 82 | 77 | 13 | 5 | 85.6% | 93.9% |
| all | 1642 | 1493 | 1320 | 322 | 173 | **80.4%** | **88.4%** |

## 3. The two findings that change what Stage 3 can gate on

**READ BOTH BEFORE DESIGNING STAGE 3's GATE.** They are in plan §7.1 too.

**1. The transduction denominator is 0 (L1) and 2 (SVA). This route cannot gate Stage 3's
transduction work on HG002, and that is a property of the data rather than a mis-set threshold.**
Attrition (job `49482568`): 197 L1 truth records → 189 align to the consensus → **20** carry a 3′
residual ≥100 bp → 10 map back at MAPQ ≥30 → 7 also align over ≥80% → 4 do not land inside a
reference L1 → **0** sit within 10 kb of a full-length (≥5000 bp) reference L1.

Two measured reasons, and the second settles it:

- `hg38_L1.fa.out` is fragment-dominated — 978,710 intervals, median **285 bp**, only 1.1% at
  ≥5000 bp — so "adjacent to a full-length reference L1" asks for a rare coincidence.
- **Most long 3′ residuals are divergent element body, not transduced sequence.** 11 of the 20 map
  *inside* a reference L1 and **every one is an older subfamily**: L1PA6/7/8 at 4-7% divergence,
  L1MB/MD/ME at 16-25%. `consensus/LINE1.fa` is 4 L1HS-centric records, so an older L1's body
  partly fails to align and surfaces as "residual". The one record that used to pass the whole
  conjunction sat inside an **L1ME3G at 25.2% divergence** — divergent body next to a full-length
  L1 by coincidence. Adding the containment test took L1 from 1 to **0**.

So do not loosen the thresholds: that admits body, not transductions. The live options are a
family-complete consensus, the paper's supplement if it enumerates transductions, or a different
sample. SVA is less affected (88-record library, no SVA residual lands inside a reference SVA,
2 survive) — but 2 is still not a gate.

**2. All 37 of BASELINE-0's L1 orphan-transduction calls are off-truth, and they are 37 of L1's 48
false positives.** An orphan transduction *site* is where transduced sequence landed **without** the
element body — not an MEI insertion, and this truth set enumerates MEI insertions, so an orphan call
cannot match whatever the caller did. The split is total rather than noisy, which is what makes it a
category error and not a quality signal: the 125 L1 calls with no `TD_SRC` are 118/125 on-truth, the
26 sibling calls 22/26, the 37 orphan calls **0/37**. The report prints both precisions as a named
stratum. Do not quote either alone.

## 4. NEXT

1. **Read job `49480961`'s result** and act per §1. If PASS, flipping `--td-gather-scan` to default
   is a one-line change plus a note in §4 Stage 1's rollback clause.
2. **Decide how transductions get scored at all** (§3 item 1) *before* Stage 3, not during it.
   The threshold sweep has been done and it is not a threshold problem: the answer has to change
   the *input*, not the cutoffs. Cheapest experiment first — rebuild `consensus/LINE1.fa` to cover
   L1PA/L1M subfamilies and re-run `derive-td`; if the residual population collapses, the
   divergent-body confound was the whole story and the surviving residuals are real candidates.
   `docs/handoff/runs/derive_truth_td.sbatch` re-runs in ~2 minutes.
3. **Then Stage 2c**, which is the remaining large piece. Suggested decomposition, informed by what
   was measured today rather than by the plan's prose:
   1. **The locus evidence table.** `scan_clip_realign.tsv` already carries
      `chrom pos family nL nR mate_in_rep` — three rows per position, i.e. per-family support side
      by side, for 8.6M positions on WGS. The table Stage 2c needs is a read of that plus the disc
      side, not a new computation.
   2. **Extract the assignment rule.** It currently lives inside `LocusTable._place` as
      priority-plus-TPRT-override over three finished call lists. Lift it to a function over
      per-family evidence so the same rule can run at either end of the pipeline, and diff the two.
   3. **Move the decision up front**, run each family's chain over only its assigned loci, and gate
      on the 22-locus population in §2c above.
   4. **Then** the manifest (§2b) and `--emit-per-type`, which is deliberately last: flipping
      per-type VCFs to opt-in breaks every gate script that `cmp`s them, so it wants to land when
      the gates are being rewritten anyway.

## 5. What NOT to redo

Plan §10 is the authority and gained seven entries this session. The ones most likely to cost
someone a day:

- **`--sva` reaches cns_filter and nothing else**, so SVA's `leftmost_cns = 5000` is dead code.
  Do not "fix" it as tidying — see §2a.
- **`consensus/*.fa` are multi-copy libraries, not single consensus sequences.** SVA 88 records
  (longest 3183 bp), LINE1 4 (6342), ALU 47 (490). A threshold against "the consensus length" must
  be against the longest *record*.
- **Plan §7.1's ALT-residual figures (60 / 51) were total soft clip, both ends** — reproduced
  exactly at 17/60/51/4619. The 3′ residual its own criterion asks for is **20 / 12**. Anything
  sized off the 60 is wrong.
- **bwa splits a long query rather than clipping it when the far end also hits the reference.**
  38 of 189 aligned L1 truth ALTs carry a supplementary hit to the consensus. "Residual = the
  primary's soft clip" counts that second consensus hit as non-consensus sequence.
- **A minus-strand insertion's ALT has its 3′ end on the LEFT.** VCF ALTs are reference-forward.
  Taking the primary's trailing soft clip — which looks right, since bwa revcomps SEQ — returns the
  5′ end on half the records.
- **`hg38_L1.fa.out` is fragment-dominated**: 978,710 intervals, median 285 bp, 1.1% at ≥5000 bp.
- **`TD_SRC`'s region delimiter is `-` and an unknown start is literally `-1`**, so
  `chr10:-1-108181069` must not be split naively. Split on the LAST `-`; the parser raises if the
  end is ever negative.
- **Never quote L1 precision on this truth set without naming the stratum**: 74.5% overall,
  92.7% excluding orphan sites.

## 6. Where things live

| | |
|---|---|
| **the agenda** | `docs/handoff/DOWNSTREAM_OVERHAUL_PLAN.md` (§4 Stage 2 and §7.1 both updated today) |
| the family bundle | `xtea/family_params.py`; tests `tests/test_family_params.py` |
| the conflict recorder | `xtea/locus_table.py` (`Conflict`, `LocusTable.conflicts`) → `merged/locus_family_conflicts.tsv` |
| the benchmark | `benchmark/benchmark_truth_set.py`; tests `tests/test_benchmark_truth_set.py` |
| the derived truth TD table | `/n/scratch/.../truth_td_2222ec6/truth_td.tsv` (job `49482568`) |
| the reusable refactor gate | `docs/handoff/runs/refactor_gate_chr1.sbatch <worktree> [ref_run]` |
| the WGS flip gate | `docs/handoff/runs/stage1_flip_wgs.sbatch`; output `stage1_flip_wgs_49480961*` |
| the truth-TD derivation | `docs/handoff/runs/derive_truth_td.sbatch` |
| BASELINE-0 | `/n/scratch/.../baseline0_hg002_{L1-Alu-SVA,L1,SVA}/BASELINE-0/` (intermediates kept) |
| truth sets | `test_data/HG002/truth/` (headerless VCFs + the paper's `x_cmp.py`) |
| worktrees added today | `/n/scratch/.../stage1_flip_wt` (`58d2203`), `/n/scratch/.../stage2a_wt` (`15081a2`) |
