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
> **Still open in 2c**, deliberately not built: the single per-locus loop that replaces the three
> subprocess chains, the shared gather manifest (§2b), and `--emit-per-type`. The chains are still
> three, and there are still two decision points.

### Stage 3 — The evidence-model rewrite (truth-set gated)

Agreement-gating stops being coherent here, because the decisions themselves change. See §5.
**The truth-set harness must be working before this stage starts.**

**Gate:** per-family sensitivity/PPV against the HG002 truth set, no regression beyond a margin
calibrated on BASELINE-0; BASELINE-0 delta reported for attribution but **not** as pass/fail.

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

**Unify sibling / orphan / regular transduction into one call class with a typed source.** These are
epistemological categories, not biologically distinct events — the evidence shape is identical (clip
polyA signature + one-sided disc source, or two-sided anchors on a common region). Model as
`source ∈ {reference_element, polymorphic_element, inferred_region}`, one decision function for all
three, distinction surviving only as an output label. This collapses `x_orphan_transduction.py` into the
main path.

**Alignment passes: target two, not three.** Today: clip→flank library, disc→flank library,
disc→consensus. The flank library already has consensus concatenated in (`x_transduction.py:186`), so the
third pass's marginal contribution is position refinement plus the cross-family decoy tally — measure
whether it changes any call before keeping it.

> **On merging the per-family flank libraries into one reference:** a design pass argued this is safe
> because the families are evolutionarily distinct. **That reasoning is wrong for SVA** — SVA is
> literally SINE-VNTR-**Alu**, carrying Alu-derived sequence, so SVA/Alu mapping competition is a real
> hazard. The existing three-way decoy tally (`x_transduction.py:1103-1124`), which explicitly counts
> wrong-family hits, is evidence the original authors hit exactly this. Treat library merging as an
> **experiment with a measured outcome**, not a refactor, and run it only after the truth-set gate exists.

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

**New:** `xtea/downstream_unified.py` (locus evidence table, per-locus loop, processor classes);
`xtea/transduction_evidence.py` (`SiteEvidence`, `decide()`, `TransductionParams`);
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
