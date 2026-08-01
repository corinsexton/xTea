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

**Progress — 2026-08-01.** Stage 0 is now partly implemented:

| Stage 0 item | status |
|---|---|
| phase timers (measure 1) | **DONE** `6ce36c7` — 13 phases + total, unconditional, `[TD-PHASE]` to stderr |
| phase breakdown (measure 2) | **DONE** — job `49437990`. **Gather 73-86%**, bwa 12-25%, SAM parse <2%. Stage 1 conditional resolves to the C++ gather; see §4 Stage 0 |
| type-set-dependence fix | **DONE** `de96727` — the BED-plan helper, + `tests/test_mate_in_rep_type_partition.py` |
| `_tprt`'s meaning | **RESOLVED** — §7.1; canonical-mechanism filter, `_all` is the denominator |
| `col5` → `mate_in_rep_clip_cnt` rename | **DONE** — `mate_in_rep_bed_plan`, CLI `--mate-in-rep-families` (+ deprecated alias) |
| BASELINE-0 (both type sets) | not started |
| chr1 slices (HG002, SMHT005-3G) | not started |
| how transductions are scored | **design settled** — §7.1: derive truth status from the ALT residual (feasibility verified), 5-way error taxonomy, two denominators. `benchmark_truth_set.py` still to implement |

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
**Rollback:** legacy remains the default; the flag selects the new gather.

> **Ordering determinism is load-bearing.** `cns_emit.hpp:798-817` documents that bwa's tie-break among
> equal-scoring placements depends on read index in the input, so reordering the FASTQ silently changes
> which repeat copy multi-mapping clips land on. Candidate-major emission order must be preserved
> exactly, or Stage 1 cannot be agreement-gated at all.

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

**Deriving truth transduction status — feasibility VERIFIED, not assumed.** Aligning all 197
`L1_all.vcf` ALT sequences against `consensus/LINE1.fa` (24,702 bp over 4 subfamily entries) with
`bwa mem` and measuring the primary alignment's soft-clip residual:

| non-L1 residual | count of 197 |
|---|---|
| median 18 bp (polyA + TSD — no transduction) | — |
| ≥ 100 bp | **60** |
| ≥ 200 bp | **51** |
| max | 4620 bp |

So ~26-30% of truth L1 insertions carry substantial non-consensus sequence, and the residuals are
long enough (200-4600 bp) to map uniquely back to the reference and yield a **source locus derived
independently of xTea**. That is the truth-side transduction annotation the VCFs lack.

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

**Still to implement** (Stage 3 prep, in `benchmark/benchmark_truth_set.py`): the residual→genome
mapping and the full-length-L1 adjacency test, the `TD_SRC` parser, and the taxonomy above. Alu is
expected to carry no transductions at all (non-autonomous, §6.2), so this applies to L1 and SVA only.

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
`xtea/mosaic_adapter.py` (§6.4); `benchmark/benchmark_truth_set.py`; a new `xtea_scan` subcommand.

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
