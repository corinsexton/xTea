# Mosaic mode: prepared plan

**Prereq:** CC.3 (case/control real-data gate) has an open, unresolved discrepancy (xtea2 vs master's
LFS 2921/L1 diff: 1/2 calls match exactly, the other diverges — chr6:82922665 (xtea2) vs
chr7:5734825 (master), not a near-miss). Investigating that is the very next thing after this doc,
not blocking mosaic planning, but worth knowing this repo has one open correctness question in flight.

**User's stated priority for this feature (their words):** "The most important thing with mosaic
besides accuracy is speed! This is prohibitively slow on 300x data. It's okay if memory scales up. We
are trying to avoid the intermediate files and exorbitant runtimes. I can give many cores." This is the
single biggest driver of the design below — more important than replicating master's algorithm
structure step-for-step (unlike case/control, where matching master's exact command shape was the
whole point).

## What two research passes pinned (so implementation is de-risked)

### 1. Mosaic is single-sample-per-run; multi-sample "joint calling" is optional and separate
Earlier project memory flagged mosaic as needing "multi-sample joint-calling orchestration xTea2's
architecture has no slot for at all" — **this was too pessimistic, confirmed by reading the code:**
- The mosaic clip entry point (`x_TEA_main.py:828`, `elif options.mosaic:`) takes a bam_list the same
  shape as every other mode — it CAN contain multiple BAMs from the same individual (the code comment
  at `x_TEA_main.py:845` says "merge the list from different bams of the same individual"), but that's
  the SAME multi-BAM-list mechanism xtea2 already supports today, not a new orchestration shape.
- Cross-*sample* comparison (`-J`/`--joint`, `x_TEA_main.py:104`, dispatched at `x_TEA_main.py:1012`,
  `x_joint_calling.py:18` `MosaicJointCalling.call_mosaic_from_multi_samples`) is a **separate, optional,
  post-hoc step**: run mosaic independently per sample, then feed a list of each sample's own result
  files into `--joint --postFmosaic`, which does a singleton filter (a call present in only one sample
  of the list is kept; calls repeated across samples are treated as germline-like and dropped). **Not a
  joint caller in the statistical sense** — no per-sample read-level comparison, just a post-hoc set
  operation on already-finished per-sample candidate lists.
- **Conclusion: single-sample mosaic support needs zero new orchestration shape in `run_xtea.py`** (same
  one-CRAM-in, one-VCF-out contract as everything else). `--joint` is out of scope for a first pass —
  it's a separate, cheap, file-level operation that can be added later without touching the core path.

### 2. What's actually different about mosaic detection (confirmed by reading the code, not guessing)
- **Clip step** (`x_TEA_main.py:828-859` dispatch → `x_TEI_locator.py:182`
  `call_TEI_candidate_sites_from_multiple_alignmts_mosaic` → `x_TEI_locator.py:654`
  `call_TEI_candidate_sites_from_clip_reads_v2_mosaic` → `clip_read.py:867`
  `run_cnt_clip_part_aligned_to_rep_by_chrm_mosaic`): same clip-collection/realignment shape as
  germline, PLUS a per-clipped-read polyA tally (`clip_read.py:967`,
  `xpolyA.is_consecutive_polyA_T(s_clip_seq_ck)`, `x_polyA.py:54` — a cheap substring-membership check,
  `"AAAAA" in seq` etc, O(read length) per read) accumulated per position (`m_polyA[pos]`, used later by
  the mosaic-specific post-filter). Automatic clip cutoffs still come from the same
  `automatic_gnrt_parameters` table germline uses (`x_TEA_main.py:853`, no tumor/purity involved — mosaic
  isn't framed as a tumor concept in master, just "normal illumina data" per the code comment at
  `x_TEA_main.py:828`).
- **Disc step**: reuses the *exact same* disc code as case-control's `--somatic` disc mode — the only
  difference is `PEAK_WINDOW` drops from 100bp to 30bp when `options.postFmosaic or options.somatic`
  (`x_TEA_main.py:440-442`). Not a new algorithm; already the same code path xtea2 runs today (and can
  already be routed through `--scan-disc` if that flag threading is extended to pass the smaller
  window — open item, see below).
- **Post-filtering** (`x_mosaic_calling.py:20` `MosaicCaller(XPostFilter)`, `run_call_mosaic` at line 30,
  per-type methods `call_mosaic_L1_from_bulk` etc at lines 56/101/139): pure in-memory record filtering
  over the already-generated candidate list — polyA-signal requirement, allele-frequency cutoff
  (hardcoded `(0.35, 0.02)` per type, `x_mosaic_calling.py:173-179`), a "polyA-dominant-on-both-sides"
  rejection (catches germline-like sites), a same-family repeat-region rejection, and a blacklist
  dedup. **Zero BAM re-reads** — unlike case-control's postF, which re-collects/realigns reads from the
  control BAM per candidate (`x_somatic_calling.py`), mosaic's postF only touches the already-loaded
  candidate-list records and a precomputed repeat interval tree. Cheap and does not scale with coverage.

### 3. Where the actual 300x cost lives — confirmed, not assumed
Cost model, germline (30x) vs mosaic (300x), by step:

| Step | What scales with depth | Verdict |
|---|---|---|
| Clip (candidate generation + realign) | Yes — same as germline clip, O(#clipped reads) | Already the exact cost `xtea-scan`'s C++ rewrite solved for germline; mosaic reuses the identical Python-loop shape (`clip_read.py`), just not yet ported |
| Clip's polyA tally | Yes, same loop, same reads | **Not itself expensive computation** (a few substring-membership checks per read, cheap even in Python) — the real cost driver is running *any* meaningful per-read Python-level work at 10x the depth, i.e. the same root cause `xtea-scan` already fixed for plain clip collection |
| Disc | Yes — O(#discordant reads genome-wide) | Identical code to today's disc step (germline/case-control); already has a C++ path (`--scan-disc`) if wired through |
| PostF (mosaic-specific) | **No** — O(#candidates), pure in-memory record filtering, no BAM access | Already cheap; not a target |

**Corrected framing vs. the first research pass:** the bottleneck isn't "the polyA check is expensive" —
`"AAAAA" in seq` is trivial. The bottleneck is doing **any** per-read Python-loop processing at 300x,
which is exactly the problem `xtea-scan` was built to solve for germline clip detection. Mosaic's clip
step just hasn't been folded into that C++ engine yet, and it happens to need one extra cheap
per-read fact (polyA membership) that the current C++ engine doesn't expose.

### 4. A related, general finding (not mosaic-specific, found while chasing this) — worth doing anyway
The **germline** (non-mosaic) clip collector, `clip_read.py:collect_clipped_parts_by_chrm` (~line
452-540), *also* does a polyA check today — a different one (`xpolyA.contain_enough_A_T`, a min-count
check, not substring membership), applied only to very short clips (length in
`[MINIMUM_POLYA_CLIP=7, BWA_REALIGN_CUTOFF=9)`, `global_values.py:26-27` — a narrow 7-8bp band). **This
check is NOT implemented anywhere in `xtea-scan.cpp` or `scan_adapter.py` today** — confirmed by direct
grep, zero hits. So xtea2's existing `--scan-clip` fast path has been silently omitting this filter for
every germline and case-control run so far. Given the narrow length band affected, this is plausibly
small enough to have been absorbed into the pipeline's already-accepted ~99.2%/98.3% WGS wobble rather
than showing up as a discrete bug — but it means **the C++ engine work below should implement polyA
support generally** (both the germline min-count filter and the mosaic substring-tally), not as a
mosaic-only add-on, since it closes a real existing gap in the already-shipped fast path too.

### 5. No real executed mosaic run found on disk (unlike case/control's free LFS 2921/L1 gate)
Searched `/home/cos689/data1/DATA/`, `/home/cos689/LFS_light_et_al/`, broader `/home/cos689/` and
`/n/data1/hms/dbmi/park/corinne/` for anything mosaic-named or mosaic-flagged (config/script content)
— nothing found, including in `mei-pediatric-cancer/individual_cases/PAUTSB/` (a pediatric cancer case
where mosaicism would be a plausible real use case, previously found during case/control research).
**No free real-data validation target exists for mosaic**, unlike case/control. Validation options:
1. **Generate our own golden**, the same way the whole xTea2 project's germline goldens were made: run
   master's own `--mosaic` pipeline end-to-end in the `xTea` conda env (the "magical env" — see
   memory `xtea-magical-env`) against some real or high-depth BAM, and diff xtea2's mosaic output
   against that. This is the established, proven pattern for this whole project — not a new risk.
2. Self-vs-self / synthetic smoke tests (same shape as CC.0-CC.2) for early plumbing gates before the
   real golden comparison.
3. **Resolved:** user will supply a real 300x CRAM next session — M.3 can proceed once it's in hand,
   no data-sourcing step needed first.

## Staged implementation (mirrors the case/control CC.0→CC.3 pattern; each gated, STOP for approval
before starting any stage — nothing below is authorized to start yet)

- **M.0 — generalize polyA support in `xtea-scan.cpp` (biggest lift, does double duty).** Add both
  polyA operations to the existing per-clipped-read hot loop (the same loop that already produces
  clip candidates for `--scan-clip`): (a) the germline short-clip min-count filter
  (`contain_enough_A_T`-equivalent, gated to the same 7-8bp band) as a bug-fix/generalization with no
  behavior-changing intent outside that narrow band; (b) the mosaic per-position polyA/non-polyA tally
  (`is_consecutive_polyA_T`-equivalent substring check), emitted as extra columns the mosaic clip
  adapter can consume. Gate: existing germline chr1/WGS goldens re-run with the new filter active,
  confirm no change outside the narrow 7-8bp band (or a small, explained, accepted change if it turns
  out to matter more than expected) — this stage is validated against data we already have, before any
  mosaic-specific code exists.
- **M.1 — mosaic clip step in xtea2, chr1 scale.** New `--mosaic` flag; a `mosaic_adapter.py` (mirrors
  `scan_adapter.py`'s shape) that reads `xtea-scan`'s per-position polyA tally output (from M.0) instead
  of shelling out to master's slow Python mosaic clip function at all — this is the key departure from
  the case/control precedent (there, we replicated master's own step chain faithfully; here, the user's
  priority is speed over structural fidelity, so we go straight to the fast path rather than building a
  slow-then-fast pair the way germline/case-control did historically). Gate: chr1, confirm the adapter's
  candidate list is sane (plausible call count, no crash) — no real target to diff against yet at this
  stage.
- **M.2 — mosaic post-filtering (`MosaicCaller`) wired in.** Since this step is already cheap and pure
  in-memory (per §2), no C++ work needed — just wire `--postFmosaic`-equivalent Python logic into
  `run_xtea.py`'s per-type chain (a new step, gated on `--mosaic`, analogous to how `--case-control`
  added two new steps). Gate: chr1 end-to-end (clip → disc(30bp window) → cns/transduction/sibling →
  postF-mosaic → gene → gVCF), sane output.
- **M.3 — real/deep-coverage validation + timing.** Resolve the open data-sourcing question above
  first. Generate a master-run golden (§5 option 1) at real depth (300x if available), diff xtea2's
  mosaic output against it for correctness, AND measure wall-clock/CPU at that depth as the actual
  bar-clearing check for this feature's stated #1 priority — this is the gate that matters most, same
  role CC.3 played for case/control (the earlier stages just prove the plumbing doesn't crash).
- **M.4 — `--joint` multi-sample singleton filtering (optional, likely deferred).** A separate, cheap,
  file-level post-hoc step per §1 — not required for core mosaic correctness or speed, low risk to add
  later once M.0-M.3 are solid. Not staged with a number ahead of M.3 on purpose: nothing here affects
  the speed story, so it shouldn't compete for attention with the parts that do.

**Sequence rationale:** M.0 first because it's the one piece of work that pays for itself twice (fixes
a real, if narrow, existing gap in the already-shipped germline/case-control fast path, AND is the
prerequisite for mosaic's own speed). M.1-M.2 get a working, fast, correctness-unvalidated pipeline
using data we already have (chr1). M.3 is the actual bar (deep-coverage correctness + the timing number
that answers the user's actual ask) and needs a data/target decision first. M.4 is explicitly
lowest-priority — optional, doesn't touch the speed goal.

## Open items to confirm before M.0 starts
- Confirm no other master mosaic-adjacent behavior depends on the *exact* Python polyA functions being
  used verbatim (vs a faithful-but-C++ reimplementation) — low risk given both are simple string
  operations, but worth a quick pass over `x_polyA.py` for any other callers before assuming this is
  a clean lift-and-port.
- Whether `--joint` (M.4) is wanted at all in this pass, or should be dropped from scope entirely rather
  than deferred — low cost either way, just a scope call.

## A related, unrelated-to-mosaic investigation also in flight — TSD is NOT a speed concern
User asked whether TSD (target site duplication) computation has the same kind of per-read efficiency
problem as polyA. Checked: no. `x_TSD.py:get_TSD_seqs` (60 lines total) is called only once, from
`x_gvcf.py`, at final gVCF-writing time, on the already-filtered final call list — one
`pysam.FastaFile.fetch()` reference lookup per surviving call, O(#final calls), zero BAM/read access.
Not a target for the C++ work above. (The `TSD=NULL` seen in the CC.3 example VCF rows is a call-type
characteristic, not a bug: `x_TSD.py:22-28` sets it explicitly whenever a record lacks valid two-sided
left/right clip positions, which the "orphan_or_sibling_transduction" subtype calls we happened to be
looking at don't have.)
