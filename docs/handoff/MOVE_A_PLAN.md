# Phase C — Move A (locus-keyed type resolution): prepared plan

**Prereq:** WGS seeded deep-fold (move B) comes back green. Auto-gate running (bg `bkvel57dx`).

## What I pinned during prep (so implementation is de-risked)

### 1. The type-resolution rule is trivial and quality-free
`x_vcf_merger.merge_single_sample_ins` + `cmp_quality_of_two_ins` (x_vcf_merger.py:629/746):
- Load order = **L1 → SVA → Alu → HERV** (= priority).
- Two calls within **±50 bp** = one locus (conflict).
- Resolution: **keep the already-placed call UNLESS the incoming one has `two_side_tprt_both`
  in its VCF line and the placed one does not.**
- **`cmp_quality_of_two_ins` never compares numeric QUAL** — the "check quality first" comment is
  dead code. Resolution depends ONLY on (TPRT-both flag, load order). This is the whole thing.

### 2. The divergence surface is 18 loci genome-wide (enumerated)
From the frozen WGS golden (`conflict_analysis.py`): 1538 per-type calls → 1520 merged → **18 conflict
loci**. Pairs: L1/Alu 13, Alu/SVA 4, L1/SVA 1. 13/18 resolved by TPRT-override, 5/18 by pure priority.
One 3-way pileup (chr1:~241980550, L1+SVA+Alu → SVA). These 18 are the *entire* set of loci move A
could flip; every other locus is single-type and cannot diverge. Full list saved in
`wgs_seeded_gate_RESULTS.txt` sibling analysis.

### 3. TPRT-both is known at the post-filter stage
`two_side_tprt_both` is set in `x_post_filter.py` (:264/307/782). So the type decision needs only
per-type output **through postF** — NOT gene/genotype/VCF. Type can be assigned immediately after postF.

### 4. cns_filter is ~type-agnostic (design §6.2, already mapped)
Only per-type diffs: consensus ref (scan already realigns all 3) + SVA's `IS_CALL_SVA` relaxed
consistency. So a single per-locus processor with a "consensus slot + one `if sva`" covers all types.

## Consequence: move A is LOWER risk than the design doc feared
Design §5 worried borderline multi-type loci flip vs golden. But a flip can happen ONLY if either the
per-type **call set** changes or the **TPRT-both flag** changes — both already covered by the existing
call-set concordance gate. Given the same per-type calls + same TPRT flags, the 18 resolutions are
deterministic and reproduce the merger exactly. So move A's risk ≈ "did we preserve per-type calls +
TPRT flags," which the harness already tests.

## Staged implementation (each harness-gated; STOP for approval before starting)

- **A.1 — Re-express C.1 as an explicit locus table (zero behavior change).** Replace the
  `SampleVCFMerger` call in run_xtea.py (:564–570) with our own `LocusTable` that implements the pinned
  rule (load L1→SVA→Alu, ±50bp cluster, TPRT-override). Gate: merged VCF call-identical to current C.1
  golden (chr1 32 / WGS 1520). Deliverable: merger deleted, `Locus`/`assign_type()` are first-class,
  validated. Lowest risk; foundation for A.2.
- **A.2 — Run the type-agnostic tail ONCE.** After the per-type loop produces postF call sets, build
  the locus table (assign type via A.1), then run **gene → genotype → VCF once** over the unified set
  instead of 3×. Gate: merged VCF unchanged; per-type recall preserved. Removes 3× tail invocation +
  the merge step. (This is design's move C.2 + A combined, now cheap because type assignment is
  quality-free.)
- **A.3 — Assign type at cns/postF; one downstream loop.** Collapse the two-phase per-type loop into a
  single per-locus loop: cns evidence per type side-by-side (move B already collects once) →
  transduction branched by assigned type (L1/SVA only, move D) → postF → single tail. Highest-value,
  highest-risk; gate hard on WGS AND explicitly re-check the 18 conflict loci resolve identically.

**Sequence rationale:** A.1 proves the rule in isolation; A.2 collapses the tail with the rule proven;
A.3 reorders compute-vs-resolve last. Each step keeps the concordance harness green (100% golden recall,
watch PPV) and preserves `--keep-intermediates` for debugging (design §5).

## Open items to confirm at A.1 start (cheap)
- Exact VCF INFO field carrying `two_side_tprt_both` (merger matches substring in the whole line) —
  confirm it survives into the final per-type VCF unchanged so the locus table reads it the same way.
- HERV is in the merger's order but we only run L1/Alu/SVA — locus table can ignore HERV.
- Seeded vs frozen golden: run the 18-conflict analysis on the *seeded* merged golden once regenerated
  (seed may shift borderline calls, but conflict COUNT should be ~stable); use it as the A.1 gate target.
