# Phase C — Unified single-pipeline downstream (design)

**Status:** design approved 2026-07-18 (user: "push forward with exactly this plan"). Supersedes the
earlier "parallelize the three per-type pipelines" framing of Phase C in the plan file.

**Goal.** Replace the three separate per-type downstream pipelines (each a serial chain of master
`x_TEA_main.py` subprocess steps, run L1→Alu→SVA) with **one pipeline keyed by genomic locus**, so the
tool is the salmon/bwa-mem shape end to end: **BAM → one process → one VCF**, no per-type forks, no
final VCF merge/dedup.

---

## 1. Why this is natural (grounding facts)

1. **The three pipelines process the same loci.** On NA12878 chr1, of ~46,232 candidate clip positions,
   **46,224 are identical across all three types** (99.98% overlap; L1 has 1 unique, Alu 7). The
   "three pipelines" differ only in (a) which consensus reads are realigned to, and (b) type-specific
   downstream (transduction). Keying by **locus** instead of **type** removes triplicated work, not
   meaningful independence.
2. **Master already unifies — but only at the very end.** `x_vcf_merger.merge_single_sample_ins`
   (`x_vcf_merger.py:629`, conflict logic ~722–736) takes the three type-VCFs and, when two types call
   the same site, resolves by quality then a fixed type-priority order. So "one insertion → one type"
   already happens; master just pays to compute all three type-calls independently first.
3. **Realignment is done twice.** Two distinct realignment stages exist today:
   - **Scan (A1, `xtea-scan`):** one BAM pass; extracts **clipped parts**, realigns them to each
     family's consensus → per-family `nL/nR` clip support (+col5). Purpose: **candidate detection**.
   - **cns_filter (`-N`, still master Python):** **re-opens the BAM**, re-collects clip **and** disc
     reads per candidate, realigns them (clip + disc, to consensus+flank) → the actual **call**
     evidence (LDISC/RDISC/LCLIP/RCLIP, breakpoint, TSD, polyA). The clip realignment is redundant
     with the scan, and the whole BAM re-read is redundant with the scan pass.

---

## 2. Target architecture

```
                 ┌─────────────────────────────────────────────────────────────┐
  BAM/CRAM  ───► │ xtea-scan  (ONE pass)                                        │
                 │   clip + disc detection                                      │
                 │   realign clip parts AND disc mate-reads to all 3 consensi   │
                 │   emit per-LOCUS × per-TYPE evidence table                   │
                 └─────────────────────────────────────────────────────────────┘
                                        │  one locus table (per-type support side by side)
                                        ▼
                 ┌─────────────────────────────────────────────────────────────┐
                 │ assign TYPE per locus (best consensus support; merger rules) │
                 │   → one call per locus (or drop)                             │
                 └─────────────────────────────────────────────────────────────┘
                                        │
                                        ▼
                 ┌─────────────────────────────────────────────────────────────┐
                 │ type-BRANCHED analysis:  transduction (L1/SVA only), sibling │
                 └─────────────────────────────────────────────────────────────┘
                                        │
                                        ▼
                 ┌─────────────────────────────────────────────────────────────┐
                 │ type-AGNOSTIC finish (ONCE):                                 │
                 │   depth/coverage filter → gene annotation → genotype → VCF   │
                 └─────────────────────────────────────────────────────────────┘
                                        │
                                        ▼
                                   ONE sorted VCF   (no per-type VCFs, no merge step)
```

**Core data model — the locus table** (one row per candidate locus):

```
Locus {
  chrom, pos,
  clip_L, clip_R, col5,            # from scan (detection)
  disc_L, disc_R,                  # from scan disc stream
  support: { L1:(clip_aln, disc_aln, breakpoint, tsd, polyA, …),
             Alu:(…),
             SVA:(…) },            # per-type call evidence  (fold cns realignment here)
  assigned_type,                   # decided from support
  transduction, gene, genotype …   # filled downstream
}
```

## 3. The four moves (as approved)

- **A — Type as an attribute of a locus, not a pipeline.** One locus list; per-type support side by
  side; assign type up front (merger's quality/priority rules); downstream runs once per locus → one
  VCF, merge step deleted.
- **B — Fold ALL realignment into the scan; cns_filter becomes a pure in-memory filter.** Scan already
  holds every read during its pass — extend it to also realign **disc** mate-reads (to consensus+flank)
  and emit the per-locus×per-type evidence table. `cns_filter` then never touches the BAM or bwa; it
  just thresholds the table. Kills the redundant clip realignment AND the 3× BAM re-read.
- **C — Collapse the type-agnostic tail.** Gene annotation, genotyping, VCF emission are type-agnostic
  → run once over the unified call set, not 3×.
- **D — Keep type-specific only where biology demands it.** Transduction (L1/SVA 3′ transduction; Alu
  none) becomes a conditional branch in the per-locus loop, not a separate pipeline.

## 4. Sequencing (low → high risk; each stays harness-gated)

1. **C (tail unification)** — lowest risk: gene/genotype/VCF once instead of 3×; still consumes the
   three per-type cns outputs. Proves the unified call set + single-VCF path with no algorithmic change.
2. **B (fold cns realignment into scan)** — remove cns_filter's BAM re-read and bwa re-run; cns becomes
   a filter over the scan's evidence table. Biggest compute win; medium risk (must reproduce cns's
   LDISC/RDISC/LCLIP/RCLIP + breakpoint/TSD/polyA exactly).
3. **A (full locus-unification + drop merger)** — key by locus, resolve type up front. Highest value,
   highest risk: changes *when/how* type is resolved → diverges from golden more than a mechanical
   change; the concordance harness must gate it hard.

## 5. Risks / watch-items

- **Transduction** is the gnarliest code to fold into a unified loop (L1 transduction, SVA
  transduction, orphan/sibling). Map it carefully before touching.
- **Type-resolution divergence:** A changes the order of "compute vs resolve type," so borderline
  multi-type loci may flip vs the golden. Gate with the concordance harness (100% golden recall on WGS;
  watch PPV).
- **cns is cross-candidate-sensitive** (known): small changes to the candidate set move cns output.
  Reproduce its per-candidate math exactly in B before optimizing.
- Keep `--keep-intermediates` / a `--debug-bundle` so the (deleted-by-default) per-type intermediates
  are recoverable for debugging.

---

## 6. First task — the cns_filter contract (linchpin for A and B)

Map exactly what `cns_filter` (`-N`, `x_clip_disc_filter.py`) computes **that the scan does not yet**.
**DONE 2026-07-18.** Full contract: **`docs/cns_filter_contract.md`** (entry chain, per-candidate
read collection, realignment, 32-col output, thresholds, type-specificity, BAM re-reads — all cited
to file:line). Distilled findings + implications:

### 6.1 What cns_filter does (entry: `x_TEA_main.py -N` → `call_MEIs_consensus`, `x_clip_disc_filter.py:1527`)

- **Per-candidate read collection** (`collect_clipped_disc_reads_by_region`, :81): ONE pysam
  `.fetch(chrm, pos-extnd, pos+extnd)` per (BAM, candidate); window `extnd = 400`, widened to
  `max_is = mean_is + 3·std` when larger. Extracts clip parts → FASTQ and disc **mate reads** (by
  name, `extract_mate_reads_by_name`) → FASTA. **This is the redundant BAM re-read** (~1 fetch per
  candidate per BAM; coverage calc may add one more).
- **Realignment (bwa mem) — to BARE CONSENSUS, not consensus+flank:**
  - clip parts → consensus (polyA cohort `-T7 -k7 -c70`, non-polyA `-T9 -k9`), and
  - **full disc mate reads → consensus** (`-T -1`).  ← the scan does NOT do this yet (it realigns
    only clip parts, and to `copies_with_flank`). **This disc-read→consensus realignment is the main
    new work move B must absorb.**
- **Evidence per call** (`_extract_detailed_info_TEI`, :1089): LCLIP/RCLIP + LDISC/RDISC aligned
  counts, clip/disc cluster spans on the consensus, peak clip positions, **TSD = |rpeak−lpeak|**,
  refined breakpoint, polyA counts, orientation (RC vs not), left/right coverage. → the **32-col**
  `candidate_disc_filtered_cns.txt` (columns fully listed in the contract §6).
- **Thresholds:** nclip (both-side total ≥ cutoff; one-side ≥ cutoff/2+1), ndisc (≥ cutoff; one-side
  ≥ cutoff/2, disc cluster ratio ≥ 0.45), clip–disc distance ≤ max_is/2, mapped-fraction ≥ 0.65,
  clip cluster ratio ≥ 0.4 within a 15 bp window, depth-ratio ≤ 0.35, coverage cap 3× average, polyA
  required unless the clip reaches the consensus boundary (cns_len−150).

### 6.2 Implications for the unified design

- **cns_filter is ~type-agnostic.** Read collection, realignment prep, clustering, TSD, breakpoint,
  polyA, thresholds are all identical across L1/Alu/SVA. The ONLY per-type variation is (a) the
  **consensus reference** — and the scan already realigns to all three consensi — and (b) SVA's
  **relaxed consistency** under `IS_CALL_SVA` (:1353) + a one-side-flank high-confidence tweak. → Move
  **A** (one locus, per-type evidence side-by-side, assign type) is well-supported; the per-locus
  processor is one code path with a consensus slot and one `if sva` relaxation.
- **Move B scope is now concrete.** The scan already reads every read and realigns clip parts. To make
  cns a pure in-memory filter, the scan must additionally, per candidate: (1) realign **disc mate
  reads** to each consensus, (2) cluster clip/disc on the consensus + compute peaks/TSD/breakpoint/
  polyA/orientation, (3) compute left/right coverage. All three use reads the scan already has in the
  pass → no extra BAM I/O.
- **Reconcile a realignment-reference mismatch (watch-item):** the scan's *detection* clip realign
  uses `copies_with_flank`; cns's *call* realign uses **bare consensus**. Move B must realign the
  call-evidence reads to bare consensus (or prove flank-vs-consensus doesn't move the 32-col output),
  or keep the two realignment purposes distinct.
- **Coverage caching gap** (contract "Unclear"): `calc_coverage_of_two_regions` may re-open the BAM
  per candidate — folding coverage into the scan pass removes that too.
