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

1. **C (tail unification)** — lowest risk. **C.1 DONE (commit 8738fb6):** emit one merged VCF via
   master's `SampleVCFMerger` after the type loop; additive, per-type VCFs untouched; golden merged
   frozen (chr1 32, WGS 1520). **C.2 (collapse gene/genotype compute) DEPRIORITIZED — low value:** the
   tail runs on each type's *filtered* call set (~27 Alu + 5 SVA + 0 L1 on chr1; ~1,538 at WGS), which
   are **partitioned** across types, NOT triplicated. Running gene/genotype "once over the merged set"
   saves only per-invocation overhead (3 model-loads/process-starts/BAM-opens → 1), not site-work. So
   the single-VCF deliverable (C.1) is the real win of move C; C.2 is deferred as cleanup.
   **→ The 3× redundancy lives in cns_filter (move B), not the tail:** cns_filter re-reads the BAM +
   realigns per *candidate*, and the candidate sets are the ~46k clip positions that are 99.98% shared
   across types — so its heavy per-candidate work runs ≈3× at the same loci. **Move B is where the
   compute win is; go there next.**
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

## 7. Move B.1 — the seam, and the concrete approach (proceeding 2026-07-18)

**The seam in `call_MEIs_consensus` (`x_clip_disc_filter.py:1527`) is clean and exact:**
- **Collect+realign (the redundant BAM read + bwa), :1536–1564:** `collect_clipped_disc_reads` fetches
  per-candidate clip parts → `candidate_sites_all_clip.fq` and disc mate reads →
  `candidate_sites_all_disc.fa` **from the BAM**, then `realign_clipped_read_with_polyA` and
  `realign_disc_reads` bwa-align both to **`sf_rep_cns`** → `temp_clip.sam`, `temp_disc.sam`.
- **Parse+cluster+filter (the math), :1570+:** consumes ONLY `temp_clip.sam` + `temp_disc.sam`
  (`parse_clip_realignment_consensus` :1571, `parse_disc_algnmt_consensus` :1619) + the candidate list
  + coverage. So the two SAMs are the entire interface between the I/O and the math.

**Reconciliation question — ANSWERED:** the call-evidence realign target is **bare consensus**
(`sf_rep_cns` = `consensus/<TYPE>.fa`, the `-r` arg), NOT the scan's detection reference
(`copies_with_flank`). Confirmed: the golden cns step ran `-r consensus/SVA.fa`; the math clusters in
consensus coordinates (`lc_cns_start` …). → B.1's scan must realign the call-evidence reads to the
**bare consensus**, as a target distinct from (in addition to) the detection realign.

**Ordering constraint (key design fact):** candidates are only known AFTER the scan (the disc adapter
produces `candidate_list_from_disc.txt` post-scan). So the scan cannot collect "per-candidate" reads at
scan time. Resolution: the scan realigns **all** detected clip parts + disc reads to the consensi
(candidate-agnostic, done once, concurrently — it already does this for clip parts vs
copies_with_flank), and a thin adapter filters those alignments to the candidate positions to build
each type's `temp_clip.sam`/`temp_disc.sam`. No BAM re-read; realign happens once in the scan.

**Two BAM re-reads to remove** (both in cns_filter): (1) `collect_clipped_disc_reads` (reads for
realign), (2) `calc_coverage_of_two_regions` (~:1679, per-candidate depth). B.1 targets (1) first;
coverage (2) folds in after.

**B.1 implementation sequence:**
- **B.1.i** — scan additionally realigns (a) clip parts and (b) **disc mate reads** to the **bare
  consensus** per family, emitting alignments in the exact SAM shape `parse_clip_realignment_consensus`
  / `parse_disc_algnmt_consensus` expect (read-name metadata `@chrm~pos~L/R~is_rc~ins_pos~ncnt~sample`).
- **B.1.ii** — adapter filters the scan's consensus alignments to each type's candidate positions →
  `temp_clip.sam` / `temp_disc.sam`.
- **B.1.iii** — cns_filter gains a "SAMs provided" path that skips collect+realign (:1536–1564) and
  enters at parse (:1570). Gate: `candidate_disc_filtered_cns.txt` matches the per-type golden.
- **B.1.iv** — fold per-candidate coverage into the scan pass; remove the second BAM re-read.

**SAM contract — MAPPED (full: `docs/cns_realign_sam_contract.md`).** Confirmed against the golden's
preserved `SVA/tmp/cns/temp_clip.sam` (11,080 rows) + `temp_disc.sam` (260 rows):
- **Reference = `consensus/<TYPE>.fa`, a MULTI-SUBFAMILY set** (SVA 88 seqs incl. SVA_A…SVA_s39, ALU
  47, LINE1 4). RNAME = the subfamily the read best matched; POS is within that subfamily. (Earlier
  "bare single consensus" was wrong — it's per-subfamily; supplementary hits appear, FLAG 2048.)
- **temp_clip.sam** qname (7 fields, `~`-sep): `chrm~map_pos(0-based)~L|R~is_rc(0/1)~ins_pos~n_cnt~
  sample_id`. SEQ = the clipped part (left=cigar[0] bases, right=cigar[-1]), no RC. bwa: two cohorts —
  non-polyA (≥9 bp) `bwa mem -T9 -k9`, polyA (<9 bp) `-T7 -k7 -c70`, merged. Parser gate: mapped ≥
  **0.65**·(len−dels).
- **temp_disc.sam** qname (10 fields, parser uses NEGATIVE indices [-8..-1] — field count is load-
  bearing): `read_name~is_first~b_lclip~b_rclip~anchor_rc~mate_rc~anchor_map_pos~ins_chrm~ins_pos~
  sample_id`. SEQ = the mate read (full, no RC). bwa: `bwa mem` (standard). Parser: `pos_in_consensus
  = ref_start + len(seq)/2`; gate mapped ≥ **0.70**. 
- **sample_id** is assigned by BAM order at collection time (0 for our single BAM) — inject it, it's
  not in the alignment.

**Ground-truth validation strategy (de-risks the whole seam):** the golden preserved the temp SAMs, so
each B.1 sub-step diffs against them. **Refined sequence — do B.1.iii FIRST (cheapest, proves the
seam):**
- **B.1.iii first** — add a "SAMs-provided" bypass to `call_MEIs_consensus` (skip collect+realign
  :1536–1564 when temp_clip.sam/temp_disc.sam already exist), feed it the golden's own temp SAMs, and
  confirm it reproduces the golden `candidate_disc_filtered_cns.txt`. Proves the seam is exactly the two
  SAMs before any scan-side work. Additive/opt-in — default path unchanged.
- **B.1.i** — scan realigns clip parts + disc reads to `consensus/<TYPE>.fa`, emits the two SAM
  formats; diff vs golden temp SAMs (or compare parsed evidence).
- **B.1.ii** — adapter filters scan alignments to candidate positions.
- **B.1.iv** — fold per-candidate coverage into the scan.
