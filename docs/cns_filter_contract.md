# xTea `cns_filter` Step Contract

## 1. Entry Chain

**Handler in x_TEA_main.py** (line 477–534):
- `-N` / `--filter_csn` option triggers `elif options.filter_csn:` block
- Line 531: Instantiates `XClipDiscFilter(sf_bam_list, s_working_folder, n_jobs, sf_ref)`
- Line 532–534: Calls **`XClipDiscFilter.call_MEIs_consensus()`** with these parameters:
  - `sf_candidate_list` — input candidates
  - `sf_raw_disc` — raw discordant reads file
  - `iextnd` = 400 (line 487)
  - `bin_size` = 50,000,000 (line 488)
  - `sf_cns` — consensus reference (line 489, `-r` / `--reference`)
  - `sf_flank` — flanking region (line 495, `--fflank`)
  - `i_flank_lenth` — flank length (line 496)
  - `bmapped_cutoff` = `MIN_CLIP_MAPPED_RATIO` = 0.65 (line 490; global_values.py:118)
  - `i_concord_dist` = 550, auto-adjusted to `max_is` if needed (lines 492, 509–513)
  - `f_concord_ratio` = 0.45 (line 493)
  - `n_clip_cutoff` (line 521, auto-computed or user-set via `options.cliprep`)
  - `n_disc_cutoff` (line 522, auto-computed or user-set via `options.ndisc`)
  - `sf_output` — output file (line 494)

**Key method chain in x_clip_disc_filter.py:**
- **`call_MEIs_consensus()`** (line 1527) orchestrates the full pipeline

---

## 2. Per-Candidate Read Collection

**Window extension values:**
- Default `iextnd` = **400 bp** (x_TEA_main.py:487)
- Adjusted upward if `max_is` (mean_is + 3×std_var) exceeds it (x_TEA_main.py:509–513)
- Typical range: [400, max_is where max_is ≈ mean_insert_size + 3×std]

**Read collection method:**
1. **`collect_clipped_disc_reads()`** (x_clip_disc_filter.py:3306) iterates over each BAM in `sf_bam_list`
2. For each BAM, creates `XClipDisc` instance (line 3317)
3. Calls **`XClipDisc.collect_clipped_disc_reads_of_given_list()`** (line 294):
   - Uses **multiprocessing Pool** with `n_jobs` workers (line 302)
   - Each worker calls **`collect_clipped_disc_reads_by_region()`** (line 81) per candidate
   - **Opens BAM via `pysam.AlignmentFile()`** (x_clip_disc_filter.py:104)
   - **Fetches reads using `.fetch(chrm_in_bam, start_pos, end_pos)`** (line 109)
   - Writes clipped reads to FASTQ, discordant read anchors to separate files
   - Returns low-mapq clip counts

**BAM access per candidate:** ONE fetch per (BAM, candidate) pair using pysam range-based fetch.

---

## 3. Realignment

**Clipped reads realignment:**
- **Reference:** Consensus only (`sf_rep_cns` / `sf_cns`)
- **Tool:** BWA mem (bwa_align.py:31–45)
- **Method:** `realign_clipped_read_with_polyA()`
  - Splits reads by length: clipped parts <9 bp → polyA cohort; ≥9 bp → standard
  - Realigns non-polyA reads to consensus with `-T 9 -k 9` (line 22)
  - Realigns polyA reads to consensus with `-T 7 -k 7 -c 70` (line 68–71)
  - Merges both SAM outputs (lines 48–60)

**Discordant mate reads realignment:**
- **Reference:** Consensus only (`sf_rep_cns`)
- **Tool:** BWA mem (bwa_align.py:120–123)
- **Method:** `realign_disc_reads()`
  - Standard BWA mem with `-T -1` (no threshold cutoff)
  - Input is full mate-pair sequences (not clipped parts)

**Key distinction:**
- UPSTREAM scan (clip and disc steps) collects clipped-part alignments to consensus
- **cns_filter ADDITIONALLY realigns:**
  - **CLIPPED PARTS** (already extracted) to consensus → detect breakpoint in consensus
  - **FULL DISCORDANT MATE READS** to consensus → detect where mate pairs land in repeat

**NO realignment to consensus+flank happens here** (flank is passed but not used for realignment; used only for annotation).

---

## 4. Computed Evidence Per Call

**Per-candidate evidence compiled in `_extract_detailed_info_TEI()`** (x_clip_disc_filter.py:1089–1226):

**Clip evidence:**
- `n_total_lclip` — total left-clipped reads (line 1217)
- `n_total_rclip` — total right-clipped reads (line 1217)
- `nlpolyA` — left polyA count (line 1217, from realignment analysis)
- `nrpolyA` — right polyA count (line 1217)
- `lc_cns_start, lc_cns_end` — left clip cluster start/end on consensus (line 1219)
- `rc_cns_start, rc_cns_end` — right clip cluster start/end on consensus (line 1219)
- `lpeak_pos, rpeak_pos` — left/right peak clip positions on consensus (lines 1106–1112)

**Discordant evidence:**
- `n_total_ldisc` — total left-discordant reads (line 1217)
- `n_total_rdisc` — total right-discordant reads (line 1217)
- `ldisc_cns_start, ldisc_cns_end` — left disc cluster start/end on consensus (line 1219)
- `rdisc_cns_start, rdisc_cns_end` — right disc cluster start/end on consensus (line 1219)
- `rc_rcd = (n_l_rc, n_l_not_rc, n_r_rc, n_r_not_rc)` — orientation counts (line 1220)
- `anchor_dir_rcd` — additional orientation tracking (line 1220)

**Breakpoint refinement:**
- `refined_pos` — best clip position or fallback to disc start (lines 1192–1194)
- `TSD` = abs(rpeak_pos − lpeak_pos) if both peaks exist (lines 1188–1190)

**Coverage:**
- `flcov` — left local coverage (from `calc_coverage_of_two_regions()`, x_clip_disc_filter.py:1679)
- `frcov` — right local coverage (same)

**Extraction source:**
- `parse_clip_realignment_consensus()` (x_clip_disc_filter.py:2946) → clip positions, polyA flags
- `parse_disc_algnmt_consensus()` (x_clip_disc_filter.py:3209) → disc positions, orientation, sample info
- `check_clip_alone_consistency2()` (x_clip_disc_filter.py:2330) → clip clustering
- `check_seprt_disc_consistency()` (x_clip_disc_filter.py:2512) → disc clustering

---

## 5. Filter Thresholds

**Clip-side cutoffs:**
- Minimum clipped reads: `nclip_cutoff` (typically 3–5, set via `-cliprep` or auto-computed)
  - If both sides present: total must reach `nclip_cutoff` (line 2349)
  - If one side only: minimum `(nclip_cutoff/2)+1` (line 2363)
- Low-mapq clip ratio: max 65% (global_values.py:36, `MAX_LOWQ_CLIP_RATIO=0.65`)

**Disc-side cutoffs:**
- Minimum discordant reads per side: `ndisc_cutoff` (typically 4–6, set via `-ndisc` or auto-computed)
  - If both sides cluster: total must ≥ `ndisc_cutoff` (line 2598)
  - If one side only (transduction): ≥ `ndisc_cutoff/2` (lines 2601–2604)
- Half-cutoff used for partial events: `ndisc_cutoff/2` (line 1625)

**Positional consistency:**
- Clip-disc distance: both clip and disc must be within insert-size window (x_clip_disc_filter.py:1231–1241)
  - `idist = max_is / 2 + 1` (line 1232)
- Disc cluster ratio on consensus: ≥ `f_concord_ratio` = 0.45 (x_clip_disc_filter.py:2591)
- Disc clustering distance on consensus: `icorcord` = 550, auto-adjusted to max_is if needed

**Clip clustering:**
- Search window: `CLIP_SEARCH_WINDOW` = 15 bp (global_values.py:112)
- Consistency ratio: `CLIP_CONSISTENT_RATIO` = 0.4 (global_values.py:113) — ≥40% of reads must cluster

**Coverage-based filtering:**
- Maximum coverage threshold: `MAX_COV_TIMES * ave_coverage` = 3 × ave_cov (global_values.py:114)
  - Sites exceed this are tagged as high-coverage islands, can be filtered (x_clip_disc_filter.py:1690)
- Depth ratio: 0.35 (x_clip_disc_filter.py:1716) — left–right depth can differ by ≤35% for validity

**Consensus boundary:**
- Must reach consensus end: min_l1_end = cns_lth − 150 (x_clip_disc_filter.py:1715)

**PolyA requirement (RNA-mediated, default ON):**
- If both sides clipped but zero polyA: filtered unless hits consensus boundary (lines 1381–1387)

---

## 6. Output Columns

**`dump_TEI_info()` format** (x_clip_disc_filter.py:735–778, also `dump_TEI_info_with_ori_pos()` at 844):

| Column | Header | Meaning |
|--------|--------|---------|
| 1 | #chrm | Chromosome |
| 2 | refined-pos | Refined insertion position (best clip or disc start) |
| 3 | lclip-pos | Left peak clip position on consensus |
| 4 | rclip-pos | Right peak clip position on consensus |
| 5 | TSD | Target-site duplication length (rclip − lclip) |
| 6 | nalclip | Total left-clipped reads |
| 7 | narclip | Total right-clipped reads |
| 8 | naldisc | Total left-discordant reads |
| 9 | nardisc | Total right-discordant reads |
| 10 | nalpolyA | Total left-polyA reads |
| 11 | narpolyA | Total right-polyA reads |
| 12 | lcov | Left local read depth |
| 13 | rcov | Right local read depth |
| 14 | nlclip | Number of (unique) left-clipped reads in peak cluster |
| 15 | nrclip | Number of (unique) right-clipped reads in peak cluster |
| 16 | nldisc | Number of left-discordant reads in cluster |
| 17 | nrdisc | Number of right-discordant reads in cluster |
| 18 | nlpolyA | Left-polyA read count (from cluster) |
| 19 | nrpolyA | Right-polyA read count (from cluster) |
| 20 | lclip-cns-start:end | Left-clip cluster region on consensus |
| 21 | rclip-cns-start:end | Right-clip cluster region on consensus |
| 22 | ldisc-cns-start:end | Left-disc cluster region on consensus |
| 23 | rdisc-cns-start:end | Right-disc cluster region on consensus |
| 24 | Transduction-info | not_transduction \| transduction type label |
| 25 | ntsd-clip | Transduction clipped reads count |
| 26 | ntsd-disc | Transduction discordant reads count |
| 27 | ntsd-polyA | Transduction polyA count |
| 28–31 | ldisc-same, ldisc-diff, rdisc-same, rdisc-diff | Orientation counts (RC vs non-RC) for left and right disc |
| 32 | 3mer-inversion | NOT-5prime-inversion \| 5-prime-inversion label |
| (additional columns in final output) | — | Sample info, confidence class, genotype features, insertion length |

**Output file:** `candidate_disc_filtered_cns.txt` (specified by `-o` / `--output`)

**Secondary output:** `candidate_disc_filtered_cns.txt.high_confident` (high-confidence subset from `_second_stage_classify_dump()`, x_clip_disc_filter.py:843–1040)

---

## 7. Type-Specificity

**Type-agnostic operations (apply to all L1, Alu, SVA, etc.):**
- Read collection from BAM (window, fetch logic) — **same for all types**
- Realignment to consensus — **same reference preparation for all types**
- Clip/disc clustering detection — **same thresholds**
- TSD/breakpoint refinement — **same algorithm**

**Type-dependent operations:**

| Operation | L1 | Alu | SVA | Mechanism |
|-----------|----|----|-----|-----------|
| **Consistency check** | Strict: both clip+disc must be distance-consistent | Strict (same) | **Relaxed if `IS_CALL_SVA==True`** | global_values.py:99–102; x_clip_disc_filter.py:1353 |
| **PolyA requirement** | Required unless consensus boundary hit | Required (same) | Required (same) | global_values.py:7–10 (GLOBAL_RNA_MEDIATED); x_clip_disc_filter.py:1381 |
| **Annotation source** | Per-type `.out` file (`-a` / `--annotation`) | Per-type (same) | Per-type, may have SVA_ANNOTATION_EXTND=200 (global_values.py:103) | x_TEA_main.py:491 sets `sf_annotation` |
| **Flank reference** | Provided separately (`--fflank`) | Provided separately | Provided separately | **Same usage: only for annotation, NOT realignment** |
| **Per-type consensus** | `-r` / `--reference` specifies consensus | Same parameter | Same parameter | **Different consensus per run, not per-type within run** |

**Flag controls:**
- `--sva` / `IS_CALL_SVA` (line 100–102):
  - If TRUE: SVA runs relax consistency check (line 1353)
  - Also changes high-confidence classification for one-side flanking cases (lines 1025–1026)

**Practical implication:**
- A **single run** filters one type at a time (separate invocations for L1, Alu, SVA)
- Each run receives:
  - Per-type candidates (from prior `-C` clip step on that type)
  - Per-type consensus (FASTA, aggregated or per-copy)
  - Per-type annotation (`.out` repeat mask)
  - Per-type flank regions (if applicable)
- **Thresholds (nclip, ndisc, etc.) are shared** across types within a run but can differ between runs

---

## 8. BAM Re-reads: Full Enumeration

**Every location the step opens/fetches the BAM:**

1. **Line 3314–3324 (`collect_clipped_disc_reads()` wrapper):**
   - Opens BAM list file, iterates each BAM
   - For each BAM, calls `collect_clipped_disc_reads_of_given_list()` (line 3322)

2. **Line 303 (multiprocessed collection):**
   - **`unwrap_self_collect_clip_disc_reads()`** dispatches to `collect_clipped_disc_reads_by_region()`
   - **Line 104: `pysam.AlignmentFile(sf_bam, "rb")`** — opens BAM
   - **Line 109: `samfile.fetch(chrm_in_bam, start_pos, end_pos)`** — range fetch per candidate
   - Called once per (BAM, candidate) pair in multiprocessing pool

3. **NO additional BAM reads for realignment:**
   - Clipped reads and disc reads are EXTRACTED and written to temp FASTQ/FASTA (lines 99–102, 101–102)
   - Realignment uses only the extracted sequences, not BAM re-access
   - Disc mate pairs are RETRIEVED using `extract_mate_reads_by_name()` on read names saved from step 2 (line 338)

4. **Coverage calculation (optional, if not pre-computed):**
   - `calc_coverage_of_two_regions()` (x_coverage.py) may re-open BAM per region if coverage not cached
   - Called at line 1679 for left/right depth per candidate
   - **Likely** re-opens BAM for each site in a search window (default COV_SEARCH_WINDOW=1000 bp)

**Summary:** 
- **Primary BAM access:** ONE pysam fetch per (BAM, candidate) pair (line 109)
- **Secondary (coverage):** Potentially ONE additional re-open per BAM per candidate at line 1679, depending on caching
- **Total:** Expect 1–2 BAM file handles per (BAM, candidate) pair
- **Scale:** With ~5000 candidates and 10 BAMs, ~50,000–100,000 BAM fetch operations

---

## Unclear/Gaps

- **Coverage caching:** Whether `calc_coverage_of_two_regions()` caches across candidates or re-opens BAM each time (affects BAM re-read count)
- **Flank realignment:** Whether flank+consensus realignment happens in a secondary pass not shown in `call_MEIs_consensus()` (comments in bwa_align.py:149–150 suggest possible two-stage, but not invoked in cns_filter)
- **Sample-specific thresholds:** When CHECK_BY_SAMPLE==True (global_values.py:94), code enforces per-sample disc support (line 1629), but default is False; affects which candidates pass

---

## Design Doc Refs

- **Entry point:** x_TEA_main.py:477–534 (handler block)
- **Main orchestrator:** x_clip_disc_filter.py:1527 (`call_MEIs_consensus`)
- **Read collection:** x_clip_disc_filter.py:294, 3306 (entry functions)
- **Realignment:** bwa_align.py:31, 120 (BWA wrapper methods)
- **Evidence extraction:** x_clip_disc_filter.py:1089, 2946, 3209
- **Filtering logic:** x_clip_disc_filter.py:1261 (`filter_by_clip_disc_polyA_consistency`)
- **Output:** x_clip_disc_filter.py:735, 844 (dump methods)
- **Constants:** global_values.py:74–159 (cns_filter-related)
