# xTea2 — Clip/Disc step contract specification

**Purpose.** This is the **detection-semantics contract** for reimplementing xTea's **clip** (`-C`)
and **discordant / disc** (`-D`) read detection. It captures exactly *which reads and which candidate
sites master finds*, and the numeric thresholds that decide that — because **that is what determines
the final results**. Every value here is extracted verbatim from the `master`/`xTea2` source at commit
`46520ab` (`/home/cos689/data1/corinne/claude_sandbox/xTea/xtea/`). Line numbers refer to that tree.

> ### What is and isn't a constraint (per user direction, 2026-07-17)
> - **CONSTRAINT — results:** the final VCF call-set must match master (concordance gate: ±25 bp +
>   same MEI type), and at the read level we must **find the same clipped/discordant reads** at each
>   candidate site. The unit tests `test_clip_reads.py` / `test_disc_reads.py` assert *set equality of
>   the reads/sites found* — a **semantic** invariant, independent of file layout.
> - **NOT a constraint — intermediate file formats/layout.** We are explicitly free to change the
>   result layout, **merge clip+disc into one pass**, integrate with downstream, drop the
>   per-chromosome FASTQ/pos file explosion, and pick whatever data model is best for
>   speed/memory/disk/file-count. Master's current formats are recorded below **for reference only**
>   (to build the golden fixtures and to understand the semantics) — they are NOT targets to reproduce.
> - **Secondary goal, nearly as important:** the "best" implementation — optimized for speed, memory,
>   disk space, and number of files — using software-engineering best practices. See §3 (Architecture).
> - **Ongoing:** flag logical bugs and remove dead code as we go. See §4 (Bug / dead-code log).

> Scope: **short reads only** (Illumina). All `LRD_*` (long-read) constants are listed for
> completeness but are OUT OF SCOPE for xTea2. 10X/barcode paths are dropped.

---

## 0. Constants (from `global_values.py`)

### MAPQ / mapping-quality cutoffs
| Name | Value | Line | Note |
|---|---|---|---|
| `MINIMUM_CLIP_MAPQ` | **12** | 31 | min MAPQ for a read to count as a clip (comment says 10, code is 12) |
| `CLIP_LOW_MAPQ` | 5 | 35 | clip read with MAPQ ≤ 5 → counted as background |
| `MINIMUM_DISC_MAPQ` | **20** | 39 | min MAPQ when collecting discordant reads |
| `MAX_BKGRND_LOW_MAPQ` | 10 | 44 | below this MAPQ → background low-mapq read |
| `MAX_BKGRND_LOW_MAPQ_RATIO` | 0.1 | 45 | |
| `UNIQUE_MAP_CUTOFF` | 60 | 213 | "unique map" threshold |

### Clip detection: lengths, counts, windows
| Name | Value | Line | Note |
|---|---|---|---|
| `INITIAL_MIN_CLIP_CUTOFF_ILLUMINA` | 2 | 19 | |
| `INITIAL_MIN_CLIP_CUTOFF` | 2 | 21 | min total clips at an exact position to record it |
| `MIN_CLIP_FOR_CANDIDATE` | 2 | 75 | filters (1,1)/(1,0)/(0,1) L/R patterns |
| `NEARBY_REGION` | **50** | 76 | site-merge / nearby window |
| `NEARBY_CLIP` | 50 | 111 | |
| `CLIP_SEARCH_WINDOW` | 15 | 112 | |
| `CLIP_FREQ` | 10 | 77 | |
| `TRIM_CLIP_FREQ` | 2 | 78 | |
| `PEAK_WINDOW` | 100 | 79 | |
| `TWO_SIDE_CLIP_MIN_LEN` | 8 | 120 | two-side clip needs clip len > 7 |
| `FLAG_LEFT_CLIP` | `"L"` | 55 | left-clip flag char in FASTQ header |
| `FLAG_RIGHT_CLIP` | `"R"` | 56 | right-clip flag char |

### Base-quality filter
| Name | Value | Line | Note |
|---|---|---|---|
| `CLIP_PHRED_SCORE_CUTOFF` | **15** | 29 | per-base Phred cutoff |
| `MAX_LOWQ_CLIP_RATIO` | 0.65 | 36 | if > this frac of clips are low-q → FP |
| `CK_POLYA_CLIP_WIN` | 25 | 61 | |
| `CK_POLYA_SEQ_MAX` | 20 | 62 | check at most 20 bases for polyA (comment wrongly says 10) |

### Poly-A / poly-T
| Name | Value | Line | Note |
|---|---|---|---|
| `MINIMUM_POLYA_CLIP` | 7 | 27 | min clip len to apply polyA test |
| `POLYA_RATIO` | **0.4** | 63 | ≥ 40% A or T |
| `N_MIN_A_T` | 5 | 110 | min consecutive A/T |
| `MAX_POLYA_RATIO` | 0.85 | 153 | |
| `MIN_POLYA_CLIP_RATIO` | 0.1 | 187 | (comment wrongly says 20%) |
| `DOMINANT_POLYA_T_MIN_RATIO` | 0.75 | 188 | |

### Discordant criteria
| Name | Value | Line | Note |
|---|---|---|---|
| `DISC_THRESHOLD` | 2000 | 106 | |
| `MAX_NORMAL_INSERT_SIZE` | 2000 | 304 | |
| `DFT_IS` | 550 | 292 | default insert size |
| `MIN_ORPHAN_DISC_PAIR_INSERT_SIZE` | 100000 | 185 | 0.1 Mb |
| `TWO_SITES_CLOSE_DIST` | 100 | 178 | |
| `TSD_CUTOFF` | 100 | 107 | |
| `MIN_RAW_DISC_CLUSTER_RATIO` | 0.3 | 68 | |
| `MIN_CLUSTER_RC_RATIO` | 0.6 | 186 | ≥ this frac of cluster reads reverse-comp |

### BWA / read length
| Name | Value | Line | Note |
|---|---|---|---|
| `BWA_REALIGN_CUTOFF` | 9 | 26 | max clip len handled by the polyA short path |
| `BWA_HALF_READ_MIN_SCORE` | 45 | 289 | |
| `READ_LENGTH` | 100 | 297 | assumed read length |

### Coverage / filtering
| Name | Value | Line | Note |
|---|---|---|---|
| `MIN_CLIP_MAPPED_RATIO` | 0.65 | 118 | |
| `MIN_DISC_MAPPED_RATIO` | 0.7 | 119 | |
| `MAX_COV_TIMES` | 3 | 114 | |
| `AVE_COVERAGE` | 30 | 282 | |
| `max_cov_cutoff` | **15 × --cov** (default 600) | x_TEA_main.py | computed at runtime: `int(15*options.cov)` |

### Tool paths (resolved from `$PATH`)
`BWA_PATH="bwa"` (47), `SAMTOOLS_PATH="samtools"` (48), `BCFTOOLS_PATH="bcftools"` (233),
`TABIX_PATH="tabix"` (234).

### Output suffixes / delimiters (exact filenames the C tool must emit)
| Name | Value | Line |
|---|---|---|
| `SEPERATOR` / `S_DELIM` | `"~"` | 57 / 192 |
| `CLIP_FQ_SUFFIX` | `".clipped.fq"` | 49 |
| `CLIP_POS_SUFFIX` | `".clip_pos"` | 52 |
| `CLIP_RE_ALIGN_POS_SUFFIX` | `".clip_realign_pos"` | 53 |
| `ALL_DISC_SUFFIX` | `".initial.all.disc"` | 58 | **written but never read — DELETE (free disk win)** |
| `DISC_SUFFIX` | `".discord_pos.txt"` | 84 |
| `DISC_POS_SUFFIX` | `".disc_pos"` | 123 |
| `DISC_NAME_SUFFIX` | `".disc_names"` | 122 |
| `RAW_DISC_TMP_SUFFIX` | `".clip_sites_raw_disc.txt"` | 90 |
| `CLIP_FOLDER` / `DISC_FOLDER` | `"clip"` / `"disc"` | 82 / 83 |

### Control flags
`GLOBAL_RNA_MEDIATED=True` (7) — non-RNA-mediated cases don't require polyA support at filtering.

> **Known source quirks to preserve (do NOT "fix"):** the `SEPERATOR` misspelling; several inline
> comments contradict the code (line 31 says MAPQ 10 but is 12; line 62 says 10 bases but is 20;
> line 187 says 20% but is 0.1). The **code value** is authoritative — reproduce the code, not the comment.

---

## 1. CLIP step (`-C`)

### 1.1 Entry / dispatch (`x_TEA_main.py` 372–422)
`options.clip` → builds `TE_Multi_Locator(sf_bam_list, s_working_folder, n_jobs, sf_ref)` and calls
`call_TEI_candidate_sites_from_multiple_alignmts(sf_annotation, sf_rep_cns, sf_rep, b_se,
cutoff_left_clip, cutoff_right_clip, cutoff_clip_mate_in_rep, b_mosaic, wfolder_pub_clip, b_force,
max_cov_cutoff, sf_out)`.

CLI args and defaults: `--lc`/`cutoff_left_clip` (def 3, we pass **4**), `--rc`/`cutoff_right_clip`
(def 3, we pass **4**), `--cr`/`cutoff_clip_mate_in_rep` (def 1, we pass **2**), `--sc`/`siteclip`
(def 2 → `set_initial_min_clip_cutoff`), `-n`, `--cp` public clip folder, `-r` repeat copies,
`-a` annotation, `--cns` consensus, `--ref` genome, `-p` workdir, `-i` bam list, `-o` output,
`--cov` (→ `max_cov_cutoff = int(15*cov)`, default 600).

### 1.2 Read scan — per-alignment "is this a clip?"
Core loop: `clip_read.py:collect_clip_positions_by_chrm()` (295–399); quality/polyA gating in
`x_clip_disc_filter.py` (~144–233).

**Flag filters (skip the read if):**
- `is_duplicate` (clip_read.py 103–104)
- `is_unmapped` (108–109)
- `is_supplementary or is_secondary` (144–145)

**MAPQ filter:** `if mapping_quality < MINIMUM_CLIP_MAPQ (12): continue` (333–334).

**CIGAR clip detection** (pysam op codes: 0=M,1=I,2=D,4=S,5=H):
- **Left clip** iff `l_cigar[0][0] in {4,5}` (353).
- **Right clip** iff `l_cigar[-1][0] in {4,5}` (370).

**Base-quality gate** (`_is_qualified_clip`, clip_read.py 224–233): over the **clipped bases only**,
count bases with Phred `< CLIP_PHRED_SCORE_CUTOFF (15)`; **reject if count > len/2** (i.e. reject
if > 50% of clipped bases are sub-Q15). Left clip uses `query_quality[:l_cigar[0][1]]`
(x_clip_disc_filter 178–180); right clip uses `query_quality[start_pos:]` (223–225).

**Short-clip / poly-A gate** (x_clip_disc_filter 171–175 left, 217–219 right):
- if `MINIMUM_POLYA_CLIP(7) ≤ len_clip_seq ≤ BWA_REALIGN_CUTOFF(9)`: keep only if
  `xpolyA.contain_enough_A_T(clipped_seq, len_clip_seq-1)` is True (needs ≥ `POLYA_RATIO`=40% A/T).
- `elif len_clip_seq < BWA_REALIGN_CUTOFF(9)`: **skip** (too-short non-polyA clip).
- (longer clips proceed to bwa realignment against the consensus.)

### 1.3 Breakpoint position math
- `map_pos = algnmt.reference_start` (338, 0-based).
- **Left clip position = `reference_start`.** Clipped seq = `query_seq[:l_cigar[0][1]]`.
- **Right clip position:** starting from `reference_start`, add `lenth` for every CIGAR op in
  `l_cigar[:-1]` whose type is **not** in {4,5,1} (i.e. count ref-consuming ops M/D/…; skip S/H/I)
  (371–376). Clipped seq = `query_seq[-l_cigar[-1][1]:]`.

### 1.4 Outputs — **reference only, NOT a reproduction target** (used to build golden fixtures)
**(a) Per-chromosome clipped-part FASTQ** `{chrm}.clipped.fq` (`CLIP_FQ_SUFFIX`), 4 lines/record
(x_clip_disc_filter 183–192):
```
@{chrm}~{map_pos}~{L|R}~{is_rc}~{insertion_pos}~{n_cnt_clip}
{clipped_seq}            # clipped portion ONLY, not the whole read
+
{clipped_quality}
```
`is_rc` = read reverse-complemented (0/1); `insertion_pos`/`n_cnt_clip` per the record builder.

**(b) Per-chromosome initial position table** `{chrm}.clip_pos` (`CLIP_POS_SUFFIX`),
tab-separated, written for positions with `left+right ≥ i_clip_cutoff` (clip_read.py 393–398):
```
pos <TAB> left_cnt <TAB> right_cnt <TAB> mate_in_rep_cnt
```

**(c) Merged candidate list** `candidate_list_from_clip.txt` (`-o`) — join of `.clip_pos` with
`.clip_realign_pos` counts, keeping a position iff `NOT (left < cutoff_left AND right < cutoff_right)`
and `left+right ≤ max_cov_cutoff` (clip_read.py 1043–1092). Columns:
```
chrm <TAB> pos <TAB> left_cnt <TAB> right_cnt <TAB> mate_in_rep_cnt <TAB> realign_left <TAB> realign_right
```
(missing realign entry → `0\t0`).

### 1.5 Determinism (must sort before writing)
`for pos in m_clip_pos` (clip_read.py 393–398, 1074–1092) and `for pos in m_sites_chrm`
(x_clip_disc_filter 841–846) iterate dicts. **C tool must emit rows sorted by (chrm, pos)** so output
is byte-stable regardless of insertion order / Python version.

---

## 2. DISC step (`-D`)  — detection semantics

### 2.1 Entry / dispatch (`x_TEA_main.py` 424–473)
`options.discordant` → `filter_candidate_sites_by_discordant_pairs_multi_alignmts(m_sites_clip_peak,
iextend, i_is, f_dev, n_disc_cutoff, sf_annotation, sf_tmp, sf_raw_disc, b_tumor)`.
- `m_sites_clip_peak`: candidate clip-peak sites from the clip step (disc only looks near these).
- `iextend`: search half-window around a site ≈ `mean_is + 3*std + read_length` (from isize stats).
- `i_is`: large-insert constant (100000) used for orphan/cluster logic.
- `n_disc_cutoff`: `--nd` (we pass **5**).

### 2.2 Discordant-read criteria (`x_alignments.py` 107–123)
Skip the read if any of: `is_duplicate`, `is_supplementary`, `is_secondary`, `is_unmapped`,
`mate_is_unmapped`. MAPQ gate: `if mapping_quality < MINIMUM_DISC_MAPQ (20): continue` (113).
Decoy/non-primary contigs for the mate are filtered (122–123).

**Discordant test** (`x_clip_disc_filter.py` 278–287, `is_discordant`):
```
if chrm != mate_chrm:            return True    # mate on a different chromosome
elif abs(mate_pos - map_pos) > is_threshold:  return True    # same chrom, large insert
else:                            return False
```
`is_threshold = DISC_THRESHOLD = 2000`. (Distinct from `iextend`, which is the *site* search window.)

### 2.3 Association to a site + clustering
For each clip-peak site, count discordant pairs on the **left** window `[pos-iextend, pos]` and the
**right** window `[pos+1, pos+iextend]` separately (`x_TEI_locator.py` 1108–1111,
`cnt_discordant_pairs`). Mate positions are clustered (`disc_cluster.py:form_one_side_cluster`,
14–64): sort mate positions per chrom, break clusters on gaps `> i_is`, take the largest cluster;
return "consistent" iff `max_count/all_count > MIN_RAW_DISC_CLUSTER_RATIO` (**0.3**, line 68) plus the
cluster's chrom + midpoint. (Note: one extraction reported 0.6 here — that is `MIN_CLUSTER_RC_RATIO`;
**verify against code**, see §4-V1.)

### 2.4 Site acceptance (`--nd`) — results-determining (`x_TEI_locator.py` 452–464)
- **non-tumor (our path):** keep site iff
  `((n_left > cutoff/2) AND (n_right > cutoff/2)) OR (n_left > cutoff) OR (n_right > cutoff)`.
- tumor: keep iff `(n_left + n_right) >= cutoff`.

### 2.4b Disc realignment uses the MATE read (the TE-side read) — NOT the anchor
Confirmed in master (`x_clip_disc_filter.py` ~250-338): the scan records, per discordant anchor, a
`disc_pos` line with the anchor position + **mate coords + read name**; then
**`extract_mate_reads_by_name`** fetches the **mate reads' sequences** into
`candidate_sites_all_disc.fa`, which is what gets realigned to the consensus. So:
- The anchor (MAPQ≥20, maps uniquely near the breakpoint) only **localizes** the site; it is never
  realigned to the consensus.
- The **mate** (the read sitting in the inserted element) is the TE evidence — its sequence is
  realigned to the L1/Alu/SVA consensus (single-read alignment; pairing was resolved at genome-align
  time). This yields LDISC/RDISC ("aligned to copies").
- **xtea2 scan implication:** disc detection must ultimately provide **mate sequences**. Chosen
  approach = master-style: pass-1 records anchor+mate-coord+name; a second, **tileable**, binned sweep
  gathers the mate sequences by name (bounded memory, no per-read random seeks). The current prototype
  implements pass-1 (anchor + mate coords) only; mate-sequence gathering + disc realignment is TODO.

### 2.4c Realignment is single-read to the consensus — and that loses no power (rationale)
Paired bwa alignment to the consensus is **not applicable**: clip fragments have no mate; disc mate
reads' partners are the anchors, which are unique *genomic* sequence absent from the TE consensus, so
there is nothing to pair against. (Master realigns `candidate_sites_all_disc.fa` single-end for the
same reason.) Paired-end **power is preserved where it lives**: (1) candidate *detection* — the
discordant relationship + mate coords + **orientations** (`is_rc`/`is_mate_rc`); (2) **per-site
evidence aggregation/voting** across all clip+disc reads (where L1/Alu/SVA discrimination happens).
The real sensitivity lever is the **reference choice, not pairing**: realign against
`*_copies_with_flank.fa` (full TE copies WITH flanking genomic sequence) so junction-spanning reads
align (flank→genome, body→TE). **Design decisions:** (a) capture `is_mate_reverse` in disc records;
(b) realign clip fragments + disc mate reads against the with-flank references (bare consensus is a
fallback); (c) clip+disc candidate reads can share **one FASTQ per family → 3 bwa calls total**
(header tags clip|disc + site + side).

### 2.5 Outputs — **reference only, NOT a reproduction target**
For building golden fixtures / understanding only:
- `candidate_list_from_disc.txt`: `chrm \t pos \t n_left \t n_right`.
- `…​.clip_sites_raw_disc.txt` (`RAW_DISC_TMP_SUFFIX`): `chrm \t pos \t n_raw_left \t n_raw_right \t
  n_left_in_rep \t n_right_in_rep` (written iff a side is raw-consistent). **This file IS consumed by
  transduction** downstream — so whatever xTea2 produces, the transduction stage must get equivalent
  raw-disc info (as data, not necessarily this file).
- `.initial.all.disc` BAM (`ALL_DISC_SUFFIX`): **written but never read → dead; drop it** (§4-D1).

### 2.6 Determinism (master bug to avoid)
`x_TEI_locator.py` 453–454 & 523–524 iterate `m_merged_sites` / `m_merged_raw_sites` dicts in
first-encountered order → row order varies across runs. xTea2 must emit sorted by (chrm, pos).

---

## 3. Proposed architecture (the "best" way — for review)

Guiding results invariant: xTea2's candidate-site set + supporting reads must be a **superset that
converges to the same final calls**; the read-level tests pin the clip/disc read sets; the e2e gate
pins the VCF. Everything else is free to change.

### 3.1 Core idea: one type-agnostic single pass replaces clip + disc + repeated per-type scans
Master scans the BAM **clip pass + disc pass, separately, once per MEI type (L1/Alu/SVA)** — i.e.
~6 full scans of the alignments. The clipped parts and discordant reads collected are *identical
regardless of MEI type* (type only matters at consensus realignment). So:

- **Stage 1 — `xtea-scan` (C/htslib), ONE streaming pass over BAM/CRAM for all types at once:**
  per read, in a single visit, evaluate clip signal (CIGAR ends, MAPQ≥12, base-qual rule, short/polyA
  gate; L/R breakpoint math) **and** disc signal (mate discordant per §2.2, MAPQ≥20, flag filters).
  Accumulate compact per-position counters (small) and buffer only the *candidate* clipped
  parts / disc mate reads (a tiny fraction of reads). Computes isize stats in the same pass (no
  separate `.mean_temp_is`). CRAM handled natively via the reference.
  - **Parallelism — CONFIRMED by user (2026-07-17): NOT per-chromosome. Partition by fixed-size
    base-pair tiles** (configurable, e.g. 5–25 Mb) laid across *all* contigs and dispatched to a
    worker pool. Per-chromosome workers straggle badly (chr1 ≈ 5× chr21; thousands of decoy/alt
    contigs add overhead) and idle cores; equal-bp tiles load-balance and saturate every core.
    - **Correctness of tiling:** each read is owned by exactly the tile containing its **POS**
      (leftmost mapped coord); a worker doing a region query for `[a,b]` skips reads with `POS < a`
      so boundary-spanning reads are counted once. Small contigs are batched many-per-tile.
    - **Disc under tiling is safe:** discordance is classified from the read's OWN record
      (mate chrom/pos via RNEXT/PNEXT + mate flags) — a worker never needs to fetch the mate read.
      Per-site disc-mate clustering (§2.3) happens in the post-aggregation merge, position-indexed.
    - Each worker emits sorted partial results → deterministic ordered merge.
    - **Future option:** BGZF/CRAM **byte-range** chunking (partition the compressed file, not
      coordinates) for near-perfect I/O balance; deferred (heavier ownership bookkeeping).
  - **Memory:** streaming aggregation, not per-read Python dicts (master's ~300–500 MB sink).
  - **Files:** no per-chrom/per-type FASTQ explosion — candidate reads go to bwa via **stdin/FIFO**
    (see Stage 2), or at most one consolidated, sorted temp per run.

- **Stage 2 — consensus realignment, shared + concurrent:** the once-collected clipped/disc parts are
  realigned to the L1/Alu/SVA consensus libraries **concurrently** (3 bwa procs, reading FASTQ from a
  pipe), with **prebuilt/cached consensus indexes** (stop rebuilding every run). SAM piped → parsed in
  memory. This is the biggest wall-clock win.

- **Stage 3 — downstream in Python (readable):** site filtering, transduction (fed the raw-disc data
  from Stage 1), sibling, post-filter, gene annotation, **ONNX** genotyping, VCF merge. Reads from the
  compact Stage-1/2 result objects, not a maze of temp files.

### 3.2 Data model / file-count
Target: **final VCF only** as the durable output, plus an optional `--debug-bundle` (one tarball or a
handful of typed files) when explicitly requested. Internally prefer in-memory hand-off / pipes;
if a spill file is needed, use **one sorted, typed table + one BAM of candidate reads**, never
per-chrom/per-type fragments. Honor `--keep-intermediates` (default off).

### 3.3 Boundary decision — **CONFIRMED (user, 2026-07-17): "Collector + realign"**
C (`xtea-scan`) owns **Stage 1 (single-pass clip+disc collection) + Stage 2 (orchestrating concurrent
bwa realignment against cached consensus indexes)**. **Python owns Stage 3**: site filtering →
transduction → sibling → post-filter → gene annotation → ONNX genotyping → VCF. Rationale: hot loop
and I/O in C; the logic-heavy, readability-sensitive downstream stays Python (original ask).
Boundary artifact = a compact in-memory/one-file hand-off of **candidate sites + per-site supporting
clip/disc reads + their consensus-realignment results**, consumed by the Python filter stage.

### 3.4 Correctness harness (reframed)
- `test_clip_reads.py` / `test_disc_reads.py`: assert **set equality of clipped/disc reads found** at
  candidate sites vs the master golden (semantic; parse master's fixtures into read-sets, compare to
  xtea2's found read-sets). Fixtures in BOTH bam + cram.
- e2e concordance gate on the NA12878 chr1 VCF (±25 bp + type).

---

## 4. Bug / dead-code log (running)

Legend: **D** = dead code to remove, **V** = to verify, **B** = suspected logic bug, **Q** = quirk,
**DC** = deliberate design change (xTea2 intentionally diverges from master; needs truth-set validation).

- **DC1 (`mate_in_rep_clip_cnt` = union mate-in-repeat, replacing master's accidental Alu-only — user-approved 2026-07-18)**
  — **NAMING (2026-08-01):** this signal's canonical name is **`mate_in_rep_clip_cnt`** — *the number of
  clipped reads at a site whose mate lands inside any selected repeat family's annotation*. It was
  previously called **`col5`** throughout the C++, the Python, the CLI and this spec; that named a
  **column position** (it is the 6th realign-out field), not a concept, and is retired. The CLI flag is
  `--mate-in-rep-families`, with `--col5-families` kept as a deprecated alias. Master's name for it was
  `mate_within_rep_cnt`.
  — Master's `mate_within_rep_cnt` counts only clipped reads whose mate lands in the **Alu**
  annotation, and does so **non-deterministically** (whichever type wins the shared-`pub_clip` race —
  see Q2). xTea2 **intentionally changes** it to count a mate landing in **ANY** MEI repeat
  (`L1 | Alu | SVA`, each family's own boundary extend baked into its BED: Alu/L1 = 0, SVA = 200),
  emitted identically for all three types. Rationale: an SVA/L1/Alu insertion's clipped-read mates
  naturally anchor in repeat elements; Alu-only discards the SVA-mate evidence (measured: at golden SVA
  call sites 111639406 and 168451690 the Alu-only count is **0** while there are **8** and **6** reads
  whose mates land in SVA). Union ≥ Alu-only always, so more marginal SVA sites clear the cns
  consistency check → expected higher SVA recall, at a possible PPV cost. **This diverges from the
  frozen golden**, so the 100%-recall concordance gate no longer strictly applies to SVA; validate by
  measuring WGS SVA recall (does it recover the ~6 `--sva`-era misses?) AND PPV (no new false SVA).
  **Scope:** `mate_in_rep_clip_cnt` is consumed ONLY by the SVA cns path (`IS_CALL_SVA` branches, `x_clip_disc_filter`
  ~1025/1353), so L1/Alu call sets are unaffected by construction (verify chr1 Alu=27, L1=0 unchanged).
  Impl: single union count in `xtea_scan.cpp` (`clip_matein` int; `mate_in_rep()` ORs the regidx
  BEDs it receives) + `realign.hpp` (same value per family). Family set is selectable via
  `run_xtea.py --mate-in-rep-families` (C++ ORs whichever BEDs are passed), and is **independent of
  `--types`** — see the type-partition note below.
  **ADOPTED DEFAULT (user, 2026-07-18): `alu,sva` (SVA ∪ Alu)** — the biologically tightest set. It
  feeds only the SVA path, and SVA is a **composite** element (Alu-like domain + VNTR + SINE-R), so a
  mate landing in Alu is genuine SVA evidence; L1 is not part of SVA's composition and contributed ~0
  at every test site, so it's excluded by default. `l1,alu,sva` (full union) and `sva` (per-family
  strict — known to undercount) remain available. Still validate SVA recall+PPV on WGS.
  **TYPE-PARTITION NOTE (fixed 2026-08-01).** The contributing family set is a property of this
  signal's *definition*, not of which families the caller requested. `run_xtea.py` used to build the
  annotation BEDs by iterating the **requested** types and skipping any outside the mate-in-rep family
  set, making the BED set `requested_types ∩ mate_in_rep_families`:
  `--types L1,Alu,SVA` → Alu ∪ SVA; `--types SVA` → SVA only; `--types L1` → **0 at every site**.
  Since union ≥ subset always and this value gates the SVA cns consistency check, **`--types SVA`
  could emit fewer SVA calls than the same sample run as `--types L1,Alu,SVA`** — a live correctness
  bug, since asking for less must not change what you get for what you did ask for. The BED set is now
  derived from `--mate-in-rep-families` alone (`run_xtea.mate_in_rep_bed_plan`, whose signature
  deliberately excludes `types`), pinned by `tests/test_mate_in_rep_type_partition.py`.
  Not-yet-included: HERV / full genome-wide RepeatMasker (no combined hg38 annotation shipped).


- **D1** — `.initial.all.disc` BAM (`ALL_DISC_SUFFIX`, `clip_read.py:75-76,198`) written, never read. Remove.
- **V1** — disc cluster ratio: `MIN_RAW_DISC_CLUSTER_RATIO`=0.3 (global_values:68) vs an extraction
  that reported 0.6 in `form_one_side_cluster`. Confirm which constant the code actually passes.
- **B1?** — determinism: unsorted dict iteration in clip (`clip_read.py:393,1074`;
  `x_clip_disc_filter.py:841`) and disc (`x_TEI_locator.py:453,523`) → nondeterministic row order.
  Not a call-set bug (downstream re-sorts?) but must be fixed for stable goldens.
- **B3 (cns_filter is non-deterministic — FLIPS ACTUAL CALLS — high severity, found 2026-07-18)** —
  Two byte-identical invocations of the cns step (`x_TEA_main.py -N` → `call_MEIs_consensus`) on the
  SAME candidate list + BAM produce **different SVA call sets**, not just different field values.
  Measured on chr1 SVA (105 candidates, `-n 4`): two normal runs disagreed on ~10 sites including the
  golden calls **11770163 / 37775710 / 64384495** (each present in one run, absent in the other), plus
  coverage values off by ±1 and consistency flags flipped. So the **frozen golden VCF is one draw from
  a distribution**, and borderline SVA recall genuinely wobbles run-to-run (this is very likely what we
  chased as the "A2 dropped-SVA" and the mate-in-rep marginal-site sensitivity). Suspected source: the per-run
  coverage / insert-size estimation (`X_BasicInfo.get_cov_is_rlth`, written to
  `tmp/cns/basic_cov_is_rlth_info.txt`) run under multiprocessing without a fixed seed/sample, which
  shifts `disc_cutoff`/`iextend` and tips borderline sites; compounded by B1's unsorted-dict ordering.
  **Implications:** (1) byte-identity is an invalid validation bar for anything downstream of cns —
  use call-set concordance and accept borderline flips, or PIN the estimation first; (2) `-n 1` does
  NOT fix it — two single-process runs still flip 11770163 / 111639421 (plus reorder several sites via
  B1's unsorted-dict iteration), so it is not merely multiprocessing. Prime suspect: **`PYTHONHASHSEED`
  randomization** (default-on for str keys since Py3.3) varying set/dict iteration → tipping borderline
  cluster/sampling decisions [under test]; (3) **xTea2 fix (aligns with the plan's determinism goal):** compute
  coverage/insert-size ONCE, deterministically (fixed sample or fold into the single scan pass), so cns
  becomes reproducible — move B's coverage-folding (B.1.iv) does exactly this as a side effect. NOTE:
  this does not implicate the B.1 collect+realign bypass, which reuses byte-identical SAMs and diverged
  LESS than two normal runs.
  **✅ FIXED 2026-07-18 (user-approved "Approach A" — seed the root RNG).** Root cause pinned NOT to
  `PYTHONHASHSEED` but to `X_BasicInfo._random_slct_site` (`x_basic_info.py:456`): it draws
  `N_RANDOM_SITES=3000` sites with the **unseeded global `random`** (`random.randint`), so every call
  to `automatic_gnrt_parameters` re-estimates `mean_is`/`std`/`cov` on a different site set →
  `max_is = mean_is + 3·std` wobbles (747/749/757/763…) → shifts `i_concord_dist` and tips borderline
  SVA/ALU. (Multiprocessing is NOT the cause: per-site results are summed / median'd = order-invariant.)
  **Fix:** `global_values.BASIC_INFO_RANDOM_SEED=42` + a **scoped** `random.Random(seed)` in
  `_random_slct_site` (env `XTEA_BASIC_INFO_SEED` overrides; `none`/non-int → master's unseeded
  fallback). Leaves the global RNG untouched. Because every step funnels through the same seeded draw,
  cns/transduction/sibling all reproduce ONE estimate with no value-threading. **Verified:** golden
  chr1 run1≡run2 and deep-fold chr1 run1≡run2 (byte-identical call sets); `max_is=750` fixed across all
  steps. **Consequence:** the seed pins one deterministic answer, so borderline calls that master
  included only by a lucky draw now consistently fall out — seeded golden chr1 SVA=3 (frozen unseeded
  was 5; the 2 dropped are the flippers 11770151/37775710). **User accepted** deterministic-borderline
  losses as the correct trade for reproducibility. Seeded goldens frozen at `master_golden_seeded/`
  (chr1) and `master_golden_wgs_seeded/` (WGS); these are the new concordance contract.
- **B4 / DC2 (master applies base-quality at clip DETECTION inconsistently; xTea2 keeps it — master
  bug + accepted divergence, 2026-07-18)** — Master's clip-position **counter**
  `clip_read.py:collect_clip_positions_by_chrm` (feeds `candidate_list_from_clip.txt`) filters on
  **mapq only** (`MINIMUM_CLIP_MAPQ=12`), NOT base quality. `_is_qualified_clip` (the Q15/half-of-bases
  test, `clip_read.py:224`) is called ONLY when extracting clip sequences for realign
  (`clip_read.py:501,546`) — so master **counts low-quality clips it would reject at realign.**
  xtea-scan (`xtea_scan.cpp:is_qualified_clip`, `kClipPhredCutoff=15`, reject if >half bases <15 —
  logic byte-identical to master's `_is_qualified_clip`) applies mapq **and** baseq at detection.
  **Effect (chr1:104688481 ALU, the one seeded-golden miss):** ground truth = 26 right-clips at the
  peak 104688489; both drop 4 mapq<12 + 2 dup → master keeps 20, xtea-scan keeps 17. The 3-read gap is
  exactly the 3 lowest-quality polyA clips (mean BQ 17.0/17.3/18.8; each with >half its bases <Q15,
  e.g. 64/109 and 56/110 — genuinely low quality). 20 clears the cns threshold, 17 does not → master
  calls it, deep-fold drops it. **This is master applying its OWN baseq rule non-uniformly** (counts at
  detection, filters at realign). **DECISION (user, 2026-07-18): keep xtea-scan STRICTER** — reject
  low-quality clips at detection; do NOT reproduce master's inconsistency. So 104688481 is a documented
  accepted divergence (chr1 gate = 29/30, PPV 1.0), the same "xTea2 is more principled than master"
  pattern as DC1/`mate_in_rep_clip_cnt`. Toggle exists (`--no-baseq` disables xtea-scan's detection baseq) if a future
  reproduce-master run is ever needed. Consumed via the clip candidate count → affects any type, but
  only tips genuinely-borderline sites. NOTE: this baseq strictness may also drop other borderline
  low-quality calls on WGS — watch WGS recall/PPV; revisit if it costs real calls.
- **Q1** — `SEPERATOR` misspelling (value "~"); comment/code mismatches (MAPQ cmt 10 vs code 12;
  polyA-check cmt 10 vs code 20; ratio cmt 20% vs code 0.1). Reproduce the **code**, not the comment.
- **B2 (misleading comments — treat as bugs, per user 2026-07-18)** — inline comments that state a
  value the code does not use are defects, not cosmetics: they actively misled this refactor. Confirmed
  offenders to fix when the file is touched: `clip_read.py:333` `if algnmt.mapping_quality <
  MINIMUM_CLIP_MAPQ:  #by default this is set to 10` — the constant is **12** (`global_values.py:31`);
  plus the polyA (cmt 10 / code 20) and ratio (cmt 20% / code 0.1) cases from Q1. Action: correct the
  comments (or delete them) in the ported code so the value in the comment always matches the constant.
- **Q2 (`mate_in_rep_clip_cnt` / master's `mate_within_rep_cnt` is shared across all types, from ONE annotation — MAJOR)** —
  the 5th column of `candidate_list_from_clip.txt` (`mate_within_rep_cnt`) is **NOT computed per
  candidate type**. All three `-C` steps (L1/Alu/SVA) share one `--cp pub_clip/` dir; the FIRST `-C`
  to run populates `pub_clip/<sample>.clipped.fq` + `pub_clip/<chrm>.clip_pos`, and every later `-C`
  **skips `collect_clip_positions`** (gate: `x_TEI_locator.py:611` `if os.path.islink(sf_all_clip_fq)
  ==False`) and reuses those shared files. So the value is the mate-in-repeat count of **whichever type
  won the pub_clip race — Alu in BOTH the chr1 and WGS goldens** (`pub_clip/0/chrN.clip_pos` →
  `Alu/tmp/clip/0/chrN.clip_pos`), computed against the **Alu annotation with `boundary_extnd=0`**
  (Alu has no `--sva`, so `IS_CALL_SVA=False`). Verified read-for-read at chr1:11770163: golden count=5
  = Alu-extend-0 mate-in-rep (mates ~400 bp upstream in a nearby Alu); SVA-extend-200 gives only 1,
  L1 gives 0. **Consequence for xTea2:** the value must be the **Alu-BED (extend 0) mate-in-rep, emitted
  identically for L1/Alu/SVA** — NOT a per-family value and NOT SVA-extend-200. This overturns the
  earlier (wrong) "the count undercounts because of mapq≥10/no-baseq" diagnosis: mapq is 12 for both and
  baseq is irrelevant here; the only variable was the **annotation**. It also means master's value is
  **order-dependent/non-deterministic** (a different pub_clip winner would change SVA calls) — xTea2
  makes it deterministic by pinning it to the Alu annotation. It only affects the **SVA cns path**
  (`IS_CALL_SVA` branches), which is why this surfaced as an SVA-only recall gap.
- **B5 (SampleVCFMerger's per-type merge silently mis-resolves conflicts — found + fixed 2026-07-20,
  Phase C move A.1)** — `x_vcf_merger.append_sites_from_vcf` (:707-743) keys every locus by its own
  exact VCF position in ONE dict shared across all types. On a conflict, `cmp_quality_of_two_ins`
  decides `b_keep_first` (True = keep the already-placed/higher-priority call), but the code only
  `continue`s on the **override** branch (`b_keep_first==False`); when `b_keep_first==True` it falls
  through to an unconditional `m_sites[chrm][ins_pos] = line` (line 743), writing the **losing**
  record anyway. Effect depends on whether the two positions are exactly equal: (a) **different
  positions** — verified via direct invocation (synthetic L1@X + non-TPRT Alu@X+10bp): BOTH records
  survive in the merged VCF, a spurious duplicate call at a near-identical locus; (b) **exact same
  position** — the fallthrough write lands on the SAME dict key, silently clobbering the correctly
  "kept" entry, so whichever type loaded LAST (always **Alu**, load order L1→SVA→Alu) wins regardless
  of priority or `two_side_tprt_both`. Verified against the real WGS seeded merged golden
  (`master_golden_wgs_seeded/NA12878/NA12878.wgs.merged.golden.vcf`, 1529 calls): **all 5 of its
  non-TPRT-override conflicts are exact-position ties, and all 5 resolve to Alu** — the opposite of
  the intended L1/SVA priority winner in each case (chr1:241980544, chr3:110694547, chr5:142000088,
  chr12:99387860, chr17:33737479). This merger is not part of master's own pipeline (master's `-f
  5907` stops at per-type VCFs) — it is code xTea2 reuses for its single-VCF deliverable — so there
  is no "reproduce master" obligation here, only a real correctness bug to not carry forward.
  **FIX (Phase C move A.1, `xtea/locus_table.py`):** the new `LocusTable` keeps each chromosome's
  placed loci in a **list**, not a position-keyed dict, so a resolved conflict always
  replaces-or-drops the loser in place — never a silent key collision, never a duplicate, regardless
  of whether the two positions happen to be equal. Validated: re-running the corrected rule over the
  WGS seeded per-type goldens reproduces the old 1529-call merged golden **byte-for-byte except at
  those exact 5 loci**, which now resolve to the intended L1/SVA priority winner (matches
  `docs/handoff/conflict_analysis.py`'s independent list-based simulation, which never had this bug).
  chr1 (30 calls, no conflicts hit the bug) is unchanged. Unit tests: `tests/test_locus_table.py`.
- **Dead code (known, from planning):** `x_cmp.py`, `union_find_set.py`, `x_insert_size.py`, duplicate
  mains (`x_TEA_main.py` extra), CWL, 10X/barcode paths, commented case-control/de-novo.
- **V2** — golden L1 = **0** high-confidence calls on NA12878 chr1 (funnel: 75046 clip → 99 disc →
  4 cns → 1 cns2 → 0 HC). Clean collapse at consensus realignment, no crash. Plausible for one
  chromosome but low; validate against a truth set (GIAB/1KG MEI) in the deferred algorithmic phase.
- **V3** — master's `candidate_list_from_clip.txt` positions are cutoff-filtered + repositioned
  (nearby-region clustering + L/R TSD offset), so they differ from raw breakpoints by up to a few
  tens of bp for low-support sites. Oracle (raw breakpoints) covers 94% of golden Alu clip sites at
  ±25bp, 97% at ±50bp. Not a bug — a convention/filtering difference. Implication: read-level tests
  compare oracle↔xtea-scan (exact), NOT oracle↔master positions (fuzzy).
- **B6 (DF21 → ONNX genotype delta is 3 sites on the chr1 Alu golden, not 2 — measured 2026-07-25)**
  — `HANDOFF_2026-07-25.md` records the DF21→ONNX migration as disagreeing at **two** chr1 sites
  (`chr1:81058430`, `chr1:232452014`), measured over 29 sites. Running the downstream chain against
  the **Alu chr1 golden's own 27 high-confidence calls** (`tests/test_downstream_from_cns.py`) finds
  **three**: those two plus **`chr1:104688489`**. All three are `1/1` (DF21) → `0/1` (ONNX) — same
  direction, pure zygosity, no call gained or lost, and every other column of the genotype table is
  byte-identical. The two measurements are over different site sets, so this is a refinement rather
  than a contradiction, but **3 is the number for the Alu chr1 golden**. Consequence: none for the
  concordance gate (`tests/utils/concordance.py:4` — genotype is explicitly not required to match)
  and none for VCF columns 1-8, which are byte-identical. The harness pins this: only column 50 of
  the genotype table, and only column 10 of the VCF, may move.
- **B7 (SVA transduction emits a different intermediate site set than master's golden — found
  2026-07-25, downstream-neutral)** — Running the post-cns chain from the golden's own cns outputs,
  the SVA transduction step (`x_TEA_main.py --transduction`) produces
  `candidate_disc_filtered_cns2.txt` with **12 rows where the golden has 11** (extra site
  `chr1:182305180`), and `candidate_disc_filtered_cns2.txt.high_confident` with **9 where the golden
  has 10** (gains `chr1:182305180`, loses `chr1:50381610` and `chr1:182510697`). **Not
  nondeterminism:** two back-to-back runs from byte-identical staged input are byte-identical to each
  other — this is a stable difference between our `x_TEA_main.py` and the one that made the golden
  (`/home/cos689/src/xtea/xtea/`). **It washes out completely:**
  `candidate_disc_filtered_cns_post_filtering.txt`, the high_confident post_filtering table and the
  final `NA12878.chr1_SVA.vcf` are all byte-identical to the golden (5 calls either way), so postF
  filters all three sites regardless. Root cause not yet attributed — candidates are the xTea2 edits
  to `x_clip_disc_filter.py` (84 lines changed vs master) and the py3.8→3.10 interpreter change.
  Tolerated explicitly in `KNOWN_INTERMEDIATE_DIVERGENCES` so everything else stays a byte-exact
  gate. **Worth attributing before trimming transduction**, since a trim there would otherwise be
  credited or blamed for this pre-existing difference.
- **P1 (the interval tree is 96% of postF — measured 2026-07-25, the single biggest optimization)**
  — `cProfile` of one `--postF` invocation on chr1 Alu (~30 sites): **104.8 s total, of which
  `x_post_filter.construct_interval_tree` is 100.9 s (96%)**. Inside it,
  `x_annotation.index_rmsk_annotation_interval_tree` is **77.8 s** inserting **955,824** intervals
  one at a time into the vendored pure-Python AVL tree (13.4M `node.add` calls, 14.3M `rotate`), and
  `load_rmsk_annotation_with_extnd_div_with_lenth_cutoff` is **23.1 s** parsing the annotation. The
  actual post-filtering work is **~4 s**. This tree is rebuilt from the same annotation file in
  **~12 places per 3-type germline run** (postF ×2 per type, transduction ×1 per type, disc_adapter
  ×1 per type, plus gene/blacklist) — and the cost is *fixed*, independent of sample size, so it is
  the same ~25 s each on WGS as on chr1. The annotations are **static**: built once, queried
  read-only (17 `is_within_repeat_region_interval_tree` call sites, no mutation after construction).
  A sorted array + `np.searchsorted`, built once and reused, replaces `intervaltree.py` (1,142) +
  `interval.py` (305) + `node.py` (593) = **2,040 lines** with roughly 50. Biggest deletion and
  biggest speedup in the codebase are the same change.
- **D2 (mechanical dead code inventory — measured 2026-07-25 by AST survey)** — across `xtea/`:
  **13 modules / 7,940 lines unreachable** from the real entry points (`run_xtea` + `x_TEA_main`);
  **739 commented-out code lines**; **184 dead local variables** (assigned, never read); **26 unused
  imports**. The dead locals cluster exactly in the downstream: `x_orphan_transduction.py` **55**,
  `x_post_filter.py` **24**, `x_clip_disc_filter.py` **19**, `x_transduction.py` **10**,
  `x_somatic_calling.py` **10**. Unreachable modules, largest first: `gnrt_pipeline_local_chm13.py`
  1349, `gnrt_pipeline_local_v38.py` 1347, `gnrt_pipeline_local.py` 1251, `x_TEA_main_germline_cloud.py`
  1104, `x_vcf_merger.py` 1006, `gnrt_pipeline_germline_cloud_v38.py` 828, `x_cmp.py` 520,
  `x_insert_size.py` 262, `x_sv.py` 100, `x_population_AF.py` 72, `union_find_set.py` 65,
  `x_ref_TE_deletion.py` 18, `x_somatic_panel_normal.py` 18. Verified by hand that `x_vcf_merger`,
  `x_insert_size` and `union_find_set` appear only in comments/docstrings, not imports — the four
  `gnrt_pipeline_*` generators are what `run_xtea.py` replaced, and `x_vcf_merger.py` is what
  `locus_table.py` replaced (see B5). **Caveat before deleting:** a "dead" local whose right-hand
  side is a call with side effects is not dead — check each RHS, do not bulk-delete on the name.
- **D2-EXECUTED (2026-07-25 evening → step 1 of `handoff/TRIM_PLAN.md` is DONE)** — two of D2's four
  counts were re-measured during execution and **both were wrong in the safe direction**:
  - **Modules: 13 / 7,940 confirmed exactly.** Re-derived independently by an AST import-graph walk
    from `run_xtea` + `x_TEA_main` (45 of 58 modules reachable). Two methods agreeing is why this
    went in as one commit. Deleted in `120e207`.
  - **Unused imports: 23, not 26.** The other 3 lived inside the modules deleted in the same step.
    Removed in `bceff69`. Only one carried side-effect risk (`x_parameter.py`'s unused
    `import global_values`) and it was cleared by hand: `global_values` is constants plus mutators,
    stays a `sys.modules` singleton for its real importers, so dropping one unused binding is inert.
  - **Commented-out code: 438, not 739.** D2's "verified by regex as code, not prose" was too loose.
    Removing comments cannot change behavior, so the classifier — not the edit — was the work here.
    **The rule that held up:** a comment is disabled *code* only if executing it would DO something
    (a non-`Expr` statement, or an `Expr` containing a `Call`/`Await`/`Yield`). A commented-out
    *pure* expression has no effect as a statement, so it documents a convention instead. Removed in
    `52bb453`.
    **Four false-positive classes were caught only because the candidates were read**, and each
    would have destroyed something load-bearing:
    (a) `# ##INFO=<ID=...>` VCF header specs in `x_gvcf.py` — stripping one `#` leaves a comment,
    which `ast.parse` accepts as an *empty* module, and empty-parse was scoring as code;
    (b) `#coverage is abnormal` — English prose that parses as a valid `is` comparison;
    (c) **`interval.py:46,48` — `begin <= C < end`, the half-open interval convention, which is
    exactly what P1's differential test must reproduce**; (d) `run_xtea.py:821`, which records
    master's `"Ave coverage is <cov_eff>: automatic parameters (clip, disc, clip-disc)"` log format,
    i.e. the assertion-1 contract.
  - **A codebase convention worth knowing:** of 441 final candidates, exactly 3 had **two or more
    leading `#`**, and all 3 were documentation (two tuple layouts in `x_TEI_locator.py`, one
    `intervaltree` invariant). `##`+ comments are annotations here; `#` is disabled code. All 438
    removed lines are single-`#`; 772 comments were spared on this basis alone.
  - **Gate for the whole step:** `compileall` exit 0; 45/45 flat import sweep (note: `xtea/` is not
    a package and the modules use flat `import global_values` style, so a sweep must put `xtea/` on
    `sys.path` — importing them as `xtea.foo` fails on *every* module and is a broken test, not a
    regression); `--dry-run` 24 steps exit 0 with the **full step plan byte-identical** to the
    pre-trim baseline; fast suite 69 passed / 1 skipped throughout. **`xtea/`: 32,973 → 24,584
    lines (−8,389, −25%), 58 → 45 modules.** Dead locals (P3) deliberately deferred to last.
- **CC.3 job `48774776` (LFS 2921, both-fix, 3 types) — PASS, read 2026-07-25 evening.** The
  `b_tumor` plumbing fix works end-to-end. **Assertion 1b clears first**, which is what makes the
  rest meaningful: raw-disc rows are L1 7476 pass / 123,321 raw vs the pre-fix 3517 / 21,605 —
  completely different, so `b_tumor` did reach `disc_adapter`. Then: assertion 1 disc_cutoff=4
  matches master; **Alu somatic target `chr12:121357375` now survives every downstream stage into
  the final VCF** (pre-fix it was lost after the disc step); germline Alu **exact-position recall
  96.6% → 99.4%** (overlap 1104/1111, master-only 7, xtea2-only 26); final somatic VCF 6 calls vs 3.
  **Two gaps remain, both pre-existing rather than caused by the fix:**
  (i) **Alu `chr7:64580547` is still not recovered** — present in `candidate_list_from_clip.txt`,
  `..._from_disc.txt` and `.clip_sites_raw_disc.txt`, but **dropped at `candidate_disc_filtered_cns.txt`**,
  so the loss is inside cns filtering, not scan. It is now the *sole* remaining missed somatic target
  and the natural next attribution target after the trim.
  (ii) an **extra SVA call `chr14:39063995`** against an expected "SVA none". L1 `chr6:82922665`
  remains the already-known baseline extra (master's own clip step emits it too).
- **B8 (which repeat copy is "the" hit is decided by Python set iteration order — and it changes
  whether calls are kept — measured 2026-07-25 evening, blocks P1 as specified)**
  `x_annotation.is_within_repeat_region_interval_tree` (`x_annotation.py:602`) does:
  ```python
  set_rslt = tmp_tree[pos]          # set of ALL intervals overlapping pos
  for rcd in set_rslt:
      return True, rcd[0]           # <-- first element of a SET iteration
  ```
  When several extended repeat copies overlap the query point it returns **an arbitrary one** —
  whichever CPython's set happens to yield first. This is not a rare tie at an exact boundary; the
  annotations are loaded with a boundary extension (`_boundary_extnd=100` for Alu/SVA,
  `L1_boundary_extnd=400` for L1, `x_post_filter.py:53,59`), and since Alu copies sit a few kb apart
  genome-wide, extension makes neighbouring copies overlap constantly. Measured on the real hg38
  annotations at the real parameters, querying every interval start:

  | | intervals | multi-overlap | set-first ≠ leftmost |
  |---|---|---|---|
  | Alu, extnd 100, minlen 200 | 978,832 | 173,030 (17.7%) | 89,899 (**9.2% of all queries**) |
  | L1, extnd 400, minlen 300 | 471,881 | 213,440 (45.2%) | 114,414 (**24.3% of all queries**) |

  > **⚠ CORRECTION (2026-07-25, late night): the L1 row above was measured at the wrong
  > parameter.** `min_copy_len=300` is `REP_LINE_MIN_LEN`, but the L1 branch of
  > `x_post_filter.py:415` *overrides* it with `L1_min_ref_copy_len=125`, and 125 is what postF
  > actually runs. Re-measured at the real parameter, L1 is materially worse:
  >
  > | | intervals | multi-overlap | set-first ≠ min-div | **filter verdict flips** |
  > |---|---|---|---|---|
  > | Alu, extnd 100, minlen 200 | 978,832 | 173,030 (17.7%) | 86,997 (8.9%) | **28,238 (2.88%)** |
  > | L1, extnd 400, minlen 300 *(wrong param)* | 471,881 | 213,440 (45.2%) | 111,297 (23.6%) | 27,429 (5.81%) |
  > | **L1, extnd 400, minlen 125 *(real)*** | **769,133** | **429,799 (55.9%)** | **235,048 (30.6%)** | **53,856 (7.00%)** |
  >
  > The re-measuring harness reproduces the original row exactly at minlen 300 (45.2%, and
  > 24.2% vs the recorded 24.3%), which is what validates it. Note the last column here is
  > **set-first vs. min-divergence** — the delta the chosen fix actually introduces — whereas the
  > original 2.95%/5.79% figures below are set-first vs. *leftmost*. Both are real; the
  > min-divergence column is the one that predicts the re-baselining.
  >
  > **Consequence:** the re-baselined goldens will move more than recorded — **7.00%** of L1
  > positions' worth of opportunity, not 5.79%.

  **This is results-determining, not cosmetic.** The returned position is fed straight to
  `get_div_subfamily(ins_chrm, i_pos)`, and the divergence rate of *that* copy gates a filter that
  drops the call — `x_post_filter.py:90,97` (and :153,160): `if (b_in_rep is True) and
  (div_rate < global_values.REP_DIVERGENT_CUTOFF): continue`, cutoff **15%**. Two overlapping copies
  routinely fall on opposite sides of 15% (a young AluY at ~3% next to an old AluJ at ~18%), so the
  arbitrary pick decides call retention. Measured share of positions where the set's choice and the
  leftmost choice give **opposite filter verdicts**: **Alu 28,910 (2.95%), L1 27,299 (5.79%)**.

  Note the boolean is *not* affected — "does anything overlap" is order-independent, so the call
  sites that discard the position (`clip_read.py:346`, `x_alignments.py:139`, `disc_adapter.py:119`,
  which passes `_`) are safe under any rule. Only the position-consuming sites are exposed.

  **Consequence for P1:** the plan's "replace the tree with sorted arrays + `np.searchsorted`" cannot
  be byte-exact, because the natural sorted-array answer is leftmost-start and the reference answer
  is arbitrary. The differential test the plan mandates would *correctly fail*. P1's risk note
  ("tie behavior at exact boundaries") understated this by an order of magnitude.

  **The design that resolves it — and still gets the speedup.** The 78 s being optimized is interval
  *construction*, not querying (querying is ~4 s of 104.8 s). So: use sorted arrays + `bisect` to
  find candidates fast, and when there is more than one candidate, reproduce the original choice by
  putting the survivors in a real `set()` and taking `next(iter(...))`. That is cheap (≤5 elements)
  and is exactly bit-compatible, because `Interval.__hash__` is `hash((begin, end))`
  (`interval.py:128-134`) — a plain tuple of ints, so it is stable across runs and unaffected by
  `PYTHONHASHSEED`. The one residual risk is that set *insertion* order can matter on a hash
  collision, and the original inserts in tree-traversal order while the replacement inserts in
  sorted order; that is precisely what the differential test must check, over all ~1.45M positions.

  **DECISION (user, 2026-07-25 evening): FIX IT — this is a bug, not an artifact to preserve.**
  Step 2 therefore stops being a pure optimization and becomes an optimization plus a deliberate,
  signed-off behavior change. Consequences accepted: calls move at ~3-6% of positions' worth of
  opportunity, **the goldens must be re-baselined**, and B6/B7's deltas must be re-measured against
  the new baseline (they were measured against the old one).

  **The rule chosen, and why it is the only one the code justifies: the MINIMUM-divergence
  overlapping copy** (tie-broken by leftmost start for total determinism). The filter exists because
  a site falling in a *young* copy of the same TE type is suspect — reads from a near-consensus
  reference copy mismap and mimic a novel insertion, whereas a diverged copy does not cause that
  confusion. So if **any** overlapping copy is young, that false-positive risk is present. Taking the
  minimum divergence is the reading that matches the filter's purpose, and it makes the predicate
  monotone and order-independent. "Leftmost" would only make the arbitrariness *deterministic*, not
  correct — it asserts nothing biological.

  **Direction of the change, stated up front:** minimum-divergence is strictly more conservative than
  the status quo. In every straddling case the old code sometimes picked the *old* copy and kept a
  call that the new rule drops, so xTea2 will **drop somewhat more calls** than master here —
  trading sensitivity for specificity. That is the intended direction for a false-positive filter,
  but it is a real trade and the re-baselined goldens should be inspected with it in mind.

- **B8b (same defect, live, output-affecting but not call-affecting) — `x_gene_annotation.py:169`
  `query_by_position` builds its hit list by iterating a set.** `l_hits.append(start_pos)` inside
  `for rcd in set_rslt:` means that when a position overlaps **multiple genes**, the *order* of the
  returned genes is set-iteration order. `get_gene_id_type` preserves that order into `l_rslt`, and
  `annotate_results` (`:217-223`) joins them with `"/"`. So the VCF `GENE_INFO` field can read
  `intron:A/intron:B` or `intron:B/intron:A` arbitrarily. This does **not** move calls (VCF columns
  1-8 are unaffected), so the current gate cannot see it — but it is a genuine reproducibility wart
  in a user-facing annotation. **Fix alongside B8: sort the hits (by start position) before joining.**

- **B8c (same pattern, LATENT — do not "clean up" into a live bug) — `x_black_list.py:52`
  `fall_in_region`.** Identical `for rcd in set_rslt: return True, start_pos` shape, and the author
  annotated it `#here only check hit or not, doesn't matter how many records in all`. Verified that
  the claim holds: all five callers (`x_post_filter.py:74,134,216`,
  `x_orphan_transduction.py:249,681`) bind the position as `tmp_pos` / `tmp_pos2` and **never read
  it** — `x_post_filter.py:236` reassigns `tmp_pos` from `ins_pos` for an unrelated purpose. Since
  the boolean is order-independent, this one is harmless *today*. It is recorded because those dead
  bindings are in D2's dead-local list and will be touched in step 5: **if anyone ever starts using
  that returned position, it inherits B8's bug.** Either fix the selection here too, or make the
  function return only the boolean so the trap cannot be re-armed.

  Common root cause across all three: `intervaltree`'s point query returns an unordered `set`, and
  three separate call sites treat "some element of that set" as if it were "the" answer. `set` is the
  wrong container to hand to code that then picks one element.

- **B9 (the hand-rolled binary search is py3-broken and still on a live code path — found
  2026-07-25 late night while evaluating P1)** `x_annotation.py:582`, inside
  `is_within_repeat_region`, does `mid = (lo + hi) / 2`. Under py2 that is integer division;
  under py3 it is a **float**, and the next line indexes a list with it. Confirmed by direct
  execution against a two-interval annotation:

  ```
  TypeError: list indices must be integers or slices, not float
  ```

  Every call reaching it raises. The three call sites are `clip_read.py:128` (dispatched by
  `pool.map(unwrap_self_cnt_clip, …)` at `clip_read.py:249`), `clip_read.py:629` (via
  `x_TEI_locator.call_TEI_candidate_sites_from_clip_reads`), and `x_alignments.py:805-806` (via
  `x_TEI_locator:871`). All three sit in the **master pysam clip/disc fallback** — the ~2,700
  lines P5b deletes in step 4 — which is why xTea2's default flow (C++ scanner) never trips it.

  **This changes a premise in the plan.** P5b's stated "accepted cost" is that the fallback stays
  available as the controlled variable for attribution runs (`a2_validate_chr1_cnsemit.sh` varies
  the scan binary while holding everything else fixed). But the fallback **cannot run under
  py3.10 at all** — so that capability is already lost, and has been since the py3 port. The
  argument for sequencing P4 before P5b is correspondingly weaker: attributing B7 against
  master's untouched clip/disc path is a git-checkout-plus-fix exercise either way, not a flag
  flip. Do not spend a run discovering this.

  Not fixed here deliberately: the function is deleted in step 4, and fixing it now would mean
  fixing a *wrong algorithm* (see the evaluation below — plain binary search on starts is
  incorrect under overlap, independent of the float bug).

- **B10 (master double-counts disc support: the mate FASTA emits one record per mate
  ALIGNMENT, not per mate READ — found + fixed 2026-07-28, MOVES CALLS)**
  `x_alignments.py:extract_mate_reads_of_region` builds the disc evidence FASTA by walking
  `bamfile.fetch(chrm, bin_start, bin_end)` and writing one record per `(s_read_id, insertion)`
  it matches, where `s_read_id = qname~is_read1`. It applies **no secondary/supplementary
  filter**, so every additional alignment of the same mate in the bin writes another copy of
  that mate — same header, same sequence. Downstream (`parse_disc_algnmt_consensus`) counts
  records, so one read becomes two or three units of disc support at that site.

  Two distinct paths in: (i) a **soft-clipped supplementary** alignment of the mate (`is_read1`
  is set on it too, and soft-clipped supplementaries carry the full `query_sequence`, so the
  emitted record is byte-identical to the primary's); (ii) the bins **overlap** — `fetch`
  returns reads overlapping `[bin_start, bin_end)` while `_parse_read_name_pos_of_region`
  (`x_alignments.py:343`) admits mates up to 100 bp beyond either edge, so a mate straddling a
  boundary is emitted by both adjacent bins.

  **Measured on NA12878 WGS (job `48842713`, Alu):** 28,173 of 583,489 distinct disc records —
  **4.8%** — appear more than once (28,071 twice, 102 three times). Duplicates concentrate in
  pericentromeric chr16 (10,170) and chr17 (12,360), i.e. exactly where supplementary alignments
  cluster. Worked example: `A00217:76:HFLT3DSXX:1:1101:10366:17190`, anchor chr16:34588134,
  mate read1 primary at chr16:46394458 (flag 81) *plus* a read1 supplementary at chr16:34902509
  (flag 2113, `98S52M`) — both in chr16's first 50 Mb bin, so both candidates (34587860 /
  34588350) got two identical records each.

  **Master's own clip collection does not do this** — `collect_clipped_disc_reads_by_region`
  explicitly `continue`s on `is_supplementary or is_secondary` (`x_clip_disc_filter.py:161,201`).
  Neither does xTea2's C++ mate gather (`cns_emit.hpp:emit_cns_from_spill`), which skips
  supp/sec and stores one seq per `(qname, is_read1)`. This path was the outlier, and that is
  what made it visible: see the sensitivity analysis in
  `handoff/HANDOFF_2026-07-26_wgs_run.md` — the two cns collection paths turned out to collect
  **the same reads** (clip: set-identical for all three types; disc: identical 583,489-key set),
  and the *only* substantive difference was this duplication.

  **Fix:** skip secondary/supplementary, and own each alignment by its start
  (`reference_start >= bin_start`) so the bins become a true partition. Owning by start is
  lossless: the bins tile each contig from 0, and the `mate_pos` recorded in the disc-pos file
  *is* the mate's `reference_start`, so every mate has exactly one owning bin. Regression tests:
  `tests/test_disc_mate_extract.py` (synthetic BAM, no external data — supplementary, secondary,
  boundary-straddle, and a not-lost-across-the-boundary case).

  **Consequence — this changes the call set.** Removing inflated disc support is strictly more
  conservative at sites that only cleared the cns cutoffs on double-counted reads. Every golden
  on file was produced with the double count, so the goldens must be re-baselined (see the
  handoff's re-baselining item, which is deliberately sequenced after this decision).

- **P1-EVALUATED (why a tree was there, and what replaces it — measured 2026-07-25 late night,
  supersedes P1's design and its cost figures)**

  **Why the tree exists.** `intervaltree.py` is a vendored copy of the PyPI `intervaltree`
  package: a general AVL structure with dynamic insert/delete, range overlap, merging, chopping
  and slicing. xTea uses exactly **one** operation — `tree[point]` — after a build-once phase.
  Nothing is inserted or deleted after construction and no range query is ever issued.

  The tree *replaced* something, and that something is still in the file: `index_rmsk_annotation`
  + `is_within_repeat_region` (`x_annotation.py:565-591`), a sorted-array binary search. It is
  **incorrect under overlap** — it binary-searches on start positions and returns at the first
  midpoint hit, so when copies overlap it can walk past a genuine hit. So the tree was almost
  certainly adopted as a **correctness** fix, not for its data-structure properties. The right
  answer to that correctness problem is not a tree either: it is sorted arrays plus a
  prefix-max-end scan, which is correct under overlap *and* ~100× cheaper to build.

  **The plan's profile numbers are cProfile-inflated ~3×.** Measured unprofiled, warm cache:

  | | Alu (extnd 100, minlen 200) | L1 (extnd 400, minlen 125) |
  |---|---|---|
  | parse | 2.75 s | 2.13 s |
  | **AVL build** | **23.6 s** | **18.0 s** |
  | sorted-array build | 0.22 s (**105×**) | 0.24 s (**75×**) |
  | merged-boolean build | 0.21 s (113×) | 0.14 s (129×) |
  | AVL point query /100k | 0.75 s | 1.00 s |
  | sorted-array query /100k | 0.06 s (12×) | 0.08 s (13×) |
  | numpy vectorized /100k | 0.02 s (37×) | 0.02 s (55×) |

  All four structures returned identical answers on 100k random probes. P1's "77.8 s insert /
  23.1 s parse" were under cProfile; the 96%-of-postF *ratio* survives (~90%), and P1's own
  "postF ~29 s" prediction matches the unprofiled measurement, so the conclusion stands and only
  the magnitudes move.

  **Why the sorted-array scan is cheap here:** intervals are short and tightly bounded — Alu
  median 498 / p99 532 / max 1,013; L1 median 1,181 / p99 6,774 / max 9,304. The backward scan
  from `bisect_right(starts, p)` therefore terminates within a handful of entries. This is a
  property of the data, so record it: if an annotation with long intervals is ever introduced,
  re-check it.

  **Where the run-level time goes:** ~12 constructions per 3-type germline run × ~22 s ≈ **270 s**.
  Sorted arrays alone → ~33 s; adding the build-once cache → ~14 s. So the structure swap is the
  dominant win and the cache is a worthwhile second (~19 s) — *not* the other way round.

  **Options evaluated and rejected**, so they are not re-proposed:

  | option | verdict |
  |---|---|
  | pandas/numpy column parse | **slower** — 3.5 s vs the Python loop's 2.75 s |
  | numpy vectorized batch | 37–55×, but the hot sites query one position per alignment inside pysam loops where decode dominates; batching means restructuring call sites for a non-bottleneck |
  | NCLS / cgranges / pyranges | a new dependency to beat 0.22 s |
  | on-disk precomputed index | parse is 2.5 s and the cache covers repeats; not worth the invalidation surface |
  | genome-wide bitset | ~390 MB for hg38 |

  **The design change vs. P1 as written: split the API by what callers need.** The consumers
  divide cleanly, verified site by site:

  | class | sites | needs | volume |
  |---|---|---|---|
  | hot | `clip_read:346`, `x_alignments:139`, `disc_adapter:119`, `x_orphan_transduction:1119`, `x_mosaic_calling:60/96/132` | **boolean only** — position is bound and never read | per-alignment (millions) |
  | cold | `x_post_filter:78/138`, `x_transduction:80` | boolean + position + divergence | per-candidate (hundreds) |

  Giving the hot class a boolean-only entry point makes **B8 and B8c non-bugs by construction** —
  there is no position to mis-pick, so the trap cannot be re-armed by a later edit. That is what
  B8c asks for, achieved structurally rather than by a comment.

  **One hazard the split must respect:** the min-divergence rule can only be applied where
  divergence was actually loaded. `get_div_subfamily` blindly reads field `[0][1]`, but only the
  four `..._divgnt...` / `..._div_...` loaders put `mismatch_rate` there. The non-div loader
  `load_rmsk_annotation_with_extnd_with_lenth_cutoff` — used by `clip_read`, `disc_adapter` and
  `x_TEI_locator` — puts **`b_rc`, a bool**, in that slot. It is harmless today only because
  those three are exactly the boolean-only consumers. Gate the rule on a flag set by the
  div-bearing loaders; fall back to leftmost-start otherwise.

- **P1-EXECUTED (step 2 outcome — 2026-07-25 late night, commits `4a9e824` evidence / `65fca02` code)**

  The AVL tree is deleted (`intervaltree.py` + `interval.py` + `node.py`, **1,996 lines**) and
  replaced by `xtea/interval_index.py` (**147 lines**, stdlib-only; a further +110 net across the
  11 edited modules is mostly the comments explaining the change, so **−1,739 net**). `xtea/` is **22,845 lines / 43
  modules**, from 24,584 / 45. `sortedcontainers` dropped from `envs/xtea2.yml` — it backed only
  the deleted tree.

  **The differential, run against the tree before deleting it**, over every interval start in the
  real hg38 annotations at the real parameters — Alu (100/200), L1 (400/125) and the SVA parameter
  set (200/300), **2,188,624 positions**:

  | | boolean diffs | hit-set diffs | position changed | verdict flipped | newly kept |
  |---|---|---|---|---|---|
  | Alu | **0** | **0** | 86,997 (8.9%) | 28,238 (2.88%) | **0** |
  | L1 | **0** | **0** | 235,048 (30.6%) | 53,856 (7.00%) | **0** |
  | SVA params | **0** | **0** | 28,988 (6.6%) | 7,381 (1.67%) | **0** |

  Boolean and hit-set identity everywhere is the correctness result. **Zero newly-kept** is the
  interesting one: every verdict flip is in the newly-dropped direction, which is the empirical
  confirmation that min-divergence is monotone — predicted in B8, now measured, not assumed.

  **Downstream gate after the change: PASSES, 8m25s (was ~11 min).** No final call moved on the
  chr1 golden. That is *not* evidence the fix is inert — the percentages above are over
  genome-wide annotation positions, while the golden has ~30 candidate sites, none of which
  happened to sit on a straddling copy. Four L1 transduction intermediates did move
  (`...unique_trsdct_half_clip`, `...clip_sites_raw_disc.txt.slct`, `.slct2`, `.slct2.tmp`),
  visible only because 5c kept intermediate comparison as diagnostics. **They are in the code
  path P4/step 3 attributes, so read them before attributing B7.**

  **Why the gate did not hit the predicted "~2-3 min":** the tree was 96% of *postF*, not of the
  whole downstream chain, so the saving was bounded by postF's share. Construction itself did
  improve ~100× (23.6 s → 0.22 s) and is now cached. The remaining time is other work — profile it
  before optimizing, which is how P1 was found to begin with.

  **Still owed:** re-baseline the goldens, and re-measure B6 and B7 against the new baseline.

- _(append findings as the refactor proceeds)_

---

## 5. Golden fixtures & the L1 coverage gap
From the master baseline (`master_golden/NA12878/{L1,Alu,SVA}/`, all 3 jobs COMPLETED exit 0):
parse the clip/disc outputs into read-sets/site-sets to seed the semantic unit tests.

**Frozen golden call counts (NA12878 chr1, `-f 5907`):** Alu **27**, SVA **5**, L1 **0** (final HC/VCF).

**Coverage gap → RESOLVED by the whole-genome golden (2026-07-17).** `master_golden_wgs/` (NA12878
whole genome, primary contigs) gives **L1 203, Alu 1262, SVA (running)** calls genome-wide.

**IMPORTANT methodology finding (V4):** a single-chromosome subset BAM is **NOT a faithful input for
L1 calling**. chr1-only golden = **0** L1; the same chr1 inside the WGS run = **9** L1. Mechanism:
xTea L1 evidence (`LDISC`/`RDISC` = discordant reads "aligned to copies") relies on **mate reads that
map into L1 copies elsewhere in the genome**; subsetting to chr1 drops those mates, so chr1 L1 support
collapses. Alu/SVA are clip/TSD-local and survive subsetting (27, 5). **Consequence:** use the
**whole-genome golden as the primary L1/Alu/SVA concordance baseline**; keep chr1-only only for fast
Alu/SVA unit tests. The read-level clip/disc oracle is unaffected (local per-read detection).
