# xTea2 — Clip/Disc step contract specification

**Purpose.** This is the **detection-semantics contract** for reimplementing xTea's **clip** (`-C`)
and **discordant / disc** (`-D`) read detection. It captures exactly *which reads and which candidate
sites master finds*, and the numeric thresholds that decide that — because **that is what determines
the final results**. Every value here is extracted verbatim from the `master`/`xTea2` source at commit
`46520ab` (`/home/cos689/src/xtea/xtea/`). Line numbers refer to that tree.

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

Legend: **D** = dead code to remove, **V** = to verify, **B** = suspected logic bug, **Q** = quirk.

- **D1** — `.initial.all.disc` BAM (`ALL_DISC_SUFFIX`, `clip_read.py:75-76,198`) written, never read. Remove.
- **V1** — disc cluster ratio: `MIN_RAW_DISC_CLUSTER_RATIO`=0.3 (global_values:68) vs an extraction
  that reported 0.6 in `form_one_side_cluster`. Confirm which constant the code actually passes.
- **B1?** — determinism: unsorted dict iteration in clip (`clip_read.py:393,1074`;
  `x_clip_disc_filter.py:841`) and disc (`x_TEI_locator.py:453,523`) → nondeterministic row order.
  Not a call-set bug (downstream re-sorts?) but must be fixed for stable goldens.
- **Q1** — `SEPERATOR` misspelling (value "~"); comment/code mismatches (MAPQ cmt 10 vs code 12;
  polyA-check cmt 10 vs code 20; ratio cmt 20% vs code 0.1). Reproduce the **code**, not the comment.
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
