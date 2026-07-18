# xTea SAM Contract Analysis: temp_clip.sam and temp_disc.sam

## A. CLIP READ FASTQ / SAM Contract

### A.1 Read-name (query-name) Format in candidate_sites_all_clip.fq

**Initial format (before sample_id addition):**
```
@chrm~map_pos~FLAG~is_rc~insertion_pos~n_cnt_clip
```

**Final format (after sample_id addition in collect_clipped_disc_reads):**
```
@chrm~map_pos~FLAG~is_rc~insertion_pos~n_cnt_clip~sample_id
```

**Field breakdown (0-indexed):**
- Field 0: `chrm` - chromosome name from candidate sites list (string)
- Field 1: `map_pos` - original mapping position on genome where read was clipped (integer, 0-based genomic coordinate)
- Field 2: `FLAG` - clipping direction: `L` (global_values.FLAG_LEFT_CLIP) or `R` (global_values.FLAG_RIGHT_CLIP)
- Field 3: `is_rc` - read orientation: `0` (forward) or `1` (reverse complement)
- Field 4: `insertion_pos` - candidate insertion position (integer, 0-based)
- Field 5: `n_cnt_clip` - counter for reads clipped at same position (integer, incremented per clipped read at that position)
- Field 6: `sample_id` - sample index, added during merge (integer, represents BAM file index in sample list)

**Separator:** `~` (global_values.SEPERATOR, line 57 of global_values.py)

**Read name generation locations:**
- Left clip: x_clip_disc_filter.py:183-186 (collect_clipped_disc_reads_by_region)
- Right clip: x_clip_disc_filter.py:208-211 (collect_clipped_disc_reads_by_region)
- Sample ID addition: x_clip_disc_filter.py:3346 (collect_clipped_disc_reads, FASTQ header line at n_cnt % 4 == 0)

### A.2 Sequence Content
- **What is written:** Only the clipped portion of the read
- **Left clip sequence:** query_seq[:l_cigar[0][1]] (first l_cigar[0][1] bases), x_clip_disc_filter.py:168
- **Right clip sequence:** query_seq[start_pos:] where start_pos = -1 * l_cigar[-1][1], x_clip_disc_filter.py:215
- **No reverse-complementing:** The clipped sequence is written as-is (same orientation as in BAM)
- **Quality:** Converted to ASCII via _cvt_to_Ascii_quality (x_clip_disc_filter.py:182, 227)
- **Length constraints:**
  - Minimum: >= global_values.MINIMUM_POLYA_CLIP (7 bp) if appears to be polyA/T (global_values.MINIMUM_POLYA_CLIP=7, line 27)
  - OR >= global_values.BWA_REALIGN_CUTOFF (9 bp) for non-polyA (global_values.BWA_REALIGN_CUTOFF=9, line 26)
  - Maximum: <= global_values.BWA_REALIGN_CUTOFF (9 bp) enforced at line 171 for polyA check

### A.3 The bwa Command for Clipped Reads

Function: `realign_clipped_read_with_polyA` (bwa_align.py:31-60)

**Process:**
1. **Split reads by length:** Lines 35-41 separate into polyA cohort (length < 9) and non-polyA cohort (length >= 9)
2. **Non-polyA realignment:** bwa_align.py:43
   ```
   bwa mem -t {n_jobs} -T 9 -k 9 -o {sf_sam1} {sf_ref} {sf_reads1}
   ```
   - `-T 9`: global_values.BWA_REALIGN_CUTOFF = 9 (minimum seed length)
   - `-k 9`: global_values.BWA_REALIGN_CUTOFF = 9 (minimum interval size)

3. **PolyA realignment:** bwa_align.py:45, 67-74
   ```
   bwa mem -t {n_jobs} -T 7 -k 7 -c 70 {sf_ref} {sf_reads2}
   ```
   - `-T 7`: global_values.MINIMUM_POLYA_CLIP = 7 (minimum seed length for polyA)
   - `-k 7`: global_values.MINIMUM_POLYA_CLIP = 7 (minimum interval size for polyA)
   - `-c 70`: self.BWA_SEED_MEDIUM_FREQ = 70 (seed occurrence threshold for polyA)

4. **Merge SAMs:** bwa_align.py:48-56 merge both SAM files, dropping headers from polyA SAM

### A.4 Parsing: parse_clip_realignment_consensus

Location: x_clip_disc_filter.py:2946-3141

**Query-name fields read (from parsed ~-delimited query_name):**
- Field [0]: ori_chrm (chromosome, line 2967)
- Field [1]: ref_mpos (original mapping position, line 2968)
- Field [2]: ori_clip_flag (FLAG_LEFT_CLIP or FLAG_RIGHT_CLIP, line 2969)
- Field [3]: ori_is_rc (reverse complement flag, line 2970)
- Field [4]: ori_insertion_pos (insertion position, line 2971)
- Field [6]: sample_id (sample index, line 2972)

**SAM fields read per alignment:**
- FLAG: Implicitly checked for unmapped status (line 2983)
- CIGAR: Used via algnmt.cigar to compute qualified alignment fraction (line 2986)
- POS (reference_start): Map position on repeat consensus (line 3052)
- MAPQ: Not explicitly gated (commented out at line 2960), but collected for polyA detection
- Query sequence: algnmt.query_sequence for polyA detection (line 2992)
- Reference name (mapped copy ID): algnmt.reference_name, format "chrm~start~end~subfamily" (line 3052 context)

**Qualification check (bmapped_cutoff):**
- Function: is_clipped_part_qualified_algnmt (x_clip_disc_filter.py:2731-2764)
- Logic: n_map >= (n_total * ratio_cutoff) where:
  - n_map: number of matched bases (CIGAR type 0, x_clip_disc_filter.py:2757)
  - n_total: total query length excluding deletions (CIGAR type != 2, x_clip_disc_filter.py:2760)
  - ratio_cutoff (bmapped_cutoff): default 0.65 (global_values.MIN_CLIP_MAPPED_RATIO=0.65, line 118)
  - Hard requirement: at least 65% of the clipped read must map
- Skip if both-end clipped and both > 8 bp (global_values.MAX_CLIP_CLIP_LEN=8, line 2749)

**Returns:**
- m_ins_pos_left: dict[chrm][insertion_pos][ref_mpos] = [pos_in_consensus, ...]
  - ref_mpos: original mapping position on genome where clipped part originates
  - pos_in_consensus: position computed on repeat consensus (line 3052-3056)
  - Computation (line 3052-3056):
    - Start with algnmt.reference_start (map position on consensus)
    - If (b_clip_part_rc AND b_left) OR (NOT b_clip_part_rc AND NOT b_left), add n_map_bases (number of matched bases)
    - This accounts for whether alignment is reverse-complemented

- m_ins_pos_right: same structure for right-clipped reads

- m_ins_pos_sample: dict[chrm][insertion_pos][sample_id] = [left_count, right_count]
  - Tracks which samples support each insertion and on which side

- m_polyA: dict[chrm][insertion_pos] = [left_polya_count, right_polya_count, left_polya_false, right_polya_false]
  - Detects polyA/T poly-tail signals (line 3032-3049)

- m_transduction: dict[chrm][insertion_pos][transduct_pos] = [pos_in_consensus, ...]
  - For clipped parts mapped to flank regions (potential transductions)
  - Only if: map_pos_copy < (i_flank_lenth - 5) OR map_pos_copy > (copy_end - copy_start + i_flank_lenth - 5)
  - AND mapq >= global_values.MINIMAL_TRANSDUCT_MAPQ (10) and not polyA

---

## B. DISC READ FASTA / SAM Contract

### B.1 Read-name Format in candidate_sites_all_disc.fa

**Initial format (before sample_id addition):**
```
>read_name~is_first~b_lclip~b_rclip~s_anchor_rc~s_anchor_mate_rc~anchor_map_pos~ins_chrm~ins_pos
```

**Final format (after sample_id addition in collect_clipped_disc_reads):**
```
>read_name~is_first~b_lclip~b_rclip~s_anchor_rc~s_anchor_mate_rc~anchor_map_pos~ins_chrm~ins_pos~sample_id
```

**Field breakdown (0-indexed):**
- Field 0: `read_name` - original read pair name from BAM (string)
- Field 1: `is_first` - is_read1 flag: `1` if read1, `0` if read2
- Field 2: `b_lclip` - anchor read has left clip: `0` or `1`
- Field 3: `b_rclip` - anchor read has right clip: `0` or `1`
- Field 4: `s_anchor_rc` - anchor read reverse complement: `0` (forward) or `1` (reverse)
- Field 5: `s_anchor_mate_rc` - mate read reverse complement: `0` (forward) or `1` (reverse)
- Field 6: `anchor_map_pos` - left-most mapped base position of anchor read (integer, 0-based)
- Field 7: `ins_chrm` - chromosome where anchor read (with clipped read) maps
- Field 8: `ins_pos` - insertion position (candidate site, integer 0-based)
- Field 9: `sample_id` - sample index, added during merge (integer, represents BAM file index in sample list)

**Separator:** `~` (global_values.SEPERATOR)

**Read name generation locations:**
- Generated during collection: x_clip_disc_filter.py:259-262 (collect_clipped_disc_reads_by_region) writes to disc_pos file
- Extracted and transformed: x_alignments.py:393-394 (extract_mate_reads_of_region)
- Sample ID addition: x_clip_disc_filter.py:3339 (collect_clipped_disc_reads, FASTA header line at n_cnt % 2 == 0)

### B.2 Sequence Content

**Which read of the pair:** The mate read (the one NOT anchored to the clipped read)
- The anchor read maps near the clipped read
- The mate read (extracted) is the one that maps into a repeat region or elsewhere
- Extracted from BAM via x_alignments.py:383: read_seq = alignmt.query_sequence

**Sequence format:** Full read sequence, no truncation
- Written as-is from BAM (same orientation as stored in BAM)
- Full length (depends on original read length, typically 100bp for Illumina)

**Read selection logic:**
- Anchor read must map to region near insertion site (within anchor_chrm:insertion_pos ± 100bp)
- Mate must map to different chromosome OR same chromosome but > disc_threshold distance
  - global_values.DISC_THRESHOLD = 2000 (line 106)
  - Checked in is_discordant() x_clip_disc_filter.py:278-287
- Mate must have MAPQ >= global_values.MINIMUM_DISC_MAPQ (20, line 39)

### B.3 The bwa Command for Disc Reads

Function: `realign_disc_reads` (bwa_align.py:120-123)

```
bwa mem -t {n_jobs} -o {sf_out_sam} {sf_ref} {sf_reads}
```
- No special flags
- Default bwa mem behavior
- n_jobs: thread count
- sf_ref: repeat consensus reference (or repeat copy library with flanks)
- sf_reads: FASTA file of mate reads

### B.4 Parsing: parse_disc_algnmt_consensus

Location: x_clip_disc_filter.py:3209-3295

**Query-name fields read (from parsed ~-delimited query_name):**
- Field [-8]: s_anchor_lclip (anchor left clip: 0 or 1, line 3231)
- Field [-7]: s_anchor_rclip (anchor right clip: 0 or 1, line 3232)
- Field [-6]: s_anchor_rc (anchor read RC: 0 or 1, line 3233)
- Field [-5]: s_anchor_mate_rc (mate read RC: 0 or 1, line 3234)
- Field [-4]: anchor_map_pos (left-most mapped base of anchor, line 3235)
- Field [-3]: ins_chrm (insertion chromosome, line 3236)
- Field [-2]: ins_pos (insertion position, line 3237)
- Field [-1]: sample_id (sample index, line 3238)

**SAM fields read per alignment:**
- FLAG: Checked for unmapped status (line 3221)
- CIGAR: Used via algnmt.cigar to compute qualified alignment (line 3224)
- POS (reference_start): Base position on repeat consensus (line 3251)
- MAPQ: Not explicitly gated
- Query sequence: algnmt.query_sequence for polyA detection (line 3239)
- is_reverse: Alignment orientation (line 3261)

**Qualification check (bmapped_cutoff):**
- Same as clip reads: is_clipped_part_qualified_algnmt (x_clip_disc_filter.py:2731-2764)
- Default bmapped_cutoff: 0.7 (global_values.MIN_DISC_MAPPED_RATIO=0.7, line 119)
- At least 70% of disc read must map

**Position computation (pos_in_consensus):**
- Line 3251: `pos_in_consensus = algnmt.reference_start + len(read_seq) / 2`
- Middle map position: start of alignment + half the read length
- This approximates the "middle" position where the read maps on consensus

**Returns:**
- m_disc_pos: dict[ins_chrm][ins_pos] = [(anchor_map_pos, b_anchor_lclip, b_anchor_rclip, pos_in_consensus, b_same_dir_xor, b_anchor_rc), ...]
  - b_same_dir_xor: whether alignment orientation differs from anchor orientation (line 3282)
  - Used to detect inversions

- m_disc_sample: dict[ins_chrm][ins_pos][sample_id] = count
  - Counts how many disc reads per sample per insertion site

- m_polyA: dict[ins_chrm][ins_pos] = count
  - Counts polyA/T signals in disc reads (line 3241)
  - Uses is_consecutive_polyA_T() check

---

## C. SAM Contract Summary

### temp_clip.sam Required Fields

**Per-record SAM fields REQUIRED by parser:**
1. QNAME (query name): Must preserve all 7 fields separated by ~
   - [0] chrm, [1] map_pos, [2] FLAG, [3] is_rc, [4] insertion_pos, [5] n_cnt_clip, [6] sample_id
2. FLAG: Must indicate mapped (not unmapped flag)
3. RNAME (reference name): Repeat consensus ID in format: chrm~start~end~subfamily
4. POS (1-based): Position on repeat consensus (used as reference_start in pysam, 0-based)
5. MAPQ: Not strictly gated but affects polyA detection
6. CIGAR: Critical - must accurately represent alignment to compute mapped base count
   - Parsed to verify >= 65% (MIN_CLIP_MAPPED_RATIO) of clipped read is mapped
7. SEQ: Query sequence for polyA detection
8. QUAL: Quality scores (though not explicitly used by parser)

**Hard gates:**
- Unmapped reads (FLAG bit 4) skipped
- Alignment must satisfy: mapped_bases >= 0.65 * total_length
- If both-end clipped and both > 8bp: skipped

**Metadata preservation:**
- Query name fields encode all necessary collection metadata:
  - Original clip position (map_pos, FLAG, is_rc)
  - Insertion candidate site (insertion_pos)
  - Read counter within position (n_cnt_clip)
  - Sample identity (sample_id) - **critical for per-sample support tracking**
- Sample ID provenance: Added sequentially during BAM-list iteration (0 = first BAM, 1 = second, etc.)

---

### temp_disc.sam Required Fields

**Per-record SAM fields REQUIRED by parser:**
1. QNAME (query name): Must preserve all 9 fields separated by ~
   - [0] read_name, [1] is_first, [2] b_lclip, [3] b_rclip, [4] s_anchor_rc, [5] s_anchor_mate_rc, [6] anchor_map_pos, [7] ins_chrm, [8] ins_pos, [9] sample_id
   - Note: Fields accessed by negative indexing from end: [-8] to [-1]
2. FLAG: Must indicate mapped (not unmapped flag)
3. RNAME: Repeat consensus ID
4. POS: Position on repeat consensus
5. MAPQ: Not strictly gated
6. CIGAR: Critical - must verify >= 70% (MIN_DISC_MAPPED_RATIO) of read is mapped
7. SEQ: Query sequence for polyA detection and pos_in_consensus calculation (used for middle position)
8. QUAL: Quality scores

**Hard gates:**
- Unmapped reads skipped
- Alignment must satisfy: mapped_bases >= 0.70 * total_length

**Metadata preservation:**
- Query name encodes anchor read information:
  - Original read pair ID (read_name) - used to connect back to paired reads
  - Read orientation in pair (is_first)
  - Anchor read clipping status (b_lclip, b_rclip)
  - Anchor read orientations (s_anchor_rc, s_anchor_mate_rc) - used for inversion detection
  - Anchor position on genome (anchor_map_pos)
  - Insertion candidate site (ins_chrm, ins_pos)
  - Sample identity (sample_id) - **critical for per-sample support tracking**

---

## D. Key Implementation Challenges

### For C++ Scanner Emitting Byte-Compatible Alignments:

1. **Exact read-name format preservation**: Query names MUST be preserved exactly with correct field count and ~ separators, as negative-index parsing is used for disc reads
   
2. **CIGAR accuracy**: Must accurately represent clipped vs matched vs inserted/deleted bases
   - For clipped reads: typically CIGAR like "4M" (full match) to "2M2I5M" (with insertions)
   - Fraction of mapped bases is computed as: matched_bases / (all_bases_except_deletions)
   - Off-by-one or incorrect CIGAR will cause false filter-outs

3. **Reference name format**: For temp_clip.sam, RNAME must be in format "chrm~start~end~subfamily" to enable repeat copy position tracking

4. **Sequence accuracy**: Full sequence must be present (not truncated) especially for disc reads where length/2 is used to compute consensus position

5. **Sample ID tracking**: Sample IDs are assigned sequentially (0, 1, 2, ...) based on BAM file order in input list and MUST be preserved through the pipeline

6. **Position computation edge cases**:
   - Clipped read consensus position depends on: whether read is RC AND whether left/right clipped (line 3055)
   - Disc read middle position: start + len(seq)/2 uses floating point division
   - Both require accurate reference_start (POS field) and sequence length

7. **PolyA detection requirements**: 
   - Clipped reads: checks consecutive 5+ A or T (is_consecutive_polyA function)
   - Disc reads: checks full sequence for consecutive polyA/T
   - Quality of MAPQ may affect this

8. **Unmappable reads**: Records with FLAG bit 4 set are completely skipped - ensure properly mapped reads don't have this bit set

### Global Constants Referenced (from global_values.py):

| Constant | Value | Purpose |
|----------|-------|---------|
| BWA_REALIGN_CUTOFF | 9 | Min seed length (-T/-k) for non-polyA clipped reads |
| MINIMUM_POLYA_CLIP | 7 | Min seed length for polyA clipped reads, min clip length threshold |
| BWA_SEED_MEDIUM_FREQ | 70 | -c parameter for polyA clipped reads |
| MIN_CLIP_MAPPED_RATIO | 0.65 | Minimum mapped fraction for clipped reads to pass filter |
| MIN_DISC_MAPPED_RATIO | 0.7 | Minimum mapped fraction for disc reads to pass filter |
| MAX_CLIP_CLIP_LEN | 8 | If both ends clipped > this, skip |
| SEPERATOR | ~ | Field delimiter in query names |
| FLAG_LEFT_CLIP | L | Constant for left-clipped reads |
| FLAG_RIGHT_CLIP | R | Constant for right-clipped reads |
| DISC_THRESHOLD | 2000 | Minimum insert size to call discordant |
| MINIMUM_DISC_MAPQ | 20 | Minimum MAPQ for collecting disc reads (collection phase, not parsing) |
| MINIMAL_TRANSDUCT_MAPQ | 10 | MAPQ threshold for transduction detection |
| N_MIN_A_T | 5 | Minimum consecutive A/T for polyA detection |

