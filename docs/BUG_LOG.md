# xTea2 bug log

This log covers bugs found in xTea (`master`) that change calls, allele fractions or genotypes,
and how xTea2 (branch `xTea2`) fixes them. Each entry gives the mechanism, the measured effect,
the fix commit, the test that pins it and where the full evidence lives. Open issues are listed
at the end.

IDs B1–B10 match the running engineering log in `docs/clipdisc_contract_spec.md` §4, which also
holds design changes (DC*), dead code (D*) and comment fixes (B2, Q*). IDs from B11 are new here.
All bugs below were checked against `master` and are present there. File paths and line
numbers point at `master` as of `46520ab`; view them with `git show master:<path>`.

## Changelog

Newest first.

- **2026-11-13** (branch `lowq-filter`): B16 fixed: the low-MAPQ clip filter counts primary
  alignments only (`2bc904d`). 1KG: +56 TP, +20 FP, -1 FP; HG002 and the mosaic HapMap BAMs
  unchanged. Exact right-clip end and alt-aware counting were measured and not adopted (no
  effect on final calls). Open issue resolved.
- **2026-11-13** (branch `vaf-site-v2`):
  - Per-site VAF model (`x_vaf_site`) rebased onto xTea2. Its PE yields now match per-pair
    counting (reference term `E[n_C]`, not twice it); the spanning-pair term was dropped; defaults
    refit on HG002 NovaSeq 30x (e_PE 1.28, l_PE 50; was 1.33, 137) (`d091488`). Germline call
    sets, GT and GQ are unchanged on all 8 datasets; only `VAF`, `VAF_CI` and the
    `##xtea2_vaf_params` header change.
  - Mosaic AF_PE ceiling set per sample: 0.60 x the median AF_PE of the sample's germline het
    calls of that family, falling back to 0.35 with fewer than 20 of them (`91059b3`).
  - Alu-spanning pairs: modeled, no effect on Alu VAF bias, left out.
- **2026-11-13** (branch `gq-calib`): GQ calibrated with a per-call error floor (0.001, GQ ≤ 30)
  fit on HG002 NovaSeq 30x (GIAB Q100 labels) and HG00438 (`9f288cc`). GT unchanged. Most
  high-GQ errors against the xTea-paper HG002 labels were label errors (open issue resolved).
- **2026-11-13** (branch `counting-v2`):
  - B15 fixed: the pair checks on right-clipped reads use the alignment start (`6551542`).
  - Genotyper defaults refit on the B15 counts; only the PE biases move, by 0.001 (`320f610`).
  - Supplementary junction reads in AF_SR: measured, implemented and reverted, not merged
    (new open issue).
- **2026-11-12:** B14 (pairs counted once for AF_PE), with AF_PE cutoffs converted.

## Summary

| ID | Bug | Master code | Effect | Fix |
|---|---|---|---|---|
| B3 | Unseeded insert-size sampling | `xtea/x_basic_info.py:472,478` | 1–2.5% of calls change between identical runs | `216b84d` |
| B5 | Cross-family merge keyed by exact position | `xtea/x_vcf_merger.py:740,743` | Wrong family or duplicate call at conflicts | `381f41f` |
| B8 | Arbitrary repeat copy picked from an unordered set | `xtea/x_annotation.py:610-615` | Young-repeat filter keeps or drops calls arbitrarily | `65fca02` |
| B10 | Mate FASTA written per alignment, not per read | `xtea/x_alignments.py:381-396` | 4.8% of Alu discordant support double counted | `ee0742d` |
| B11 | Multi-BAM feature merge reads a stale loop variable | `xtea/x_genotype_feature.py:106-114` | One unrelated site's counts added to every site's AF | `dca7363` |
| B12 | Both clip sides' reference reads subtracted | `xtea/x_genotype_feature.py:327,331` | Split-read AF inflated at right-dominant sites | `b25879c` |
| B13 | Mosaic writer appends a second repeat column | `xtea/x_mosaic_calling.py:89,127,165` | Mosaic VCF writes the gene annotation as the sample's GT | `d71a72d` |
| B14 | Concordant support counted per read, not per pair | `xtea/x_genotype_feature.py:182-185,305-324` | Read-pair AF (AF_PE) biased low, by library | `ec1b6f0` |
| B15 | Right-clipped reads use the clip position in the pair checks | `xtea/x_genotype_feature.py:261-265,305,312,316` | Pairs misclassified near the site; pair counts shift slightly; one final call changes over the 8 validation sets | `6551542` |
| B16 | Low-MAPQ clip filter counts supplementary alignments | `xtea/x_clip_disc_filter.py:148-157` | With `bwa mem -Y` input, true insertions dropped at the rate of false ones | `2bc904d` |

## Fixed bugs

### B3. Unseeded insert-size sampling makes calls nondeterministic

- **Master code:** `xtea/x_basic_info.py:456-480` (`_random_slct_site`). The unseeded draws
  are `random.randint` at lines 472 (chromosome) and 478 (position); callers at 119, 149
  and 209.
- **Mechanism:** 3000 random sites are drawn without a seed to estimate mean insert size, SD and
  coverage. `max_is` (mean + 3 SD) therefore moves a few bp from run to run and tips borderline
  SVA and Alu calls in the consensus filter.
- **Effect:** identical input gave different call sets. Master runs shared only 97.5–98.9% of
  call positions, and 1 L1 and 8 Alu truth calls flickered in and out.
- **Fix:** a scoped `random.Random(42)` for the draw (`BASIC_INFO_RANDOM_SEED`; the
  `XTEA_BASIC_INFO_SEED` environment variable overrides it, and `none` restores master's
  behavior). The global RNG is untouched.
- **Evidence:** three xTea2 replicates give byte-identical call sets.
- **Detail:** contract spec §4, B3.

### B5. Cross-family merge keyed by exact position

- **Master code:** `xtea/x_vcf_merger.py`, class `SampleVCFMerger` (line 605).
  `load_in_sites_from_vcf` stores L1 calls by position (line 706); `append_sites_from_vcf`
  (lines 707-743) adds SVA, then Alu, then HERV. Line 740 overwrites by position when the new
  call wins; line 743 also stores the new call when the old one wins (`b_keep_first` is True),
  so the loser is written too. The quality comparison is `cmp_quality_of_two_ins` (line 746).
- **Mechanism:** loci were stored in a dict keyed by exact VCF position. A conflicting call at
  the same position overwrote the winner, so the family loaded last (always Alu) won. A losing
  call at a different position was written as well, giving a duplicate.
- **Effect:** 5 HG002 WGS loci resolved to the wrong family.
- **Fix:** `LocusTable` (`xtea/locus_table.py`), with list-based storage per chromosome.
  Reproduces the old 1529-call golden byte for byte except at those 5 loci.
- **Test:** `tests/test_locus_table.py`.
- **Later design change (not a bug fix):** `11bdca6` replaced master's priority rule
  (L1 > SVA > Alu) with an evidence score and records the losing calls (`OTHER_FAMILIES` INFO,
  `merged/all_calls.tsv`). It recovered 2 true HG002 SVA insertions.

### B8. Arbitrary repeat copy decides the young-repeat filter

- **Master code:** `xtea/x_annotation.py:605-615` (`is_within_repeat_region_interval_tree`):
  line 610 queries the interval tree for a set of hits, and lines 613-615 return the start of
  whichever comes out of the set first. The consumer is `xtea/x_post_filter.py`: lines 78-79
  and 139-140 take that copy's divergence (`get_div_subfamily`), and lines 90, 97, 155 and 163
  drop the call when it is in a repeat with divergence below `REP_DIVERGENT_CUTOFF`. B8b is the
  same pattern in `xtea/x_gene_annotation.py:206-209` (`query_by_position`), which sets the
  order of the gene hits reported by `get_gene_id_type` (line 213).
- **Mechanism:** when several (boundary-extended) repeat copies overlap a position, the lookup
  returned the first element of a Python set, which is arbitrary. That copy's divergence gates
  the filter that drops calls inside young repeats (divergence < 15%), and neighboring copies
  often straddle 15% (a young AluY beside an old AluJ).
- **Effect:** the verdict depends on which copy is picked at 2.95% of Alu and 5.79% of L1
  positions.
- **Fix:** the minimum-divergence overlapping copy, ties broken leftmost, inside a new static
  sorted-array interval index that replaced the vendored interval tree (about 100x faster to
  build). Callers that only need "is this in a repeat" get a boolean, so they cannot misuse a
  position (B8c). B8b (gene annotation order) was fixed with it.
- **Test:** `tests/test_interval_index.py`.
- **Detail:** contract spec §4, B8, B8b, B8c and P1-EXECUTED.

### B10. Discordant support double counted

- **Master code:** `xtea/x_alignments.py:368-397` (`extract_mate_reads_of_region`). The loop
  at lines 381-396 writes a record for every alignment of a wanted mate, with no check for
  supplementaries or names already written. Bins admit mates up to 100 bp past either edge
  (`_parse_read_name_pos_of_region`, line 351), so neighboring bins can both emit a mate. The
  records are counted by `parse_disc_algnmt_consensus` (`xtea/x_clip_disc_filter.py:3209`).
- **Mechanism:** the mate FASTA got one record per mate alignment, not per mate read. A
  soft-clipped supplementary alignment of the mate produced an identical second record, and
  overlapping bins emitted mates near a bin edge twice.
- **Effect:** 4.8% of Alu discordant records were duplicates, concentrated in pericentromeric
  chr16 and chr17.
- **Fix:** one record per mate read.
- **Test:** `tests/test_disc_mate_extract.py`.
- **Detail:** contract spec §4, B10.

### B11. Multi-BAM feature merge adds one unrelated site to every site

- **Master code:** `xtea/x_genotype_feature.py:51` (`call_genotype`). The merge loop at lines
  101-114 iterates `ins_chrm`/`ins_pos`; line 105 uses them, but lines 106-114 read
  `m_info[chrm][pos]`, where `chrm` and `pos` are left over from the candidate-list loop at
  lines 69-70 (the last candidate read).
- **Mechanism:** field 0 was indexed by the loop's site; fields 1–9 read `m_info[chrm][pos]`,
  where `chrm` and `pos` were leftovers from an earlier loop. Every site seen in 2+ BAMs got
  one fixed unrelated site's counts added to its full-map, discordant and concordant counts,
  which are the AF denominators.
- **Effect:** wrong AF for every multi-BAM site. On a 2-BAM HapMap mixture, the SVA clip AF read
  about 0.005 instead of 0.14, so all 59 SVA truth records that reached the mosaic filter were
  rejected (SVA 0/179 against truth).
- **Fix:** index every field by the loop's own site.
- **Test:** `tests/test_gntp_feature_merge.py`.
- **Follow-on:** an Alu-specific mosaic ceiling fitted to the corrupted numbers was withdrawn in
  `0aede86` (see the `MOSAIC_AF_BOUNDS` comment in `x_mosaic_calling.py`).

### B12. Split-read AF inflated at right-dominant sites

- **Master code:** `xtea/x_genotype_feature.py:326-332` (end of `collect_features_one_site`).
  Line 327 subtracts `n_l_full_map` unconditionally. When the right clip dominates, line 331
  also subtracts `n_r_full_map` from the already-reduced count, so both sides are removed.
- **Mechanism:** the split-read AF uses the dominant clip side only, so only that side's
  breakpoint-adjacent reference reads should be removed. The right-dominant branch subtracted
  the right side's reads on top of the left side's already removed.
- **Effect:** on HG002, het right-dominant calls read SR AF 0.556 against 0.466 for
  left-dominant calls at equal clip depth.
- **Fix:** subtract exactly one side. The genotyper defaults were refit afterwards (`b9249b8`).

### B13. Mosaic VCF writes the gene annotation as the sample's GT

- **Master code:** `xtea/x_mosaic_calling.py:89, 127, 165` (in `call_mosaic_L1_from_bulk`,
  `call_mosaic_Alu_from_bulk`, `call_mosaic_SVA_from_bulk`): `rcd[-1].rstrip() + "\t" + s_in_rep`.
  The VCF writer reads the genotype at `xtea/x_gvcf.py:289` (`s_gntp=l_fields[50]`, in
  `_cvt_to_vcf_rcd`, line 238).
- **Mechanism:** the mosaic writers appended a second repeat-copy label to rows that already
  end with one (field 48). That pushed the gene annotation to field 50, where the VCF writer
  reads the genotype.
- **Effect:** mosaic VCFs carried the repeat label as `GENE_INFO` and the gene string as GT.
- **Fix:** write rows unchanged. The VCF writer (`x_gvcf._cvt_to_vcf_rcd`) now accepts only the
  known row widths and checks that field 50 looks like a genotype, so a future extra column
  fails loudly.
- **Test:** `tests/test_gvcf_row_width.py` (8 of its 12 tests fail on the old code).

### B14. Concordant support counted per read, not per pair (AF_PE)

- **Master code:** `xtea/x_genotype_feature.py:141-346` (`collect_features_one_site`), which
  feeds `AF_DISC` / `AF_CONCORDNT` (fields 39 and 40), so AF_PE = disc / (disc + concordant).
  - Lines 183-185: the secondary/supplementary skip is commented out, so every record of a
    read is seen.
  - Lines 312-313: `n_disc_pairs += 1` per record.
  - Line 316: every non-discordant record goes to `l_check_concord`; lines 319-324 count each
    one, so both mates of a spanning pair are counted.
  - Line 320: the only exclusion is `m_clip_qname`, filled at lines 242 and 276 from
    soft-clipped primaries (CIGAR op 4) only, so hard-clipped supplementaries pass.
  - Cutoffs on this scale: `xtea/x_post_filter.py:1078-1084` (germline 0.075),
    `:1093` (transduction 0.08) and `xtea/x_mosaic_calling.py:173-179` (mosaic 0.35).
- **Mechanism:**
  - A reference pair spanning the site has both mates inside the fetch window, so it was
    counted twice.
  - A hard-clipped supplementary piece of an insertion read (primary on a TE copy elsewhere,
    mate fields pointing at the anchor) also passed as concordant, so some alt templates were
    counted as reference. The clip exclusion only recognized soft clips.
  - Discordant counts were inflated slightly (about 1.05x).
- **Effect:** AF_PE biased low by a library-dependent amount. HG002 NovaSeq het/hom medians
  were 0.38/0.75 (0.57/0.92 after the fix); hom medians varied 0.19 across libraries (0.07
  after). Per-read counting also overdispersed the counts, making GQ overconfident.
- **Fix:** per-pair sets keyed by query name; secondary records skipped; concordant evidence
  from non-supplementary records only; a template counted as discordant is never also
  concordant. Supplementaries still count on the discordant side, because in repeats a junction
  read's only record at the site can be a supplementary.
- **Downstream changes in the same commit:**
  - Genotyper defaults refit: PE het/hom 0.380/0.761 to 0.569/0.913.
  - Germline AF_PE floor 0.075 to 0.14 and transduction floor 0.08 to 0.15: master's values
    converted with 2f/(1+f) (concordant counts halve). A sweep found nothing better.
  - Mosaic AF_PE ceiling kept at 0.35 as an interim value (see open issues; replaced by a
    per-sample ceiling on 2026-11-13).
- **Validation (8 datasets against the previous commit):** identical call sets; GT errors 128
  to 120 (all 5 HG002 sets better; the three 1KG samples 0.1–0.2 points worse).
- **Test:** `tests/test_gntp_pair_counts.py`.
- **Detail:** `docs/supplemental/README_NOTES.md`, "Read-pair counting for AF_PE"; sweep data,
  figures and scripts in `docs/validation/2026-11-12_afpe_pair_counting/`.

### B15. Right-clipped reads use the clip position in the pair checks

- **Master code:** `xtea/x_genotype_feature.py:261-265` (`collect_features_one_site`). For a
  right-clipped read the loop adds the aligned lengths to `map_pos` to get the clip position.
  The same, now shifted, `map_pos` is then reused for the pair checks: `is_discordant` at line
  305, the `DFT_IS` window at line 312 and the concordant record appended at line 316.
- **Mechanism:** for a right-clipped read the pair checks measured from the clip position (the
  alignment end), not the alignment start. That shifted the insert-size test by the aligned
  length (about one read length), the `DFT_IS` window for discordant reads, and the span used to
  decide whether a concordant pair covers the site. A read ending just before the site could
  count as a spanning concordant pair, and a pair whose read was clipped after the site could be
  missed.
- **Effect:** small. Only the pair counts (fields 39 and 40, `AF_DISC` / `AF_CONCORDNT`)
  change; AF_SR is identical.
- **Downstream changes:** genotyper defaults refit, PE het/hom 0.569/0.913 to 0.570/0.914; SR
  and hom prior unchanged (`320f610`).
- **Validation (8 datasets against `0553edb`):**
  - TP unchanged (9823), FP 1986 to 1987.
  - The one gained call is an FP in HG00438 (Alu, chr6:32572129, HLA region). It loses 10
    discordant reads, so its `f_ef_disc` crosses the 0.075 cutoff.
  - In all, 7 records flip `f_ef_disc`, 5 the AF_PE gate and 6 the transduction gate; only
    that one changes a final call. The AF_PE cutoffs were not retuned.
  - GT errors 120 to 118. GQ calibration unchanged.
  - Mosaic HapMap mixtures at disc ceiling 0.35 (TP / FP, B14 to B15): BCM 1393/623 to
    1393/623, UW 189/3636 to 187/3517, WashU 688/204 to 687/205. The mean concordant count per
    mosaic record falls by 0.2-1.2.
- **Fix:** the clip position goes in its own variable (`i_rclip_pos`) and `map_pos` stays the
  alignment start.
- **Test:** `tests/test_gntp_pair_counts.py` (four tests that fail on the old code; one old test
  that asserted the buggy outcome was replaced).
- **Detail:** `docs/supplemental/README_NOTES.md`, "Right-clip pair checks"; tables and scripts
  in `docs/validation/2026-11-13_counting_v2/`.

### B16. Low-MAPQ clip filter counts supplementary alignments

- **Master code:** `xtea/x_clip_disc_filter.py:148-157` (`XClipDisc`, the per-candidate read
  loop). Reads with MAPQ < 5 and a soft clip at either end within 35 bp of the candidate
  (`NEARBY_LOW_Q_CLIP`) are counted as low-quality background. The count runs before the
  supplementary/secondary skip at line 164, so supplementary and secondary alignments count.
  The site is then dropped in `check_clip_alone_consistency2` when the low-MAPQ share of its
  clips exceeds `MAX_LOWQ_CLIP_RATIO` (0.65; `_is_high_lowq_clip`, line 2314).
- **Mechanism:** a junction read whose primary alignment lies in the inserted element leaves a
  supplementary alignment at the site, clipped exactly at the junction, usually with MAPQ 0.
  That read is evidence for the insertion, but it was counted as background against it.
  Whether it is seen depends on the aligner: `bwa mem -Y` soft-clips supplementaries, so they
  count; plain `bwa mem` hard-clips them, and a hard-clipped read has no op-4 clip and was
  never counted. The filter therefore behaved differently on the same sample aligned with and
  without `-Y`.
- **Effect:** on the 1KG CRAMs (`-Y`), supplementaries were 40-44% of the counted (read,
  candidate) pairs; on HG002 (no `-Y`) 0%. On the three 1KG sets the filter removed 60 TP and
  25 FP final calls (1.5% of TPs, 1.9% of FPs), so it did not discriminate. Counting primary
  alignments only, it removes 4 TP and 6 FP there (0.1% vs 0.45%). Against the old filter:
  HG00438 +12 TP / +4 FP, HG02257 +20 TP / +8 FP / -1 FP, HG02622 +24 TP / +8 FP. HG002 is
  identical by construction; the mosaic HapMap BAMs have no soft-clipped supplementaries or
  secondaries (one BAM per center checked), so they are unchanged too.
- **Not changed:** the right-clip position estimate (start + read length - clip length) differs
  from the exact alignment end in 0.5-1% of counted pairs (indels, hard clips) and changed no
  final call in any dataset where it was run. Counting only reads with an alternative hit in
  the primary assembly (alt-contig suspect) changed one call (+1 TP, HG02622). Both kept as in
  master.
- **Fix:** `count_as_lowq_clip` in `src/csrc/cns_emit.hpp` (both `emit_cns` and
  `emit_cns_multi`) and `low_mapq_clip_anchor` in `x_clip_disc_filter.py` skip flag 0x900.
- **Test:** `tests/test_lowq_clip_count.py` (scanner single-file and manifest modes, and the
  Python port; supplementary/secondary reads at the junction must not count). The scanner tests
  fail on the old code (5 counted, 2 expected).
- **Detail:** `docs/supplemental/README_NOTES.md`, "Low-MAPQ clip filter"; tables, figure and
  scripts in `docs/validation/2026-11-13_lowq_filter/`.

## Open issues

| Issue | Status |
|---|---|
| Mosaic AF_PE ceiling is a fixed 0.35 | **Resolved 2026-11-13** (`91059b3`). The ceiling is now 0.60 x the sample's germline het AF_PE median per family (0.28-0.31 on the HapMap mixtures), with 0.35 as the logged fallback below 20 het calls. Left-one-run-out F1 against 0.35: BCM 0.597 vs 0.592, UW 0.065 vs 0.059, WashU 0.379 vs 0.383; calls at HG005 germline sites 23 to 3 in total. Most of the gain is from the lower level, not the per-sample adaptation (a fixed 0.28-0.30 does about as well). Detail: `docs/supplemental/README_NOTES.md`, "Per-sample mosaic ceiling"; `docs/validation/2026-11-13_mosaic_v2/`. |
| GQ overconfident above 30 | **Resolved 2026-11-13** (`9f288cc`). 35 of the 37 HG002 errors at GQ ≥ 30 were errors in the xTea-paper truth labels (GIAB Q100 v5.0q agrees with xTea). GQ now has a per-call error floor of 0.001 (GQ ≤ 30). Held out, Q100 labels: top bin 0.11% observed vs 0.13% predicted; GT unchanged. Remaining top-bin errors are complex or nested sites, not modeled separately. Detail: `docs/supplemental/README_NOTES.md`, "GQ calibration"; `docs/validation/2026-11-13_gq_calibration/`. |
| Supplementary junction reads not counted in AF_SR | **Measured, not merged.** `x_genotype_feature.py` counts only soft-clipped primaries as split reads: CIGAR op 4 only (master lines 233, 259), supplementaries skipped (234, 267-268). In repeats a junction read's primary often maps to a TE copy elsewhere, so the piece at the site is a supplementary (hard-clipped with bwa defaults, soft-clipped with `bwa -Y`). Measured before any fix on truth-matched calls: 92% of HG002 NovaSeq 30x calls (1158 of 1260) have such reads, 5.6 per call on average. Counting them would move het AF_SR medians by family: Alu 0.433 to 0.542, L1 0.467 to 0.519, SVA 0.491 to 0.560; hom by 0-0.02. Other HG002 sets: +0.08 to +0.10. The 1KG CRAMs have only soft-clipped supplementaries and show the same shift (HG00438 Alu het 0.467 to 0.536, HG02622 0.458 to 0.543; pooled +0.07 to +0.08). A fix plus refit was tried (`6325b4d`, reverted). It was not merged for three reasons. The Alu polyA-dominant gate (`is_polyA_dominant_one_side`) treats the extra clips as non-polyA and dropped 44 TP and 13 FP Alu one_half calls. GT errors went from 120 to 123, so no gain. Het medians overshot 0.5 (0.46-0.60). Next step: double-check that the supplementary counting doesn't overcount templates before revisiting. Patch, tests and measurements are in `docs/validation/2026-11-13_counting_v2/` (`supp_reads_not_merged/`, `measure/`, `supp_reads_arm/`). |
| Soft-clipped supplementaries skipped for pair counts | In `collect_features_one_site` a soft-clipped (`bwa -Y`) supplementary hits the supplementary `continue` in the clip branch before the pair logic. A hard-clipped one has no op-4 clip, so it reaches the pair logic and can count as discordant or concordant. AF_PE therefore depends slightly on the aligner's supplementary style. Not measured. |
| Low-MAPQ clip filter | **Resolved 2026-11-13** (`2bc904d`, B16). The TP loss on the 1KG CRAMs came from soft-clipped supplementary alignments (`bwa mem -Y`), not from alt contigs or the right-clip estimate. Only primary alignments count now. Exact end: no final call changed. Alt-aware counting: +1 TP. The 65% cutoff was not tuned. Detail: `docs/supplemental/README_NOTES.md`, "Low-MAPQ clip filter". |
| Per-site VAF model (`x_vaf_site`) | **Resolved 2026-11-13** (`d091488`). Reference term is per pair, refit on the B15 counts. Germline het VAF/truth 1.00-1.03, 95% CI coverage 0.93-0.98. HapMap mixtures, 5-6%: 0.60-1.10 of truth (UW lowest; 0.71 before), coverage 0.91-0.93. Low-VAF mixture calls (1-3%) remain biased low on short-insert, high-background libraries (UW). Detail: `docs/supplemental/README_NOTES.md`, "Per-site VAF model". |
| Alu-spanning reference pairs | **Checked 2026-11-13, not added.** Modeled as an expected alt-concordant count `2 E[max(0, F - 2R - S + 2c)]`; it adds 0-2 pairs for a 300-bp Alu at R = 151 and changed median Alu VAF/truth by at most 0.02. Germline Alu hets show no length trend. The low Alu VAF in the mixtures is not explained by it. Detail: README_NOTES, "Alu-spanning pairs". |
