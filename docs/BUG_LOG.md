# xTea2 bug log

This log covers bugs found in xTea (`master`) that change calls, allele fractions or genotypes,
and how xTea2 (branch `xTea2`) fixes them. Each entry gives the mechanism, the measured effect,
the fix commit, the test that pins it and where the full evidence lives. Open issues are listed
at the end.

IDs B1–B10 match the running engineering log in `docs/clipdisc_contract_spec.md` §4, which also
holds design changes (DC*), dead code (D*) and comment fixes (B2, Q*). IDs from B11 are new here.
All bugs below were checked against `master` and are present there. File paths and line
numbers point at `master` as of `46520ab`; view them with `git show master:<path>`.

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
  - Mosaic AF_PE ceiling kept at 0.35 as an interim value (see open issues).
- **Validation (8 datasets against the previous commit):** identical call sets; GT errors 128
  to 120 (all 5 HG002 sets better; the three 1KG samples 0.1–0.2 points worse).
- **Test:** `tests/test_gntp_pair_counts.py`.
- **Detail:** `docs/supplemental/README_NOTES.md`, "Read-pair counting for AF_PE"; sweep data,
  figures and scripts in `docs/validation/2026-11-12_afpe_pair_counting/`.

## Open issues

| Issue | Status |
|---|---|
| Mosaic AF_PE ceiling is a fixed 0.35 | **Interim.** On the B14 counts it costs 12% of true calls on the short-insert UW library. Replace with a per-sample ceiling, e.g. a fraction of the sample's own germline het AF_PE median. |
| GQ overconfident above 30 | Observed GT error 0.4–0.6% against predicted near 0. B14 narrowed the gap; reads at a site are overdispersed relative to the binomial model. |
| Hard-clipped supplementary junction reads not counted as split-read evidence | Master `xtea/x_genotype_feature.py:233,259` test CIGAR op 4 (soft clip) only, and 234/267 skip supplementaries. Clip counts only see soft clips, so AF_SR is likely undercounted on the alt side, worst at hom Alu sites. Measure on HG002 before changing; needs a genotyper refit. |
| Right-clipped reads use the clip position in the pair checks | Master `xtea/x_genotype_feature.py:261-265` shift `map_pos` to the clip position for right-clipped reads; it is then reused at 305, 312 and 316 for the discordant/concordant checks. Still present on xTea2. Effect not measured. |
| Low-MAPQ clip filter | Which reads count (exact end vs estimate) and the 65% cutoff on alt references need review. |
| Per-site VAF model (`x_vaf_site`) | Its reference-pair term assumed per-read counting (2 x expected pairs); rebase and refit on the B14 counts. |
| Alu-spanning reference pairs | Pairs that span a short Alu insertion count as reference under any counting scheme; a geometric correction could account for it. |
