# xTea2 — supplemental notes

Background, rationale and developer-facing options moved out of the main [README](../../README.md)
to keep it focused on installing and running xTea2. Nothing here is needed for a standard run.

---

## Install details

### Genotyping model

Genotyping (`xtea/x_genotype_lik.py`) is a closed-form two-genotype (0/1, 1/1) likelihood:
a beta-binomial (overdispersion `rho` = 0.05) on split reads (`AF_CLIP` of `AF_CLIP+AF_FMAP`)
plus one on read pairs (`AF_DISC` of `AF_DISC+AF_CONCORDNT`). The expected alt fraction of each
read type under each genotype, and the hom fraction, are fitted per sample by EM over all of the
sample's calls, pooled across families. It runs once in the shared tail, over the merged
all-family call set. On HG002, fitting each family separately gave the same accuracy as the
pooled fit, so pooled is the default. If there are fewer than 50 informative calls, or EM gives a
degenerate fit (het bias not at least 0.10 below hom bias, a het bias or hom prior on its bound),
the HG002 fit is used instead: SR 0.441/0.958, PE 0.570/0.914, hom prior 0.308 (fit on counts
with the AF_FMAP fix, each read pair counted once on both sides of AF_PE, and the right-clip
pair checks fixed). `GQ` is the phred-scaled posterior error with a per-call error floor of
0.001 (`GQ_ERR_FLOOR`), so it tops out at 30; see "GQ calibration" below. The parameters
used, including `rho` and `gq_err_floor`, are recorded in the VCF header line
`##xtea2_genotype_params`.

It replaced an ONNX random forest (75% GT accuracy on HG002, 10% hom recall) and with it the
onnxruntime, scikit-learn and pandas dependencies. Python 3.10 remains the floor.

**Validation** (end-to-end runs of this genotyper; GT accuracy over truth-matched calls with a
het/hom label; the old ONNX model's accuracy on the v2.0.0 runs of the same inputs for reference):

| Dataset | Matched calls | GT accuracy (old) | Het correct | Hom correct |
|---|---|---|---|---|
| HG002 NovaSeq 30x | 1284 | **98.9%** (75.3%) | 99.1% | 98.3% |
| HG002 HiSeq 30x | 1272 | **98.3%** (74.8%) | 98.6% | 97.7% |
| HG002 downsampled 20x | 1210 | **98.3%** (72.8%) | 98.6% | 97.4% |
| HG002 downsampled 15x | 1109 | **97.4%** (71.5%) | 97.8% | 96.4% |
| HG002 downsampled 10x | 830 | **95.7%** (64.9%) | 96.0% | 95.1% |
| HG00438 (1KG 30x CRAM)\* | 1065 | **99.7%** (86.8%) | 99.9% | 99.4% |
| HG02257 (1KG 30x CRAM)\* | 1364 | **99.8%** (93.7%) | 99.8% | 99.5% |
| HG02622 (1KG 30x CRAM)\* | 1375 | **99.9%** (93.2%) | 99.9% | 99.6% |

HG002 truth: the xTea paper's HG002 set (H1/H2 = het, H1,H2 = hom). \*Truth from HPRC release 2
PAV and dipcall assembly calls, scored on tier1 sites where PAV and dipcall agree (tier3-only sites
lean hom, especially SVA, and are less reliable). At 10x the limit is detection (52.7% recall),
not genotyping.

### Per-site VAF model

`VAF` and `VAF_CI` come from `xtea/x_vaf_site.py` (numpy and scipy only). The module docstring
holds the full derivation; in short:

- **Read yields.** For one site, the expected numbers of alt and ref reads of each type, per unit
  fragment-start density, follow from the read length R, the insert-size distribution (normal,
  from `tmp_basic/basic_cov_is_rlth_info.txt`), the insertion length S (`SVLEN`; the family
  default when missing: Alu 300, L1 6000, SVA 2000) and the TSD length t. With K = 70 (longest
  clip a flank-primary read carries), c = 2 and a = R - K:
  - PE counts are per pair (one count per query name, since B14); SR counts are per read.
  - discordant: `D_alt = 2 E[max(0, min(F - R, S + 2a - R))]` (one anchor read per pair, summed
    over the two junctions)
  - concordant at the alt allele: `C_alt = l_PE + t`
  - concordant at the ref allele: `C_ref = E[F]`, folded at mean + 3 sd
  - split reads: `A_SR = 2(K - c)` alt, `B_SR = 2(R - 44)` ref, `l_SR` fully mapped at alt
  - Before B14 the two reads of a concordant pair were counted separately, and the model doubled
    `C_ref` and the TSD term to match. The refit on per-pair counts moved the HG002 defaults from
    e_PE 1.33, l_PE 137 to e_PE 1.28, l_PE 50.
  - Pairs that span the whole insertion (`n_S = max(0, F - 2R - S + 2c)`) are left out. See
    "Alu-spanning pairs" below.
- **Alt fractions** at allele fraction v, with efficiencies e and background odds b:
  - `p_SR(v) = (v e_SR A + b_SR B) / (v e_SR A + b_SR B + (1 - v) B + v l_SR)`
  - `p_PE(v) = (v e_PE D_alt + b_PE C_ref) / (v e_PE D_alt + b_PE C_ref + v C_alt + (1 - v) C_ref)`
- **Estimate.** VAF maximizes the beta-binomial likelihood of `AF_CLIP | AF_CLIP + AF_FMAP` and
  `AF_DISC | AF_DISC + AF_CONCORDNT` over a grid (log-spaced below 0.1). `VAF_CI` is the 95%
  profile interval (log L within 1.92 of the maximum); its lower end is 0 when the likelihood
  stays flat down to v = 1e-4.
- **Calibration.** e_SR, e_PE, l_SR, l_PE and the hom fraction are fitted per sample by maximum
  likelihood over its calls, as an unlabeled mixture of v = 0.5 and v = 1 (rho 0.05). Fewer than
  50 informative calls, or a fit outside plausible bounds (e 0.05-4, l_SR 0-60, l_PE 0-300,
  hom fraction 0.02-0.95), falls back to the HG002 NovaSeq fit (e_SR 1.22, e_PE 1.28, l_SR 8.6,
  l_PE 50). Germline runs set b = 0. Across 7 genomes the fitted l_PE is 33-65. Libraries with few
  discordant pairs fit a low e_PE (HG002 2x250: 0.11), which is what keeps their het calls at 0.5.
- **Mosaic mode** (`--mosaic` in the pipeline, `--mosaic` in `xtea2-genotype`): no genotype,
  rho 0.02, calibration on the sample's germline calls (each type's germline post-filter list;
  the mosaic call set has none), and a background b measured per sample. `x_vaf_background.py`
  counts the genotyping features at 1000 random autosomal clip peaks of the BAM, fits b at v = 0,
  and multiplies by a context factor (SR 1.87, PE 1.36): MEI call sites carry more background
  than random positions. The factor was fitted on the 300x BCM HapMap mixture (refit on per-pair
  counts 2026-11-13). Without a BAM the BCM value is used (b_SR 0.017, b_PE 0.065). Background is
  library-specific: measured this way it is b_SR 0.017-0.057 and b_PE 0.064-0.16 across the four
  HapMap mixture BAMs (UW highest on both).
- **Header.** `##xtea2_vaf_params=<model=site_yield_betabinom,mode,source,reason,n_calib,
  read_len,mean_is,sd_is,K,c,rho,e_SR,e_PE,l_SR,l_PE,b_SR,b_PE,bg_source,hom_frac>`.
  `xtea2-genotype` rebuilds the read length and insert size from it.
- **Options.** `x_TEA_main.py --gntp_classify` takes `--basic-info`, `--vaf-mosaic`,
  `--vaf-calib FILE:RTYPE,...` and `--vaf-background B_SR,B_PE|FILE|auto`. `run_xtea.py` sets
  them. `xtea2-genotype` takes `--basic-info` or `--read-length/--insert-mean/--insert-sd`,
  `--mosaic` and `--vaf-background`.

Validation (2026-11-13, per-pair counts), on counts at truth sites computed with the pipeline's
own feature code. Tables, figures and scripts: `docs/validation/2026-11-13_mosaic_v2/`.

| data | truth | median VAF / truth | 95% CI coverage | before (per-read model) |
|---|---|---|---|---|
| 7 germline genomes (HG002 Nova 30x, HiSeq 30x, 10x, 2x250; HG00438, HG02257, HG02622) | het | 1.00-1.03 | 0.93-0.98 | 1.00-1.03, 0.92-0.98 |
| same | hom | 0.99-1.0 | 0.90-0.99 | 0.99-1.0, 0.90-0.98 |
| HapMap mixture BCM 300x (background factor fitted here) | 5-6% | 0.89 | 0.93 | 0.91, 0.94 |
| HapMap mixture WashU 100x | 5-6% | 1.10 | 0.92 | 1.10, 0.93 |
| HapMap mixture NYGC 100x | 5-6% | 0.80 | 0.93 | 0.78, 0.95 |
| HapMap mixture UW 100x (insert size 249) | 5-6% | 0.60 (33% estimated 0) | 0.91 | 0.71, 0.92 |
| HapMap mixture, 4 BAMs | 7-20% | 0.89-1.21 | 0.91-0.95 | |
| HapMap mixture, 4 BAMs | 2-3% | 0.0-0.95 (8-55% estimated 0) | 0.91-0.93 | 0.0-1.02, 0.92-0.95 |
| HapMap mixture, 4 BAMs | 1% | 0.01-0.62 (39-61% estimated 0) | 0.89-0.93 | 0.0-0.74, 0.91-0.94 |
| HapMap mixture, 4 BAMs | 0.25-0.5% | 0.02-0.17, NYGC 0.5% 2.09 (38-68% estimated 0) | 0.77-0.91 | |

Mixture rows use the HG002 default calibration and each BAM's measured background. In the 5-6%
bin, by family across the four BAMs: Alu 0.59-1.07, L1 0.73-1.23 and SVA 0.58-1.15. UW is the
worst for every family: its short inserts and high background push a third of its 5-6% calls to
an estimate of 0. The per-pair model is no better than the per-read one on the mixtures (UW 5-6%
moved from 0.71 to 0.60, the others by 0.02 or less); its gain is that the model now matches the
counts xTea2 writes.

With one constant background (the BCM value) instead, UW and NYGC 5-6% calls come out at 1.91x
and 1.18x. The old VAF (the AF along a per-sample bias line) put HG002 2x250 hets at 0.16, and
short Alu hets (< 200 bp) at 0.34 on all genomes; the site model gives about 0.5 on both. At 1%
and below the counts carry little beyond "below a few percent": the interval typically runs from
0 to 5-10%.

Calibration in a mixture: the mosaic step calibrates on the sample's germline calls, assuming they
sit at 0.5 and 1. In the HapMap mixture they sit at 0.42 and 0.84 (83.5% HG005), so the pipeline's
own estimates there are about 1.3-1.5x higher than the table above (which uses the HG002
calibration). HG005 germline hets in the mixtures come out at a median 0.43-0.46 (truth 0.4175).
In a tissue sample, germline calls are at 0.5 and 1 and the issue does not arise.

End to end (the 8 gtval germline datasets, genotype and gVCF steps rerun from the same inputs with
and without the model): the call set, GT, GQ and every field other than `VAF`, `VAF_CI` and the
`##xtea2_vaf_params` header are identical. `VAF` changes at 737-1812 calls per dataset. Median VAF
of 0/1 calls is 0.499-0.502 (1/1: 1.0), and the median 95% interval width is 0.39-0.52 for 0/1
and 0.13-0.25 for 1/1, about twice as wide as before. The model does not feed back into
genotyping.

**Alu-spanning pairs.** A pair with one read anchored on each side of a short insertion is
counted as concordant. Its expected count, `2 E[max(0, F - 2R - S + 2c)]`, was tried as an extra
alt-concordant term (l_PE refit). It needs F > S + 2R - 2c, so for a 300-bp Alu at R = 151 and
inserts of 250-420 bp it adds 0-2 pairs against l_PE 35-50 and C_ref 250-420. Adding it changed
the median Alu VAF/truth by at most 0.01 in the mixtures (bins of 1% and up), 0.006 by germline
sample and 0.022 by Alu length bin, and coverage by at most 0.024. An empirical insert-size
distribution instead of the normal one did no better. Germline Alu hets show no length trend
(0.94-1.04 by SVLEN bin). The term is left out: it does not explain the low Alu VAF in the
mixtures.

### Read length in the cns scan step

`xtea_scan` placed a low-MAPQ right clip at
`start + 100 - clip length` whatever the read length. It now gets the measured read length
(`--cns-read-length`), as master's Python step does. On the 5 HG002 runs (hg38 no-alt reference)
the calls are unchanged, except for 1 extra call on HiSeq 30x. On the 3 HPRC genomes (GRCh38 with
alts, decoys and HLA, which has more low-MAPQ reads) it removes 12-26 of 1582-2015 calls. Recall
drops by 0.5-0.9 points and PPV is unchanged (74.6-76.5%).

### Scan engine: how the binary is found, and when to rebuild

`xtea2` locates `xtea_scan` in this order: **`--scan-bin`** → **`$XTEA2_SCAN_BIN`** → the in-repo
build at `src/csrc/xtea_scan` → **`xtea_scan` on `PATH`**. If you installed the package outside the
repo, either export `XTEA2_SCAN_BIN` or copy `xtea_scan` somewhere on `PATH`.

The build verifies htslib at compile time and embeds an rpath, so the binary finds `libhts` at
runtime without `LD_LIBRARY_PATH`.

Rebuild after editing **any** file under `src/csrc/`, headers included: `realign.hpp`,
`cns_emit.hpp` and `td_gather.hpp` carry most of the engine and are `#include`d rather than compiled
separately. The Makefile lists them as prerequisites, so a plain `make` picks up a header-only edit.

### Test suite and the scan binary

About 30 tests cover the C++ engine and resolve it from `src/csrc/xtea_scan`, so without a build they
skip and a green run means less than it looks. A binary **older than its own sources** is a hard
failure rather than a skip, because that is the state in which a passing suite is testing code that
is not in your tree. See [tests/README.md](../../tests/README.md) for the environment variables that
point the suite at your own data.

---

## Removed seeding flags (`--nclip`, `--cr`, `--nd`)

There is no `--nclip`, `--cr` or `--nd`, and passing one is a hard error. The seeding cutoffs those
names refer to (`n_clip_cutoff`, the clipped-in-repeat minimum, and `n_disc_cutoff`) are **derived
from the input BAM's own coverage** at runtime, and scaled by `--purity` under `--tumor`. They were
effectively inert on legacy xTea's command line too, which overrode them from coverage unless
`--user` was passed. Some diagnostic banners still say `--nclip` when they mean the derived
`n_clip_cutoff`; that is a wart in the log text, not a flag you can pass.

---

## Transduction detection options

These control how transductions (a TE carrying flanking sequence from its source locus) are
detected. **Two are ON by default and change calls relative to legacy xTea.** The `--no-*` forms are
the rollback path.

```
  --td-decoy-structural      Tell a decoy consensus record from a flank record by the NAME'S
  --no-td-decoy-structural   STRUCTURE (a flank id is chrom~start~end~subfamily~rc, with
                             integer coordinates; a consensus id is bare) rather than by a
                             substring match on LINE1/SVA/ALU/HERV.
                             ** ON BY DEFAULT. **
                             Why: the substring test classifies 100% of hg38's 5,454 SVA flank
                             records as the SVA consensus decoy, so NO reference SVA copy can
                             be a transduction source.

  --td-decoy-concat          Concatenate the 45-record decoy consensus block into the
  --no-td-decoy-concat       transduction bwa reference when the packaged flank library lacks
                             it. ** ON BY DEFAULT. **
                             Why: hg38's SVA flank library carries NONE of the 45 records (both
                             L1 libraries carry all 45), so without this SVA has no decoy
                             channel on either the clip or the disc side. The block is rebuilt
                             from <rep_lib>/consensus/ and is byte-identical to the packaged
                             one; appending is BY NAME, so it is a structural no-op on L1 and
                             affects SVA only.

  --td-polya-cutoff N        Gate [A]'s polyA-signature threshold, decoupled from the derived
                             clip cutoff. Default: derive it as n_clip_cutoff/2-1. Supply it
                             only when a clip-support sweep needs to be single-variable.

  --td-gather-scan           Gather the transduction step's reads natively with xtea-scan.
  --no-td-gather-scan        ** ON BY DEFAULT. ** The native gather is 13-14x faster at WGS
                             scale and byte-identical to the Python path on real data, so
                             --no-td-gather-scan is a rollback switch, not a tuning knob.
                             NOTE: the Python path has a known bug with N>1 input BAMs (it
                             truncates one shared all_disc_pos per BAM, so only the last BAM's
                             positions survive). Do not combine it with multi-file input.

  --td-gather-driver         Gather every family's transduction reads once in the driver
                             instead of once per family. Implies --locus-universe. Structural;
                             not a measured speed win.
```

**A default run does not reproduce legacy xTea's output, and that is deliberate.** Both decoy flips
exist because without them **no reference-copy SVA transduction is detectable at all** on hg38. The
driver says so on every run:

```
[decoy] the transduction decoy test is STRUCTURAL (DEFAULT ...): this run does NOT reproduce BASELINE-0
[decoy-concat] the transduction bwa reference GAINS the decoy consensus block (DEFAULT): ...
```

To reproduce a call set frozen before them, pass **both** `--no-td-decoy-structural` and
`--no-td-decoy-concat`.

---

## Scan-path options

The C++ engine owns the clip and discordant passes outright: a single BAM pass collects both, and
there is no Python fallback for either. The cns_filter inputs are also collected by the engine by
default.

```
  --cns-emit-scan            Collect the cns_filter step's SAM inputs natively with xtea-scan
  --no-cns-emit-scan         (merged windows, threaded) instead of one indexed fetch and one
                             temp file per candidate in Python. ** ON BY DEFAULT. **
                             --no-cns-emit-scan is the A/B rollback; expect it to be far
                             slower at mosaic candidate counts.
  --scan-clip, --scan-disc   Accepted for compatibility; no-ops. The scan pass is always on.
  --scan-bin                 Path to the xtea_scan binary (see the resolution order above).
  --bwa                      Path to bwa. Default: bwa on PATH
  --scan-min-clip            Min in-repeat clipped reads to emit a clip candidate. Default: 2
  --scan-disc-threshold      Same-chromosome discordant distance threshold (bp) for emission.
                             Default: 2000 (deliberately below the effective cutoff so the
                             stream is a superset; the adapter re-applies the exact threshold)
  --scan-realign-ref         Reference for realigning clipped parts.
                             copies_with_flank (default, more sensitive) | consensus
```

---

## Advanced and diagnostic options

```
  --mate-in-rep-families     Which families' repeat annotations contribute to
                             mate_in_rep_clip_cnt (the # clipped reads at a site whose mate
                             lands in ANY of them). Default: alu,sva — the biologically
                             tightest set, since the value is consumed only by the SVA cns
                             path and SVA carries an Alu-derived domain. Independent of
                             --types by design.
  --col5-families            DEPRECATED alias for --mate-in-rep-families.
  --family-assign-early      Decide each locus's family once, right after cns_filter, instead
                             of resolving at the postF_hc merge. Same rule either way; no-op
                             for a single --types entry.
  --locus-universe           Also enumerate the union locus universe once in the driver.
                             Diagnostic: each family still computes its own, and the two must
                             be byte-identical.
  --global-barrier           Run cns_filter -> transduction as one global join instead of the
                             default staged form. Same output, slower. For A/B measurement.
  --adopt-existing           With --resume: accept scan_pass/cns_emit output from a run that
                             predates resume markers. It trusts files it cannot verify (a
                             truncated output looks complete), so confirm in that run's own
                             log that the phase finished before using it.
```

### `--region-mask` details

The mask is applied to **candidate positions only**, right after the scan pass. The discordant-read
evidence around surviving candidates is deliberately kept whole: dropping discordant reads by their
own position would strip flank support from real insertions just outside a masked region. See
`xtea/region_mask.py` for the measurements behind this choice.

The built-in names (`simple_repeats`, `segdup`, `centromeres`) currently resolve to GRCh38 files on
the Park lab's O2 cluster. Off O2, pass your own BED paths or `--region-mask none`.

---

## Multi-file pooling details

The scan engine enforces that all input files share one reference: before scanning it compares every
file's header against the first and fails if the sequence count, the `@SQ` names or their order, or
any sequence length differs.

Pooling is a flat sum: coverage adds, read length and insert size are averaged, and clip/disc counts
are summed, with no per-file weighting. Duplicate reads cannot be detected *across* files, since the
duplicate flag is per-file, so a fragment independently sequenced at two centers is counted twice.

---

## Read-pair counting for AF_PE (2026-11-12)

This is B14 in `docs/BUG_LOG.md`.

### The bug

AF_PE is the share of read pairs spanning the breakpoint that are discordant:
`AF_DISC / (AF_DISC + AF_CONCORDNT)`, fields 39 and 40 of the post-filter record and the
`AF_DISC` and `AF_CONCORDNT` INFO keys. Before this fix, `collect_features_one_site` in
`x_genotype_feature.py` counted reads rather than pairs on the concordant side:

- A reference pair spanning the site has both mates inside the fetch window, so it was counted
  twice.
- A hard-clipped supplementary piece of an insertion read (its primary on a TE copy elsewhere,
  its mate fields pointing at the anchor about 300 bp away) also passed as concordant, so some
  alt templates were counted as reference. The clip exclusion only recognized soft clips.
- Discordant counts were inflated slightly (about 1.05x).

AF_PE was therefore biased low at every site, and by a library-dependent amount. Now each pair
is counted once on both sides, keyed by query name: secondary records are skipped, concordant
evidence comes from non-supplementary records only, and a template that is discordant is never
also counted as concordant. Supplementary records are still kept for the discordant side,
because in repeats a junction read's only record at the site can be a supplementary. The
full-map count (field 36) is unchanged.

### Effect on the counts

Concordant counts fell by a factor of 2.0 to 2.3, and the discordant counts barely moved. The
median AF_PE of the final calls changed as follows:

| dataset | before | after |
|---|---|---|
| HG002 NovaSeq 30x | 0.420 | 0.611 |
| HG002 HiSeq 30x | 0.326 | 0.507 |
| HG002 20x / 15x / 10x | 0.436 / 0.447 / 0.526 | 0.630 / 0.640 / 0.714 |
| HG00438 | 0.432 | 0.604 |
| HG02257 | 0.405 | 0.576 |
| HG02622 | 0.412 | 0.582 |

If the concordant count halves, a ratio f becomes 2f/(1+f). That is the conversion used for
every AF_PE cutoff below.

### Genotyper refit

The PE half of the genotype likelihood is fitted on AF_PE. The HG002 default parameters, used
when a sample has too few informative calls for its own EM fit, were refit on the new counts
(NovaSeq 30x, 1494 calls):

| parameter | before | after |
|---|---|---|
| SR het / hom | 0.441 / 0.956 | 0.441 / 0.958 |
| PE het / hom | 0.380 / 0.761 | 0.569 / 0.913 |
| hom prior | 0.311 | 0.308 |

The SR parameters did not move, which confirms that the change is confined to the PE counts.
GT accuracy improved on all five HG002 datasets (+0.1 to +0.7 percentage points) and fell
slightly on the three 1KG samples (−0.1 to −0.2 points, 99.6% after). GQ calibration, pooled
over the 8 datasets (observed vs predicted GT error):

| GQ bin | observed | predicted |
|---|---|---|
| 0–10 | 25.9% | 29.2% |
| 10–20 | 4.4% | 4.0% |
| 20–30 | 0.5% | 0.4% |

Above GQ 30 the observed error is 0.4–0.6% against a predicted error near zero, so GQ is
still overconfident there. The fix narrows this gap (fewer sites reach GQ 60+: 4731 to 4131)
but does not close it. Part of the remaining gap is overdispersion: reads at a site are not
independent binomial draws.

### AF_PE cutoffs

Three post-filter cutoffs read AF_PE. None has a recorded derivation in master, and all three
were set on the double-counted scale. Each was swept by rerunning the real post-filter, then
merge, gene, genotype and gVCF, and scoring the final VCF as `score_final.sh` does.

| cutoff | master | converted 2f/(1+f) | swept | chosen |
|---|---|---|---|---|
| germline AF_PE floor (`GERMLINE_AF_PE_FLOOR`) | 0.075 | 0.1395 | 0–0.30, off | **0.14** |
| transduction floor (`TD_AF_PE_FLOOR`) | 0.08 | 0.148 | 0–0.30 | **0.15** |
| mosaic ceiling (`MOSAIC_AF_BOUNDS[2]`) | 0.35 | 0.52 | 0.04–1.0 | 0.35, then **per sample** (below) |

- **Germline floor.** It was chosen on HG002 NovaSeq 30x plus HG00438, where F1 is flat from 0
  to 0.30. HiSeq 30x, the 10x/15x/20x downsamples, HG02257 and HG02622 were held out. On those,
  the plateau middle (0.15) differs from the converted value by one FP on HiSeq, so the converted
  value is kept.
  - Leaving 0.075 on the new scale lets in 5 extra calls, all FPs.
  - Above about 0.18 the floor starts removing TPs on HiSeq (−38 TPs at 0.30).
- **Transduction floor.** It applies to 48–89 L1 transduction records per dataset. On all 8
  datasets the final calls do not change anywhere from 0 to 0.17. It is kept at the converted
  value, and folding it into the germline floor would give the same calls.
- **Mosaic ceiling.** The ceiling rejects calls whose pair fraction looks germline.
  - **Data.** HapMap mixtures BCM 300x, UW 100x and WashU 100x. Their post-filter inputs were
    moved onto the new counts per record, and calls were matched one-to-one to the mixture's
    somatic sites.
  - **Converted value.** 0.52 reproduces the old call sets to within 1.3% of calls.
  - **Sweep.** Leave-one-run-out picks 0.26–0.30, and all three runs together give an F1
    plateau at 0.28–0.30.
  - **First choice: 0.35, as an interim value** (B14). Same number as master, but on the new
    counts much stricter than master's (whose like-for-like value is 0.52). Per run, against
    0.52:

    | ceiling | BCM TP / FP | UW TP / FP | WashU TP / FP |
    |---|---|---|---|
    | 0.52 | 1402 / 1223 | 215 / 5789 | 692 / 637 |
    | 0.35 | 1393 / 623 | 189 / 3636 | 688 / 204 |
    | 0.30 | 1390 / 578 | 167 / 2512 | 678 / 182 |

  - **Replaced by a per-sample ceiling** (2026-11-13). See "Per-sample mosaic ceiling" below.

The B14 sweep data, figures and scripts are in the bug-2 review folder (`validation/tune/`).

### Per-sample mosaic ceiling (2026-11-13)

Where a true mosaic call's pair fraction lands depends on insert size relative to read length, so
a fixed ceiling favors some libraries. The mosaic step now sets the ceiling from the sample's own
germline calls:

- **Rule.** `disc_upper = 0.60 x median AF_PE` of the sample's germline het calls for that family,
  read from the germline high-confidence list (`postF_hc`) that the pipeline already writes. A
  germline het call is one with AF_SR ≥ 0.2 and 0.2 ≤ AF_PE < 0.8. The other two bounds
  (`MOSAIC_AF_BOUNDS[0:2]`) are unchanged.
- **Fallback.** With fewer than 20 germline het calls, or no `--germline-hc` list (direct
  `x_TEA_main.py --mosaic` calls), the fixed 0.35 is used. The mosaic step logs which value it
  used and why: `[mosaic] disc_upper (AF_PE ceiling) = 0.290 (per-sample: 0.6 x median AF_PE ...)`
  or `(fixed fallback: ...)`.
- **Data.** The same three HapMap mixture runs (BCM 300x, UW 100x, WashU 100x), post-filter
  outputs scored one-to-one against the mixture's somatic sites (±200 bp, same family). Calls at
  known HG005 germline sites are counted separately. Germline het median AF_PE is 0.46-0.52 in
  every run and family (48-4910 het calls).
- **Fit and test.** k was fitted on BCM + WashU and tested on UW, and also left one run out:

  | held out | k fitted on the other two | TP / FP / HG005 / F1 | best fixed ceiling on the other two | TP / FP / HG005 / F1 |
  |---|---|---|---|---|
  | BCM | 0.60 | 1388 / 571 / 1 / 0.597 | 0.28 | 1385 / 564 / 1 / 0.597 |
  | UW | 0.63 | 162 / 2147 / 0 / 0.065 | 0.30 | 165 / 2392 / 2 / 0.063 |
  | WashU | 0.57 | 670 / 175 / 1 / 0.379 | 0.28 | 670 / 175 / 1 / 0.379 |

- **Per run, fixed 0.35 vs the chosen rule** (k = 0.60, fitted on all three):

  | run | 0.35: TP / FP / HG005 / F1 | 0.60 x het median: TP / FP / HG005 / F1 | ceilings Alu / L1 / SVA |
  |---|---|---|---|
  | BCM 300x | 1393 / 623 / 6 / 0.592 | 1388 / 571 / 1 / 0.597 | 0.290 / 0.294 / 0.293 |
  | UW 100x | 187 / 3517 / 8 / 0.059 | 161 / 1953 / 0 / 0.067 | 0.277 / 0.305 / 0.309 |
  | WashU 100x | 687 / 205 / 9 / 0.383 | 676 / 183 / 2 / 0.381 | 0.300 / 0.289 / 0.304 |

- **What the gain is.** Most of it comes from lowering the ceiling from 0.35 to about 0.29, not
  from adapting it per sample: a fixed 0.28-0.30 does as well on BCM and WashU and within 0.002 F1
  on UW. The per-sample rule was kept because it costs nothing, removes nearly all calls at HG005
  germline sites (23 → 3 summed over the three runs) and follows libraries whose pair fractions
  sit elsewhere, which three runs cannot test. The per-run F1 plateaus are 0.25-0.32 (BCM),
  0.25-0.28 (UW) and 0.30-0.37 (WashU).
- **Variants tried.** The median over all germline-like calls (k 0.58) and a median pooled over
  families (k 0.60) give the same results to within 0.003 F1. A rule based on the VAF model
  (reject a call whose estimated VAF exceeds v_max; leave-one-out v_max 0.22-0.25) did worse on
  UW (F1 0.051 against 0.065 on the same offline scorer), whose germline calibration is distorted
  by a contaminated germline list (e_PE 2.07, l_PE 139 on Alu). It was not used.
- **Caveat.** The scoring above is per family before the cross-family merge, from the post-filter
  outputs rerun at each ceiling. The end-to-end check is below.

## Right-clip pair checks (2026-11-13)

This is B15 in `docs/BUG_LOG.md`, in `collect_features_one_site` (`x_genotype_feature.py`). For
right-clipped reads the pair checks (discordant test, `DFT_IS` window, concordant span) used the
clip position, which is the alignment end, instead of the alignment start. The clip position now
has its own variable. Only the pair counts (fields 39/40, `AF_DISC` / `AF_CONCORDNT`) change.
Split-read counts and AF_SR do not: the AF_SR medians of truth-matched calls are identical in
every dataset, family and truth GT (`sr_af_medians.tsv`).

A second counting change, supplementary junction reads as split-read evidence, was measured and
tried on the same branch but not merged. See the open issues in `docs/BUG_LOG.md`.

### Genotyper refit

Same method as B14: HG002 NovaSeq 30x rerun with the B15 counts and the old defaults, EM fit
from the VCF header, rounded (1494 informative calls):

| parameter | before | after |
|---|---|---|
| SR het / hom | 0.441 / 0.958 | 0.441 / 0.958 |
| PE het / hom | 0.569 / 0.913 | 0.570 / 0.914 |
| hom prior | 0.308 | 0.308 |

All 8 validation datasets use their own EM fit, so the defaults only set the EM start there.

### Validation (8 datasets, control `0553edb` vs `320f610`)

Pooled over families, scored as `score_final.sh` does:

| dataset | TP | FP | GT accuracy | GT errors |
|---|---|---|---|---|
| HG002 NovaSeq 30x | 1299 → 1299 | 168 → 168 | 98.99% → 98.99% | 13 → 13 |
| HG002 HiSeq 30x | 1289 → 1289 | 212 → 212 | 98.66% → 98.66% | 17 → 17 |
| HG002 20x | 1225 → 1225 | 132 → 132 | 98.51% → 98.51% | 18 → 18 |
| HG002 15x | 1122 → 1122 | 107 → 107 | 97.57% → 97.57% | 27 → 27 |
| HG002 10x | 841 → 841 | 61 → 61 | 96.39% → 96.39% | 30 → 30 |
| HG00438 | 1143 → 1143 | 389 → 390 | 99.62% → 99.62% | 4 → 4 |
| HG02257 | 1441 → 1441 | 443 → 443 | 99.56% → 99.63% | 6 → 5 |
| HG02622 | 1463 → 1463 | 474 → 474 | 99.63% → 99.71% | 5 → 4 |
| total | 9823 → 9823 | 1986 → 1987 | | 120 → 118 |

- **Call changes.** One gained FP, none lost. The gained call is HG00438 Alu chr6:32572129
  (HLA region). Ten fewer reads count as discordant there, so its `f_ef_disc` rises from 0.068
  to 0.075 and crosses the 0.075 cutoff.
- **Gate flips** (records in both arms' post-filter input, current cutoffs): across all 8 sets,
  7 records flip `f_ef_disc`, 5 the AF_PE gate and 6 the transduction gate. No other gate flips.
  Only the record above changes a final call.
- **GT.** Two fewer GT errors, both in 1KG samples; within noise.
- **GQ calibration** (pooled, observed vs predicted GT error, before → after):

  | GQ bin | before | after |
  |---|---|---|
  | 0–10 | 25.9% / 29.2% | 25.0% / 29.0% |
  | 10–20 | 4.4% / 4.0% | 4.4% / 4.0% |
  | 20–30 | 0.46% / 0.40% | 0.47% / 0.40% |
  | 30–50 | 0.58% / 0.02% | 0.58% / 0.02% |
  | 50–60 | 0.36% / 0.00% | 0.42% / 0.00% |
  | 60+ | 0.46% / 0.00% | 0.43% / 0.00% |

  Unchanged; GQ above 30 is still overconfident.
- **AF_PE cutoffs.** At `GERMLINE_AF_PE_FLOOR` 0.14 and `TD_AF_PE_FLOOR` 0.15, 0-3 records per
  dataset flip the AF_PE gate and 0-3 the transduction gate, in both directions. None changes a
  final call, so the cutoffs were not retuned.
- **Mosaic** (HapMap mixtures at `MOSAIC_AF_BOUNDS` 0.35; post-filter inputs moved per record by
  the recount difference between the two code versions, as for B14). Per run and family, B15
  changes field 39 in 6-29% of records and field 40 in 15-62%, and lowers the mean concordant
  count per record by 0.2-1.2:

  | run | TP / FP before (B14) | TP / FP after (B15) |
  |---|---|---|
  | BCM 300x | 1393 / 623 | 1393 / 623 |
  | UW 100x | 189 / 3636 | 187 / 3517 |
  | WashU 100x | 688 / 204 | 687 / 205 |

  The B14 inputs reproduce the B14 numbers exactly with the new code copy.

Tables and scripts are in `docs/validation/2026-11-13_counting_v2/`.

## GQ calibration (2026-11-13)

Before this change GQ predicted essentially no GT errors above 30, while the observed error there
was 0.4–0.6% (pooled over the 8 validation datasets, xTea-paper HG002 labels). GQ is now
`-10 log10(eps + (1 - 2 eps) perr)` with `eps` = 0.001, capped at 30. The floor is monotone and
does not change GT or `GTL`; `rho` stays 0.05 and the EM fit is unchanged.

**Where the high-GQ errors came from.**
- Of the 39 errors at GQ ≥ 30, 37 were in HG002, and they were the same ~10 sites in every HG002
  set, the independent HiSeq library included.
- The GIAB Q100 v5.0q HG002 benchmark agrees with xTea's call at 35 of those 37. Most of the
  gap is therefore label error in the xTea-paper HG002 set, which disagrees with Q100 at about
  0.8% of labeled sites.
- Under Q100 labels the remaining high-GQ errors are complex or nested sites, not count noise.
  Examples: two Q100 insertions (849 and 1912 bp) at one call; SVLEN 75 against 455 bp;
  duplicate calls at one insertion; SVLEN 8.

**Why a floor and not more overdispersion.**
- A per-sample `rho` fit is not stable. The marginal-likelihood MLE is about 0.02 at 30x and
  falls to 0–0.005 at 10x.
- A labeled joint fit on the training data puts `rho` between 0.02 and 0.06. 0.05 is within
  about 5 log-likelihood units of the optimum, so `rho` was kept.
- Raising `rho` does not remove the top-bin errors. At 0.1, 5 held-out errors stay at GQ ≥ 30
  with predicted error ≈ 0, and held-out GT errors rise from 84 to 126 (446 at 0.2).
- The error floor was fit by maximum likelihood on HG002 NovaSeq 30x (Q100 labels) plus HG00438
  (HPRC tier1): 0.00084, 95% CI 0.00013–0.0025, rounded to 0.001. With the xTea-paper labels the
  fit gives 0.0045, which reflects their label errors, so it was not used.

**Calibration, held-out sets** (HiSeq 30x, 20x/15x/10x downsamples, HG02257, HG02622). Q100
labels for HG002, HPRC tier1 for 1KG. Observed / predicted GT error:

| GQ bin | before | after |
|---|---|---|
| 0–10 | 24.7% / 28.7% | 24.6% / 28.7% |
| 10–20 | 4.9% / 3.8% | 4.6% / 3.8% |
| 20–30 | 0.52% / 0.41% | 20–25: 0.29% / 0.64%; 25–29: 0.59% / 0.21% (4/682) |
| 30–50 | 0.16% / 0.02% | — |
| 50–60 | 0 / 0.00% | — |
| 60 | 0.16% / 0.00% | — |
| 29–30 | — | 0.11% / 0.13% (6/5681) |

- Pooled over all 8 datasets, the top bin (29–30) is 0.10% observed against 0.13% predicted
  (8 of 7796).
- With the xTea-paper HG002 labels the top bin is 0.47% observed (26/5575). The excess is the
  label errors above.

**GT.** No genotype changes in any dataset. GT errors (xTea-paper / Q100 labels):

| Dataset | Errors |
|---|---|
| NovaSeq 30x | 13 / 7 |
| HiSeq 30x | 17 / 14 |
| 20x | 18 / 13 |
| 15x | 27 / 20 |
| 10x | 30 / 28 |
| HG00438 | 4 / 4 |
| HG02257 | 5 / 5 |
| HG02622 | 4 / 4 |

**End to end.** `--gntp_classify` and `--gVCF` were rerun with the new code on HG002 NovaSeq 30x
and HG02257. The records are identical to the B15 VCFs except for GQ. GT and GQ match the
offline evaluation at every record.

**Not measured.**
- Q100 labels cover HG002 only. Two HG02257 errors at GQ 31 and 38 (hom-looking reads, het
  label) are probably truth errors, but this was not checked.
- The floor does not model complex or nested sites separately. They make up most of the
  remaining top-bin errors.

Tables, the reliability plot and scripts are in `docs/validation/2026-11-13_gq_calibration/`.

## Benchmark methodology (HG002)

All runs use HG002 30× WGS, the same BAM, `-c 24`, `--mem=64G`, the whole pipeline, and default
xTea2 settings (`--parallel-downstream`). Three sets of numbers are reported:

- **xTea master**, three replicates (Aug 2026).
- **xTea2, Aug 2026 build** (`ea08f45`), three replicates, run concurrently with master's.
- **xTea2 v2.0.0** (`ffc1b46`), one run (job 59228997, 2026-11-10).

**Speed and memory**

| | xTea (n=3) | xTea2 Aug build (n=3) | xTea2 v2.0.0 (n=1) |
|---|---|---|---|
| wall clock | 1:03:18 – 1:21:18 | 0:37:52 – 0:39:45 | **0:39:10** |
| peak RSS | 38.12 – 38.19 GB | 17.83 – 18.21 GB | **21.7 GB** |

The Aug-build speed comparison (1.67×) was taken on a node that ran *both* tools concurrently (same
hardware, same instant, same page cache), because a cross-node ratio confounds the comparison with
how many jobs share each node's page cache for the BAM. The v2.0.0 run was a single run on its own
node, so its ~1.6× is a cross-run comparison. Peak memory rose from about 18 GB to 21.7 GB between
the two builds. On memory generally: xTea master holds **one** element family resident at a time (it
loops L1/Alu/SVA serially), while xTea2 holds **all three**, so part of the wall-clock win is
concurrency master does not have, not per-family speed.

**Sensitivity and precision** against an orthogonal HG002 truth set (positions, ±200 bp). Precision
here excludes the orphan/sibling *sites* the caller creates, which are not MEI insertions and which
the truth set does not enumerate, so they can never match:

| family | xTea sens / prec | xTea2 Aug build | **xTea2 v2.0.0** |
|---|---|---|---|
| L1 (197 truth) | 71.4% / 92.5% | 71.1% / 92.7% | 71.1% / 92.7% |
| Alu (1355) | 81.5% / 89.9% | 81.4% / 90.2% | 81.4% / 90.2% |
| SVA (90) | 88.9% / 90.2% | 85.6% / 91.7% | **88.9% / 93.0%** |
| **pooled (1642)** | 80.67% / 90.15% | 80.39% / 90.53% | **80.57% / 90.62%** |

Counting every call, including created orphan/sibling sites (the convention behind the README's
headline table), pooled precision is 87.6% for xTea, 88.3% for the Aug build and **88.6%** for
v2.0.0. The difference is entirely in L1, where those sites lower v2.0.0's precision to 75.7%.

**SVA.** The Aug build reproducibly missed 3 SVA truth records (`chr1:88507126`, `chr14:20707961`,
`chr5:56152167`). v2.0.0 recovers all three, bringing SVA sensitivity level with xTea's.

**Reproducibility** (Aug build). Across three runs, xTea2 produced byte-identical call sets (100% of
call positions shared, all three families). xTea master did not: 22 of its calls appear in some runs
and not others. So while master recovered slightly more truth records *per run* (1324.7 vs 1320.0 on
average), the two recovered **exactly the same number reproducibly** (1320 each); master's extra
hits were a different handful each time. v2.0.0 recovers 1323. Reproducibility has not been
re-measured on v2.0.0.

These numbers are **not comparable to the xTea paper's**: this BAM is hs38DH-aligned rather than
`no_alt`, and MAPQ filtering behaves differently in exactly the repeat-homologous regions that
matter.

### Other genome builds

Other builds run if the repeat library carries them, but are not benchmarked. CHM13's L1 library
ships **zero** decoy records, so the transduction decoy behavior is untested there.

### Aligner

The benchmark BAMs are bwa-mem alignments, and bwa-compatible aligners (Sentieon, DRAGEN) behave
the same way. Discordant-pair evidence comes from anchors whose mate bwa placed on some reference
copy of the element. An anchor whose mate is **unmapped** is skipped: the scan's discordant pass
requires a mapped mate (`src/csrc/xtea_scan.cpp`, the `BAM_FMUNMAP` test in the DISC block), and
legacy xTea skips these reads too. Aligners that leave repeat-derived mates unmapped instead of
placing them (e.g. Novoalign) therefore lose almost all discordant-pair support, and sensitivity
collapses. HG002 GIAB 2x250 Novoalign gave 25 calls. This is documented, not fixed: realign such
data with `bwa mem` first.
