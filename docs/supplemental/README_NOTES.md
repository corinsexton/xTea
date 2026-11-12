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
the HG002 fit is used instead: SR 0.441/0.958, PE 0.569/0.913, hom prior 0.308 (fit on counts
with the AF_FMAP fix and with each read pair counted once on both sides of AF_PE). The parameters
used are recorded in the VCF header line `##xtea2_genotype_params`.

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

This is bug 2 in `docs/BUG_LOG.md`.

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
| mosaic ceiling (`MOSAIC_AF_BOUNDS[2]`) | 0.35 | 0.52 | 0.04–1.0 | **0.35** (interim) |

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
  - **Chosen: 0.35, as an interim value.** Same number as master, but on the new counts it is
    much stricter than master's (whose like-for-like value is 0.52). Per run, against 0.52:

    | ceiling | BCM TP / FP | UW TP / FP | WashU TP / FP |
    |---|---|---|---|
    | 0.52 | 1402 / 1223 | 215 / 5789 | 692 / 637 |
    | 0.35 | 1393 / 623 | 189 / 3636 | 688 / 204 |
    | 0.30 | 1390 / 578 | 167 / 2512 | 678 / 182 |

    0.35 keeps most of the FP cut of the F1 pick (0.30), including nearly all calls at known
    HG005 germline sites, but halves the sensitivity cost on UW (−12% TPs vs −22%). UW has the
    shortest inserts (430 bp), which push true mosaic pair fractions up.
  - **Why interim.** Where a true mosaic call's pair fraction lands depends on insert size
    relative to read length, so any fixed ceiling favors some libraries. The plan is a
    per-sample ceiling (for example a fraction of the sample's own germline het AF_PE median),
    as part of the per-site VAF-correction work.

The sweep data, figures and scripts are in the bug-2 review folder (`validation/tune/`).

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
