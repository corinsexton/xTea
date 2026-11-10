# xTea2 — supplemental notes

Background, rationale and developer-facing options moved out of the main [README](../../README.md)
to keep it focused on installing and running xTea2. Nothing here is needed for a standard run.

---

## Install details

### Why Python ≥ 3.10 and onnxruntime ≥ 1.19

Python 3.10 is a hard floor: the legacy DF21 genotyping model published only cp36–cp39 wheels with
no sdist, which is why genotyping moved to ONNX.

The onnxruntime floor matters because onnxruntime 1.16.x is compiled against numpy 1.x and, under
numpy 2, raises `AttributeError: _ARRAY_API not found` while dumping an ~880 MB core file.
Genotyping is on the critical path to the final VCF, so an environment without a working
onnxruntime cannot complete a run.

### Why `PYTHONNOUSERSITE=1` is required

`conda env create` does not carry environment variables from the YAML, and there is no YAML key that
sets them. Without `conda env config vars set PYTHONNOUSERSITE=1 -n xtea2` you get an environment
that passes its own test suite and then **SIGSEGVs at the genotype step**, because a user-site
`~/.local` onnxruntime compiled against numpy 1.x shadows the environment's copy on `sys.path`.

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
  --model                    Genotyping model: an ONNX FILE, not a directory.
                             Default: the bundled genotype_model_6_7_2024.onnx.
                             Requires onnxruntime >= 1.19.
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

## Benchmark methodology (HG002)

Measured on HG002 30× WGS, same BAM, `-c 24`, `--mem=64G`, **three replicates per arm** against
xTea master, whole pipeline both sides.

| | xTea (n=3) | xTea2 (n=3) |
|---|---|---|
| wall clock | 1:03:18 – 1:21:18 | **0:37:52 – 0:39:45** |
| peak RSS | 38.12 – 38.19 GB | **17.83 – 18.21 GB** |

The speed figure is taken from a node that ran *both* arms concurrently (same hardware, same instant,
same page cache), because a cross-node ratio would confound the comparison with how many jobs were
sharing each node's page cache for the BAM. On memory: xTea master holds **one** element family
resident at a time (it loops L1/Alu/SVA serially), while xTea2 holds **all three**. Correspondingly,
part of the wall-clock win is concurrency master does not have, not per-family speed.

**Sensitivity and precision** against an orthogonal HG002 truth set (positions, ±200 bp):

| family | xTea sens / prec | xTea2 sens / prec |
|---|---|---|
| L1 (197 truth) | 71.4% / 92.5% | 71.1% / 92.7% |
| Alu (1355) | 81.5% / 89.9% | 81.4% / **90.2%** |
| SVA (90) | **88.9%** / 90.2% | 85.6% / **91.7%** |
| **pooled (1642)** | 80.67% / 90.15% | 80.39% / **90.53%** |

Precision here excludes created orphan/sibling *sites*, which are not MEI insertions and which the
truth set does not enumerate, so they can never match. Including them depresses L1 precision on both
tools equally, to 72.1% and 74.5%.

**Reproducibility.** Across three runs, xTea2 produced byte-identical call sets (100% of call
positions shared, all three families). xTea master did not: 22 of its calls appear in some runs and
not others. So while master recovers slightly more truth records *per run* (1324.7 vs 1320.0 on
average), the two recover **exactly the same number reproducibly** (1320 each); master's extra hits
are a different handful each time.

The one systematic sensitivity difference is SVA: xTea2 reproducibly misses 3 truth records
(`chr1:88507126`, `chr14:20707961`, `chr5:56152167`), offset by +3 on Alu.

These numbers are **not comparable to the xTea paper's**: this BAM is hs38DH-aligned rather than
`no_alt`, and MAPQ filtering behaves differently in exactly the repeat-homologous regions that
matter.

### Other genome builds

Other builds run if the repeat library carries them, but are not benchmarked. CHM13's L1 library
ships **zero** decoy records, so the transduction decoy behavior is untested there.
