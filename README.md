# xTea2

**Transposable element insertion calling from short-read WGS, in a single process.**

xTea2 identifies mobile element insertions (**L1**, **Alu**, **SVA**) from paired-end Illumina
BAM/CRAM. It is a reimplementation of [xTea](https://github.com/parklab/xTea): where the original
generates a shell script that you then submit, xTea2 runs start to finish in one invocation, and a
**C++ scan engine** collects clipped and discordant reads in a single BAM pass instead of one pass
per element type.

Against the original on the same HG002 30× WGS BAM, same 24 cores and 64 GB:

| | xTea | **xTea2** |
|---|---|---|
| wall clock | 1:03 | **0:38** (**1.67×** faster) |
| peak memory | 38.1 GB | **18.0 GB** (**47%**) |
| run-to-run call sets | vary | **byte-identical** |
| sensitivity / precision | 80.7% / 87.6% | 80.4% / **88.3%** |

See [Performance and accuracy](#performance-and-accuracy) for how these were measured and what the
numbers do and do not support.

---

## Requirements

**External tools**, on `PATH`:

- **bwa** ≥ 0.7.17 (the `-o` option is required)
- **samtools** ≥ 1.0
- **htslib ≥ 1.22**, discoverable via `pkg-config` — needed to *build* the C++ scan engine
- **pkg-config**, to locate htslib at build time

**Python ≥ 3.10.** This is a hard floor with nothing below it: the legacy DF21 genotyping model
published only cp36–cp39 wheels with no sdist, which is precisely why genotyping moved to ONNX.

**Python packages** (installed for you by either route below): `numpy`, `pysam`, `pyyaml`, `scipy`,
`pandas`, and **`onnxruntime` ≥ 1.19**. That floor matters — onnxruntime 1.16.x is compiled against
numpy 1.x and, under numpy 2, raises `AttributeError: _ARRAY_API not found` while dumping an ~880 MB
core file. Genotyping is on the critical path to the final VCF, so an environment without a working
onnxruntime cannot complete a run.

## Install

### 1. Get the code and the repeat library

```bash
git clone https://github.com/corinsexton/xTea.git
cd xTea
git checkout xTea2
```

The pre-processed repeat library is a data dependency shared with upstream xTea:

```bash
wget https://github.com/parklab/xTea/raw/master/rep_lib_annotation.tar.gz
tar -xzf rep_lib_annotation.tar.gz
```

You also need a gene annotation GFF3, decompressed, from
[GENCODE](https://www.gencodegenes.org/human/) — the GRCh38 release for `--genome hg38`, the
GRCh37-lift release for `hg19`.

### 2. Install the Python package

Either a conda environment (recommended — it pins the verified combination):

```bash
mamba env create -f envs/xtea2.yml
conda env config vars set PYTHONNOUSERSITE=1 -n xtea2
conda activate xtea2
pip install -e .
```

**The `conda env config vars` line is not optional.** `conda env create` does not carry environment
variables from the YAML and there is no YAML key that sets them. Skip it and you get an environment
that passes its own test suite and then **SIGSEGVs at the genotype step**, because a user-site
`~/.local` onnxruntime compiled against numpy 1.x shadows the environment's copy on `sys.path`. With
the variable registered, `conda activate xtea2` is the entire setup.

Or plain pip into any Python ≥ 3.10:

```bash
pip install .
```

Either way you get an **`xtea2`** command on `PATH`. (Editable mode, `pip install -e .`, is worth it
if you intend to modify the code.)

### 3. Build the C++ scan engine

The binary `src/csrc/xtea_scan` is **not tracked in git** — only its sources are — so every fresh
clone must build it:

```bash
make -C src/csrc
```

If htslib is not on the default pkg-config path (a local or conda install), point at it:

```bash
PKG_CONFIG_PATH=/path/to/htslib/lib/pkgconfig make -C src/csrc
```

The build verifies htslib at compile time and embeds an rpath, so the binary finds `libhts` at
runtime without `LD_LIBRARY_PATH`.

`xtea2` locates the binary in this order: **`--scan-bin`** → **`$XTEA2_SCAN_BIN`** → the in-repo
build at `src/csrc/xtea_scan` → **`xtea_scan` on `PATH`**. If you installed the package outside the
repo, either export `XTEA2_SCAN_BIN` or copy `xtea_scan` somewhere on `PATH`.

Rebuild after editing **any** file under `src/csrc/`, headers included: `realign.hpp`,
`cns_emit.hpp` and `td_gather.hpp` carry most of the engine and are `#include`d rather than compiled
separately. The Makefile lists them as prerequisites, so a plain `make` picks up a header-only edit.

### 4. Verify the install

```bash
python3 -m pytest -m "not slow"
```

Build `xtea_scan` first: about 30 tests cover the C++ engine and resolve it from
`src/csrc/xtea_scan`, so without a build they skip and a green run means less than it looks. A binary
**older than its own sources** is a hard failure rather than a skip — that is the state in which a
passing suite is testing code that is not in your tree. See **[tests/README.md](tests/README.md)**
for the environment variables that point the suite at your own data.

## Quick start

```bash
xtea2 \
  -i sample.cram \
  -s my_sample \
  -o /path/to/output \
  --reference /path/to/reference.fa \
  --rep-lib-dir /path/to/rep_lib_annotation \
  --gene-gff3 /path/to/genes.gff3 \
  --types L1,Alu,SVA \
  --genome hg38 \
  -n 24
```

Produces a merged VCF at `<output>/my_sample/my_sample.vcf`.

Add `--dry-run` to print the per-step plan and exit without executing — the quickest way to see
which defaults resolved on for your input.

> Working inside a clone without installing? `python xtea/run_xtea.py …` is equivalent — the
> `xtea2` command is a thin entry point onto the same `main()`, so the two cannot drift apart.

## Options

### Required

```
  -i, --input                Input BAM or CRAM (sorted and indexed).
                             Pass multiple space-delimited paths in ONE argument to
                             pool libraries: -i "a.cram b.cram c.cram"
  -s, --sample               Sample identifier (used in output file names)
  -o, --output-dir           Output directory
  --reference                Reference genome FASTA (also used to decode CRAM)
  --rep-lib-dir              Path to the decompressed rep_lib_annotation/ directory
  --gene-gff3                Gene annotation in GFF3 format
```

### Types and resources

```
  --types                    Comma-separated MEI types. Default: L1,Alu,SVA
                             Choices: L1, Alu, SVA (aliases: line1, alu, sva)
  --genome                   Genome build subdirectory under --rep-lib-dir. Default: hg38
                             (hg19, CHM13 etc. work if the rep-lib carries them; see
                             Scope and limitations)
  -n, --cores                Worker cores. Default: 8
  --parallel-downstream      Run the per-type downstream chains concurrently rather than
                             serially. Pure wall-clock win, results unchanged; each type
                             gets cores//ntypes workers.
  --keep-intermediates       Retain intermediate files and per-step logs on success.
                             By default only the merged VCF is kept. Failed runs are
                             NEVER cleaned, regardless of this flag.
```

### Thresholds

```
  --nfclip                   Min clipped reads at the consensus/final filter. Default: 4
  --nfdisc                   Min discordant reads at the consensus/final filter. Default: 5
  --flklen                   Flank length for realignment and transduction (bp). Default: 3000
  --blacklist                Blacklist BED (centromeres, segmental dups). Default: none
  --dry-run                  Print the per-step commands and exit without executing
```

**There is no `--nclip`, `--cr` or `--nd`**, and passing one is a hard error. The seeding cutoffs
those names refer to — `n_clip_cutoff`, the clipped-in-repeat minimum, and `n_disc_cutoff` — are
**derived from the input BAM's own coverage** at runtime, and scaled by `--purity` under `--tumor`.
(They were effectively inert on legacy xTea's command line too, which overrode them from coverage
unless `--user` was passed.) Some diagnostic banners still say `--nclip` when they mean the derived
`n_clip_cutoff`; that is a wart in the log text, not a flag you can pass.

### Tumor samples

```
  --tumor                    Input is a tumor sample: scale the coverage-derived clip/disc
                             cutoffs by --purity. Does not by itself require case-control.
  --purity                   Tumor purity, 0..1. Default: 0.45. Used only with --tumor.
```

### Case-control (somatic) mode

```
  --case-control             Somatic filtering against a matched normal. Requires --normal-bam.
                             OUTPUT IS SOMATIC CALLS ONLY — no germline VCF is written.
  --normal-bam               Matched normal BAM/CRAM. Accepts the same space-delimited
                             multi-file convention as -i.
  --ignore-rg                Skip the @RG SM-tag consistency check across multi-file input.
                             By default differing SM tags are a hard error, to prevent
                             accidental sample mixing. Pass this when pooling libraries
                             that legitimately carry different SM tags.
```

### Transduction detection

These control how transductions — a TE carrying flanking sequence from its source locus — are
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

> **A default run does not reproduce legacy xTea's output, and that is deliberate.** Both flips
> above exist because without them **no reference-copy SVA transduction is detectable at all** on
> hg38. The driver says so on every run:
>
> ```
> [decoy] the transduction decoy test is STRUCTURAL (DEFAULT ...): this run does NOT reproduce BASELINE-0
> [decoy-concat] the transduction bwa reference GAINS the decoy consensus block (DEFAULT): ...
> ```
>
> To reproduce a call set frozen before them, pass **both** `--no-td-decoy-structural` and
> `--no-td-decoy-concat`.

### Scan-path options

The C++ engine combines passes that were separate. These are opt-in; leaving one unset keeps the
corresponding Python step.

```
  --scan-clip                Use xtea-scan for clip detection (one pass per input file)
  --scan-disc                Take the discordant stream from the same pass. Implies --scan-clip
  --cns-emit-scan            Collect the cns_filter step's SAM inputs in one pass and reuse
                             them, instead of re-reading the BAM per candidate. Requires
                             --scan-disc.
  --scan-bin                 Path to the xtea_scan binary. See Install step 3 for the
                             resolution order when this is not given.
  --bwa                      Path to bwa. Default: bwa on PATH
  --scan-min-clip            Min in-repeat clipped reads to emit a clip candidate. Default: 2
  --scan-disc-threshold      Same-chromosome discordant distance threshold (bp) for emission.
                             Default: 2000 (deliberately below the effective cutoff so the
                             stream is a superset; the adapter re-applies the exact threshold)
  --scan-realign-ref         Reference for realigning clipped parts.
                             copies_with_flank (default, more sensitive) | consensus
```

### Advanced and diagnostic

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
```

## Examples

**Germline, scan-optimised** — the configuration the performance numbers were measured on:

```bash
xtea2 \
  -i sample.cram -s my_sample -o /path/to/output \
  --reference /path/to/reference.fa \
  --rep-lib-dir /path/to/rep_lib_annotation \
  --gene-gff3 /path/to/genes.gff3 \
  --types L1,Alu,SVA \
  --scan-clip --scan-disc --cns-emit-scan \
  --parallel-downstream \
  -n 24
```

**Case-control (somatic):**

```bash
xtea2 \
  -i tumor.cram -s patient_01 -o /path/to/output \
  --reference /path/to/reference.fa \
  --rep-lib-dir /path/to/rep_lib_annotation \
  --gene-gff3 /path/to/genes.gff3 \
  --case-control --normal-bam normal.cram \
  --tumor --purity 0.72 \
  -n 24
```

Output is **somatic calls only**; no germline VCF is written in this mode.

## Multi-file input

`-i` and `--normal-bam` accept space-delimited paths in a single argument, to pool several libraries
of one sample:

```bash
xtea2 -i "lib1.cram lib2.cram lib3.cram" -s sample_01 ...
```

By default all files must carry the same `@RG SM` tag; pass `--ignore-rg` when they legitimately
differ (e.g. different sequencing centers).

All input files must be aligned against the same reference, and the scan engine **enforces** this
rather than trusting it: before scanning it compares every file's header against the first and fails
loudly if the sequence count, the `@SQ` names or their order, or any sequence length differs.

Pooling is a flat sum — coverage adds, read length and insert size are averaged, clip/disc counts
are summed, with no per-file weighting. Note that duplicate reads cannot be detected *across* files,
since the duplicate flag is per-file: a fragment independently sequenced at two centers is counted
twice.

## Output

- **Merged VCF** at `<output>/<sample>/<sample>.vcf`, unified across types and deduplicated. This is
  the only output retained by default.
- `--keep-intermediates` additionally retains per-type VCFs, per-step logs and all intermediates.
- **Failed or incomplete runs are never cleaned**, regardless of options.

Reserve several GB of scratch per run: on a whole-genome benchmark the working tree reached about
7 GB before cleanup, against a final VCF of about 1 MB.

## Performance and accuracy

Measured on HG002 30× WGS, same BAM, `-c 24`, `--mem=64G`, **three replicates per arm** against
xTea master, whole pipeline both sides.

**Speed and memory**

| | xTea (n=3) | xTea2 (n=3) |
|---|---|---|
| wall clock | 1:03:18 – 1:21:18 | **0:37:52 – 0:39:45** |
| peak RSS | 38.12 – 38.19 GB | **17.83 – 18.21 GB** |

**1.67× faster on 47% of the memory.** The speed figure is taken from a node that ran *both* arms
concurrently — same hardware, same instant, same page cache — because a cross-node ratio would
confound the comparison with how many jobs were sharing each node's page cache for the BAM. Read the
memory number the strong way: xTea master holds **one** element family resident at a time (it loops
L1/Alu/SVA serially), while xTea2 holds **all three**. Correspondingly, part of the wall-clock win is
concurrency master does not have, not per-family speed.

**Sensitivity and precision**, against an orthogonal HG002 truth set (positions, ±200 bp):

| family | xTea sens / prec | xTea2 sens / prec |
|---|---|---|
| L1 (197 truth) | 71.4% / 92.5% | 71.1% / 92.7% |
| Alu (1355) | 81.5% / 89.9% | 81.4% / **90.2%** |
| SVA (90) | **88.9%** / 90.2% | 85.6% / **91.7%** |
| **pooled (1642)** | 80.67% / 90.15% | 80.39% / **90.53%** |

(Precision here excludes created orphan/sibling *sites*, which are not MEI insertions and which the
truth set does not enumerate, so they can never match. Including them depresses L1 precision on both
tools equally, to 72.1% and 74.5%.)

**Reproducibility is the real difference.** Across three runs, **xTea2 produced byte-identical call
sets** — 100% of call positions shared, all three families. xTea master did not: 22 of its calls
appear in some runs and not others. Consequently, while master recovers slightly more truth records
*per run* (1324.7 vs 1320.0 on average), the two recover **exactly the same number reproducibly**
(1320 each) — master's extra hits are a different handful each time. For a one-off call set that
lottery is a small gain; for a cohort, or any comparison across samples, it is a liability.

The one systematic sensitivity difference is **SVA: xTea2 reproducibly misses 3 truth records**
(`chr1:88507126`, `chr14:20707961`, `chr5:56152167`) — offset by **+3 on Alu**.

## Scope and limitations

What is validated, and what is not:

- **Validated**: hg38, short-read Illumina WGS, germline, L1/Alu/SVA. This is the configuration all
  the numbers above come from.
- **Other genome builds** run if the repeat library carries them, but are not benchmarked here.
  Note that CHM13's L1 library ships **zero** decoy records, so the transduction decoy behaviour is
  untested there.
- **Case-control/somatic mode is implemented but has no regression gate.** Treat its output as
  unvalidated.
- **Not carried over from xTea**: long reads (PacBio/Nanopore), 10X linked reads, mosaic and de-novo
  trio calling, reference TE deletion, population allele frequency, and the `gnrt_pipeline_*` script
  generators. Use [upstream xTea](https://github.com/parklab/xTea) and its branches for those.
- These numbers are **not comparable to the xTea paper's**: this BAM is hs38DH-aligned rather than
  `no_alt`, and MAPQ filtering behaves differently in exactly the repeat-homologous regions that
  matter.

## Citation

If you use xTea, please cite:

> Chu, Chong, Rebeca Borges-Monroy, Vinayak V. Viswanadham, Soohyun Lee, Heng Li, Eunjung Alice Lee,
> and Peter J. Park. "Comprehensive identification of transposable element insertions using multiple
> sequencing technologies." *Nature Communications* 12, no. 1 (2021): 3836.
