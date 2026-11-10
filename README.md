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

See [Performance](#performance) for details.

---

## Requirements

- **bwa** ≥ 0.7.17 and **samtools** ≥ 1.0, on `PATH`
- **htslib** ≥ 1.22 and **pkg-config**, to build the C++ scan engine
- **Python** ≥ 3.10, with `numpy`, `pysam`, `pyyaml`, `scipy`, `pandas` and **`onnxruntime` ≥ 1.19**
  (installed for you by either route below)

## Install

### 1. Get the code and the repeat library

```bash
git clone https://github.com/corinsexton/xTea.git
cd xTea
git checkout xTea2

wget https://github.com/parklab/xTea/raw/master/rep_lib_annotation.tar.gz
tar -xzf rep_lib_annotation.tar.gz
```

You also need a decompressed gene annotation GFF3 from
[GENCODE](https://www.gencodegenes.org/human/): the GRCh38 release for `--genome hg38`, the
GRCh37-lift release for `hg19`.

### 2. Install the Python package

Conda (recommended; it pins a verified combination):

```bash
mamba env create -f envs/xtea2.yml
conda env config vars set PYTHONNOUSERSITE=1 -n xtea2   # required, see note below
conda activate xtea2
pip install -e .
```

**Do not skip the `PYTHONNOUSERSITE=1` line.** Without it, an onnxruntime in `~/.local` can shadow
the environment's copy and crash the genotyping step.

Or plain pip into any Python ≥ 3.10:

```bash
pip install .
```

Either way you get an **`xtea2`** command on `PATH`.

### 3. Build the C++ scan engine

The binary is not tracked in git, so every fresh clone must build it:

```bash
make -C src/csrc
# if htslib is not on the default pkg-config path:
PKG_CONFIG_PATH=/path/to/htslib/lib/pkgconfig make -C src/csrc
```

If you installed the package outside the repo, export `XTEA2_SCAN_BIN=/path/to/xtea_scan` or put
`xtea_scan` on `PATH`.

### 4. Verify the install

```bash
python3 -m pytest -m "not slow"
```

Build `xtea_scan` first, or the scan-engine tests are skipped. See
[tests/README.md](tests/README.md) to point the suite at your own data.

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
  --parallel-downstream \
  -n 24
```

This writes a merged VCF to `<output>/my_sample/my_sample.vcf`. Add `--dry-run` to print the
per-step plan and exit without running anything.

## Options

### Required

```
  -i, --input                Input BAM or CRAM (sorted and indexed). Pass several
                             space-delimited paths in ONE argument to pool libraries
                             of one sample: -i "a.cram b.cram c.cram"
  -s, --sample               Sample identifier (used in output file names)
  -o, --output-dir           Output directory
  --reference                Reference genome FASTA (also used to decode CRAM)
  --rep-lib-dir              Path to the decompressed rep_lib_annotation/ directory
  --gene-gff3                Gene annotation in GFF3 format
```

### Common

```
  --types                    Comma-separated MEI types. Default: L1,Alu,SVA
  --genome                   Genome build subdirectory under --rep-lib-dir. Default: hg38
  -n, --cores                Worker cores. Default: 8
  --parallel-downstream      Run the per-type steps concurrently. Faster; results unchanged.
  --nfclip                   Min clipped reads at the final filter. Default: 4
  --nfdisc                   Min discordant reads at the final filter. Default: 5
  --flklen                   Flank length for realignment and transduction (bp). Default: 3000
  --blacklist                Blacklist BED (centromeres, segmental dups). Default: none
  --region-mask SPEC         Drop candidates in these regions: BED/BED.gz paths, built-in
                             names (simple_repeats, segdup, centromeres), 'default' or
                             'none'. Unset: 'default' under --mosaic, 'none' otherwise.
  --ignore-rg                Allow multi-file input whose @RG SM tags differ (otherwise
                             a hard error, to prevent accidental sample mixing)
  --keep-intermediates       Keep per-type VCFs, logs and intermediates on success
  --resume                   Reuse completed scan output from a previous run when its inputs
                             and the scan binary are unchanged
  --dry-run                  Print the per-step plan and exit
```

Clip and discordant read cutoffs used to seed candidates are derived from each BAM's own coverage;
there is no `--nclip`, `--cr` or `--nd`.

### Tumor and case-control (somatic)

```
  --tumor                    Tumor sample: scale the coverage-derived cutoffs by --purity
  --purity                   Tumor purity, 0..1. Default: 0.45
  --case-control             Somatic calls against a matched normal. Requires --normal-bam.
                             Output is SOMATIC CALLS ONLY; no germline VCF is written.
  --normal-bam               Matched normal BAM/CRAM (same multi-file convention as -i)
```

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

### Mosaic (experimental)

`--mosaic` calls low-allele-fraction insertions from high-coverage (200×+) bulk data from one
sample. **It is implemented but not yet validated; treat its output as experimental.**

```
  --mosaic                   Low clip/disc cutoffs, polyA tally, mosaic post-filter and
                             --region-mask default
  --mosaic-min-clip N        Clip-support cutoff. Default: 2
  --mosaic-min-disc N        Discordant-support cutoff. Default: 2
```

```bash
xtea2 \
  -i "rep1.cram rep2.cram" -s tissue_01 -o /path/to/output \
  --reference /path/to/reference.fa \
  --rep-lib-dir /path/to/rep_lib_annotation \
  --gene-gff3 /path/to/genes.gff3 \
  --mosaic --region-mask simple_repeats.bed,segdup.bed.gz,centromeres.bed.gz \
  -n 24
```

> The built-in region masks currently point to files on the Park lab's O2 cluster and are
> GRCh38-only. Elsewhere, pass `--region-mask` with your own GRCh38 BEDs (simple repeats,
> segmental duplications, centromeres) or `--region-mask none`, or `--mosaic` will stop with
> `region mask file(s) not found`.

Transduction, scan-engine and diagnostic options are documented in
[docs/supplemental/README_NOTES.md](docs/supplemental/README_NOTES.md). Run `xtea2 --help` for the
full list.

## Output

- **Merged VCF** at `<output>/<sample>/<sample>.vcf`, unified across types and deduplicated. By
  default this is the only file kept.
- Failed or incomplete runs are never cleaned up.
- Reserve several GB of scratch per run: a whole-genome run reached about 7 GB before cleanup.

## Performance

On HG002 30× WGS (24 cores, 64 GB, three replicates per tool), xTea2 ran **1.67× faster than xTea
on 47% of the memory** and produced **byte-identical call sets** across runs, while xTea's varied.
Pooled sensitivity was essentially unchanged (80.4% vs 80.7%) and precision was slightly higher.
Methodology, per-family results and caveats are in
[docs/supplemental/README_NOTES.md](docs/supplemental/README_NOTES.md#benchmark-methodology-hg002).

## Scope and limitations

- **Validated:** hg38, short-read Illumina WGS, germline, L1/Alu/SVA.
- **Implemented but not validated:** other genome builds, case-control/somatic mode, mosaic mode.
- **A default run does not reproduce legacy xTea's output exactly:** transduction detection was
  fixed so that SVA transductions from reference copies can be found on hg38 (see the
  [notes](docs/supplemental/README_NOTES.md#transduction-detection-options)).
- **Not carried over from xTea:** long reads, 10X linked reads, de-novo trio calling, reference TE
  deletion and population allele frequency. Use [upstream xTea](https://github.com/parklab/xTea) for
  those.

## Citation

If you use xTea, please cite:

> Chu, Chong, Rebeca Borges-Monroy, Vinayak V. Viswanadham, Soohyun Lee, Heng Li, Eunjung Alice Lee,
> and Peter J. Park. "Comprehensive identification of transposable element insertions using multiple
> sequencing technologies." *Nature Communications* 12, no. 1 (2021): 3836.
