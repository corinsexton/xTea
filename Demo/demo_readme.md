## Note: Here we use `NA12878` downloaded from Illumina (Platinum Genomes) as the demo input bam. Please replace it with your own bam if you already have a processed one.

> **Updated for xTea2.** This demo previously ran upstream's two-stage workflow —
> `run_gnrt_pipeline.sh` invoked `gnrt_pipeline_local.py` to *generate* a `submit_jobs.sh`, which
> you then submitted. xTea2 replaced that with a single self-contained entry point, so the
> generator script and its `sample_id.txt` / `sample_bam.txt` input lists have been removed.
> The reference output under `output/` is unchanged.

## 1. Download and preprocess the bam file
### 1.1 Download
From `https://github.com/Illumina/PlatinumGenomes` (`https://www.ebi.ac.uk/ena/data/view/PRJEB3381`), specifically with wget:
```
wget ftp://ftp.sra.ebi.ac.uk/vol1/run/ERR194/ERR194147/NA12878_S1.bam
```
### 1.2 Preprocess (make sure the bam is sorted and indexed)
```
samtools index NA12878_S1.bam
```

## 2. Run xTea2

There is no script-generation step. Run the pipeline start-to-finish in one invocation:

```
python xtea/run_xtea.py \
    -i NA12878_S1.bam \
    -s NA12878 \
    -o demo_out/ \
    --reference <ref.fa> \
    --rep-lib-dir <rep_lib_annotation/> \
    --gene-gff3 <genes.gff3>
```

Add `--dry-run` first to print the step plan without executing it. See the **Usage: run_xtea.py**
section of the top-level [`README.md`](../README.md) for the full argument list and for worked
examples (germline, tumor/somatic, and case-control).

## 3. Demo output
The demo output are 3 VCF files (for L1, Alu and SVA respectively) under the `output` folder.
