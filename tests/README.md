# xTea2 test harness

Reproducibility + read-level tests guarding the refactor. **Results are the contract, not
intermediate file layout** — see `docs/clipdisc_contract_spec.md`.

## Layout
- `utils/vcf_io.py` — parse xTea MEI VCFs into normalized `MEICall`s.
- `utils/concordance.py` — call-set concordance (±25 bp + MEI type) → sensitivity / PPV. **The gate.**
- `utils/clipdisc_io.py` — golden clip/disc site parsers + a **pysam oracle** that re-derives
  clipped/discordant reads from a BAM/CRAM per the spec (the read-level ground truth).
- `test_concordance_comparator.py` — comparator unit tests (synthetic; always run).
- `test_e2e_concordance.py` — golden self-concordance now; **xtea2-vs-golden gate** activates when
  `XTEA2_OUTPUT_DIR` is set.
- `test_clip_reads.py` / `test_disc_reads.py` — read-level "same reads found" invariant, incl.
  **BAM↔CRAM equivalence**; oracle-vs-`xtea-scan` activates when `XTEA2_SCAN_BIN` is set.
- `test_rg_consistency.py` — the multi-file `@RG` SM-tag guard: differing SM tags must fail loud
  before any pipeline step, `--ignore-rg` must opt out, on both `-i` and `--normal-bam`
  (synthetic BAMs + `--dry-run`; always runs, no external data).

## Running
```bash
conda activate xtea2         # create it with: mamba env create -f envs/xtea2.yml
cd /path/to/xTea             # repo root
pytest                       # all
pytest -m "not slow"         # fast only
pytest tests/test_concordance_comparator.py   # comparator units (no data needed)
```

## Inputs (override via environment variables)
All paths are overridable via environment variables so the suite is portable. When an env var is
not set, the harness falls back to a `$HOME`-relative location — no absolute per-user path is
committed. Tests requiring a missing input are skipped (not failed) and the skip message names
the env var that would fix it, so the pure-Python units always run regardless of external data
availability.

| Env var | Meaning | Fallback when unset |
|---|---|---|
| `XTEA2_TESTDATA` | root of test-data directory | `$HOME/scratch/xtea2_testdata` |
| `XTEA2_GOLDEN_DIR` | master golden VCFs (per-type subdirs) | `$XTEA2_TESTDATA/master_golden/NA12878` |
| `XTEA2_REFERENCE` | reference FASTA for CRAM decoding | `$HOME/hg38_no_alt.fa` |
| `XTEA2_TEST_BAM` / `XTEA2_TEST_CRAM` | BAM and CRAM read-level fixtures | `$XTEA2_TESTDATA/NA12878.chr1.bam` and `.../NA12878.final.cram` |
| `XTEA2_OUTPUT_DIR` | pipeline-under-test output VCFs | unset (activates real e2e gate when set) |
| `XTEA2_SCAN_BIN` | the C `xtea-scan` binary | unset (activates oracle-vs-scan tests when set) |

**To run on your machine:** set the env vars to point to your data before running pytest:
```bash
export XTEA2_TESTDATA=/path/to/your/testdata
export XTEA2_REFERENCE=/path/to/your/hg38_no_alt.fa
pytest
```

## Golden call counts (NA12878 chr1, `-f 5907`)
Alu **27**, SVA **5**, L1 **0**. The L1=0 means the L1 e2e path is untested on chr1 — a whole-genome
NA12878 golden (planned) will exercise it.
