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

## Running
```bash
conda activate xtea2         # python + pysam + pytest + numpy
cd /home/cos689/data1/corinne/claude_sandbox/xTea
pytest                       # all
pytest -m "not slow"         # fast only
pytest tests/test_concordance_comparator.py   # comparator units (no data needed)
```

## Inputs (override via env; default to the O2 NA12878 chr1 golden)
| Env var | Meaning | Default |
|---|---|---|
| `XTEA2_TESTDATA` | test-data dir | `/home/cos689/scratch/xtea2_testdata` |
| `XTEA2_GOLDEN_DIR` | master golden (per-type subdirs) | `$XTEA2_TESTDATA/master_golden/NA12878` |
| `XTEA2_REFERENCE` | reference FASTA (CRAM decode) | `/home/cos689/hg38_no_alt.fa` |
| `XTEA2_TEST_BAM` / `XTEA2_TEST_CRAM` | read-level fixtures | chr1 BAM / full CRAM |
| `XTEA2_OUTPUT_DIR` | pipeline-under-test VCFs → activates the real e2e gate | unset |
| `XTEA2_SCAN_BIN` | the C `xtea-scan` binary → activates oracle-vs-scan | unset |

Tests skip (not fail) when an input is missing, so the pure-Python units always run.

## Golden call counts (NA12878 chr1, `-f 5907`)
Alu **27**, SVA **5**, L1 **0**. The L1=0 means the L1 e2e path is untested on chr1 — a whole-genome
NA12878 golden (planned) will exercise it.
