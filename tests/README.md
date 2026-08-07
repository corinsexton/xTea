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
- `test_downstream_from_cns.py` — **the fast downstream gate** (`slow` marker). Stages the golden's
  own post-cns outputs and runs only the tail (transduction → sibling → postF ×2 → gene → genotype
  → gVCF), ~11 min for all three types vs ~18 for the full chr1 gate. **Gate = VCF columns 1-8**;
  every intermediate is still compared byte-for-byte and printed as diagnostics. Use
  `pytest -k Alu tests/test_downstream_from_cns.py -s` while iterating. Needs `XTEA2_REP_LIB_DIR`,
  `XTEA2_GENE_GFF3`, `XTEA2_BLACKLIST` (plus the golden dir, BAM and reference).
- `test_scan_multifile.py` — multi-file scanning in the C++ engine: pooling a name-split pair of
  files reproduces the unsplit scan exactly, and each way `@SQ` can differ (count, name/order,
  length) fails loud naming the bad file. Pooling tests need `XTEA2_TEST_BAM` + `XTEA2_SCAN_BIN`;
  the `@SQ` tests need only the binary.

## Running
```bash
conda activate xtea2         # the whole env setup; build it per docs/environments.md
cd /path/to/xTea             # repo root
make -C src/csrc             # BUILD xtea-scan FIRST -- see below; ~30 skips without it
pytest                       # all
pytest -m "not slow"         # fast only
pytest tests/test_concordance_comparator.py   # comparator units (no data needed)
```

### Build `xtea-scan` before you trust a green run
About 30 tests cover the C++ engine that owns `scan_pass` — the single largest step in a WGS run.
They resolve the binary from **`src/csrc/xtea_scan`**, the shipped Makefile's own output, so a plain
`make -C src/csrc` is all that is needed (`PKG_CONFIG_PATH=<htslib-prefix>/lib/pkgconfig` if
htslib's `.pc` is not already on the path). `$XTEA2_SCAN_BIN` still overrides.

Three outcomes, and the difference matters:

| binary | outcome | why |
|---|---|---|
| built, newer than its sources | **runs** | the only state in which a pass means anything |
| absent | **skips**, naming the `make` line | not-built is legitimate; the fast units must still run |
| **older than any of `xtea_scan.cpp`, `realign.hpp`, `cns_emit.hpp`, `td_gather.hpp`** | **FAILS** | a stale binary is the one state that *lies* — the suite goes green against code that is not in the tree |

The staleness check is not hypothetical housekeeping. The in-tree binary was once found five days
older than `xtea_scan.cpp`, and several harnesses under `docs/handoff/runs/` default
`XTEA2_SCAN_BIN` to a **different commit's** worktree build — so inheriting someone else's default
is enough to test a binary that never matched these sources. It watches the headers as well as the
`.cpp` because they carry most of the code and are `#include`d rather than compiled separately,
which is the same reason the Makefile lists them as prerequisites.

Set `XTEA2_ALLOW_STALE_SCAN_BIN=1` to downgrade the failure to a warning, for the deliberate case
of testing a frozen binary against a moved tree.

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
| `XTEA2_SCAN_BIN` | the C `xtea-scan` binary | **`src/csrc/xtea_scan`, the in-tree build** — absent skips, stale FAILS (see above) |
| `XTEA2_ALLOW_STALE_SCAN_BIN` | downgrade a stale binary from failure to warning | unset (stale is a failure) |
| `XTEA2_REP_LIB_DIR` | `rep_lib_annotation` dir (downstream gate) | `$HOME/xTea/rep_lib_annotation` |
| `XTEA2_GENE_GFF3` | gene GFF3 for the annotation step | `$HOME/data1/corinne/ref/gencode.v45.basic.annotation.gff3` |
| `XTEA2_BLACKLIST` | blacklist BED | `$XTEA2_REP_LIB_DIR/blacklist/hg38/centromere.bed` |

**To run on your machine:** set the env vars to point to your data before running pytest:
```bash
export XTEA2_TESTDATA=/path/to/your/testdata
export XTEA2_REFERENCE=/path/to/your/hg38_no_alt.fa
pytest
```

## Golden call counts (NA12878 chr1, `-f 5907`)
Alu **27**, SVA **5**, L1 **0**. The L1=0 means the L1 e2e path is untested on chr1 — a whole-genome
NA12878 golden (planned) will exercise it.
