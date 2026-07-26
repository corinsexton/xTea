# xTea2 handoff — 2026-07-25 (late evening session)

One-line status: **The environment is now a single lever (`conda activate xtea2`), a fast
byte-exact downstream gate exists (`tests/test_downstream_from_cns.py`), and the Python trim is
planned and prioritized in [`TRIM_PLAN.md`](TRIM_PLAN.md) — with the biggest win identified by
profiling, not guesswork.** Working tree has uncommitted changes (listed below).

Read [`TRIM_PLAN.md`](TRIM_PLAN.md) first — it is the work. This file is what happened and what
state things are in.

## Start here

```bash
conda activate xtea2        # that is the entire env setup now; nothing else to export
cd /n/data1/hms/dbmi/park/corinne/claude_sandbox/xTea
XTEA2_SCAN_BIN=$PWD/src/csrc/xtea_scan pytest -m "not slow"    # 69 passed, 1 skipped, ~92s
pytest -k Alu tests/test_downstream_from_cns.py -s             # the new downstream gate, one type
```

## 1. Environment: two levers → one (DONE)

`conda activate xtea2` is the whole setup. `PYTHONNOUSERSITE=1` now lives in the env's own config
(`conda env config vars set PYTHONNOUSERSITE=1 -n xtea2`), exported on activate and unset on
deactivate. The old two-lever rule is gone; `docs/environments.md` was rewritten around this.

- **22 scripts converted** under `/n/scratch/users/c/cos689/xtea2_testdata/` to the standard
  preamble. Backups in `.env_conversion_backup_2026-07-25/` there. Conversion was idempotent.
- The legacy py3.8 `xTea` env is referenced by nothing in the repo. Not deleted, just unused.
- **The one sharp edge:** `conda env create` does not carry env vars from YAML, so a rebuild must
  re-run the `config vars set` command. Documented in `envs/xtea2.yml`, the README, and
  `docs/environments.md`.
- **The other way to lose it:** calling `/n/.../envs/xtea2/bin/python` by absolute path instead of
  activating. Always activate.

Verified in a clean `env -i` non-login shell under `set -euo pipefail`: py3.10.20,
`PYTHONNOUSERSITE=1`, onnxruntime 1.22.0, samtools 1.22 / bwa 0.7.17 still from
`/home/cos689/software/bin`, all 58 modules import, `--dry-run` → 24 steps exit 0.

## 2. The fast downstream gate (NEW — this is the tool for the trim)

`tests/test_downstream_from_cns.py`. Stages the golden's own post-cns outputs and runs only the
tail (transduction → sibling → postF ×2 → gene → genotype → gVCF), comparing byte-for-byte.

- **The gate is VCF columns 1-8** (decision 3 below). Intermediate byte-differences are printed as
  diagnostics, not asserted — read the `NOTE — intermediate moved` lines.
- **59 files byte-identical** at the time of writing: L1 31, Alu 7, SVA 21.
- **~11 min** for all three types. Use `-k Alu` while iterating.
- Argv comes from `run_xtea.build_type_steps(..., start_from="transduction")` — a new parameter
  added this session, symmetric with the existing `stop_after` — so the harness cannot drift from
  the real chain.
- Comparison is **self-maintaining**: every file present on both sides is compared, so a new
  downstream output starts being gated automatically.
- Two tolerated exception classes, both enumerated in code. Everything else is byte-exact.

New conftest fixtures: `rep_lib_dir`, `gene_gff3`, `blacklist` (env vars `XTEA2_REP_LIB_DIR`,
`XTEA2_GENE_GFF3`, `XTEA2_BLACKLIST`), skip-if-missing per existing convention.

**Why it is ~11 min and not 2 min:** 96% of postF is interval-tree construction (P1 below). Fix
that and this gate becomes fast — which is the argument for doing P1 first.

## 3. Findings, all recorded in `docs/clipdisc_contract_spec.md` §4

- **P1 — the interval tree is 96% of postF.** Profiled: 104.8 s total, `construct_interval_tree`
  100.9 s, of which 77.8 s is inserting 955,824 intervals one at a time into the vendored
  pure-Python AVL tree. Actual filtering work: ~4 s. Rebuilt ~12× per 3-type run; cost is fixed
  (same on WGS as chr1). The annotations are static and queried read-only. **Deleting 2,040 lines
  and the biggest speedup available are the same change.** This is priority 1.
- **B6 — the DF21→ONNX genotype delta is 3 sites on the Alu chr1 golden, not 2.** The earlier
  handoff records two (`chr1:81058430`, `chr1:232452014`) over a 29-site set; over the Alu golden's
  own 27 calls it is those plus **`chr1:104688489`**. All `1/1`→`0/1`, no call gained or lost, VCF
  cols 1-8 byte-identical. Different site sets, so a refinement rather than a contradiction.
- **B7 — SVA transduction differs from the golden at the intermediate level, and it washes out.**
  `cns2.txt` 12 rows vs golden's 11 (extra `chr1:182305180`); `.high_confident` 9 vs 10 (gains that
  site, loses `chr1:50381610` and `chr1:182510697`). **Deterministic** — two back-to-back runs from
  identical staged input are byte-identical to each other, so this is a stable difference between
  our `x_TEA_main.py` and the golden's, not B3-style nondeterminism. postF filters all three: the
  final SVA VCF is byte-identical, 5 calls. **Attribute this before trimming transduction.**
- **D2 — mechanical dead-code inventory.** 13 modules / 7,940 lines unreachable; 739 commented-out
  code lines; 184 dead locals (clustered in the downstream: `x_orphan_transduction` 55,
  `x_post_filter` 24, `x_clip_disc_filter` 19); 26 unused imports. Caveat recorded: a dead local
  whose RHS is a call with side effects is not dead.

## 4. Uncommitted changes in the working tree

```
M  README.md                            env section -> one lever
M  docs/environments.md                 rewritten (235 lines changed)
M  docs/handoff/HANDOFF_2026-07-25.md   item 7 resolved + env-transition section appended
M  envs/xtea2.yml                       two-command build documented
M  tests/README.md                      one-line env fix
M  tests/conftest.py                    rep_lib_dir / gene_gff3 / blacklist fixtures
M  xtea/run_xtea.py                     start_from param on build_type_steps
M  docs/clipdisc_contract_spec.md       §4: B6, B7, P1, D2
A  tests/test_downstream_from_cns.py    the new gate
A  docs/handoff/TRIM_PLAN.md            the plan
A  docs/handoff/HANDOFF_2026-07-25_evening.md
```

Nothing has been committed — deliberately, pending review.

## 5. Decisions — ALL MADE (user, 2026-07-25). Nothing is blocked.

1. **Legacy entry points: DELETE all 5,879 lines** — four `gnrt_pipeline_*` + `x_TEA_main_germline_cloud.py`.
2. **Scan-only: YES** (~2,700 lines — `x_TEI_locator.py`, much of `clip_read.py`, the `-C`/`-D`
   branches). Do this *after* attributing B7, since B7's attribution wants master's untouched
   clip/disc path still present.
3. **Gate scope: call-set concordance only** (VCF cols 1-8). Already implemented. The byte-level
   intermediate comparison still runs as *diagnostics* — kept deliberately, because that is how B7
   was found. A green gate no longer means nothing changed; read the NOTE lines in the log.
4. **Package layout: LATER** — keep `xtea/`, move to `src/xtea/` + pyproject after the trim.

**Execution order (from `TRIM_PLAN.md`):**

```
1. pure deletions   -8,679   risk-free; do first, makes everything after readable
2. interval tree    -1,990   biggest speedup; gate drops to ~2-3 min after this
3. attribute B7        --    before scan-only, while master's path still exists
4. scan-only        -2,700
5. dead locals      small    check every RHS for side effects
                    ------
         32,963 ->  ~19,600  (40% smaller)
```

## 6. Loose ends

- **`src/csrc/xtea_scan` is stale** relative to `xtea_scan.cpp` (source 17:56, binary 12:09).
  Usage-string-only change, deferred earlier in the day while `48774776` was running. Rebuild:
  `PKG_CONFIG_PATH=/home/cos689/software/lib/pkgconfig make` in `src/csrc`.
- `48774776` (CC.3, LFS 2921) — **COMPLETED**, 04:52:43, exit 0:0. Result not yet read; the log is
  `/n/scratch/users/c/cos689/xtea2_testdata/cc3_lfs_2921_all3_scan_bothfix_*.out`. Read assertion 1b
  in `HANDOFF_2026-07-25.md` BEFORE interpreting the Alu result.
