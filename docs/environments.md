# xTea2 runtime environment

**The rule, in one line:**

```bash
conda activate xtea2
```

That is the entire environment setup. There is nothing else to export, no interpreter path to
hard-code, no flag to pass. Everything — the pipeline, the gates, the test suite, golden
regeneration — runs in this one environment.

## What activation gives you

| | |
|---|---|
| python | 3.10.20 |
| numpy | 2.2.6 |
| onnxruntime | 1.22.0 (the genotype step) |
| pysam | 0.24.0 |
| `PYTHONNOUSERSITE` | `1`, set by the env itself |
| pinned spec | `envs/xtea2.yml` |

External tools are deliberately **not** conda-managed and are expected on `PATH` from
`/home/cos689/software/bin`: `samtools` (1.22), `bwa` (0.7.17). These are the same binaries the
golden runs used, so external-tool resolution is unchanged from the goldens' provenance.

A batch script therefore needs exactly this preamble:

```bash
set -euo pipefail
source /n/app/conda/miniforge3/24.11.3-0/etc/profile.d/conda.sh
conda activate xtea2
export PATH=/home/cos689/software/bin:$PATH
```

The `source` line is required because `sbatch` scripts are not login shells, so `conda` is not
otherwise a defined function. Verified to work under `set -u` in a fully clean environment
(`env -i`).

## Why `PYTHONNOUSERSITE=1` is baked into the env

This is the one non-obvious thing, and it matters when the env is ever rebuilt.

Python automatically prepends the *user site-packages* directory to `sys.path`, ahead of the
environment's own. `~/.local/lib/python3.10/site-packages` on this system carries **onnxruntime
1.16.3** — compiled against numpy 1.x, and it segfaults under this env's numpy 2.2.6. Left alone it
shadows the env's working 1.22.0. Twelve other `~/.local` packages shadow the env the same way; the
environment is not self-contained for *anything* without this.

The variable is registered into the environment itself:

```bash
conda env config vars set PYTHONNOUSERSITE=1 -n xtea2
```

`conda activate` exports it, `conda deactivate` unsets it. This is why activation alone is
sufficient, and it is the reason the old two-lever rule is gone.

**The one way to lose it:** invoking the env's python by absolute path
(`/n/.../envs/xtea2/bin/python`) instead of activating. That bypasses `conda activate`, so the
variable is never set and the broken onnxruntime comes back. Always activate; never hard-code the
interpreter path.

An environment variable rather than `python -s` because `run_xtea.py:186` runs every step as a
fresh subprocess (`[sys.executable, XTEA_MAIN] + argv`) with no interpreter flags carried over. A
`-s` on the outer invocation would protect only the parent and be silently lost in all 24 step
subprocesses — including genotype, the one that matters. `subprocess.run` inherits the parent's
environment, so an exported variable reaches every step.

## Recognising failure

| symptom | cause |
|---|---|
| `ERROR: step 'shared:genotype' failed (exit -11)` — SIGSEGV, possibly a 100 MB–900 MB `core.*` dump (gitignored) | `PYTHONNOUSERSITE` was lost — almost always because the interpreter was invoked by absolute path instead of activated. Interactively the same cause shows as `AttributeError: _ARRAY_API not found`. |
| `ModuleNotFoundError: No module named 'onnxruntime'` at the genotype step | running under the legacy py3.8 `xTea` env. The import is lazy (inside `predict_for_site_onnx`), so it fails at *predict* time — the run gets all the way through detection and downstream first. |

A germline run that reaches `shared:genotype` and dies there is almost always one of these, not a
pipeline bug: every preceding step will have already succeeded.

## Rebuilding

Two commands — the second is not optional:

```bash
mamba env create -f envs/xtea2.yml
conda env config vars set PYTHONNOUSERSITE=1 -n xtea2
```

`conda env create` does not carry env vars from the YAML, so a rebuild that skips the second command
produces an environment that passes its test suite and then segfaults at genotype. Verify with:

```bash
conda activate xtea2
python -c "import sys; print([p for p in sys.path if 'site-packages' in p])"   # env path only
python -c "import onnxruntime; print(onnxruntime.__version__)"                  # 1.22.0
```

**When validating any environment claim, confirm `~/.local` is off `sys.path` first.** A suite that
passes with user-site enabled proves nothing about whether the environment is buildable anywhere
else — this is how the env looked healthy for months while silently borrowing `sortedcontainers`,
`scipy`, `scikit-learn` and `pandas` from `~/.local`.

## The legacy `xTea` env

`/n/data1/hms/dbmi/park/corinne/.conda/envs/xTea` (python 3.8.15, numpy 1.19.5) is still on disk and
is **no longer required by anything in this repo**. It is not deleted, but nothing should reach for
it.

It cannot run xTea2 at all: every onnxruntime release requires `numpy>=1.21`, so it can never have
one, and genotyping is on the critical path to the final VCF (`genotype` → `gVCF`).

It was previously kept as the provenance of the goldens. That is no longer a reason to use it:
an attribution run varying *only* the interpreter measured **zero** differences in VCF columns 1-8
on chr1 across py3.8/numpy 1.19 → py3.10/numpy 2.2. The interpreter change is output-neutral
outside genotyping, so goldens can be regenerated under `xtea2`.

`deepforest` — the old DF21 genotyping model's dependency — ships only cp36–cp39 wheels with no
sdist, so it is py3.8-only and cannot exist here. The model and `x_genotype_classify.py` were
deleted (recoverable from git history) when genotyping moved to ONNX in `95b3ff9`. Not a live need.

## Script conversion — done 2026-07-25

All 22 pipeline-invoking scripts under `/n/scratch/users/c/cos689/xtea2_testdata/` were converted to
the preamble above: the legacy `export PATH=.../envs/xTea/bin:...` line, the `PY=` absolute-path
assignments, the hand-rolled `export PYTHONNOUSERSITE=1`, and `conda activate <path-to-xTea>` are
all gone. Backups of the pre-conversion versions are in
`.env_conversion_backup_2026-07-25/` in that directory.

Verified after conversion: all 27 scripts parse (`bash -n`), no residual references to the legacy
env or to an absolute interpreter path, and a `--dry-run` through the converted preamble emits 24
steps and exits 0.

`fetch_chr1.sh`, `fetch_full.sh`, the `monitor_*.sh` scripts and `setup_seeded_golden.sh` were left
alone — none of them invokes a python interpreter. `setup_seeded_golden.sh` rewrites master's own
generated `run_xTEA_pipeline.sh` golden scripts, which are frozen provenance and are not converted.
