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
| scipy | the genotype step (`x_genotype_lik`, numpy/scipy only since onnxruntime was dropped) |
| pysam | 0.24.0 |
| `PYTHONNOUSERSITE` | `1`, set by the env itself |
| pinned spec | `envs/xtea2.yml` |

External tools are deliberately **not** conda-managed and are expected on `PATH`: `samtools`
(1.22 verified), `bwa` (0.7.17 verified). These are the same binaries the golden runs used, so
external-tool resolution is unchanged from the goldens' provenance — which is the reason to keep
them outside the env rather than an oversight.

A batch script therefore needs this preamble:

```bash
set -euo pipefail
source <conda-install>/etc/profile.d/conda.sh
conda activate xtea2
export PATH=<dir-holding-samtools-and-bwa>:$PATH
```

Concrete values for this group's cluster are in `docs/o2_setup.md`, which is deliberately
untracked — site paths do not belong in a repo that gets pushed.

The `source` line is required because `sbatch` scripts are not login shells, so `conda` is not
otherwise a defined function. Verified to work under `set -u` in a fully clean environment
(`env -i`).

## Why `PYTHONNOUSERSITE=1` is baked into the env

This is the one non-obvious thing, and it matters when the env is ever rebuilt.

Python automatically prepends the *user site-packages* directory to `sys.path`, ahead of the
environment's own. `~/.local/lib/python3.10/site-packages` on this system carries **onnxruntime
1.16.3** — compiled against numpy 1.x, and it segfaulted the old ONNX genotype step under this env's
numpy 2.2.6. Genotyping no longer imports onnxruntime (`x_genotype_lik` replaced the ONNX model),
so that crash is gone, but twelve other `~/.local` packages shadow the env the same way; the
environment is not self-contained for *anything* without this.

The variable is registered into the environment itself:

```bash
conda env config vars set PYTHONNOUSERSITE=1 -n xtea2
```

`conda activate` exports it, `conda deactivate` unsets it. This is why activation alone is
sufficient, and it is the reason the old two-lever rule is gone.

**The one way to lose it:** invoking the env's python by absolute path
(`/n/.../envs/xtea2/bin/python`) instead of activating. That bypasses `conda activate`, so the
variable is never set and `~/.local` packages shadow the env again. Always activate; never
hard-code the interpreter path.

An environment variable rather than `python -s` because `run_xtea.py:186` runs every step as a
fresh subprocess (`[sys.executable, XTEA_MAIN] + argv`) with no interpreter flags carried over. A
`-s` on the outer invocation would protect only the parent and be silently lost in all 24 step
subprocesses — including genotype, the one that matters. `subprocess.run` inherits the parent's
environment, so an exported variable reaches every step.

## Recognising failure

| symptom | cause |
|---|---|
| `ERROR: step 'shared:genotype' failed (exit -11)` — SIGSEGV | an xTea2 build from before `x_genotype_lik` (ONNX genotyper) with `PYTHONNOUSERSITE` lost. Current builds do not import onnxruntime. |

## Rebuilding

Two commands — the second keeps `~/.local` off `sys.path` (see above):

```bash
mamba env create -f envs/xtea2.yml
conda env config vars set PYTHONNOUSERSITE=1 -n xtea2
```

`conda env create` does not carry env vars from the YAML, so a rebuild that skips the second command
produces an environment that silently borrows `~/.local` packages. Verify with:

```bash
conda activate xtea2
python -c "import sys; print([p for p in sys.path if 'site-packages' in p])"   # env path only
```

**When validating any environment claim, confirm `~/.local` is off `sys.path` first.** A suite that
passes with user-site enabled proves nothing about whether the environment is buildable anywhere
else — this is how the env looked healthy for months while silently borrowing `sortedcontainers`,
`scipy`, `scikit-learn` and `pandas` from `~/.local`.

## The legacy `xTea` env

The legacy py3.8 env (python 3.8.15, numpy 1.19.5) may still exist on a given machine and is
**no longer required by anything in this repo**. It is not deleted, but nothing should reach for
it.

It is not supported for xTea2 (Python 3.10 is the floor).

It was previously kept as the provenance of the goldens. That is no longer a reason to use it:
an attribution run varying *only* the interpreter measured **zero** differences in VCF columns 1-8
on chr1 across py3.8/numpy 1.19 → py3.10/numpy 2.2. The interpreter change is output-neutral
outside genotyping, so goldens can be regenerated under `xtea2`.

`deepforest` — the old DF21 genotyping model's dependency — ships only cp36–cp39 wheels with no
sdist, so it is py3.8-only and cannot exist here. The model and `x_genotype_classify.py` were
deleted (recoverable from git history) when genotyping moved to ONNX in `95b3ff9`; the ONNX model
and `x_genotype_classify_sklearn.py` were in turn replaced by `x_genotype_lik.py`. Not a live need.

## Script conversion — done 2026-07-25

All 22 pipeline-invoking scripts under the scratch test-data tree were converted to
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
