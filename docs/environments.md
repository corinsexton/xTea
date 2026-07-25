# xTea2 runtime environments

**The rule, in one line:** run everything with the `xtea2` python **and** `PYTHONNOUSERSITE=1`. Reach
for the legacy `xTea` env only when the baseline itself is the thing you are reproducing.

Both halves are load-bearing. Neither one alone works, and the reason is not obvious — the second
requirement is *created* by the first. This document exists because getting it wrong costs a whole
run: it has already cost one 18-minute chr1 gate (job `48777874`, SIGSEGV) and it will silently break
any script that produces a final VCF.

## The two environments

| | `xTea` (legacy) | `xtea2` |
|---|---|---|
| python | 3.8.15 | 3.10.20 |
| numpy | 1.19.5 | 2.2.6 |
| onnxruntime | **none, and cannot have one** | 1.22.0 |
| detection + downstream | works | works |
| **genotyping** | **impossible** | works |
| pinned spec | none | `envs/xtea2.yml` |
| role | frozen provenance of every golden | the pipeline environment |

`xTea` cannot simply be given onnxruntime: every onnxruntime release requires `numpy>=1.21`, so
installing one means dragging that env off numpy 1.19.5 — the exact interpreter that produced every
golden baseline. It is deliberately frozen instead. The consequence is that **`xTea` can run
everything except genotyping**, and since genotyping is on the critical path to the final VCF
(`genotype` → `gVCF`), it can no longer complete a germline run at all.

External tools are not conda-managed and are expected on `PATH` from `/home/cos689/software/bin`:
`samtools` (1.22), `bwa` (0.7.17). Neither conda env ships them, so switching interpreters does not
change external-tool resolution.

## Why two levers are needed

They answer two different questions.

**The interpreter answers "is there an onnxruntime at all?"** `xTea`'s py3.8 has none. So you must
switch to `xtea2`.

**`PYTHONNOUSERSITE=1` answers "which onnxruntime wins?"** Python automatically prepends the
*user site-packages* directory to `sys.path`, ahead of the environment's own. Measured:

```
$ python -c "import sys; [print(p) for p in sys.path if 'site-packages' in p]"
/home/cos689/.local/lib/python3.10/site-packages          <-- FIRST
  ... 12 more ~/.local entries (intervaltree, TEToolkit, bamsnap, ...) ...
/n/data1/.../envs/xtea2/lib/python3.10/site-packages      <-- LAST

$ PYTHONNOUSERSITE=1 python -c "..."
/n/data1/.../envs/xtea2/lib/python3.10/site-packages      <-- only entry
```

`~/.local` carries **onnxruntime 1.16.3** (the `cli` branch's pin), compiled against numpy 1.x. It
therefore shadows the env's working 1.22.0, and dies under numpy 2.2.6.

**The two are linked by the version number.** `~/.local/lib/`**`python3.10`** is keyed to the
interpreter's version — and 3.10 is exactly `xtea2`'s version. Switching to the `xtea2` interpreter
is *what makes that directory visible in the first place*; a py3.8 interpreter looks at
`~/.local/lib/python3.8`, which does not contain the bad wheel. So the fix for the first problem is
what exposes you to the second.

## Why the environment variable, and not `python -s`

`-s` and `PYTHONNOUSERSITE=1` both disable user-site for the process they are given to. Only one of
them survives the pipeline's structure:

`run_xtea.py:186` runs every step as a fresh subprocess, `[sys.executable, XTEA_MAIN] + argv`, with
no interpreter flags carried over. A `-s` on the outer invocation would protect only the parent and
be silently lost in all 24 step subprocesses — including the genotype step, which is the one that
matters. `subprocess.run` inherits the parent's environment by default, so an environment variable
reaches every step. Export it once; do not pass a flag.

(This also means the variable must be exported in the *submitting script*, not typed on the
`run_xtea.py` command line.)

## It is not only onnxruntime

Thirteen `~/.local` entries precede the env's site-packages, including `intervaltree`, which the
pipeline imports via `disc_adapter` → `x_annotation`. Without the variable the environment is not
self-contained for *anything*.

This is the mirror image of a problem found earlier the same day: `~/.local` had been *hiding missing*
packages, making `xtea2` look complete while it silently borrowed `sortedcontainers`, `scipy`,
`scikit-learn` and `pandas` from there. The same mechanism can equally well *override correct*
packages with broken ones. Hence the standing rule below.

**When validating any environment claim, set `PYTHONNOUSERSITE=1`.** A suite that passes without it
proves nothing about whether the environment is buildable anywhere else. (Reference point: the suite
gives identical results either way — this is how the env looked healthy for months while being
unbuildable.)

## Recognising the two failure modes

Both are the same root cause — the wrong onnxruntime, or none — but they look nothing alike:

| setup | symptom |
|---|---|
| `xTea` py3.8 (no onnxruntime) | `ModuleNotFoundError: No module named 'onnxruntime'` at the genotype step. The import is lazy (inside `predict_for_site_onnx`), so it fails at *predict* time, not import time — the run gets all the way through detection and downstream first. |
| `xtea2` py3.10 **without** the variable | `ERROR: step 'shared:genotype' failed (exit -11)` — SIGSEGV. May leave a 100 MB–900 MB `core.*` dump (gitignored via `core.*`). Imported interactively it instead shows `AttributeError: _ARRAY_API not found`. |
| `xtea2` py3.10 **with** the variable | works — onnxruntime 1.22.0, genotype and gVCF clean. |

A germline run that reaches `shared:genotype` and dies there is almost always this, not a pipeline
bug: everything before that step will have already succeeded.

## Which environment for which task

- **All pipeline runs, gates, and tests** → `xtea2` + `PYTHONNOUSERSITE=1`. This is the default and
  covers essentially all real work.
- **Reproducing or regenerating a golden** → the legacy `xTea` env. The goldens came from master's
  generated pipeline under py3.8/numpy 1.19.5, and that env is frozen precisely so they stay
  reproducible. Running a golden under `xtea2` would conflate "did my change break this" with "did
  the interpreter change this."
- **Attribution / bisection runs** → pin the interpreter deliberately as a controlled variable rather
  than defaulting. Worked example: attributing a one-call VCF difference required a run that varied
  *only* the scan binary, which is why `a2_validate_chr1_cnsemit.sh` takes the scan binary as a
  parameter. Empirical result from that run, worth knowing: the py3.8→3.10 / numpy 1.19→2.2 change
  produced **zero** differences in VCF columns 1-8 on chr1 — it is output-neutral outside genotyping.
- **DF21 genotypes specifically** → not possible in `xtea2`. `deepforest` ships only cp36–cp39
  wheels with no sdist, so it is py3.8-only; the model and `x_genotype_classify.py` were deleted and
  are recoverable from git history. Not a live need.

## Script audit — 2026-07-25

Scripts under `/n/scratch/users/c/cos689/xtea2_testdata/` invoke the pipeline with a bare `python`
resolved from a `PATH` that prepends the **`xTea` py3.8** env. Any of them that reaches the genotype
step is therefore broken until converted:

| script | reaches genotype? | status |
|---|---|---|
| `a2_validate_chr1_cnsemit.sh` | yes | **converted** — `PY=` + `export PYTHONNOUSERSITE=1` |
| `a0_validate_chr1`, `a1_validate_chr1`, `a1_validate_wgs` | yes | needs conversion |
| `a2_validate_chr1`, `a2_validate_wgs`, `a2_validate_wgs_cnsemit` | yes | needs conversion |
| `a2_moveA2_wgs_timing`, `cc0_validate_chr1` | yes | needs conversion |
| `cc1_validate_chr1_selfvself`, `cc2_validate_wgs_selfvself`, `cc3_lfs_2921_*` | **no** | fine on py3.8 |

The case-control scripts are unaffected because case-control mode skips genotyping outright
(`run_xtea.py:361-365`). That is why CC.3 (`48774776`) ran correctly under py3.8 throughout.

Conversion is mechanical — see `a2_validate_chr1_cnsemit.sh` for the pattern:

```bash
PY=/n/data1/hms/dbmi/park/corinne/.conda/envs/xtea2/bin/python
export PYTHONNOUSERSITE=1
...
$PY xtea/run_xtea.py ...
```

Leave `PATH` pointing at the `xTea` env: `samtools`/`bwa` come from `/home/cos689/software/bin`
either way, so keeping it unchanged holds external-tool resolution identical to the golden runs.

## Reducing it to one lever (optional)

Removing or renaming `~/.local/lib/python3.10/site-packages/onnxruntime` (1.16.3) would make the
`xtea2` interpreter sufficient for the genotype case, since nothing would be left to shadow 1.22.0.
It is unusable under numpy 2 regardless, so it has no working use in this project.

Not done here: it lives outside the repo and other projects may want it under a different numpy.
It would also not fix the general shadowing of the other twelve `~/.local` packages, so scripts
should keep exporting the variable either way.

## Rebuilding the environment

```bash
mamba env create -f envs/xtea2.yml
```

Verified to solve to exactly the installed versions. Rebuild from that file rather than by hand.
