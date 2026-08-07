# xTea2 handoff — 2026-08-07 (second session): three replicates per arm, and the accuracy question is settled

**This is the current entry point for state.** It supersedes `HANDOFF_2026-08-07_objective_met.md`
for the speed/memory/accuracy numbers only; that document remains the entry point for the agenda
rebase and for §6's parked list. The agenda entry point is still
`docs/handoff/DOWNSTREAM_OVERHAUL_PLAN.md`.

**The one-sentence version: master and xTea2 were each run three times, and the result is that
xTea2 is bit-reproducible where master is not, its precision edge is real at every family, and its
entire systematic sensitivity deficit is three SVA insertions.**

---

## 0. State at handoff

| | |
|---|---|
| branch | `xTea2` |
| commits | `ea08f45` the C++ test gap; `e568774` the replicate harness + aggregator; this handoff |
| `xtea/` changed | **nothing.** No shipping code moved this session |
| tests | **372 passed / 1 skipped** on a plain `pytest -m "not slow"`, was 324/30 |
| jobs | array **`50062602`**, 6 tasks, all `COMPLETED 0:0` |
| what is NOT done | the 3 SVA misses in §4.3 are diagnosed to a coordinate, not fixed; `HANDOFF_2026-08-07_objective_met.md` §6 items 2–8 remain parked |

## 1. THE C++ TOOL'S 30 TESTS WERE SKIPPING IN SILENCE (`ea08f45`)

`xtea2_scan_bin` resolved `$XTEA2_SCAN_BIN` **with no default**, so a binary built in-tree by the
shipped Makefile was invisible to pytest. A plain `python3 -m pytest -m "not slow"` therefore
skipped all 30 tests of the C++ engine that owns `scan_pass` — 1,600 s and the 11.64 GB step that
first raises the WGS high-water — and reported it as *"the C xtea-scan tool is not built yet"*,
which was false. **The tool builds and every one of those tests passes.** Nothing in the 30 needed
fixing.

Now defaults to `src/csrc/xtea_scan`, repo-relative. **324 passed / 30 skipped → 372 passed / 1
skipped**, the remaining skip being `test_xtea2_vs_golden`, which wants `XTEA2_OUTPUT_DIR`.

### 1.1 And a second defect, which mattered more

Nothing compared the binary's mtime to its sources. **The in-tree binary was dated Aug 2, five days
older than `xtea_scan.cpp`**, and 22 harnesses under `docs/handoff/runs/` default
`XTEA2_SCAN_BIN` to a **different commit's** worktree build — so testing a binary that never
matched the sources under test takes only inheriting someone else's default.

The contract is now three-valued, and the middle case is the point:

| binary | outcome | why |
|---|---|---|
| built, newer than sources | **runs** | the only state in which a pass means anything |
| absent | **skips**, naming the `make` line | not-built is legitimate; the fast units must still run |
| older than any of `xtea_scan.cpp`, `realign.hpp`, `cns_emit.hpp`, `td_gather.hpp` | **FAILS** | a stale binary is the one state that *lies* |

It watches the headers, not just the `.cpp`, for the Makefile's own stated reason: they carry most
of the code and are `#include`d rather than compiled separately. `XTEA2_ALLOW_STALE_SCAN_BIN=1`
downgrades to a warning for the deliberate frozen-binary case.

Verified on the real tree, not only synthetically: `touch src/csrc/realign.hpp` turns the four
`test_scan_tool.py` tests into errors naming `realign.hpp` and the rebuild; plain `make -C
src/csrc` rebuilds via the `HDR` prerequisite and they pass again.

**19 tests, 9 mutations, all bite. Two survived first and both were the TEST's fault** — recorded
in `tests/test_scan_bin_freshness.py` because they are general traps:

- **fail→skip survived `pytest.raises(FAILED)`.** `Skipped` derives from `BaseException`, so it
  leaked past the `raises()`, pytest reported the test as *skipped*, and the run still exited 0.
  The helper now catches `BaseException` and asserts the outcome type.
- **shrinking `SCAN_BIN_SOURCES` to just the `.cpp` survived** a test that looped over
  `SCAN_BIN_SOURCES` — the expectation was derived from the constant under mutation. Hardcoded.

## 2. THE REPLICATE EXPERIMENT (`e568774`, array `50062602`)

Corinne's direction: master vs xTea2, three runs each, on precision/recall, speed and memory. Same
BAM, `-c 24`, `--mem=64G`, priopark, L1+Alu+SVA. All six `COMPLETED 0:0`, all on **the same CPU
model** (Xeon Platinum 8468).

Two design choices were hers and both go against the local default:

1. **The seed is UNPINNED.** Every prior xTea2 gate exported `XTEA_BASIC_INFO_SEED=42`; master has
   no such knob, so pinning one arm would hold a source of variation fixed on that arm only and
   xTea2 would look more reproducible for a reason that is not the pipeline. The harness unsets it
   and **asserts** it is unset. **Consequence: these runs are NOT drop-in comparable with
   `49792012`/`49879677`/`50038678`, which were pinned.**
2. **All six ran concurrently**, so both arms saw identical cluster conditions and identical
   contention on the one 44 GB BAM. The walls are therefore under-load, not solo.

The master arm **reuses `master_wgs_hg002.sbatch`** with `MASTER_OUT` repointed rather than
reimplementing master's chain — that file drives the golden's own `run_xTEA_pipeline.sh` command
lines, which are the output of master's own `gnrt_pipeline_local_v38.py`. The xTea2 arm passes
**neither decoy flag**, and check 0 asserts both banners name the **DEFAULT** as their cause;
it passed on all three replicates.

`docs/handoff/runs/probe_repro_comparison.py` reproduces everything below offline in seconds. It
imports `load_truth_family` / `read_calls_vcf` / `match_one_to_one` from `benchmark_truth_set`
rather than copying them (plan §10), and was validated against two by-hand scores before any new
data existed.

## 3. SPEED AND MEMORY

| arm | wall (n=3) | peak RSS (n=3) |
|---|---|---|
| master `46520ab` | 1:03:18 – **1:21:18**, mean 1:09:18 | 38.12 – 38.19 GB, mean **38.15** |
| xTea2 `ea08f45` | 0:37:52 – 0:39:45, mean 0:39:07 | 17.83 – 18.21 GB, mean **18.00** |

**xTea2 uses 47% of master's peak memory, and the ranges are 20 GB apart.** Read it the strong way:
master's peak is **one family resident** (it loops L1/Alu/SVA serially) and xTea2's is **all
three**.

### 3.1 QUOTE 1.67×, AND HERE IS WHY IT IS THE DEFENSIBLE NUMBER

The naive cross-node ratio spans 1.59×–2.15× because two co-tenancy effects distort wall clock:

- **The lone job was the SLOWEST.** master rep1, alone on `-177`, took **1:21:18**; master rep2/rep3
  sharing `-178` took **1:03:18**. Co-tenancy *helped* — six jobs streaming the same 44 GB BAM warm
  its page cache for each other. So a cross-node ratio mixes the arm difference with how many
  readers each node had.
- **Same-node same-arm replicates run in LOCKSTEP.** All six started at 16:13:31; master rep2 and
  rep3 both ended at 17:16:49, xTea2 rep2 and rep3 both at 16:53:16 — identical to the second,
  twice. Such a pair is **one** independent wall measurement, so the effective n for wall is **4
  (node, arm) groups over 6 runs**, not 6. Peak RSS is *not* lockstep (38.19/38.14/38.12) and keeps
  its full n=3.

**Node `-178` hosted BOTH arms** — same hardware, same instant, same page cache — and there
master 1:03:18 vs xTea2 0:37:52 = **1.67×**. That agrees with the solo n=1 estimate of 1.66–1.82×
from `HANDOFF_2026-08-07_objective_met.md` §1, which is the useful convergence.

### 3.2 A CORRECTION: the "3.4 GB noise floor" was not noise

`HANDOFF_2026-08-07_objective_met.md` §3.4 established a noise floor of "199 s and 3.41 GB" from
jobs `49879677` vs `50038678`, two runs of the same configuration, and used it to discipline every
other memory number. **With three genuine replicates the peak-RSS spread is 0.07 GB (master) and
0.38 GB (xTea2)** — an order of magnitude tighter. Peak RSS on this pipeline is nearly
deterministic, so **something else differed between those two jobs**, and calling the gap a noise
floor hid that rather than explaining it. The wall half of that finding stands and is now explained:
wall noise is real and it is **co-tenancy** (§3.1).

Do not cite 3.4 GB as a noise floor again. The `+223 s / +0.41 GB` conclusion it was invented to
refute is still correctly refuted — by §3.1's co-tenancy effect and by the ~55 s in-job null-arm
decomposition in `HANDOFF_2026-08-06_timing_resolved.md`, which remains the right instrument.

## 4. PRECISION AND RECALL — the part that needed n=3

Positions, ±200 bp, against the HG002 truth set. Per replicate:

| family | arm | calls | TP | FN | FP | sens | prec |
|---|---|---|---|---|---|---|---|
| L1 (197) | master | 195 / 194 / 196 | 140 / 141 / 141 | 57 / 56 / 56 | 55 / 53 / 55 | **71.4%** mean | **72.1%** mean |
| | xTea2 | **188 / 188 / 188** | **140 / 140 / 140** | 57 | **48 / 48 / 48** | 71.1% | **74.5%** |
| Alu (1355) | master | 1227 / 1231 / 1228 | 1103 / 1105 / 1104 | 252 / 250 / 251 | 124 / 126 / 124 | **81.5%** | 89.9% |
| | xTea2 | **1223 ×3** | **1103 ×3** | 252 | **120 ×3** | 81.4% | **90.2%** |
| SVA (90) | master | 89 / 89 / 88 | **80 ×3** | 10 | 9 / 9 / 8 | **88.9%** | 90.2% |
| | xTea2 | **84 ×3** | **77 ×3** | 13 | **7 ×3** | 85.6% | **91.7%** |
| **pooled (1642)** | master | | 1323 / 1326 / 1325 | | 188 / 188 / 187 | **80.67%** | 87.59% |
| | xTea2 | | **1320 ×3** | | **175 ×3** | 80.39% | **88.29%** |

### 4.1 xTea2 IS BIT-REPRODUCIBLE. MASTER IS NOT. That is the headline.

| arm | L1 | Alu | SVA |
|---|---|---|---|
| master, truth records recovered by all 3 / union | 140 / **141** | 1100 / **1108** | 80 / 80 |
| master, call positions shared by all 3 | 97.5% | 98.7% | 98.9% |
| **xTea2**, truth records all 3 / union | **140 / 140** | **1103 / 1103** | **77 / 77** |
| **xTea2**, call positions shared by all 3 | **100%** | **100%** | **100%** |

Master varies run to run — **1 L1 and 8 Alu truth records** are recovered by some replicates and
not others. This is the direct measurement of the thread-sensitivity job `49075238` inferred from
master's `-n 8` output differing from its own `-n 24`. xTea2 produced **byte-identical call sets
across all three unpinned runs, in all three families.**

### 4.2 So most of the apparent accuracy gap was master's own noise

The n=1 comparison said xTea2 was −4 TP and −9 FP. Judged against each arm's own spread:

| family | TP delta | verdict |
|---|---|---|
| L1 | −0.7 | **explained** by master's own range (1) |
| Alu | −1.0 | **explained** by master's own range (2) |
| **SVA** | **−3.0** | **both arms exactly reproducible → SYSTEMATIC** |
| pooled | −4.7 | outside the observed range (3) but within 2×; n=3 understates the range |

And the precision edge is **not** noise — xTea2 has fewer false positives at every family, in
every replicate: pooled **175 vs 187.7**.

### 4.3 The whole systematic deficit is three SVA insertions, named

xTea2 misses these in **all three** replicates; master recovers them in **all three**:

```
SVA  chr1:88507126
SVA  chr14:20707961
SVA  chr5:56152167
```

**xTea2 gains zero SVA truth records master misses** — on SVA it is a strict subset. These predate
both decoy flips (SVA TP has been 77 across the structural flip and the concat flip alike), so they
are not fallout from the 2026-08-05 / 2026-08-07 defaults. **This is the one real, reproducible
sensitivity regression in the port, it is three loci wide, and it is now a diagnosable target
rather than a statistical suspicion.** Not started.

## 5. SO IS IT READY FOR PRODUCTION

For **hg38 short-read germline L1/Alu/SVA**, the evidence now supports yes: 1.67× faster, 47% of
the memory, better precision at every family, bit-reproducible output, and a sensitivity deficit
localised to three named SVA loci. Bit-reproducibility is a stronger production property than the
0.28-point pooled sensitivity difference is a liability, and master cannot match it.

What still argues against a general "yes", unchanged from the previous handoff except for item 2:

1. **Validated scope is one sample, one build, one modality.** chm13's L1 library carries 0 decoy
   records, so that hazard is untested, not refuted.
2. ~~The C++ engine is untested on the default path~~ — **FIXED, §1.**
3. **The Stage 3 no-regression margin is PROVISIONAL.** The gate has never been fed a change it
   *should* fail; the PASS is real, its power is unmeasured.
4. **22 harnesses pass neither decoy flag** and inherit defaults. This already bit `baseline0_hg002.sbatch`.
5. **Somatic/case-control is implemented but ungated** (Stage 4). Mosaic untouched.
6. **Present in master, absent here**: `x_ref_TE_deletion`, `x_population_AF`, `x_vcf_merger`,
   `x_genotype_classify`, and the `gnrt_pipeline_local*` entry points. Scope call, not a measurement.

## 6. Where things live

| | |
|---|---|
| **the replicate numbers** | `docs/handoff/runs/probe_repro_comparison.py --array-job 50062602` — offline, seconds |
| the harness | `docs/handoff/runs/repro_master_vs_xtea2.sbatch` (array of 6; 1–3 master, 4–6 xTea2) |
| the six trees | `/n/scratch/users/c/cos689/xtea2_testdata/repro_{master,xtea2}_rep{1,2,3}` |
| per-run provenance | each tree's `repro_run.tsv` — arm, rep, node, **`sacct_id`**, wall |
| the C++ freshness contract | `tests/conftest.py:xtea2_scan_bin` + `tests/test_scan_bin_freshness.py`; `tests/README.md` |
| the previous state handoff | `docs/handoff/HANDOFF_2026-08-07_objective_met.md` (agenda + parked list) |
| the standing decisions | `HANDOFF_2026-08-04_td_gate_built.md` §10 — four, unchanged |

## 7. A TRAP THIS SESSION HIT, so nobody re-derives it

**`$SLURM_JOB_ID` is the wrong sacct key for an array task.** Slurm gave task 6 a `SLURM_JOB_ID`
**equal to the array job id** (`50062602`), and `sacct -j 50062602` expands to *all six* tasks — so
that replicate matched no per-task row and silently lost its MaxRSS, in a way that looks like sacct
retention rather than a key collision. Both the harness and the aggregator now use
`<array_job>_<task>`, and the sidecar records it as `sacct_id`.

Also still true, and still worth the line: **never invoke
`docs/handoff/runs/master_serial_wgs.sbatch`** — it opens `rm -rf "$OUT"` and its default `OUT` is
a frozen NA12878 tree.
