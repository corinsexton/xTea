# xTea2 handoff — 2026-08-01 (Stage 0 complete, gate passes, Stage 1 unblocked)

**This is the current entry point for state. It is NOT the entry point for agenda — the plan is.**

> **Read `docs/handoff/DOWNSTREAM_OVERHAUL_PLAN.md` first.** It still outranks this chain on agenda.
> This handoff records what state that plan is in after Stage 0 was executed end to end; it does not
> restate the plan. Everything factual below is also written into the plan itself, at the section that
> owns it — the plan is the durable copy, this is the index.

**One-line status: Stage 0 is DONE and its gate PASSES. Nothing is in flight. The tree is clean at
`1db9cc2` on `xTea2`, 89 tests passing. Stage 1 is unblocked and its scope conditional is resolved —
the read gather is 73-87% of transduction, so the C++ per-site gather subcommand is the first move,
not collapsing the bwa passes.**

## 0. State at handoff

| | |
|---|---|
| branch / HEAD | `xTea2` @ `1db9cc2` |
| working tree | clean |
| tests | **89 passed, 18 skipped** (skips are the unbuilt C scan bin, pre-existing) |
| in flight | **nothing** |
| commits this session | 16, all Stage 0 |
| the agenda | `docs/handoff/DOWNSTREAM_OVERHAUL_PLAN.md` |

## 1. Stage 0 scorecard

| item | status |
|---|---|
| phase timers | `6ce36c7` — 13 phases + total, **unconditional**, `[TD-PHASE]` to stderr |
| phase breakdown | **gather 73-87%**, bwa 9-25%, SAM parse <3% (2 samples × 2 families) |
| type-set-dependence fix | `de96727` + `tests/test_mate_in_rep_type_partition.py` |
| `col5` → `mate_in_rep_clip_cnt` | `39dd582`, CLI `--mate-in-rep-families` (+ deprecated alias) |
| `_tprt` | **resolved** — canonical-mechanism filter; `_all` is the denominator |
| transduction scoring | **design settled**, feasibility verified (`b819533`) |
| chr1 slices | both cut and verified (`6cf5927`) |
| BASELINE-0 + gate | **frozen, gate PASSES** (`5d5272c`, `1db9cc2`) |

## 2. The result that decides Stage 1

Transduction's read gather dominates. Measured on two samples and two families, so it is not a
one-dataset artifact:

| | gather | bwa (3 passes) | SAM parse |
|---|---|---|---|
| L1 / NA12878 | 86.2% | 11.6% | 1.8% |
| L1 / HG002 | **87.3%** | 9.2% | 2.9% |
| SVA / NA12878 | 72.9% | **25.2%** | 1.5% |
| SVA / HG002 | 79.6% | 17.3% | 2.3% |

Reproducibility gate met (2.4% / 0.6% spread across repeats; outputs byte-identical, so the step is
deterministic). The 13 phases cover 99.6% of the function and the function is 95-96% of the whole
`--transduction` step.

**Two of the plan's own claims are corrected by this:** SAM parsing is **not** a cost (~18s per family
to parse all ~4.4 GB — "stream the SAMs" is refuted as a Stage 1 target), and bwa is secondary but
**consistently ~2× heavier for SVA than L1** on both samples. Amdahl: removing the gather entirely
caps the step at 7.2× (L1) / 3.7× (SVA).

## 3. BASELINE-0, and what the gate actually showed

Five HG002 arms, jobs `49441918-22`. Frozen VCFs + `PROVENANCE.txt` (commit, job, BAM md5, reference,
seed, sha256 per VCF) under each run's `BASELINE-0/`.

| run | LINE1 | ALU | SVA | merged |
|---|---|---|---|---|
| post-fix `L1,Alu,SVA` | 188 | 1223 | 82 | 1493 |
| post-fix `L1` | 200 | — | — | 200 |
| post-fix `SVA` | — | — | 85 | 85 |

`bash docs/handoff/runs/stage0_gate_check.sh` reproduces the verdict: **A PASS** (three-type
byte-identical across the rename+fix, 7/7 artifacts), **B PASS output-neutral**, **C PASS** (SVA
per-type list byte-identical between `--types SVA` and `--types L1,Alu,SVA`).

**The one result that contradicts a plan prediction.** §6.1 predicted "SVA-only output changes". It
did not, on HG002. The mechanism moved enormously — `mate_in_rep` changed at **1,434,938 of 8,622,501**
positions, a 21× increase in total signal, including **38 of the 85 called SVA sites** — and the SVA
call set is byte-identical. The fix is right (**`post >= pre` at every one of the 8.6M positions, 0
violations** — union ≥ subset, measured rather than asserted, which also pins the direction: it can
only ever *add* SVA calls). It simply flipped no decision here, because the SVA cns consistency check
was not the binding constraint at any called site. **Say neither "changed nothing" nor "recovered N
calls".**

## 4. NEXT — Stage 1, and the one thing to do before it

Per the plan, and now unblocked:

1. **The new `xtea_scan` subcommand for the unified per-site gather.** Reuse the two patterns already
   proven in `src/csrc/cns_emit.hpp`: `FetchRegion` coalescing with the exact ownership rule
   (`:502-523`) and candidate-major emission (`:798-843`). **Ordering determinism is load-bearing** —
   `cns_emit.hpp:798-817` documents that bwa's tie-break among equal-scoring placements depends on read
   index, so reordering the FASTQ silently changes which repeat copy multi-mapping clips land on.
2. **Eliminate the per-site temp-file round-trip and the per-chromosome `Pool` churn**
   (`x_clip_disc_filter.py:351-358`).
3. Stage 1's gate is a chr1 diff vs BASELINE-0 with **byte-identical calls required**, plus a two-run
   determinism check. The chr1 slice exists; use it, not the WGS run.

**Before Stage 3 (not before Stage 1):** `benchmark/benchmark_truth_set.py`. The design is settled in
§7.1 — derive truth transduction status from the ALT soft-clip residual (feasibility verified: 51-60
of 197 truth L1 records carry ≥200 bp of non-consensus sequence), five-way error taxonomy, two
denominators never conflated.

## 5. What NOT to redo

The plan's §10 is the authority and has been extended this session. Highlights of what is new:

- **Do not chase `Error happen at merge clip and disc feature step: <contig> not exist`.** It is not an
  error — `x_intermediate_sites.py:514` fires when a contig has clip candidates but no *discordant*
  ones, normal on alt/unplaced contigs. NA12878 emits 4,395; HG002 emits 87,039 across 1,625 contigs
  because its BAM is alt-rich and the sample is male. No primary chromosome affected.
- **Do not re-measure the transduction breakdown.** Four measurements, two samples, two families,
  reproducible within 2.4%.
- **Do not compare across the family merge by exact position.** Use `locus_table`'s ±50 bp window.
  Exact matching reported 8 L1 calls as "dropped" that were all reassigned to Alu.
- **`f=$(ls MISSING 2>/dev/null | head -1)` is a landmine under `set -euo pipefail`** and
  `docs/handoff/runs/port_wgs.sbatch` still contains it, latent.
- **`run_xtea.py` deletes everything but the merged VCF on success.** Any per-type comparison needs
  `--keep-intermediates`.
- **The `[TD-PHASE]` timers land in `<sample>_xtea2.log`, not the driver's stdout.**
- **Per-type VCFs are named from the BAM basename, not the sample.**

## 6. Where things live

| | |
|---|---|
| **the agenda** | `docs/handoff/DOWNSTREAM_OVERHAUL_PLAN.md` |
| BASELINE-0 | `/n/scratch/.../xtea2_testdata/baseline0_hg002_{L1-Alu-SVA,L1,SVA}/BASELINE-0/` |
| pre-fix comparators | `/n/scratch/.../xtea2_testdata/prefix_hg002_{L1-Alu-SVA,SVA}/` |
| the gate | `docs/handoff/runs/stage0_gate_check.sh` |
| HG002 chr1 slice | `test_data/HG002/chr1/HG002.chr1.bam` (3.4 GB) |
| SMHT005-3G chr1 slice | `test_data/SMHT005-3G/chr1/` (5 CRAMs, 43 GB, `bam_list.txt`; needs `--ignore-rg`) |
| germline inputs + provenance | `test_data/HG002/`, esp. `PROVENANCE.md` |
| run harnesses | `docs/handoff/runs/{baseline0_hg002,chr1_slices,transduction_phase_timing}.sbatch` |
| worktrees added this session | `baseline0_wt` (`9a43b80`), `prefix_wt` (`6ce36c7`) |
| mosaic design, already done | `docs/handoff/MOSAIC_PLAN.md` |
