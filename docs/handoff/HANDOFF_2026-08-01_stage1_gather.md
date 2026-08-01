# xTea2 handoff — 2026-08-01 (Stage 1 lead item built: the native transduction gather)

**This is the current entry point for state. It is NOT the entry point for agenda — the plan is.**

> **Read `docs/handoff/DOWNSTREAM_OVERHAUL_PLAN.md` first.** It still outranks this chain on agenda.
> This handoff records what state that plan is in after Stage 1's lead item was built; it does not
> restate the plan. Everything factual below is also written into the plan itself, at the section
> that owns it — the plan is the durable copy, this is the index. It supersedes
> `HANDOFF_2026-08-01_stage0_complete.md`, whose §4 "NEXT" item 1 is what this delivers.

**One-line status: `src/csrc/td_gather.hpp` replaces both of the transduction step's Python read
gathers with one native pass and is BYTE-identical to them — clip FASTQ, `all_disc_pos` and disc
mate FASTA — on real chr1 data, on a hermetic multi-contig fixture, and end to end through the
HG002 chr1 pipeline. It is off by default (`--td-gather-scan`). PLACEHOLDER_GATE_ONELINE**

## 0. State at handoff

| | |
|---|---|
| branch / HEAD | `xTea2` @ PLACEHOLDER_HEAD |
| working tree | clean |
| tests | PLACEHOLDER_TESTS |
| in flight | **nothing** |
| the agenda | `docs/handoff/DOWNSTREAM_OVERHAUL_PLAN.md` §4 Stage 1 |

## 1. What was built

`xtea_scan --td-sites SITES [--td-extnd N] --td-clip-fastq F --td-disc-fa F [--td-disc-pos F]`,
implemented in `src/csrc/td_gather.hpp`, called from `XClipDiscFilter.native_td_gather`, selected by
`run_xtea.py --td-gather-scan` → the inherited env var `XTEA_SCAN_TD_GATHER` (the transduction step
is its own subprocess, so a flag cannot reach it). It replaces, in one pass:

- `collect_clipped_reads` → `collect_clipped_reads_of_given_list`: a fresh `Pool` **per
  chromosome**, one forked task and **two temp files per site** (~54K sites on WGS L1), which the
  parent re-reads and unlinks;
- `collect_disc_reads` → `extract_mate_reads_by_name`: a `Pool` over **50 Mb genome bins** where
  every bin re-parses the entire `all_disc_pos` file and fetches its whole range.

Both halves in one process is safe because the disc half never depended on anything computed
between them: `collect_disc_reads` ignores its own site-list argument and reads the file the clip
half wrote (`x_clip_disc_filter.py:392`). Reuses the two patterns the plan names —
`cns_emit.hpp`'s `FetchRegion` merge with its ownership rule (`:502-523`) and candidate-major
emission (`:798-843`).

## 2. The part that is easy to get wrong: three DIFFERENT orderings

bwa's tie-break among equal-scoring placements depends on a read's index in its input
(`cns_emit.hpp:798-817`, measured), so the gather's outputs are order-sensitive and the gate is
`cmp`, not "same reads". Three separate rules carry it, each with its own test:

1. **Clip FASTQ: candidate-major in the SITE FILE's order** — contigs by first appearance,
   positions in file order within a contig. `sf_raw_slct` is **not** position-sorted (cns
   candidates first, then raw-disc lines), so `emit_cns_multi`'s position-sorted fan would have
   been wrong here. This is why it is a separate collector rather than a flag on `emit_cns_multi`.
2. **`n_cnt` (field 5) keeps master's GAPS** — master builds the right-clip read name and
   increments the counter *before* that side's gates, so a failing right clip burns an index.
   Nothing reads field 5; reproducing it is what makes `cmp` a valid oracle at all.
3. **Disc FASTA: BAM-HEADER contig order** (decoy names skipped), coordinate order within a
   contig, per mate read its insertion sites in disc-pos order.

Rules 1 and 3 are opposite orders for the two files, which no single-contig test can see — hence
the hermetic fixture in `tests/test_td_gather.py` whose header lists `chr2` first while its site
file names `chr1` first, with every mate on its anchor's own contig.

## 3. Numbers

PLACEHOLDER_NUMBERS

## 4. NEXT

1. **Stage 1's remaining scope item is the manifest** — one collection serving all three types.
   Deferred on purpose: today the three types' transduction steps are separate concurrent
   subprocesses with separate site lists, so there is nothing to fan across. The collector already
   takes an arbitrary site list and nothing in it is type-specific, so this is a **caller** change
   and it belongs with Stage 2's unification, where the driver stops forking per type. §1 of the
   plan is the reason it is worth doing: the cross-family sharing win only exists if the shared
   work runs once at full width.
2. **Then Stage 2** (unify the family dimension) per the plan. Its gate is different in kind —
   calls may move, but only as family-assignment flips, each hand-explained.
3. **Before Stage 3, not before Stage 2:** `benchmark/benchmark_truth_set.py`. Design settled in
   plan §7.1; unchanged by this session.

## 5. What NOT to redo

The plan's §10 is the authority and has been extended again. New this session:

- **The legacy disc gather's cost barely depends on the site count** — 106s for 207 chr1 sites and
  106s for 3733. Do not benchmark it on a small list and conclude anything.
- **The clip gather is index-fetch bound, not fork bound**: 1.26× at equal threads, against 9.6×
  for the disc half. Do not attribute the clip half's speedup to C++; it comes from threads.
- **`[TD-PHASE]` lands in `<sample>_<TYPE>.steps.log` under `--parallel-downstream`**, not
  `<sample>_xtea2.log` as the previous handoff says (true only for the serial path). A harness that
  greps the wrong file prints an empty timing table.
- **An exact duplicate `(chrom,pos)` site line crashes the legacy gather** (`FileNotFoundError`,
  `x_clip_disc_filter.py:373`). Observed. Hand-built site lists are where it bites.
- **The legacy transduction gather is broken for multi-BAM input** — one shared `all_disc_pos`,
  truncated per BAM, so only the last BAM's positions survive. No baseline exists; `td_gather.hpp`
  fixes it and says so.
- **`local arm=$1 armout="$OUT/$arm"` dies under `set -u`** (the whole line is expanded before
  `local` assigns). Killed the first gate submission, job `49452897`.
- `docs/handoff/runs/port_wgs.sbatch`'s `ls MISSING | head -1` landmine is **fixed**, not just
  documented. Do not re-flag it.

## 6. Where things live

| | |
|---|---|
| **the agenda** | `docs/handoff/DOWNSTREAM_OVERHAUL_PLAN.md` (§4 Stage 1 has the build note) |
| the gather | `src/csrc/td_gather.hpp`; CLI in `src/csrc/xtea_scan.cpp` (`--td-*`) |
| the switch | `run_xtea.py --td-gather-scan` → `XTEA_SCAN_TD_GATHER` → `x_transduction.py` |
| tests | `tests/test_td_gather.py` (4 fast + 1 `slow`; skip without `XTEA2_SCAN_BIN`) |
| the gate | `docs/handoff/runs/stage1_gate_chr1.sbatch` (3 arms + `test_downstream_from_cns.py`) |
| gate output | PLACEHOLDER_GATEDIR |
| BASELINE-0 (Stage 2's comparator) | `/n/scratch/.../xtea2_testdata/baseline0_hg002_{L1-Alu-SVA,L1,SVA}/BASELINE-0/` |
| HG002 chr1 slice | `test_data/HG002/chr1/HG002.chr1.bam` |
| worktree used for the gate | `/n/scratch/.../xtea2_testdata/stage1_wt` |
