# xTea2 handoff — 2026-08-01 (Stage 1 lead item built: the native transduction gather)

**This is the current entry point for state. It is NOT the entry point for agenda — the plan is.**

> **Read `docs/handoff/DOWNSTREAM_OVERHAUL_PLAN.md` first.** It still outranks this chain on agenda.
> This handoff records what state that plan is in after Stage 1's lead item was built; it does not
> restate the plan. Everything factual below is also written into the plan itself, at the section
> that owns it — the plan is the durable copy, this is the index. It supersedes
> `HANDOFF_2026-08-01_stage0_complete.md`, whose §4 "NEXT" item 1 is what this delivers.

**One-line status: `src/csrc/td_gather.hpp` replaces both of the transduction step's Python read
gathers with one native pass and is BYTE-identical to them — clip FASTQ, `all_disc_pos` and disc
mate FASTA. The Stage 1 chr1 gate PASSES on all three legs (agreement 9/9 artifacts, determinism
9/9, `test_downstream_from_cns.py` against the NA12878 golden with the native gather forced on), and
so does a WGS-scale artifact oracle at 49K/54K sites over 24 contigs. The gather is **15.25× (L1) /
17.82× (SVA)** at WGS scale, ~1.8× on the whole downstream phase. Off by default:
`--td-gather-scan`.**

## 0. State at handoff

| | |
|---|---|
| branch / HEAD | `xTea2` @ `969a2d8` (+ the docs commit that adds this file) |
| working tree | clean |
| tests | **108 passed, 1 skipped** with `XTEA2_SCAN_BIN` set (the skip is `test_e2e_concordance`, needs `XTEA2_OUTPUT_DIR`); **86 passed, 23 skipped** without it. `-m slow` adds 1 (the disc-FASTA oracle, ~2 min) |
| in flight | **nothing** (both gates finished: chr1 job `49453032` PASS, WGS oracle job `49453857` PASS) |
| the agenda | `docs/handoff/DOWNSTREAM_OVERHAUL_PLAN.md` §4 Stage 1 |

## The gate

Job `49453032`, HG002 chr1, `L1,Alu,SVA`, three arms at `49154dc` from the frozen `stage1_wt`
worktree — **A** legacy gather, **B** native, **C** native again:

| leg | result |
|---|---|
| AGREEMENT A vs B | **PASS** — 9/9 byte-identical (3 per-type call lists, 2 sibling lists, 3 per-type VCFs, merged VCF) |
| DETERMINISM B vs C | **PASS** — the same 9/9 |
| `test_downstream_from_cns.py` (native forced on) | **PASS** — reproduces the NA12878 chr1 *master* golden's downstream |
| **WGS artifact oracle** (job `49453857`) | **PASS** — all three gather artifacts byte-identical in both families at 49K/54K sites over 24 contigs, and deterministic across two native runs. §3 has the numbers |

All three arms: L1 1 / Alu 14 / SVA 2, merged 17. Whole-run 419s / 344s / 356s.

Two things about that gate are worth knowing before trusting it. It ran at `49154dc`, so it does
**not** cover the two commits after it — but those are a stderr timer line and comments only, and
the chr1 artifacts were `cmp`-verified unchanged across the timer commit. And the "arms actually
ran different paths" assertion was added to the harness *after* this submission, so for this job it
was checked by hand (arm A's `[TD-PHASE]` lines carry `clip_gather`/`disc_gather`, arm B's carry
`native_gather` plus the `[td-gather]` banner). Future runs assert it automatically — without that
check, three legacy arms agreeing with each other would PASS and mean nothing, because the native
path is selected by an inherited env var.

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

**In the pipeline** (HG002 chr1, gate job `49453032`, arm A legacy vs arm B native, `-n 6` per type
under `--parallel-downstream`):

| | legacy gather | native gather | `call_candidate_transduction_v3` | the whole step |
|---|---|---|---|---|
| L1 (2815 sites, extnd 956) | 75.1s = clip 38.5 + disc 36.7 | **21.8s (3.44×)** | 83.6s → 30.5s (2.74×) | 100.7s → **46.9s (2.15×)** |
| SVA | 82.1s = clip 45.6 + disc 36.5 | **22.4s (3.67×)** | 101.6s → 41.8s (2.43×) | 115.4s → **54.8s (2.11×)** |

Three levels because they answer different questions and quoting the wrong one oversells this: the
gather itself is 3.4-3.7×, the instrumented function 2.4-2.7×, the *step* 2.1× (the step also pays
automatic-parameter generation, the annotation tree, and the re-selection, none of which moved).

Whole chr1 run 6:59 → 5:44. That end-to-end figure is small for a reason worth remembering: on
chr1 the **detection** `scan_pass` is 203.9s of the run — half of it — and the downstream plan does
not touch it.

**Where the native gather's own time goes** (chr1 L1 site list, 2 threads, from the three
`[td-gather]` numbers): collect 19.0s + fan 0.1s + mate-gather 11.1s. Against the legacy pair at
the same thread count (clip 24.0s, disc 106.5s) that is **1.26× on the clip half and 9.6× on the
disc half.** The clip half is index-fetch bound, not fork bound — its win comes from threads, not
from leaving Python. The disc half's old cost was nearly independent of site count (~106s for 207
sites *and* for 3733, because it scanned the whole genome either way).

**At WGS scale — job `49453857`, both arms at 12 threads, the run's own `.slct` lists, peak RSS
0.94 GB:**

| | sites | legacy | native | speedup |
|---|---|---|---|---|
| L1 | 49,042 over 24 contigs | **812.0s** (clip 615.0 + disc 197.0) | **53.3s** (collect 22.0 + fan 1.8 + mate 29.1) | **15.25×** |
| SVA | 53,569 over 24 contigs | **1011.3s** (clip 817.0 + disc 194.3) | **56.8s** (collect 23.2 + fan 2.1 + mate 31.0) | **17.82×** |

All three artifacts byte-identical in both families (clip FASTQ 178/181 MB, `all_disc_pos` 90/94 MB,
disc FASTA 235/245 MB; 1.46M/1.49M clip records, 1.04M/1.08M anchors, 592K/629K mates), and two
native runs per family byte-identical to each other.

**That 15-18× corrects a projection I wrote earlier in this file from the chr1 ratio (1.5-1.7×) and
which was wrong by an order of magnitude.** The reason is worth carrying forward: the legacy clip
gather's cost is per-site overhead — one fork and two temp files each — so it scales with the site
count, and chr1 has 13× fewer sites (6.4 ms/site at 3733, 12.5 ms/site at 49,042). **Do not size
this port, or the next one, from a chr1 measurement.**

Recomputed step effect, mixing thread counts so treat as ±: L1's step 1480.7s → ~676s (2.2×),
SVA's 1772.3s → ~931s (1.9×), `downstream_parallel` 1884s → ~1040s (1.8×), with the SVA chain
becoming the critical path. Of the plan's 7.2× Amdahl ceiling on the *function*, 5.2× is realised.

**Read the evidence in this order, because the chr1 call diff is the WEAKEST link in it.** HG002
chr1 yields **17 calls** (L1 1 / Alu 14 / SVA 2) — that is normal for a chr1 slice, the NA12878
chr1 master golden is 32 (L1 0 / Alu 27 / SVA 5), and the plan already notes chr1 cripples L1. So
17 identical calls is a wiring proof, not a coverage proof. The load-bearing evidence is the
byte-identity of the gather's **output files** at 49K/54K sites over 24 contigs above. Identical
bytes into bwa is what makes identical calls a consequence rather than a coincidence.

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

**One open decision, deliberately not taken here: whether `--td-gather-scan` becomes the default.**
The plan's Stage 1 rollback clause is "legacy remains the default; the flag selects the new gather",
so it stays opt-in until someone decides otherwise. The WGS oracle criterion is now **met** (§3);
what remains before flipping is one WGS run whose calls match a BASELINE-0 arm — the chr1 gate does
not substitute, being 17 calls. Note what flipping the default would also do: it makes the multi-BAM
divergence in §5 reachable, so the mosaic path would then be running the fixed behaviour with no
baseline of its own. That is a Stage 4 concern, and it is the reason not to flip casually.

## 5. What NOT to redo

The plan's §10 is the authority and has been extended again. New this session:

- **The legacy disc gather's cost barely depends on the site count** — 106s for 207 chr1 sites and
  106s for 3733. Do not benchmark it on a small list and conclude anything.
- **Do not size a gather port from a chr1 measurement.** The legacy clip gather's cost is per-site
  overhead, so it scales with the site count: the clip half is 1.26× on chr1's 3733 sites but the
  whole gather is 15-18× on the real 49K/54K-site WGS lists. The projection this file originally
  carried (1.5-1.7×, from the chr1 ratio) was wrong by an order of magnitude.
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
| tests | `tests/test_td_gather.py` (5 fast + 1 `slow`; all skip without `XTEA2_SCAN_BIN`) |
| the gate | `docs/handoff/runs/stage1_gate_chr1.sbatch` (3 arms + `test_downstream_from_cns.py`) |
| gate output | `/n/scratch/.../xtea2_testdata/stage1_gate_chr1_49453032/` (+ `stage1_gate_chr1_49453032.out`) |
| WGS oracle | `docs/handoff/runs/{stage1_gather_oracle_wgs.sbatch,td_gather_oracle.py}`; output `stage1_oracle_wgs_49453857*` |
| BASELINE-0 (Stage 2's comparator) | `/n/scratch/.../xtea2_testdata/baseline0_hg002_{L1-Alu-SVA,L1,SVA}/BASELINE-0/` |
| HG002 chr1 slice | `test_data/HG002/chr1/HG002.chr1.bam` |
| worktree used for the gate | `/n/scratch/.../xtea2_testdata/stage1_wt` |
