# xTea2 handoff — 2026-08-02 (Stage 2c's decision point moved; the transduction question is closed)

**This is the current entry point for state. It is NOT the entry point for agenda — the plan is.**

> **Read `docs/handoff/DOWNSTREAM_OVERHAUL_PLAN.md` first.** It still outranks this chain on agenda.
> Everything factual below is written into the plan at the section that owns it; this is the index.
> Supersedes `HANDOFF_2026-08-02_stage2a_and_benchmark.md`, whose §4 NEXT items are what this session
> addressed: **item 1 is settled and closed** (§1 below), **item 2's sub-items 1-3 are built and
> gated-but-not-run** (§2), sub-item 4 is deliberately still open.

**One-line status: family assignment can now happen once, right after `cns_filter`, instead of three
chains later at the postF_hc merge — `--family-assign-early`, default OFF, with a WGS gate written
and NOT YET RUN. Building it corrected the premise Stage 2c was planned on: `scan_clip_realign.tsv`
carries no per-family signal at all, family contention is 0.6–5% at every level so 2c is not a
speedup, and the predicted movement is 4 of 1493 loci — all four toward truth. Separately, HG002's
inability to gate transductions is now settled by a third method and the search is closed.**

## 0. State at handoff

| | |
|---|---|
| branch / HEAD | `xTea2` @ `34236ea` + this session's commits |
| working tree | see `git log`; all work committed |
| tests | **152 passed, 1 skipped** with `XTEA2_SCAN_BIN` set, excluding `test_downstream_from_cns.py` (the slow chr1 golden). Was 139 — +6 in `test_locus_table.py`, +7 in the new `test_family_assign_early.py`. The skip is `test_e2e_concordance` (needs `XTEA2_OUTPUT_DIR`) |
| in flight | **job `49498881`** — the Stage 2c.3 gate, `docs/handoff/runs/stage2c_family_assign_wgs.sbatch`, submitted 2026-08-02 10:00, two arms at n=24, ~90 min. Worktree `/n/scratch/.../stage2c_wt` @ `4e8f6b1`. Output `stage2c_family_assign_49498881.out` |
| the agenda | run that gate, then Stage 3 |

## 1. NEXT item 1 is CLOSED: HG002 cannot gate transduction, stop looking

The previous handoff said: go to Chu *et al.* 2021's supplement and check whether it enumerates
transductions. **The truth files answer it without the supplement.** Their complete INFO key
vocabulary, across all five VCFs, is a GIAB SV benchmark's — `ClusterIDs`, `ExactMatchIDs`,
`NumTechs`, `MendelianError`, `HG003_GT`/`HG004_GT`, `TRall`, `segdup`, `sizecat`, `REPTYPE`, … —
plus `QNAME`/`QSTART`/`QSTRAND` (NA24385 haplotype-contig provenance), plus `xTEA_specific` on SVA
only. **No transduction field, and no schema slot one could have been dropped from.**

That is a third independent line of evidence, after the ALT-vs-consensus derivation (`49482568`) and
the ALT-vs-genome cross-check (`49482779`). The decision, written into plan §7.1:

1. **Stage 3's gate is position-level per family plus the named strata** (orphan / sibling /
   no-`TD_SRC`), which `benchmark/benchmark_truth_set.py score` already reports. Transduction calls
   are reported and stratified, **not scored right or wrong**.
2. **Transduction correctness needs a different sample** — Stage 4 work, and *not* a Stage 3 blocker.
3. **Do not re-derive this from the ALT sequence.** Three methods, one answer.

One thing to carry: **`SVA_tprt_all.vcf` has 5 `xTEA_specific` records**, so 5.6% of the 90-record
SVA denominator came from xTea itself. Mildly circular; harmless at this benchmark's precision, but
an SVA sensitivity quoted to the tenth of a percent is quoting through it.

## 2. Stage 2c: the decision point moved, and the plan's premise for it was wrong

### 2.1 The correction — read this before touching Stage 2c

Plan §4's 2c decomposition says the locus evidence table is "a read of `scan_clip_realign.tsv` plus
the disc side, not a new computation," because that file has three rows per position, one per
family. The three rows exist. **The numbers in them are the same number.**

| | |
|---|---|
| positions with clip support in all three families | 8,621,063 of 8,622,517 |
| cross-family position-pairs whose `nL/nR/mate_in_rep` differ | **1,147 of 17,244,556** |
| positions with ≥5 clip support | 155,615, of which **504** differ |
| resulting clip candidate lists | 577,293 / 577,312 / 577,347 — the same positions |

**Family-discriminating evidence starts on the disc side**, where the mate is tested against each
family's own annotation: `candidate_list_from_disc.txt` is 2,244 L1 / 9,468 Alu / 1,183 SVA, **union
12,223, only 672 rows (5.2%) claimed by more than one family.** That is also `cns_filter`'s actual
input. Any evidence table for Stage 3 should start there.

**And Stage 2c is not a speedup.** Contention is 672 of 12,895 disc candidates, 16 of 2,796 cns
rows, 22 of 1,515 postF_hc calls — 0.6–5% — of a downstream that is itself only ~21% of the run
after Stage 1 (`downstream_parallel` 457s of 2,166s; `scan_pass` ~1,500s). Do 2c because family
stops being decided by whichever list a call happened to survive on. Anyone expecting a speedup will
correctly conclude it isn't worth doing.

### 2.2 How far up front the decision can go: exactly as far as `cns_filter`

`two_side_tprt_both` — which decides 13 of the 22 conflicts — is written **by `cns_filter`**
(`x_clip_disc_filter.py:1845`). Nothing earlier has it. So `candidate_disc_filtered_cns.txt` is the
earliest file the existing rule can read, and deciding there is a **relocation**; deciding at the
candidate list would be a different rule adopted without evidence for it.

### 2.3 What was built

- **2c.2, the lift.** `locus_table.load_rows` / `LocusTable.assign_rows` — the rule now takes
  pre-parsed claims, so one copy serves both ends of the pipeline. `assign_type` is *defined* as
  `load_rows + assign_rows` and a test asserts the two agree; "moved, not rewritten" is the whole
  safety argument.
- **2c.3, the relocation.** `locus_table.assign_and_filter` resolves family over the three
  `candidate_disc_filtered_cns.txt` and rewrites each to its own winners (originals kept at
  `.pre_family_assign`, which is what a gate diffs to attribute a moved call).
  `run_xtea.downstream_waves` splits the per-family chain into `cns_filter` → **barrier** →
  `transduction…postF_hc`, reusing `build_type_steps`' existing `start_from`/`stop_after` so the argv
  cannot drift from the single-wave chain. `--family-assign-early`, **default OFF**, no-op for a
  single `--types` entry.
- **The sidecars.** `<sample>/cns_family_conflicts.tsv` (barrier) alongside the existing
  `merged/locus_family_conflicts.tsv` (merge), same format, so the two decision points are directly
  diffable.

### 2.4 The prediction, made before the gate exists

`docs/handoff/runs/probe_early_assign.py` runs offline against a completed run's kept intermediates
— no BAM, no compute. On BASELINE-0: **4 of 1,493 loci move, and all four move toward truth.**

| locus | today | at cns | against truth |
|---|---|---|---|
| `chr14:20707975` | L1 | **SVA** | truth has an **SVA at 20707961** (14 bp) — +1 SVA TP, −1 L1 FP |
| `chr12:95392885` | L1 | dropped | off-truth in every family — −1 L1 FP |
| `chr9:4710380` | L1 | dropped | off-truth in every family — −1 L1 FP |
| `chr5:101161558` | L1 | Alu | off-truth either way — an FP changes family |

**The surprise, which is the opposite of the obvious prediction:** the 12 `Alu over L1
tprt_override` conflicts — the ones that motivate the override — **do not move**. At cns time the
losing L1 record does not exist yet; transduction and sibling create it (35 of the 1,493 called loci
have no cns-level counterpart at all). So the barrier reaches those same 12 answers with no override
needed.

**One approximation, and it bounds what the gate may assume:** the probe holds each family's
transduction/sibling output fixed. Those steps create sites (L1 cns 300 → cns2 315), so a real arm
may move loci the probe does not name. The gate therefore reports every moved locus rather than
asserting this list.

## 3. The gate — SUBMITTED as job `49498881`, result pending

```
sbatch docs/handoff/runs/stage2c_family_assign_wgs.sbatch     # already done: 49498881
```

Worktree `/n/scratch/users/c/cos689/xtea2_testdata/stage2c_wt` @ `4e8f6b1`, `xtea_scan` built there,
clean. Two arms at n=24, ~90 min. Output
`/n/scratch/users/c/cos689/xtea2_testdata/stage2c_family_assign_49498881.out`.

**This is the first gate in the overhaul a byte-diff cannot settle.** Stages 0 and 1 were mechanics —
any moved call was a bug. This stage is *expected* to move calls, so the gate names them in advance.

| leg | what it decides |
|---|---|
| 0 | the arms ran the paths they claim (OFF must not log the barrier; ON must, and must write a non-empty conflicts sidecar) |
| 1 | arm OFF vs frozen BASELINE-0 VCFs — **must be identical**, or ON's delta is uninterpretable |
| 2 | movement ON vs OFF, locus by locus, **bounded at 22** (the whole contended population) with every moved locus printed |
| 3 | **the decisive leg** — `benchmark_truth_set.py score` on both arms. ON must not lose sensitivity or precision in any family; predicted delta is +1 SVA TP and −3 L1 FP |

Two arms rather than one, for the Stage 1 reason: BASELINE-0 was frozen at `9a43b80` and HEAD is
later, so a lone ON arm conflates "the barrier did this" with "something else moved since".

## 4. NEXT

1. **Read job `49498881`'s result.** If leg 3 holds, flip `--family-assign-early` to default and
   record the flip the way Stage 1's was recorded. If leg 1 fails (arm OFF ≠ BASELINE-0), stop and
   find the drift first — arm ON's delta means nothing until then.
2. **Then the rest of 2c**, still open and unchanged in scope: the single per-locus loop that
   replaces the three subprocess chains, the shared gather manifest (plan §4 2b), and
   `--emit-per-type` — deliberately last, since flipping per-type VCFs to opt-in breaks every gate
   script that `cmp`s them.
3. **Then Stage 3**, whose gate is now fully specified (§1): position-level per family plus the
   orphan/sibling/no-`TD_SRC` strata, no transduction pass/fail.

## 5. What NOT to redo

Plan §10 is the authority and gained seven entries. The ones most likely to cost a day:

- **`scan_clip_realign.tsv` is not per-family evidence** — 504 of 155,615 well-supported positions
  differ across families. Do not build a locus evidence table on it.
- **`two_side_tprt_both` is written by `cns_filter`**, so that is as early as the family rule can go.
- **Stage 2c is not a speedup** — 0.6–5% contention, on ~21% of the run.
- **The 12 `Alu over L1 tprt_override` conflicts are not contended at cns time**, because the losing
  L1 record is created later by transduction/sibling.
- **The HG002 truth VCFs' INFO vocabulary is GIAB's and has no transduction field** — third
  confirmation. Do not go to the paper's supplement for it.
- **`SVA_tprt_all.vcf` has 5 `xTEA_specific` records** (5.6% of the SVA denominator).
- Everything already listed in the previous handoff's §5 still stands, in particular that `--sva`
  reaches `cns_filter` only, and that L1 precision must always be quoted with its stratum
  (74.5% overall, 92.7% excluding orphan sites).

## 6. Where things live

| | |
|---|---|
| **the agenda** | `docs/handoff/DOWNSTREAM_OVERHAUL_PLAN.md` (§4 Stage 2c, §7.1 and §10 all updated) |
| the rule, one copy | `xtea/locus_table.py` — `load_rows`, `LocusTable.assign_rows`, `assign_and_filter` |
| the barrier + waves | `xtea/run_xtea.py` — `downstream_waves`, `run_family_assign_barrier` |
| tests | `tests/test_locus_table.py` (+6), `tests/test_family_assign_early.py` (new) |
| the prediction instrument | `docs/handoff/runs/probe_early_assign.py` — offline, no BAM |
| **the unrun gate** | `docs/handoff/runs/stage2c_family_assign_wgs.sbatch` |
| the benchmark | `benchmark/benchmark_truth_set.py`; tests `tests/test_benchmark_truth_set.py` |
| BASELINE-0 | `/n/scratch/users/c/cos689/xtea2_testdata/baseline0_hg002_L1-Alu-SVA/` (intermediates kept) |
| truth sets | `/n/data1/hms/dbmi/park/corinne/claude_sandbox/test_data/HG002/truth/` |
