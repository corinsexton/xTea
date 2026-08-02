# xTea2 handoff — 2026-08-02 (Stage 2c gated at zero movement; the transduction question is closed)

**This is the current entry point for state. It is NOT the entry point for agenda — the plan is.**

> **Read `docs/handoff/DOWNSTREAM_OVERHAUL_PLAN.md` first.** It still outranks this chain on agenda.
> Everything factual below is written into the plan at the section that owns it; this is the index.
> Supersedes `HANDOFF_2026-08-02_stage2a_and_benchmark.md`, whose §4 NEXT items are what this session
> addressed: **item 1 is settled and closed** (§1 below), **item 2's sub-items 1-3 are built and
> GATED — job `49502115` PASSED** (§2, §3), sub-item 4 is deliberately still open.

**One-line status: family is now decided right after `cns_filter` as well as at the postF_hc merge
(`--family-assign-early`, default OFF), the WGS gate PASSED, and it moves ZERO calls — because the
barrier does not replace the merge, and everywhere the two decision points overlap they agree, 7 for
7. Building it corrected three things the plan assumed: `scan_clip_realign.tsv` carries no
per-family signal at all, 2c is not a speedup (0.6–5% contention), and 12 of the 22 family conflicts
are structurally unreachable before transduction — so 2c cannot collapse to one decision point until
site creation is unified, which is a Stage 3 dependency. Separately, HG002's inability to gate
transductions is settled by a third method and the search is closed.**

## 0. State at handoff

| | |
|---|---|
| branch / HEAD | `xTea2` @ `34236ea` + this session's commits |
| working tree | see `git log`; all work committed |
| tests | **155 passed, 1 skipped** with `XTEA2_SCAN_BIN` set, excluding `test_downstream_from_cns.py` (the slow chr1 golden). Was 139 — +9 in `test_locus_table.py`, +7 in the new `test_family_assign_early.py`. The skip is `test_e2e_concordance` (needs `XTEA2_OUTPUT_DIR`) |
| in flight | **nothing** — the gate is done. Job `49502115` PASSED (`49498881` before it died on the companion-file bug, §2.5). Output `/n/scratch/.../stage2c_family_assign_49502115{,.out}` |
| the agenda | the rest of 2c (still three chains, still two decision points), then Stage 3 |

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

### 2.4 The result: ZERO calls move, and the reason is the finding

Job **`49502115`** (job `49498881` before it died — §2.5).

| | OFF | ON |
|---|---|---|
| merged VCF | 1493 | **1493** |
| per family (L1 / Alu / SVA) | 188 / 1223 / 82 | **188 / 1223 / 82** |
| call set `(chrom, pos, family)` | — | **identical to OFF** |
| truth score, all families | TP 1320 / FN 322 / FP 173 — 80.4% / 88.4% | **identical** |
| vs frozen BASELINE-0 VCFs | **4/4 identical** | — |
| conflicts resolved at the barrier | — | 16 |
| conflicts resolved at the merge | 22 | **15** |

**The barrier does not replace the merge — both run.** So it changes output only where it drops a
row that would have *won* at the merge, and it never does. Of its 16 decisions, **7 are also seen by
the merge and all 7 give the identical winner, loser and rule**; the other 9 the merge never sees.
22 − 7 = 15, exactly ON's merge count.

**The 15 the merge still resolves are unreachable from cns.** The losing family has no cns row
there, because transduction and sibling create its record *after* the barrier runs — and 12 of the
15 are the `Alu over L1 tprt_override` population. **So Stage 2c cannot end with one decision point
until site creation is unified too.** That is a Stage 3 dependency nobody had written down.

**Zero movement is still worth the run:** the relocation is now *proven* answer-preserving on real
WGS data rather than argued, which is what would license deleting the merge's resolution later.

### 2.5 Two of my own predictions were wrong. Both are corrected in place.

- **The predicted 4 moved loci did not happen**, and `probe_early_assign.py` is fixed. It compared
  the two resolutions as **alternatives** when both run, and it matched each late locus to the
  nearest early locus of **any family** — so a late L1 call at a position where L1 had no cns row
  paired with a nearby Alu/SVA row and read as a flip. The barrier left L1 at **300 of 300 rows**;
  it never took anything from L1. The probe now compares the two conflict populations and reproduces
  the run exactly (7 both / 0 disagree / 9 barrier-only / 15 merge-only).
- **The first gate run crashed on a real bug the design missed**: `SVA:transduction` died with a
  bare `KeyError: 31329020`. **A step can open an input it was never handed** —
  `x_TEA_main.py:551` derives `<candidates>.high_confident` from its own `-i`, and
  `update_high_confident_callset` looks every one of its rows up in the transduction output with no
  membership test. The barrier had filtered the candidate list and left the row in the companion.
  `locus_table.CNS_COMPANION_SUFFIXES` now carries it; `.before_calling_transduction` and
  `.before_filtering` are deliberately untouched (written and consumed inside `cns_filter`).

## 3. The gate — job `49502115`, PASS

```
sbatch docs/handoff/runs/stage2c_family_assign_wgs.sbatch
XTEA2_ANALYZE_ONLY=1 XTEA2_OUT=<existing_out> bash docs/…/stage2c_family_assign_wgs.sbatch
```

The second form re-runs the checks against arms that already exist, without recomputing them — for
iterating on the script, not for skipping a real gate.

| leg | result |
|---|---|
| 0 | arms ran the paths they claim — **OK** |
| 1 | arm OFF vs frozen BASELINE-0 VCFs — **4/4 identical**, so no drift |
| 2 | movement ON vs OFF — **0 moved** |
| 2b | **the structural leg**: 7 loci resolved at both stages, **0 disagreements**; 9 barrier-only, 15 merge-only |
| 3 | truth score both arms — **identical**, which follows from 2b having found nothing to move |

**The gate script itself had two reporting bugs that made it fail a good run**, both fixed: it
grepped `<sample>_xtea2.log` for a line `log()` writes to **stdout** (driver-level lines are in
`driver.log`), and its per-type VCF glob matched nothing because those live in
`merged/<bam>_{LINE1,ALU,SVA}.vcf` — so leg 2 compared two empty sets and called it agreement.
There is now a guard that an empty call set is not a pass.

## 4. NEXT

1. **Decide whether `--family-assign-early` should become the default.** It is proven
   answer-preserving but currently buys nothing on its own — it is a second decision point, not a
   relocated one. The honest case for flipping it is that it is a prerequisite for deleting the
   merge's resolution; the honest case against is that until then it is dead weight. My read: leave
   it OFF and flip it together with the change that makes it load-bearing.
2. **Then the rest of 2c**, still open and unchanged in scope: the single per-locus loop that
   replaces the three subprocess chains, the shared gather manifest (plan §4 2b), and
   `--emit-per-type` — deliberately last, since flipping per-type VCFs to opt-in breaks every gate
   script that `cmp`s them.
3. **Then Stage 3**, whose gate is now fully specified (§1): position-level per family plus the
   orphan/sibling/no-`TD_SRC` strata, no transduction pass/fail.

## 5. What NOT to redo

Plan §10 is the authority and gained thirteen entries this session. The ones most likely to cost a day:

- **`scan_clip_realign.tsv` is not per-family evidence** — 504 of 155,615 well-supported positions
  differ across families. Do not build a locus evidence table on it.
- **`two_side_tprt_both` is written by `cns_filter`**, so that is as early as the family rule can go.
- **Stage 2c is not a speedup** — 0.6–5% contention, on ~21% of the run.
- **The barrier and the merge BOTH run** — the barrier does not replace the merge. Reasoning about
  it as though the merge were gone predicts movement that does not happen.
- **12 of the 22 family conflicts are unreachable before transduction**, because the losing family's
  record is created later by transduction/sibling. No decision point before transduction can resolve
  them.
- **Do not compare two family assignments by nearest-locus matching across families** — it pairs a
  call with a different family's nearby row and reads as a flip. Compare the conflict populations.
- **`run_xtea.log()` writes to stdout, not `<sample>_xtea2.log`**, and the per-type VCFs are at
  `merged/<bam>_{LINE1,ALU,SVA}.vcf`. Both cost a gate run.
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
| tests | `tests/test_locus_table.py` (+9), `tests/test_family_assign_early.py` (new, 7) |
| the offline comparator | `docs/handoff/runs/probe_early_assign.py` — no BAM; compares the two conflict populations |
| **the gate (PASSED)** | `docs/handoff/runs/stage2c_family_assign_wgs.sbatch`; output `/n/scratch/.../stage2c_family_assign_49502115/` |
| the benchmark | `benchmark/benchmark_truth_set.py`; tests `tests/test_benchmark_truth_set.py` |
| BASELINE-0 | `/n/scratch/users/c/cos689/xtea2_testdata/baseline0_hg002_L1-Alu-SVA/` (intermediates kept) |
| truth sets | `/n/data1/hms/dbmi/park/corinne/claude_sandbox/test_data/HG002/truth/` |
