# xTea2 handoff — 2026-07-28 (cns decision, C++ cns port, trim step 5)

**This is the current entry point.** Supersedes `HANDOFF_2026-07-26_wgs_run.md`, whose run table,
timings and master baseline are still accurate and still the reference — read it for those. What
changed is the *open questions*: the two it was built around are both answered, both answers are
different from what it predicted, and one of them (the cns path) went on to be rewritten. Sections
are in the order the work happened, so "The two answers" below is the earlier part of the session
and the C++ port supersedes its conclusion about which collection path to keep.

One-line status: **the cns question is settled — it was a master bug (B10), not a sensitivity
trade — and B10 is fixed and WGS-validated. The candidate-agnostic "fan" is DELETED and cns
collection is now a C++ port of master's per-candidate design. The realign hypothesis is refuted.
The trim is DONE, all five steps. The one thing owed is re-baselining the goldens, which is no
longer optional.**

**Read first if you are picking this up mid-flight:** WGS job `49070030` is the three-way
comparison (master vs Python-cns vs C++-cns) and may still be running. Run
`python docs/handoff/compare_wgs_three_way.py` — it needs no arguments, finds the three output
trees itself, and skips any arm that is absent. `--quick` omits the two read-level sections.

## The two answers

### 1. The cns paths collect the SAME READS. The difference was a double-count.

The 07-26 handoff framed this as a judgment call about sensitivity — "does the C++ fan collect
reads master's indexed fetch misses, or reads it should not?" — and read the 18-of-19 loss bias as
evidence the fan was more permissive. **Both halves were wrong.** Comparing the two completed WGS
runs at record level (primary records, keyed on everything the collector decides, SEQ
canonicalised for strand):

| | fan (`48825976`) | indexed (`48842713`) | fan-only | idx-only |
|---|---|---|---|---|
| clip L1 | 273,604 | 273,604 | 0 | 0 |
| clip Alu | 604,289 | 604,289 | 0 | 0 |
| clip SVA | 400,139 | 400,139 | 0 | 0 |
| disc Alu | 583,489 keys | 583,489 keys | 0 | — |

Clip is **set-identical** for all three types. Disc has an **identical key set**. Neither path
sees a read the other misses. What differs is multiplicity: in the indexed run 28,071 Alu disc
keys appear twice and 102 appear three times, 4.8% of the total. The fan emits every key exactly
once.

**Root cause = B10** (spec §4). `x_alignments.extract_mate_reads_of_region` writes one FASTA
record per matching mate **alignment**, not per mate **read**, with no supp/sec filter — so a mate
with a soft-clipped supplementary alignment in the same bin gets a byte-identical second record
(`is_read1` is set on it, and soft-clipped supplementaries carry the full `query_sequence`).
Separately the bins overlap, because `fetch` returns reads overlapping the bin while
`_parse_read_name_pos_of_region` admits mates 100 bp beyond either edge. `parse_disc_algnmt_consensus`
counts records, so one read became two units of disc support. Duplicates concentrate in
pericentromeric chr16 (10,170) and chr17 (12,360) — where supplementary alignments cluster.

Master's own clip collection already skips supp/sec (`x_clip_disc_filter.py:161,201`), and so does
`cns_emit.hpp`'s mate gather. **This path was the outlier**, which is exactly why it surfaced as a
difference between the two paths rather than as an obvious defect in either.

Fixed in `ee0742d`: skip supp/sec, and own each alignment by its start so the bins partition.
Regression tests in `tests/test_disc_mate_extract.py` (synthetic BAM, no external data, fast
suite; the three bug cases were confirmed to FAIL without the guards).

**The residual noise is not a second bug.** 602 of 148,122 Alu clip keys differ in record count by
±1 between the runs. That is bwa, not collection: identical query sequences get different
**equal-scoring** placements (verified `AS=23` for `ALU_v8:282 33M`, `ALU_v2:322 33M` and
`ALU_v2:287 23M10S` on one 33 bp query) because the two paths feed bwa the same reads in a
different order and the tie-break is order-dependent. That reshuffles which consensus copy a
multi-mapper lands on, which is what tips borderline candidates across the cns cutoffs. Expect it
to persist; it is not fixable by making the two paths agree.

### 2. The realign bottleneck is not the parent parse, and bwa scales fine

`realign.hpp` had no timers and sent bwa's stderr to `/dev/null`, so the phase was unobservable.
It now reports, per family, spawn / time-to-first-record / total / **time blocked in `getline`
vs time in the loop body** / records / mapped / bytes, plus bwa's own log to
`<realign_out>.<FAMILY>.bwa.log`. `run_bwa_sam` (the cns path) gets the same treatment.

The getline/body split is the discriminator, and it is decisive: **`getline` 100%, body 0.1 s**
against 18–52 s of bwa per family. The parent's own work cannot be a bottleneck in a 20-minute
phase. **Do not spend a run on bwa → file + parallel parse.**

The scaling benchmark (`49066537`, `runs/bwa_scaling_bench.sbatch`) then closed off the other
mechanical options. On 24 real cores bwa mem scales **linearly** on this workload — SVA 1.00 /
2.00 / 3.97 / 7.89 / 11.83 / 23.19 at `-t 1..24`, `cpu/wall` tracking the thread count exactly. So
`bwa_threads_for`'s threads/nfamilies split is already correct, and sharding buys nothing (all six
shard counts tied). The ~18–20 min is genuine bwa work at full parallelism. **Remaining levers
there are algorithmic — fewer reads to realign, or less permissive params (`-T 7 -k 7 -c 70` yields
~4.5 alignment records per read) — not mechanical.**

> ⚠ If you benchmark anything in a Claude Code session, submit it as a batch job. The session holds
> **2 cores**, which caps every configuration at ~1.9× and reads exactly like a serial section
> inside the program under test. The first pass at the above was misread that way.

The one free win it found is memory: `-K 10000000` → `1000000` cuts bwa's resident from 2.59 GB to
0.33 GB at no time cost, three concurrent, so ~7.8 GB was being held for nothing on a run whose
peak is what forced `--mem=64G`. Alignment records byte-identical, verified for all four param sets
the pipeline uses (only bwa's `@PG CL:` header differs, and the cns parsers read via pysam).

## The cns collection is now C++, and the fan is gone

**This supersedes the fan entirely.** After B10 was fixed, the two collection paths were
provably equivalent, which turned the choice into an engineering one — and the fan lost it. Rather
than keep optimising a candidate-agnostic collector, the decision (user, 2026-07-28) was to keep
**master's algorithm** and move it to C++.

**Why the algorithm, not the language, was the problem.** ~2,800 candidates at `extnd=753` is ~4 Mb
of reference. The fan streamed all 3.1 Gb because it collected during the detection scan, *before
candidates existed*, so it could not know what to skip. Every cost that made it lose followed from
that single decision:

| | fan | reason |
|---|---|---|
| disc mates gathered | 4,151,548 | to serve the ~0.5M any type reads |
| spill written | 1.10 GB | to carry the collection across a process boundary |
| reference touched | 3.1 Gb | candidate-agnostic, so nothing could be skipped |
| gather threading | single-threaded | a whole-genome pass, serially |

Knowing the candidates first removes all four at once. `emit_cns_multi` now fetches each candidate
window with the index, with two cost-only refinements over master: windows are **unioned across
types and merged** (with an exact ownership rule — a read belongs to the region that fetched it iff
its start is at or after the previous region's end, which is provably collect-once), and the disc
mate gather **fetches recorded mate positions** instead of scanning 50 Mb bins.

**Deleted**: `collect_cns_read`, `load_cns_spill`, `lowq_near_count`, `emit_cns_from_spill` (431
lines), the `--cns-collect-*`/`--cns-fan-*` options, the spill buffers and writers, and
`scan_tile`'s spill hook — which also removed its `file_idx` parameter, since tagging the spill was
all it was for. A note sits where the code was. `--cns-emit-scan` **survives and now means the fast
path**; it invokes `--cns-emit-manifest`, which existed in C++ but had never been wired into
`run_xtea.py` — only the fan had, so `emit_cns_multi` had never run end to end before this.

### Ordering is not cosmetic — this is the part worth remembering

The port initially emitted the clip FASTQ **read-major** (per read, all candidates it overlaps).
Master emits **candidate-major** (one candidate's window at a time). Same records, different order
— and **bwa's tie-break among equal-scoring consensus placements depends on a read's index in its
input**, so regrouping silently moved which consensus copy multi-mapping clip fragments landed on.
Measured on chr1 against the Python step:

| | read set | + placement |
|---|---|---|
| read-major | identical | 570 / 2021 / 1751 differ (L1/Alu/SVA), **1 extra SVA call** |
| candidate-major | identical | **identical** |

Zero, after inverting the nesting. If you touch the emission order, re-run that comparison.

`n_cnt` (field 5) is deliberately **not** reproduced: master builds the right-clip read name and
increments the counter *before* applying that side's gates (`x_clip_disc_filter.py:205-218`), so a
failing right clip burns an index and leaves a gap. `parse_clip_realignment_consensus` reads fields
0,1,2,4,6 and never 5, so the values are free — only their order matters.

### The one accepted difference

The **disc** FASTA is still emitted in qname order where master emits it mate-position-major in
header-contig order: 115 / 325 / 80 records placed differently on chr1, which is what costs one
extra SVA call there (L1 and Alu are exactly equal). **The user accepted this rather than requiring
byte-identity.** Matching it is the same shape of fix as the clip one, with one residual that is not
reachable: two reads at the same mate position are ordered by BAM order in master and by qname here.

chr1 result: 274 s vs the Python step's 303 s, extra scratch 0 vs 2.0 MB, peak RSS 1.30 GB for both
— the peak is set by downstream genotyping, so cns was never the binding constraint at that scale.
**Whether that holds at WGS is what `49070030` answers.**

## The trim is finished

Step 5 (`e4a6050`) removed all **125** dead locals — down from the 184 recorded pre-step-4;
`x_orphan_transduction` alone went 55 → 29 because step 4 deleted the modules the old count
included. Net −151 lines (171 deleted, 20 comment lines added at the four call sites).
**`xtea/` is now 20,196 lines / 42 modules, from 32,963 / 58 (−38.7%)** (20,205 at the trim commit; the port added nine net lines to the Python side).

> Small correction while re-measuring: the 07-26 handoff recorded 20,273 lines post-step-4. The
> tree at `a152aca` (immediately pre-step-5) measures **20,356** by
> `find xtea -name '*.py' -exec cat {} + | wc -l`. Use the measured figures; the reduction is
> −38.7%, not the −39.0% that the older number implies.

The plan's RHS hazard is mechanised rather than eyeballed:
[`classify_dead_locals.py`](classify_dead_locals.py) walks the AST and splits sites by RHS purity
(121 provably call-free, 4 with calls), and refuses anything whose deletion would empty a block
(15 sites) instead of discovering it as a `SyntaxError`. Reuse it — re-run it after any deletion
pass, since deletions cascade (three passes to converge here).

Four call sites and 15 block-level sites were hand-resolved; two functions
(`x_alignments.cvt_alignmt_to_fa`, `x_contig.add_to_final_rslt`) turned out to be nothing but dead
locals with no callers and went with them. The best single find:
`x_orphan_transduction`'s `xannotation = self.prep_annotation_interval_tree(...)` sat in a
`pool.map` worker invoked **once per record**, rebuilding the whole RepeatMasker index per record
and discarding it.

Gates used, worth reusing: **pyflakes finding set diffed against HEAD** (125 gone, zero new — this
is what proves nothing deleted was read elsewhere), step plan byte-identical (18 steps), import
sweep 42/42, fast suite 88 passed / 1 pre-existing skip.

## NEXT SESSION — start here

### 1. Re-baseline the goldens. This is now required, not owed.

It was already outstanding from step 2 (checklist item 8) and deferred twice, and the 07-26
handoff correctly sequenced it *after* the cns decision. That decision is made, so the dependency
is discharged and nothing else blocks it.

It is no longer just a bookkeeping debt: **every golden on file was produced with B10's double
count.** The fix is strictly more conservative — sites that only cleared the cns cutoffs on
double-counted reads now fall out — so the goldens sit *above* the current pipeline by
construction. Until they are regenerated, every concordance gate reports a loss that is a fix.

Regenerate the WGS seeded set. `tests/test_downstream_from_cns.py` and
`tests/test_e2e_concordance.py` both compare against pre-B8-fix goldens today.

**Before regenerating, check the dropped sites' disc support** — that is the one thing that would
change the plan. If sites are being lost whose *deduplicated* support still clears the cutoffs,
something other than B10 moved and the fix is not the whole story.

**Pre-worked: the blast radius is small and bounded.** `b10_blast_radius.py` (run against the
pre-fix run `48842713`) counts, per candidate, the *mapped* primary disc records the filter actually
saw versus the deduplicated set — mapped only, because an unmapped mate contributes no consensus
position and duplicating it changes nothing:

| type | candidate sites | inflated | % | median inflation | worst |
|---|---|---|---|---|---|
| L1 | 1,891 | 138 | 7.3% | 1.07× | 2.00× |
| Alu | 8,868 | 531 | 6.0% | 1.07× | 2.00× |
| SVA | 1,478 | 54 | 3.7% | 1.09× | 2.00× |

Of Alu's 531 inflated sites, **431 lost exactly one read** and only **95 end up with ≤6
deduplicated mapped disc reads** — the band where a 1–2 read change can actually cross a cutoff.
So expect the re-baselining delta to be tens of sites, not hundreds, which is consistent with the
23 positions that separated the two runs. **Those 95 are the inspection list** — if a site outside
that band moved, B10 is not what moved it.

### 2. Read the WGS three-way comparison, and decide whether the C++ cns path ships

`--cns-emit-scan`'s old question is **resolved**: the fan is deleted and the flag now selects the
C++ port. What is open is whether the port is worth keeping over the Python step it replaces, on
the bar the user set — **memory, speed, disk**.

```bash
python docs/handoff/compare_wgs_three_way.py          # all three arms, three levels
python docs/handoff/compare_wgs_three_way.py --quick   # calls + timing only
```

| arm | job | what |
|---|---|---|
| `master` | `48825977` | master `46520ab` **as published** — still carries B10's double-count |
| `python` | `49065573` | xTea2, Python cns, B10-fixed. 1:00:37, 21.99 GB (sacct) |
| `cpp` | `49070030` | xTea2, C++ cns collection |

Master is deliberately left unfixed, which makes master-vs-either a **cumulative** comparison
across everything the fork changed, not an attribution. **`python`-vs-`cpp` is the clean isolate** —
those two differ only in which code collects the cns reads.

Already known from the two completed arms, master → python:

| type | master | python | diff |
|---|---|---|---|
| L1 | 203 | 199 | +15 / −19 |
| Alu | 1262 | 1264 | +12 / −10 |
| SVA | 73 | 67 | **+0 / −6** |

SVA loses six and gains nothing, which is exactly the shape B10 predicts — removing double-counted
disc support can only push a site below a cutoff, never above one. L1 and Alu move both ways, so
they carry other fork changes (B8, the trim, col5) and are **not** attributable to B10 alone.

**How to read the result.** chr1 gave −10% wall, disk 2.0 MB → 0, and memory unchanged. Expect the
same shape, but note two things chr1 could not settle:
- **Memory is the axis the user named first and the one chr1 could not test** — its 1.30 GB peak was
  set by downstream genotyping, not cns. WGS is where the Python path's 24 forked pool workers would
  show up, if they do.
- `downstream_parallel` should shrink because `cns_filter` now runs on the bypass (parse+filter
  only), but transduction and post-filtering dominate that phase at WGS scale and are untouched. If
  it does not shrink, the wall-clock win is smaller than chr1 implied.

If it clears the bar, the Python collection path becomes the thing to delete. If it only ties, the
honest call is to keep Python and drop the port — the 31-minute win was B10, not this.

### 3. If more realign speed is wanted, the target is the read count

Per item 2 above, the mechanical options are closed. What is left:
- `-T 7 -k 7 -c 70` produces ~4.5 alignment records per read against a tiny consensus. Nothing has
  ever measured what a less permissive `-c` costs in recall.
- Per-family cost is wildly asymmetric — SVA and L1 emit ~3× the records of Alu from the *same*
  FASTQ (SVA 878k / L1 853k / Alu 299k on a 2 Mb window). An equal core split across families is
  therefore not matched to the work, even though each bwa scales linearly. Weighting
  `bwa_threads_for` by expected records would shorten the critical path.

### 4. Housekeeping

- `--scan-clip`/`--scan-disc` remain accepted no-ops (≈15 scripts under
  `/n/scratch/users/c/cos689/xtea2_testdata` pass them). Fine, but they are now the only
  vestigial flags left.
- Trim step 3 (attribute B7) is still skipped by user decision. If picked up, `TRIM_PLAN`'s master
  path is wrong: the real untouched checkout is `/home/cos689/xtea_exploration/xTea/`.
- Worktrees, all under `/n/scratch/users/c/cos689/xtea2_testdata/`: `port_wt` (`b44a85f`) is what
  `49070030` runs from — keep until item 2 is decided. `master_cmp_wt` (`46520ab`) backs the master
  arm and should be kept as long as that comparison matters. `step45_wt`, `b10_wt` and
  `xtea_scan_prerework` are spent and can go.
- `docs/handoff/runs/` has every sbatch used, and `b10_validate_49065573.out` /
  `bwa_scaling_49066537.out` are archived next to their scripts.

### What NOT to redo

- **Do not re-open the cns sensitivity question.** The two paths collect the same reads; measured
  at record level, all three types. See the table above.
- **Do not implement bwa → file + parallel parse.** The parent's parse is 0.1 s per family.
- **Do not shard the realign FASTQ for speed.** bwa already scales linearly to 24 threads; all six
  shard counts tied. (The mapped-qname set *is* shard-invariant, verified — so it stays available
  as a correctness fact, just not as a speedup.)
- **Do not re-litigate concurrency, or run types serially at `-n 24`.** Both refuted earlier; see
  the 07-26 handoff.
- **Do not re-establish the master baseline.** `48825977`, 50:20, `master_serial_wgs.sbatch`.
- **Do not benchmark in-session.** 2 cores. See the warning above.
- **Do not rebuild candidate-agnostic cns collection.** It is deleted and the note where it lived
  explains why: a collector that runs before the candidates exist cannot know what to discard, so it
  pays for everything. That is a property of the design, not of the implementation.
- **Do not emit the clip FASTQ read-major.** It is order-equivalent to nothing: bwa's tie-break is
  input-index dependent, so it moves placements and calls. Candidate-major matches master exactly.
- **Do not trust `--dry-run` to exercise `--cns-emit-scan`.** The invocation is guarded by
  `not opt.dry_run`, so a passing dry-run says nothing about that path. The cheap real check is a
  one-line manifest over ~40 candidates run directly against `xtea_scan` (0.4 s).

## This session's commits

| commit | what |
|---|---|
| `ee0742d` | B10: stop double-counting disc support -- master's mate FASTA emits per ALIGNMENT, not per read |
| `a152aca` | Make the realign phase observable, and refute the parent-parse hypothesis |
| `e4a6050` | Trim step 5: all 125 dead locals, plus the bwa -K memory win the benchmark turned up |
| `e2e977d` | Handoff: both open questions answered, and both differently than predicted |
| `7e7837a` | Handoff: correct the post-trim line count to the measured figure |
| `bf23149` | Pre-work item 1: bound B10's blast radius before re-baselining starts |
| `d326673` | Rework the cns fan's mate gather: filter before gathering, then fetch by position in parallel |
| `6c1d553` | Port master's cns collection to C++/htslib, and delete the fold-in spill path |
| `17cf30d` | Fix the cns-emit option name, and note why --dry-run could not have caught it |
| `b44a85f` | Emit the clip FASTQ candidate-major, which removes ALL clip placement divergence from master |
| `99da02c` | WGS three-way comparison: master vs Python-cns vs C++-cns |

Jobs, all NA12878, `-c 24 --mem=64G -p priopark`:

| job | what |
|---|---|
| `49065573` | B10 WGS validation — PASS, zero dup disc records, disc set matches the fan exactly |
| `49066537` | bwa thread scaling + sharding on the real realign workload |
| `49068968` | first port comparison — its Python arm succeeded; its ported arm failed on my option-name error |
| `49069598` | ported chr1, read-major — read sets identical, 570/2021/1751 clip placements differed |
| `49069881` | ported chr1, candidate-major — clip placement **identical**, one SVA call from disc order |
| `49070030` | **WGS three-way arm 3** (C++ cns), the one still to read |

Also fixed in passing, because it silently invalidates C++ work: `src/csrc/Makefile` did not list
`realign.hpp` / `cns_emit.hpp` as prerequisites. They carry most of the code but are `#include`d
rather than compiled, so after a header-only edit `make` saw `xtea_scan` newer than
`xtea_scan.cpp`, printed nothing, and left the old binary in place.

## Where things live

| | |
|---|---|
| **the WGS three-way comparison** | `compare_wgs_three_way.py` (no arguments; `--quick` for calls+timing) |
| the cns/perf findings that still stand | `PERF_2026-07-26_why_xtea2_is_slower.md` (timings, master baseline) |
| the run table and the four WGS jobs | `HANDOFF_2026-07-26_wgs_run.md` |
| the bug log (B10 is the newest entry) | `../clipdisc_contract_spec.md` §4 |
| dead-local classifier | `classify_dead_locals.py` |
| B10 blast-radius sizer | `b10_blast_radius.py` |
| star-import safety checker | `check_star_import_safety.py` |
| run scripts | `runs/` |
| per-step timings | `<out>/NA12878.wgs/step_timings.tsv` |
