# Why xTea2 is slower than master, and what to do about it

Measured 2026-07-26. Both runs: NA12878 whole genome, one job, `-c 24 --mem=64G -p priopark`,
same node class, same BAM.

| | wall | CPU-time | CPU% of 2400% | MaxRSS |
|---|---|---|---|---|
| master `46520ab`, serial, `-n 24` (`48825977`) | **50:20** | 8:41:09 | 1036% | 40.6 GB |
| xTea2 trim steps 1-4 (`48825976`) | **1:31:21** | 13:51:37 | 910% | 24.3 GB |

**xTea2 is 1.81× slower in wall clock and burns 1.60× the CPU time.**

## First: this is not a regression from the trim

xTea2's best-ever WGS run was **1:18:23** (`48110538`, pre-trim), still 1.56× slower than
master-serial's 50:20. The trim added ~13 min on top (of which some is the 192G→64G change and
node variance), but **xTea2 has been slower than the code it replaces for its whole existence.**
Nobody saw it because the only comparator on file was the golden's 38:31 — and that was three
concurrent 8-core jobs, so it was read as "xTea2 is 2× slower but the golden had free
concurrency", which sounded like a scheduling artefact rather than a real deficit. It was real.
The `-n 24` serial comparator is what made it visible.

It also kills the concurrency hypothesis. xTea2 ran at 910% CPU and master at 1036% — nearly
identical. xTea2 is not idling; **it is doing 60% more work.**

## Where the time goes

xTea2 phases (non-overlapping, sum to the run):

```
scan_pass            27.34 min  30.7%
disc_pos_load         0.07 min   0.1%
cns_fan              30.81 min  34.6%   <-- the problem
downstream_parallel  30.24 min  34.0%
shared_tail           0.48 min   0.5%
TOTAL                88.94 min
```

master, per step label, summed over the three types (they ran sequentially, so these do sum):

```
transduction  17.30    clip  15.63    disc  7.23    cns_filter  5.47
sibling        2.02    postF_hc 0.90  postF_all 0.90  gene 0.67   gvcf/genotype 0.16
TOTAL ~50.3 min
```

Aligning the two on what they actually compute:

| stage | xTea2 | master | verdict |
|---|---|---|---|
| candidate detection | ~2 min (C++ scan) | 15.63 min (clip) | **xTea2 ~8× FASTER** |
| clip realignment | ~25 min (rest of scan_pass) | (inside the 15.63) | xTea2 slower |
| disc | 0.07 min (disc_pos_load) | 7.23 min | **xTea2 ~100× FASTER** |
| **cns collection** | **30.81 + 4.15 = 34.96 min** | **5.47 min** | **xTea2 6.4× SLOWER** |
| downstream (td/sib/postF) | 30.24 min wall | ~27.3 min | roughly a wash |

The C++ rewrite did exactly what it was built to do: detection and disc are enormous wins. Those
wins are then handed back, with interest, by one stage.

## Root cause 1 — `--cns-emit-scan` reads all 39 GB to service ~2,800 point queries

`cns_fan` is `emit_cns_from_spill`'s "ONE mate-gather pass". In `src/csrc/cns_emit.hpp`:

```c
else if (idx) it = sam_itr_queryi(idx, HTS_IDX_START, 0, 0);   // = iterate the ENTIRE BAM
```

It walks every read in the genome to find mate sequences for the anchors loaded from the spill.

Master never does this. Its `cns_filter` fetches an `iextnd` window around each candidate site
through the BAM index, and picks up mates inside those windows. The number of sites it fetches for
is small — this run's own candidate counts:

| type | xTea2 cns | master cns |
|---|---|---|
| L1 | 297 | 287 |
| Alu | 2221 | 2128 |
| SVA | 273 | 268 |

**~2,800 sites total.** Reading 39 GB sequentially to answer ~2,800 point queries is the wrong
algorithm, and the measurement is unambiguous: 34.96 min vs 5.47 min for the same result.

The design note in `run_xtea.py` says the fold-in makes it "2 BAM passes total, not 6" — which is
true and was the wrong thing to optimize. Six *indexed, targeted* fetch sets are far cheaper than
two *whole-genome* passes. The comparison that justified it was against xTea2's own earlier
structure, never against master's indexed path.

## Root cause 2 — `--parallel-downstream` starves transduction

Transduction is the dominant downstream cost in both pipelines, and it is where the thread split
bites:

| | transduction, summed over types | threads per step |
|---|---|---|
| xTea2 | 52.45 min | 8 (`cores // ntypes`) |
| master | 17.30 min | 24 |

**52.45 / 17.30 = 3.03×, and 24 / 8 = 3.0×.** Transduction scales essentially linearly with `-n`,
so this is not slower code — it is the same code given a third of the machine. Running three
chains at 8 threads each *should* wash out to the same wall clock, and roughly does (30.24 vs
27.3 min), but only because the types are imbalanced: Alu is 11.28 min of master's work against
L1's 20.43 and SVA's 18.57, so the Alu worker finishes early and its 8 cores sit idle while L1 and
SVA grind on at 8 threads each instead of 12.

That idle tail is most of the gap between xTea2's 910% and the 2400% available.

## Root cause 3 — clip realignment is slower than master's whole clip step

`scan_pass` is 27.34 min, of which detection is ~2 min: the spill files (`cns_spill_clip.txt`
645 MB, `scan_disc_pos.tsv` 380 MB) were complete 2 min after the phase started, and nothing else
was written until the phase ended. So **bwa realignment is ~25 min** — more than master's entire
clip step (15.63 min), which includes both its BAM scan *and* its own realignment.

This one is diagnosed less firmly than the other two, because `realign.hpp` sends bwa's stderr to
`/dev/null` and there is no timing inside the phase. The volume is 15.1M rows in
`scan_clip_realign.tsv` (~5M distinct positions × 3 families). The likely mechanism is that
`realign_one_file` parses bwa's SAM in the **parent, single-threaded, line by line through a
pipe** (`fgets` on a `FILE*`), one reader thread per family — so three reader threads must consume
what 24 bwa threads emit. That is a classic single-consumer bottleneck and would cap throughput
regardless of how many threads bwa gets.

**This is a hypothesis, not a measurement.** Confirming it needs a timer around the bwa spawn and
the parse loop separately, which does not exist yet.

## Proposed fix, in order of value

**1. Stop using `--cns-emit-scan`; make cns collection indexed. (~29 min, the whole gap.)**
The immediate step is to drop the flag — master's per-site fetch path is still present and is what
runs when the flag is absent, so this needs no new code. Expected 1:31:21 → ~1:02.
The durable fix is to keep the C++ fan but give it the candidate sites and use
`sam_itr_querys` per site (or one merged sorted interval list) instead of `HTS_IDX_START, 0, 0`.
That keeps the "collect once for all three types" benefit while paying only for the regions that
matter. ~2,800 windows at `cns_extnd` ≈ 750 bp is a few MB of BAM, not 39 GB.

**2. Fix the downstream thread allocation. (~5-10 min.)**
Three options, cheapest first:
   - Run types **serially at `-n 24`** (drop `--parallel-downstream`). Master's shape; transduction
     gets the whole machine. Best if transduction keeps scaling to 24.
   - Keep concurrency but **weight the split** by expected cost (L1 and SVA get more than Alu)
     instead of `cores // ntypes`.
   - Better: a **shared worker pool** across types so a finished type's cores are reused, rather
     than a fixed per-type slice. Removes the idle tail entirely.

**3. Instrument and then fix the realign parse. (unknown, up to ~20 min.)**
First add timers: bwa spawn → first record, and total parse time, per family. `realign.hpp` throws
bwa's stderr away today, so its own timing report is lost too. If the parent parse is the
bottleneck, the fix is to have bwa write SAM to a file and parse it with multiple threads, or to
split the FASTQ and run more bwa processes with fewer threads each.

Do **not** start here. It is the least certain and the least valuable of the three.

## Attribution runs in flight

| job | config | isolates |
|---|---|---|
| `48842713` | no `--cns-emit-scan`, keeps `--parallel-downstream` | fix 1 alone |
| `48842714` | no `--cns-emit-scan`, no `--parallel-downstream` (`-n 24` serial) | fix 1 + fix 2 |

Read them with:

```bash
sacct -j 48842713,48842714 --format=JobID%14,JobName%20,State%12,Elapsed,MaxRSS%10
# per-phase breakdown for each:
grep ^phase /n/scratch/users/c/cos689/xtea2_testdata/fix_nofan_out/NA12878.wgs/step_timings.tsv
grep ^phase /n/scratch/users/c/cos689/xtea2_testdata/fix_nofan_serial_out/NA12878.wgs/step_timings.tsv
```

Predictions, recorded before the results so they can be wrong:
- `48842713`: `cns_fan` disappears, `cns_filter` rises to ~5-6 min, total **~1:02**.
- `48842714`: same front end; downstream drops from 30.24 min to ~22-25 min, total **~55-58 min**,
  i.e. roughly at parity with master-serial's 50:20.

If `48842714` does not land near master, then something outside these three causes is also wrong
and the phase table in its `step_timings.tsv` is where to look next.

**Caveat on both:** they are timing runs on a shifted baseline (the B8 fix), and dropping
`--cns-emit-scan` changes which code produces the cns SAMs. Their call sets should be compared to
`48825976`'s before either config is adopted — a faster wrong answer is not a fix.
