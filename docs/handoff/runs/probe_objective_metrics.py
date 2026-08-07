#!/usr/bin/env python3
"""The four axes this project was started on, measured -- offline, from artifacts already on disk.

    python docs/handoff/runs/probe_objective_metrics.py

No allocation, no BAM, no pipeline run. Seconds. Reads `sacct` and the `step_timings.tsv` files
the runs already wrote, and prints wall clock, peak memory, sensitivity and precision for xTea2
against master on the same HG002 WGS BAM.

WHY THIS EXISTS. Plan sec 1's intended outcome is a downstream that is faster and leaner and is
judged on sensitivity/precision against orthogonal truth sets. Every number needed to say whether
that happened has existed since 2026-08-05, in two places nobody looked:

  * `sacct` still holds Elapsed and MaxRSS for all three WGS jobs.
  * `run_xtea.record_timing` has been writing a `children_peak_rss_gb` column into every run's
    `step_timings.tsv` since the instrument was added. **Nothing has ever read it.** It is the
    only per-step memory attribution that exists -- sacct gives one number for the whole job.

THE ONE THING THIS PROBE MUST NOT DO IS FLATTER THE RESULT, so every caveat that weakens a number
is printed next to it rather than in a footnote:

  * master runs the three families SERIALLY (`for t in L1 Alu SVA` in master_wgs_hg002.sbatch,
    each at -n 24); xTea2 runs them CONCURRENTLY. So master's peak is ONE family resident and
    xTea2's is THREE. That makes the memory comparison stronger than it looks, and it means the
    wall-clock ratio partly reflects a parallelism master does not have -- which is an xTea2
    change and is labelled as one, not smuggled in as like-for-like scheduling.
  * sacct MaxRSS is sampled (jobacctgather polls, ~30 s), so a short spike can be missed.
  * `children_peak_rss_gb` is `getrusage(RUSAGE_CHILDREN).ru_maxrss` -- the high-water of the
    LARGEST SINGLE CHILD, monotone, never a per-step delta and never a sum over concurrent
    children. The only sound reading is *which row first raised it*, which is what is printed.
    It is a different quantity from sacct's job-level MaxRSS and the two are not expected to
    agree exactly; that they come out close says one child dominates.
  * both jobs ran under `--mem=64G` and neither was OOM-killed, so MaxRSS is an honest
    high-water rather than a number pressed against a ceiling.

SENSITIVITY AND PRECISION ARE QUOTED, NOT RECOMPUTED. They come from the truth-set scoring
already recorded in the plan and handoffs, with the job that produced each. Recomputing them here
would need the truth sets and the VCFs and would make this probe a second implementation of
`benchmark_truth_set.py` -- the private-copy-of-a-shipping-rule mistake plan sec 10 records three
times. The provenance string is printed with every figure so a reader can go check it.
"""
import os
import subprocess
import sys

SCRATCH = "/n/scratch/users/c/cos689/xtea2_testdata"

# --------------------------------------------------------------------------------------------
# The three WGS runs. Same BAM, same -c 24, same --mem=64G, same partition, whole pipeline each.
# `expect_*` are what was on disk when this probe was written: sacct rows age out of the
# accounting DB, and a probe that silently prints nothing once they do is worse than one that
# says "sacct no longer has this; here is what it said on 2026-08-07".
# --------------------------------------------------------------------------------------------
RUNS = [
    dict(job="49763841", label="master (xTea, 46520ab)",
         tree=None,
         serial=True,
         expect_elapsed="01:10:15", expect_rss_gb=38.03,
         note="three families SERIALLY, each at -n 24"),
    dict(job="49792012", label="xTea2, concat OFF",
         tree=SCRATCH + "/gate_decoy_structural_wgs_49792012/HG002.wgs",
         serial=False,
         expect_elapsed="00:38:32", expect_rss_gb=19.50,
         note="three families CONCURRENTLY, 24 cores total"),
    dict(job="49879677", label="xTea2, concat ON (explicit flag)",
         tree=None,
         serial=False,
         expect_elapsed="00:42:15", expect_rss_gb=19.91,
         note="three families CONCURRENTLY, 24 cores total"),
    dict(job="50038678", label="xTea2, concat ON (via the DEFAULT)",
         tree=SCRATCH + "/gate_decoy_concat_wgs_50038678/HG002.wgs",
         serial=False,
         expect_elapsed="00:38:56", expect_rss_gb=23.32,
         note="SAME configuration as 49879677 -- this pair is the variance control"),
]

# The two runs whose configuration is IDENTICAL. Their spread is the only estimate of
# run-to-run noise this project has, and it is what stops a two-job difference being read as a
# cost. Keyed by job id so the pair stays honest if a run is replaced.
SAME_CONFIG_PAIR = ("49879677", "50038678")

# Truth-set figures, each with the job or commit that produced it. HG002, positions not
# transductions -- sec 7.1's obstacle: the derived transduction denominator is 0 (L1) and 2 (SVA),
# so the truth set scores POSITIONS. With 90 SVA truth records one SVA call is ~1.2 precision
# points, so a small SVA delta is quoted as the calls behind it and not as the number.
TRUTH = [
    ("L1  calls",           "188",   "byte-identical Stage 0 -> Stage 2 (jobs 49441918-22, 49792012)"),
    ("Alu calls",           "1223",  "byte-identical Stage 0 -> Stage 2"),
    ("SVA calls",           "82",    "BASELINE-0; 85 with --types SVA alone"),
    ("SVA TP / FN",         "77 / 13", "UNCHANGED across the structural flip AND the concat flip"),
    ("SVA position prec., BASELINE-0",  "93.9%", "before the 2026-08-05 structural decoy flip"),
    ("SVA position prec., structural",  "90.6%", "3 new calls, all off-truth; 1 VERIFIED 2 UNDECIDABLE (b2a7f41)"),
    ("SVA position prec., + concat",    "91.7%", "SVA FP 8 -> 7; the dropped call is a truth-set FP (49879677)"),
]

# The frozen Stage 3 verifier range, which has TWO values and is the project's most-miscited
# number. Anything quoting the old one is quoting the pre-guard baseline and must say so.
VERIFIER = [
    ("frozen 3-arm range, PRE-guard",  "59.1%-81.0%", "50 VERIFIED; td_gate_baseline_hg002_preguard/"),
    ("frozen 3-arm range, current",    "57.8%-80.4%", "48 VERIFIED, after decision 6's self-window guard"),
    ("concat OFF arm vs ON arm",       "54.0%-76.2%", "IDENTICAL in both arms (49995743)"),
]


def sacct(job):
    """(elapsed, max_rss_gb) from the job's .batch step, or (None, None) if sacct has aged out."""
    try:
        out = subprocess.run(
            ["sacct", "-j", job, "--format=JobID%20,Elapsed,MaxRSS", "--units=G", "-P", "-n"],
            capture_output=True, text=True, timeout=60).stdout
    except (OSError, subprocess.SubprocessError):
        return None, None
    for line in out.splitlines():
        f = line.split("|")
        if len(f) >= 3 and f[0].strip().endswith(".batch"):
            rss = f[2].strip().rstrip("G")
            try:
                return f[1].strip(), float(rss)
            except ValueError:
                return f[1].strip(), None
    return None, None


def high_water_rows(path):
    """The rows that RAISED the children-peak high-water, in order.

    This is the only sound reading of a monotone counter: not "this step used 4 GB" but "the
    high-water was still X when this row finished, and this row is where it moved". A per-step
    delta would be an invention.
    """
    rows, peak = [], -1.0
    with open(path) as fh:
        hdr = fh.readline().rstrip("\n").split("\t")
        if "children_peak_rss_gb" not in hdr:
            return None, hdr
        i_rss = hdr.index("children_peak_rss_gb")
        i_kind, i_name, i_wall = hdr.index("kind"), hdr.index("name"), hdr.index("wall_s")
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if len(f) <= i_rss:
                continue
            try:
                v = float(f[i_rss])
            except ValueError:
                continue
            if v > peak:
                peak = v
                rows.append((f[i_kind], f[i_name], f[i_wall], v))
    return rows, hdr


def hms_to_s(t):
    p = [int(x) for x in t.split(":")]
    while len(p) < 3:
        p.insert(0, 0)
    return p[0] * 3600 + p[1] * 60 + p[2]


def main():
    print("=" * 86)
    print("THE OBJECTIVE, MEASURED  --  HG002 WGS, same BAM, -c 24, --mem=64G, whole pipeline")
    print("=" * 86)
    print()
    print("--- 1. wall clock and peak memory, from sacct ---")
    print()
    print("  {0:<50} {1:>9}  {2:>10}  {3}".format("run", "wall", "peak RSS", "source"))
    print("  " + "-" * 90)
    measured = {}
    for r in RUNS:
        el, rss = sacct(r["job"])
        src = "sacct"
        if el is None:
            el, rss, src = r["expect_elapsed"], r["expect_rss_gb"], "RECORDED 2026-08-07 (sacct aged out)"
        elif rss is None:
            rss, src = r["expect_rss_gb"], "sacct (RSS recorded)"
        measured[r["job"]] = (el, rss)
        print("  {0:<50} {1:>9}  {2:>8.2f} G  {3}".format(
            r["label"] + " [" + r["job"] + "]", el, rss, src))
        print("  {0:<50} {1}".format("", r["note"]))
    print()

    m_el, m_rss = measured["49763841"]
    a_el, a_rss = measured[SAME_CONFIG_PAIR[0]]
    b_el, b_rss = measured[SAME_CONFIG_PAIR[1]]
    d_wall = abs(hms_to_s(a_el) - hms_to_s(b_el))
    d_rss = abs(a_rss - b_rss)
    x_runs = [measured[j][1] for j in ("49792012",) + SAME_CONFIG_PAIR]
    lo, hi = min(x_runs), max(x_runs)
    fast_lo = hms_to_s(m_el) / max(hms_to_s(measured[j][0]) for j in
                                   ("49792012",) + SAME_CONFIG_PAIR)
    fast_hi = hms_to_s(m_el) / min(hms_to_s(measured[j][0]) for j in
                                   ("49792012",) + SAME_CONFIG_PAIR)
    print("  RUN-TO-RUN NOISE FIRST, because it bounds every other number here.")
    print("  Jobs {0} and {1} are the SAME configuration on the same BAM at the same core "
          "count.".format(*SAME_CONFIG_PAIR))
    print("  They differ by {0}s of wall and {1:.2f} G of peak RSS.".format(d_wall, d_rss))
    print("    {0}: {1}  {2:.2f} G".format(SAME_CONFIG_PAIR[0], a_el, a_rss))
    print("    {0}: {1}  {2:.2f} G".format(SAME_CONFIG_PAIR[1], b_el, b_rss))
    print()
    print("  So DO NOT read a two-job wall or RSS difference as a cost. In particular the")
    print("  concat default's cost is NOT (concat-ON wall - concat-OFF wall): that difference")
    print("  is 223s and the noise floor just measured is {0}s. The flag's real cost is ~55s,".format(d_wall))
    print("  decomposed in HANDOFF_2026-08-06_timing_resolved.md from x_transduction.td_phase's")
    print("  own logging with L1 as a same-node same-job NULL ARM -- which is the only way to")
    print("  separate the flag from the node.")
    print()
    print("  xTea2 across {0} WGS runs:  wall {1}-{2},  peak RSS {3:.2f}-{4:.2f} G".format(
        len(x_runs), min(measured[j][0] for j in ("49792012",) + SAME_CONFIG_PAIR),
        max(measured[j][0] for j in ("49792012",) + SAME_CONFIG_PAIR), lo, hi))
    print("  master, n=1:                 wall {0},  peak RSS {1:.2f} G".format(m_el, m_rss))
    print("  => xTea2 is {0:.2f}-{1:.2f}x faster on {2:.0f}-{3:.0f}% of master's peak memory."
          .format(fast_lo, fast_hi, 100.0 * lo / m_rss, 100.0 * hi / m_rss))
    print("  MASTER IS n=1. Its own run-to-run spread is unmeasured, so the ratio's precision is")
    print("  bounded by a variance nobody has characterised. The DIRECTION is not in doubt --")
    print("  the gap is 15-18 G and the noise floor is {0:.1f} G -- but quote the range.".format(d_rss))
    print()
    print("  READ THE MEMORY NUMBER THE RIGHT WAY, because the honest version is the stronger")
    print("  one: master's peak is ONE family resident at a time (it loops L1, Alu, SVA")
    print("  serially); xTea2's is ALL THREE at once. So xTea2 holds three families in a")
    print("  little over half the memory master needs for one.")
    print()
    print("  AND THE WALL-CLOCK RATIO IS NOT PURE PER-FAMILY SPEED. Part of it is that xTea2")
    print("  runs the families concurrently and master cannot. That IS an xTea2 change -- it is")
    print("  not a benchmark artifact -- but it is not the C++ gather, which is separately")
    print("  measured at 15.25x (L1) / 17.82x (SVA) on the gather itself (job 49453857).")
    print()

    print("--- 2. WHICH STEP RAISES THE HIGH-WATER -- the attribution sacct cannot give ---")
    print()
    for r in RUNS:
        if not r["tree"]:
            continue
        path = os.path.join(r["tree"], "step_timings.tsv")
        if not os.path.isfile(path):
            print("  {0}: no step_timings.tsv at {1}".format(r["label"], path))
            continue
        rows, hdr = high_water_rows(path)
        if rows is None:
            print("  {0}: step_timings.tsv predates the RSS column (header: {1})".format(
                r["label"], "\t".join(hdr)))
            continue
        print("  {0}".format(r["label"]))
        print("    {0:<28} {1:>10}  {2:>12}".format("row that raised it", "wall_s", "new peak"))
        for kind, name, wall, v in rows:
            print("    {0:<28} {1:>10}  {2:>10.2f} G".format(kind + " " + name, wall, v))
        print("    {0} row(s) of the whole run ever moved it.".format(len(rows)))
        print()
    print("  `children_peak_rss_gb` is getrusage(RUSAGE_CHILDREN).ru_maxrss: the high-water of")
    print("  the LARGEST SINGLE CHILD, monotone. It is NOT a per-step delta and NOT a sum over")
    print("  concurrent children, so it does not equal sacct's job-level MaxRSS and is not")
    print("  expected to. That the two land close says one child dominates the job's footprint.")
    print()

    print("--- 3. sensitivity and precision (quoted, with provenance) ---")
    print()
    for name, val, prov in TRUTH:
        print("  {0:<34} {1:>9}   {2}".format(name, val, prov))
    print()
    print("  The truth set scores POSITIONS, not transductions -- sec 7.1: HG002's derived")
    print("  transduction denominator is 0 (L1) and 2 (SVA). With 90 SVA truth records one call")
    print("  is ~1.2 precision points, so quote small SVA deltas as the calls behind them.")
    print()

    print("--- 4. the Stage 3 verifier's own instrument ---")
    print()
    for name, val, prov in VERIFIER:
        print("  {0:<34} {1:>12}   {2}".format(name, val, prov))
    print()
    print("  ANY TEXT QUOTING 59.1%-81.0% IS QUOTING THE PRE-GUARD BASELINE and must say so.")
    print("  The two entries decision 6 removed were the pool's two HIGHEST-support ones, so the")
    print("  range does not simply shift down.")
    print()
    print("=" * 86)
    return 0


if __name__ == "__main__":
    sys.exit(main())
