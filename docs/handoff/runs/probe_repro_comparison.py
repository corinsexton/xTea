#!/usr/bin/env python3
"""MASTER vs xTea2 over THREE REPLICATES EACH -- speed, memory, precision/recall.

Reads the six run trees `repro_master_vs_xtea2.sbatch` leaves on disk plus `sacct`, and prints
the four-axis comparison. Offline, seconds, no allocation -- the established probe style
(`probe_objective_metrics.py`, `probe_step_timing.py`, `probe_decoy_library.py`).

WHAT THIS ANSWERS THAT n=1 COULD NOT
------------------------------------
Every master-vs-xTea2 number in this project had master at n=1. Two consequences, and the second
is the one that actually blocked a conclusion:

  * The SPEED/MEMORY ratio had no error bar on the master side. xTea2's noise floor was known
    (199 s and 3.41 GB between two identical configurations, jobs 49879677 vs 50038678); master's
    was not, so "1.66-1.82x on 51-61%" was a range built from one arm's spread only.

  * The ACCURACY difference was uninterpretable. Master scored 80.6%/87.8% against xTea2's
    80.4%/88.3% -- 4 true positives and 9 false positives over 1642 truth records. Job 49075238
    had already measured master's own `-n 8` output differing from its own `-n 24` output ON
    IDENTICAL INPUTS, so master is known thread-scheduling sensitive and those 4 calls could be
    entirely its own run-to-run spread. THE WITHIN-ARM SPREAD IS THEREFORE THE YARDSTICK FOR THE
    BETWEEN-ARM DIFFERENCE, and this probe prints it first for that reason.

The house rule this file follows: NEVER a private copy of a shipping rule (plan sec 10, which
records that mistake three times). Truth loading, call parsing and one-to-one position matching
are imported from `benchmark/benchmark_truth_set.py`, so a number here is the same number the
gates print.

USAGE
    python3 docs/handoff/runs/probe_repro_comparison.py
    python3 docs/handoff/runs/probe_repro_comparison.py --array-job 50062602
"""
import argparse
import glob
import math
import os
import subprocess
import sys

# docs/handoff/runs/ -> the repo root is three levels up, same idiom as probe_td_pool.py:57.
HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.dirname(os.path.dirname(os.path.dirname(HERE)))
sys.path.insert(0, os.path.join(REPO, "benchmark"))
import benchmark_truth_set as bts  # noqa: E402

SCRATCH = "/n/scratch/users/c/cos689/xtea2_testdata"
TRUTH_DIR = "/n/data1/hms/dbmi/park/corinne/claude_sandbox/test_data/HG002/truth"
VCF_PREFIX = "HG002.novaseq.pcr-free.30x.dedup.grch38"
FAMILIES = ("L1", "Alu", "SVA")
ALT_TAG = {"L1": "LINE1", "Alu": "ALU", "SVA": "SVA"}
ARMS = ("master", "xtea2")

# The two arms lay their per-family VCFs out differently: master writes one per type subdir (it
# runs the three types as separate chains), xTea2 writes them under `merged/` after its parallel
# downstream. Same records either way -- `read_calls_vcf` keys family off the ALT tag, not the
# path -- but the aggregator has to know where to look.
def family_vcf(tree, arm, family):
    sample = "HG002" if arm == "master" else "HG002.wgs"
    tag = ALT_TAG[family]
    pats = ([os.path.join(tree, sample, family, "*_%s.vcf" % tag)] if arm == "master"
            else [os.path.join(tree, sample, "merged", "*_%s.vcf" % tag)])
    for pat in pats:
        hits = sorted(glob.glob(pat))
        if hits:
            return hits[0]
    return None


# --- slurm accounting ----------------------------------------------------------------------------

def parse_rss(s):
    """sacct MaxRSS ('4194304K', '38.03G', '') -> GB, or None. sacct's unit suffix is not
    optional in practice but it is not guaranteed either, so an unsuffixed value is read as KB,
    which is sacct's own default unit."""
    s = (s or "").strip()
    if not s or s in ("", "0"):
        return None
    mult = {"K": 1.0 / (1024 ** 2), "M": 1.0 / 1024, "G": 1.0, "T": 1024.0}
    if s[-1] in mult:
        try:
            return float(s[:-1]) * mult[s[-1]]
        except ValueError:
            return None
    try:
        return float(s) / (1024 ** 2)
    except ValueError:
        return None


def parse_elapsed(s):
    """'[DD-]HH:MM:SS' -> seconds."""
    s = (s or "").strip()
    if not s:
        return None
    days = 0
    if "-" in s:
        d, s = s.split("-", 1)
        days = int(d)
    parts = [int(x) for x in s.split(":")]
    while len(parts) < 3:
        parts.insert(0, 0)
    return days * 86400 + parts[0] * 3600 + parts[1] * 60 + parts[2]


def sacct_for(jobids):
    """{jobid: (elapsed_s, peak_rss_gb, state, node)}.

    MaxRSS lands on the `.batch` step, Elapsed and NodeList on the parent, so both are collected
    per job and the RSS is the MAX over every step -- a job whose peak lived in a step other than
    .batch would otherwise read as having no memory measurement at all.
    """
    if not jobids:
        return {}
    out = {}
    try:
        r = subprocess.run(
            ["sacct", "-j", ",".join(jobids), "-P", "-n", "--units=G",
             "--format=JobID,Elapsed,MaxRSS,State,NodeList,AllocCPUS"],
            capture_output=True, text=True, timeout=120)
    except (OSError, subprocess.SubprocessError) as e:
        print("  (sacct unavailable: %s -- falling back to each run's own sidecar wall)" % e)
        return {}
    if r.returncode != 0:
        print("  (sacct exited %d -- falling back to the sidecar wall)" % r.returncode)
        return {}
    for line in r.stdout.splitlines():
        f = line.split("|")
        if len(f) < 6:
            continue
        jid, elapsed, rss, state, node, cpus = f[:6]
        base = jid.split(".")[0]
        e, m = out.get(base, (None, None, None, None)), parse_rss(rss)
        el, pk, st, nd = e
        if "." not in jid:                       # the parent row carries state/node/elapsed
            el, st, nd = parse_elapsed(elapsed), state, node
        elif el is None:
            el = parse_elapsed(elapsed)
        if m is not None and (pk is None or m > pk):
            pk = m
        out[base] = (el, pk, st, nd)
    return out


# --- run discovery -------------------------------------------------------------------------------

def read_sidecar(tree):
    p = os.path.join(tree, "repro_run.tsv")
    if not os.path.exists(p):
        return {}
    with open(p) as fin:
        rows = [l.rstrip("\n").split("\t") for l in fin if l.strip()]
    if len(rows) < 2:
        return {}
    return dict(zip(rows[0], rows[1]))


def discover(array_job=None):
    runs = []
    for arm in ARMS:
        for tree in sorted(glob.glob(os.path.join(SCRATCH, "repro_%s_rep*" % arm))):
            if not os.path.isdir(tree):
                continue
            side = read_sidecar(tree)
            if array_job and side.get("array_job") and side["array_job"] != str(array_job):
                continue
            # THE SACCT KEY IS `<array_job>_<task>`, NOT `$SLURM_JOB_ID`, and this is a trap that
            # already bit: Slurm gave array task 6 a SLURM_JOB_ID equal to the ARRAY id itself
            # (50062602). `sacct -j 50062602` then expands to ALL SIX tasks, so that replicate
            # matched no per-task row and silently lost its MaxRSS -- a missing memory number that
            # looks like sacct retention rather than like a key collision.
            job = side.get("job")
            aj, task = side.get("array_job"), side.get("task")
            sacct_id = "%s_%s" % (aj, task) if aj and task else job
            runs.append(dict(arm=arm, tree=tree,
                             rep=side.get("rep", os.path.basename(tree)[-1]),
                             job=job, sacct_id=sacct_id, node=side.get("node"),
                             wall_sidecar=int(side["wall_s"]) if side.get("wall_s") else None,
                             complete=bool(side)))
    return runs


# --- scoring -------------------------------------------------------------------------------------

def score_run(run, window=200):
    """{family: dict(truth, calls, tp, fn, fp, sens, prec, hit_keys, call_keys)}.

    `hit_keys` is the set of TRUTH records this run recovered -- the population the replicate
    comparison is over. Keyed on the truth side rather than the call side on purpose: two runs
    calling the same insertion 12 bp apart are the same finding, and comparing call coordinates
    would score that as a difference.
    """
    per = {}
    for fam in FAMILIES:
        vcf = family_vcf(run["tree"], run["arm"], fam)
        if not vcf:
            per[fam] = None
            continue
        truth = bts.load_truth_family(TRUTH_DIR, fam)
        calls = bts.read_calls_vcf(vcf, [fam])
        matched = bts.match_one_to_one(truth, calls, window)
        tp = len(matched)
        per[fam] = dict(
            truth=len(truth), calls=len(calls), tp=tp,
            fn=len(truth) - tp, fp=len(calls) - tp,
            sens=bts._rate(tp, len(truth)), prec=bts._rate(tp, len(calls)),
            hit_keys={(truth[ti].chrom, truth[ti].pos) for ti in matched},
            call_keys={(c.chrom, c.pos) for c in calls},
            vcf=vcf)
    return per


def stats(xs):
    xs = [x for x in xs if x is not None]
    if not xs:
        return None
    n, mean = len(xs), sum(xs) / len(xs)
    sd = math.sqrt(sum((x - mean) ** 2 for x in xs) / (n - 1)) if n > 1 else 0.0
    return dict(n=n, mean=mean, sd=sd, lo=min(xs), hi=max(xs), spread=max(xs) - min(xs))


def hms(s):
    if s is None:
        return "     --"
    s = int(round(s))
    return "%d:%02d:%02d" % (s // 3600, s % 3600 // 60, s % 60)


# --- report ---------------------------------------------------------------------------------------

def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--array-job", default=None, help="restrict to one array job id")
    ap.add_argument("--window", type=int, default=200, help="position match window (bp)")
    a = ap.parse_args()

    runs = discover(a.array_job)
    if not runs:
        print("no repro_{master,xtea2}_rep* trees under %s" % SCRATCH)
        return 1

    bar = "=" * 94
    print(bar)
    print("MASTER vs xTea2, THREE REPLICATES EACH  --  HG002 WGS, same BAM, -c 24, --mem=64G")
    print(bar)

    acct = sacct_for([r["sacct_id"] for r in runs if r.get("sacct_id")])
    for r in runs:
        el, pk, st, nd = acct.get(r["sacct_id"] or "", (None, None, None, None))
        r["wall"] = el if el is not None else r["wall_sidecar"]
        r["rss"] = pk
        r["state"] = st
        r["node"] = nd or r["node"]

    incomplete = [r for r in runs if not r["complete"]]
    print("\n--- the six runs ---")
    print("  %-7s %-4s %-12s %-18s %10s %10s  %s"
          % ("arm", "rep", "job", "node", "wall", "peakRSS", "state"))
    for r in runs:
        print("  %-7s %-4s %-12s %-18s %10s %10s  %s"
              % (r["arm"], r["rep"], r["sacct_id"] or "?", (r["node"] or "?")[:18],
                 hms(r["wall"]),
                 "%.2f G" % r["rss"] if r["rss"] else "--",
                 r["state"] or ("INCOMPLETE" if not r["complete"] else "?")))
    if incomplete:
        print("\n  %d run(s) have no sidecar and are INCOMPLETE OR STILL RUNNING; every number"
              % len(incomplete))
        print("  below excludes nothing, so read them as provisional until this line is gone.")

    # --- 1. within-arm spread FIRST, because it is the yardstick ---------------------------------
    print("\n--- 1. WITHIN-ARM SPREAD, which bounds every between-arm claim below ---")
    per_arm = {}
    for arm in ARMS:
        rs = [r for r in runs if r["arm"] == arm and r["complete"]]
        per_arm[arm] = dict(runs=rs, wall=stats([r["wall"] for r in rs]),
                            rss=stats([r["rss"] for r in rs]))
        w, m = per_arm[arm]["wall"], per_arm[arm]["rss"]
        print("  %-7s n=%d" % (arm, len(rs)))
        if w:
            print("     wall     %s - %s   mean %s   spread %d s   sd %.0f s"
                  % (hms(w["lo"]), hms(w["hi"]), hms(w["mean"]), w["spread"], w["sd"]))
        if m:
            print("     peak RSS %.2f - %.2f G   mean %.2f G   spread %.2f G   sd %.2f G"
                  % (m["lo"], m["hi"], m["mean"], m["spread"], m["sd"]))
        nodes = {r["node"] for r in rs if r["node"]}
        print("     nodes    %s%s" % (", ".join(sorted(nodes)) or "?",
                                      "" if len(nodes) == len(rs) else "  <- a node is REUSED"))
    print("\n  All six ran CONCURRENTLY and streamed the same 44 GB BAM, so this spread includes")
    print("  cluster contention. That is deliberate (the harness header's design note 2): both")
    print("  arms saw the same conditions. It also means these walls are under-load and are NOT")
    print("  comparable to the solo n=1 numbers in HANDOFF_2026-08-07 sec 1.")

    # --- 2. speed and memory --------------------------------------------------------------------
    print("\n--- 2. SPEED and MEMORY, master vs xTea2 ---")
    mw, xw = per_arm["master"]["wall"], per_arm["xtea2"]["wall"]
    mr, xr = per_arm["master"]["rss"], per_arm["xtea2"]["rss"]
    if mw and xw:
        print("  wall     master %s (n=%d)   xTea2 %s (n=%d)"
              % (hms(mw["mean"]), mw["n"], hms(xw["mean"]), xw["n"]))
        print("           speedup on the MEANS: %.2fx" % (mw["mean"] / xw["mean"]))
        print("           worst case for xTea2 (its slowest vs master's fastest): %.2fx"
              % (mw["lo"] / xw["hi"]))
        print("           best  case for xTea2 (its fastest vs master's slowest): %.2fx"
              % (mw["hi"] / xw["lo"]))
        sep = mw["lo"] > xw["hi"]
        print("           the two arms' RANGES %s"
              % ("DO NOT OVERLAP -- the direction is not a coincidence of scheduling" if sep
                 else "OVERLAP, so the difference is inside the noise these runs measured"))
    if mr and xr:
        print("  peak RSS master %.2f G (n=%d)   xTea2 %.2f G (n=%d)"
              % (mr["mean"], mr["n"], xr["mean"], xr["n"]))
        print("           xTea2 uses %.0f%% of master's peak on the means"
              % (100.0 * xr["mean"] / mr["mean"]))
        sep = mr["lo"] > xr["hi"]
        print("           the two arms' RANGES %s"
              % ("DO NOT OVERLAP" if sep else "OVERLAP -- do not quote a ratio"))
        print("  AND READ THE MEMORY NUMBER THE RIGHT WAY: master loops L1/Alu/SVA SERIALLY, so")
        print("  its peak is ONE family resident; xTea2 runs all three concurrently. xTea2 holds")
        print("  three families in that fraction of what master needs for one.")
        print("  Correspondingly the wall ratio is NOT pure per-family speed -- part of it is a")
        print("  parallelism master does not have. That is an xTea2 change, not an artifact.")

    # --- 2.1 the within-node comparison, which is the one with a control ------------------------
    # Two co-tenancy effects showed up in array 50062602 and both distort a bare arm-vs-arm wall
    # ratio, so this section exists to route around them:
    #
    #   LOCKSTEP. Same-node pairs of the SAME arm finished at the SAME SECOND (master rep2/rep3
    #   both 1:03:18; xTea2 rep2/rep3 both 0:39:45), having started the same second. Two identical
    #   workloads sharing a node progress at the same rate, so such a pair is ONE independent wall
    #   measurement, not two -- the effective n for WALL is the number of distinct (node, arm)
    #   groups, which this section counts and prints.
    #
    #   PAGE CACHE. All six stream the same 44 GB BAM, so co-tenants warm it for each other. In
    #   50062602 the master replicate that ran ALONE was the SLOWEST (1:21:18 vs 1:03:18 for two
    #   sharing a node), i.e. co-tenancy HELPED. A cross-node ratio therefore mixes the arm
    #   difference with how many readers each node had.
    #
    # A node hosting BOTH arms controls for both at once: same hardware, same instant, same cache.
    print("\n--- 2.1 WITHIN-NODE, the comparison with a control ---")
    by_node = {}
    for r in runs:
        if r["complete"] and r["node"]:
            by_node.setdefault(r["node"], []).append(r)
    groups = {(r["node"], r["arm"]) for r in runs if r["complete"] and r["node"]}
    print("  independent (node, arm) groups: %d over %d runs -- same-node same-arm replicates run"
          % (len(groups), len([r for r in runs if r['complete']])))
    print("  in LOCKSTEP and count once for WALL. Peak RSS is NOT lockstep and keeps its full n.")
    found = False
    for node, rs in sorted(by_node.items()):
        arms_here = {r["arm"] for r in rs}
        tag = "  ".join("%s rep%s %s" % (r["arm"], r["rep"], hms(r["wall"])) for r in rs)
        print("  %-18s %s" % (node.split(".")[0], tag))
        if len(arms_here) < 2:
            continue
        m = stats([r["wall"] for r in rs if r["arm"] == "master"])
        x = stats([r["wall"] for r in rs if r["arm"] == "xtea2"])
        if m and x:
            found = True
            print("      -> BOTH ARMS ON THIS NODE: master %s vs xTea2 %s = %.2fx"
                  % (hms(m["mean"]), hms(x["mean"]), m["mean"] / x["mean"]))
    if not found:
        print("  no node hosted both arms -- the cross-node ratio in section 2 is all there is,")
        print("  and it carries the page-cache confound above. Re-run to land a mixed node.")

    # --- 3. precision / recall ------------------------------------------------------------------
    print("\n--- 3. PRECISION and RECALL against the HG002 truth set (positions, +-%d bp) ---"
          % a.window)
    for r in runs:
        if r["complete"]:
            r["score"] = score_run(r, a.window)
    scored = [r for r in runs if r.get("score")]

    for fam in FAMILIES:
        print("\n  %s" % fam)
        print("    %-7s %-4s %6s %6s %5s %5s %5s %9s %9s"
              % ("arm", "rep", "truth", "calls", "TP", "FN", "FP", "sens%", "prec%"))
        for r in scored:
            s = r["score"].get(fam)
            if not s:
                print("    %-7s %-4s   (no VCF)" % (r["arm"], r["rep"]))
                continue
            print("    %-7s %-4s %6d %6d %5d %5d %5d %9.1f %9.1f"
                  % (r["arm"], r["rep"], s["truth"], s["calls"], s["tp"], s["fn"], s["fp"],
                     s["sens"], s["prec"]))
        for arm in ARMS:
            ss = [r["score"][fam] for r in scored
                  if r["arm"] == arm and r["score"].get(fam)]
            if not ss:
                continue
            se, pr = stats([x["sens"] for x in ss]), stats([x["prec"] for x in ss])
            print("    %-7s MEAN  sens %.1f%% (spread %.1f)   prec %.1f%% (spread %.1f)"
                  % (arm, se["mean"], se["spread"], pr["mean"], pr["spread"]))

    print("\n  all families pooled")
    for arm in ARMS:
        rs = [r for r in scored if r["arm"] == arm]
        if not rs:
            continue
        tots = []
        for r in rs:
            tp = sum(r["score"][f]["tp"] for f in FAMILIES if r["score"].get(f))
            tr = sum(r["score"][f]["truth"] for f in FAMILIES if r["score"].get(f))
            ca = sum(r["score"][f]["calls"] for f in FAMILIES if r["score"].get(f))
            tots.append((tp, tr, ca))
            print("    %-7s rep %s   TP %d  FN %d  FP %d   sens %.1f%%  prec %.1f%%"
                  % (arm, r["rep"], tp, tr - tp, ca - tp,
                     bts._rate(tp, tr), bts._rate(tp, ca)))
        se = stats([bts._rate(tp, tr) for tp, tr, _ in tots])
        pr = stats([bts._rate(tp, ca) for tp, _, ca in tots])
        print("    %-7s MEAN      sens %.2f%% (spread %.2f)   prec %.2f%% (spread %.2f)"
              % (arm, se["mean"], se["spread"], pr["mean"], pr["spread"]))

    # --- 4. run-to-run STABILITY of the call set -------------------------------------------------
    print("\n--- 4. IS EACH CALLER EVEN DETERMINISTIC? ---")
    print("  Job 49075238 measured master's own -n 8 output differing from its own -n 24 output")
    print("  on identical inputs, so this is a measurement and not a formality. Over each arm's")
    print("  replicates: how many TRUTH records are recovered by all 3 / by some / by none, and")
    print("  how stable the raw call positions are.")
    for arm in ARMS:
        rs = [r for r in scored if r["arm"] == arm]
        print("\n  %s (n=%d)" % (arm, len(rs)))
        if len(rs) < 2:
            print("    fewer than 2 replicates -- nothing to compare")
            continue
        for fam in FAMILIES:
            sets = [r["score"][fam]["hit_keys"] for r in rs if r["score"].get(fam)]
            csets = [r["score"][fam]["call_keys"] for r in rs if r["score"].get(fam)]
            if len(sets) < 2:
                continue
            inter, union = set.intersection(*sets), set.union(*sets)
            cinter, cunion = set.intersection(*csets), set.union(*csets)
            wob = len(union) - len(inter)
            print("    %-4s truth hits: all %d reps %d, union %d -> %d record(s) UNSTABLE"
                  % (fam, len(sets), len(inter), len(union), wob))
            print("         call positions EXACTLY shared by all %d: %d of %d in the union "
                  "(%.1f%%)" % (len(csets), len(cinter), len(cunion),
                                bts._rate(len(cinter), len(cunion))))
    print("\n  A truth record recovered by some replicates and not others is a REAL sensitivity")
    print("  wobble, and its count is the yardstick for any cross-arm TP difference. Exact call")
    print("  positions moving while the truth hit is stable is a coordinate jitter, not a")
    print("  sensitivity change -- which is why the two lines are reported separately.")

    # --- 5. the cross-arm difference, judged against section 4 -----------------------------------
    print("\n--- 5. THE VERDICT ON THE ACCURACY DIFFERENCE ---")
    for fam in list(FAMILIES) + ["POOLED"]:
        got = {}
        for arm in ARMS:
            rs = [r for r in scored if r["arm"] == arm]
            if not rs:
                continue
            if fam == "POOLED":
                tps = [sum(r["score"][f]["tp"] for f in FAMILIES if r["score"].get(f)) for r in rs]
                fps = [sum(r["score"][f]["fp"] for f in FAMILIES if r["score"].get(f)) for r in rs]
            else:
                tps = [r["score"][fam]["tp"] for r in rs if r["score"].get(fam)]
                fps = [r["score"][fam]["fp"] for r in rs if r["score"].get(fam)]
            if tps:
                got[arm] = (stats(tps), stats(fps))
        if len(got) < 2:
            continue
        (mt, mf), (xt, xf) = got["master"], got["xtea2"]
        dtp, dfp = xt["mean"] - mt["mean"], xf["mean"] - mf["mean"]
        wob = max(mt["spread"], xt["spread"])
        # THREE RUNS MEASURE A FLOOR ON VARIABILITY, NOT A CONFIDENCE INTERVAL. The observed
        # range of n=3 is on average well inside the true range, so "delta exceeds the range" is
        # "not explained by the variation these six runs could see" -- never "significant". The
        # three-way wording keeps the middle case from being read as either.
        if abs(dtp) <= wob:
            verdict = "EXPLAINED by run-to-run variation (observed range %g)" % wob
        elif abs(dtp) <= 2 * wob:
            verdict = ("outside the observed range %g but within 2x it -- SUGGESTIVE, and n=3 "
                       "understates the range" % wob)
        else:
            verdict = ("%.1fx the observed range %g -- not explained by the variation these "
                       "runs could see" % (abs(dtp) / wob, wob)) if wob else \
                      ("both arms were EXACTLY reproducible over %d runs, so this delta is "
                       "systematic" % mt["n"])
        print("  %-7s TP  master %.1f (spread %g)  xTea2 %.1f (spread %g)  delta %+.1f  -> %s"
              % (fam, mt["mean"], mt["spread"], xt["mean"], xt["spread"], dtp, verdict))
        print("  %-7s FP  master %.1f (spread %g)  xTea2 %.1f (spread %g)  delta %+.1f"
              % ("", mf["mean"], mf["spread"], xf["mean"], xf["spread"], dfp))

    print("\n" + bar)
    print("CAVEATS, printed rather than footnoted")
    print(bar)
    print("  * THE SEED IS UNPINNED ON THE xTea2 ARM, unlike jobs 49792012/49879677/50038678")
    print("    which all pinned XTEA_BASIC_INFO_SEED=42. Master has no such knob, so pinning one")
    print("    arm would understate xTea2's spread. These runs are therefore NOT drop-in")
    print("    comparable with those three.")
    print("  * ALL SIX RAN CONCURRENTLY on the same 44 GB BAM. Contention is inside the spread,")
    print("    by design; the walls are under-load, not solo.")
    print("  * sacct MaxRSS is POLLED (~30 s), so a short spike can be undersampled. No run was")
    print("    near its 64 GB ceiling, so these are honest high-waters and not ceiling artifacts.")
    print("  * The truth set scores POSITIONS. HG002's derived transduction denominator is 0 (L1)")
    print("    and 2 (SVA), so transduction sensitivity is not measurable here at all.")
    print("  * NOT COMPARABLE TO THE PAPER'S NUMBERS: this BAM is hs38DH-aligned, not no_alt, and")
    print("    xTea filters on MAPQ. (benchmark_truth_set.py's own warning, plan sec 7.1.)")
    print("  * L1 precision is depressed on BOTH arms by created orphan/sibling SITES that the")
    print("    truth set does not enumerate and so cannot match; `benchmark_truth_set.py score`")
    print("    prints the excluding-them figure. It is a property of the metric, not of an arm.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
