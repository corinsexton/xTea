#!/usr/bin/env python3
"""Truth-set recall at every stage of the mosaic pipeline, per arm and family.

    python3 benchmark/mosaic_funnel.py --arm broad
    python3 benchmark/mosaic_funnel.py --arm bcm --families L1
    python3 benchmark/mosaic_funnel.py --arm broad bcm washu --stage-only clip

WHY A FUNNEL AND NOT A FINAL NUMBER. A single sensitivity figure says a run is bad without
saying where it went bad, and on this pipeline the answer has been different every time: the
mosaic AF gate costs ~15%, postF_all's germline logic costs ~50%, and the scan's own detection
costs ~35% before any filter runs. Only the first of those is a threshold you can tune, so
knowing which stage owns the loss decides whether tuning is even the right lever.

READS FRESHNESS, BECAUSE THE PATHS ARE REUSED. A stage file existing does not mean the run you
are looking at produced it -- the previous attempt's output sits at the same path until
overwritten, and stages BELOW the currently-running step are routinely last-run's. Every row is
labelled FRESH or stale against the arm's running job, and the "lost here" column is only
computed between consecutive FRESH rows, since differencing a fresh row against a stale one
above it is meaningless.
"""
import argparse
import bisect
import datetime
import os
import subprocess
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from benchmark_hapmap_mosaic import load_truth, DEFAULT_TRUTH_DIR  # noqa: E402

OUT_BASE = "/n/scratch/users/c/cos689/xtea2_testdata"
STAGES = [("clip candidates", "candidate_list_from_clip.txt"),
          ("disc candidates", "candidate_list_from_disc.txt"),
          ("cns_filter", "candidate_disc_filtered_cns.txt"),
          ("+transduction (cns2)", "candidate_disc_filtered_cns2.txt"),
          ("postF_all", "candidate_disc_filtered_cns_post_filtering.txt"),
          ("postF_mosaic (FINAL)", "candidate_disc_filtered_cns_post_filtering_mosaic.txt")]


def resolve_pos(l, r):
    b, e = (l, r) if l <= r else (r, l)
    return e if b <= 0 else b


def job_start_epoch(arm):
    try:
        out = subprocess.run(["squeue", "-h", "-u", os.environ.get("USER", ""),
                              "-o", "%j|%S"], capture_output=True, text=True, timeout=15).stdout
    except Exception:
        return None
    for line in out.splitlines():
        name, _, start = line.partition("|")
        if name.strip() == arm:
            try:
                return datetime.datetime.strptime(start.strip(), "%Y-%m-%dT%H:%M:%S").timestamp()
            except ValueError:
                return None
    return None


def index_positions(path):
    idx = {}
    n = 0
    with open(path) as fin:
        for line in fin:
            f = line.rstrip("\n").split("\t")
            try:
                pos = int(f[1]) if len(f) <= 7 else resolve_pos(int(f[2]), int(f[3]))
            except (ValueError, IndexError):
                continue
            idx.setdefault(f[0], []).append(pos)
            n += 1
    for c in idx:
        idx[c].sort()
    return idx, n


def recall(idx, truth, window):
    hit = 0
    for t in truth:
        arr = idx.get(t.chrom)
        if not arr:
            continue
        i = bisect.bisect_left(arr, t.pos - window)
        if i < len(arr) and arr[i] <= t.pos + window:
            hit += 1
    return hit


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("--arm", nargs="+", required=True)
    ap.add_argument("--families", nargs="+", default=["L1", "Alu", "SVA"])
    ap.add_argument("--truth-dir", default=DEFAULT_TRUTH_DIR)
    ap.add_argument("--window", type=int, default=200)
    ap.add_argument("--stage-only", help="print just this stage, one row per arm/family "
                                         "(substring match on the stage label)")
    a = ap.parse_args(argv)

    rows = []
    for arm in a.arm:
        root = f"{OUT_BASE}/mosaic_hapmap_wgs_{arm}/hapmap_{arm}_wgs"
        js = job_start_epoch(arm)
        for fam in a.families:
            truth = load_truth(a.truth_dir, fam, tier=1)
            prev = None
            for label, fn in STAGES:
                p = os.path.join(root, fam, fn)
                if not os.path.exists(p):
                    rows.append((arm, fam, label, None, None, len(truth), "", "not written yet"))
                    continue
                mt = os.path.getmtime(p)
                fresh = js is None or mt > js
                idx, n = index_positions(p)
                hit = recall(idx, truth, a.window)
                lost = f"-{prev - hit}" if (fresh and prev is not None) else ""
                state = ("FRESH" if fresh else "stale") + \
                        f" [{datetime.datetime.fromtimestamp(mt):%H:%M}]"
                rows.append((arm, fam, label, n, hit, len(truth), lost, state))
                if fresh:
                    prev = hit

    if a.stage_only:
        print(f"{'arm':7s}{'fam':5s}{'calls':>10}{'truth found':>13}{'recall':>9}   stage/state")
        for arm, fam, label, n, hit, tot, lost, state in rows:
            if a.stage_only.lower() not in label.lower() or n is None:
                continue
            print(f"{arm:7s}{fam:5s}{n:>10}{hit:>7}/{tot:<5}{100.0*hit/tot:>8.1f}%   {label} | {state}")
        return 0

    for arm in a.arm:
        for fam in a.families:
            sel = [r for r in rows if r[0] == arm and r[1] == fam]
            if not sel:
                continue
            print(f"\n=== {arm} / {fam} -- truth {sel[0][5]} tier1 somatic loci ===")
            print(f"  {'stage':24s}{'calls':>10}{'truth':>8}{'recall':>9}{'lost':>7}   state")
            for _, _, label, n, hit, tot, lost, state in sel:
                if n is None:
                    print(f"  {label:24s}{'-':>10}{'-':>8}{'-':>9}{'':>7}   {state}")
                else:
                    print(f"  {label:24s}{n:>10}{hit:>8}{100.0*hit/tot:>8.1f}%{lost:>7}   {state}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
