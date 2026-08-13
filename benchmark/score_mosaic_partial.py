#!/usr/bin/env python3
"""Score whichever per-type mosaic call sets have landed, without waiting for the merged VCF.

    python3 benchmark/score_mosaic_partial.py --arm broad
    python3 benchmark/score_mosaic_partial.py --out-dir /path/to/sample_root

WHY THIS EXISTS. benchmark_hapmap_mosaic.py scores the merged VCF, which the shared tail
writes only after ALL THREE families finish. On these runs SVA is the long pole -- broad's
SVA chain is ~2h behind its L1 chain -- so waiting for the VCF means waiting hours past the
point where L1 could already have told you whether a change worked. This reads each family's
`candidate_disc_filtered_cns_post_filtering_mosaic.txt` as soon as it appears and scores it
with the SAME load_truth/score functions the real scorer uses, so the numbers are directly
comparable to the final report rather than a second opinion.

POSITION RESOLUTION IS THE ONE THING NOT SHARED, and it has to be replicated exactly.
The per-type file is not a VCF; its POS is derived the way x_gvcf.py:365-368 derives it:

    i_begin = min(lclip_pos, rclip_pos); if i_begin <= 0: i_begin = the other one

That `<= 0` fallback is not cosmetic. lclip_pos is -1 (a "missing" sentinel) on a real
fraction of calls -- 45 of 216 on one broad Alu run -- and keying on lclip_pos alone sends
every one of them to position -1, where they match nothing. That silently cost 15 true
positives and inflated the FP count by the same 15 before it was noticed.
"""
import argparse
import datetime
import os
import subprocess
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from benchmark_truth_set import CallRecord                      # noqa: E402
from benchmark_hapmap_mosaic import (load_truth, load_exclusions,  # noqa: E402
                                     score, DEFAULT_TRUTH_DIR)

OUT_BASE = "/n/scratch/users/c/cos689/xtea2_testdata"
FAMILIES = ("L1", "Alu", "SVA")
MOSAIC_NAME = "candidate_disc_filtered_cns_post_filtering_mosaic.txt"


def resolve_pos(lclip_pos, rclip_pos):
    """x_gvcf.py's own POS rule -- see the module docstring for why the fallback matters."""
    i_begin, i_end = (lclip_pos, rclip_pos) if lclip_pos <= rclip_pos else (rclip_pos, lclip_pos)
    return i_end if i_begin <= 0 else i_begin


def _job_start_epoch(arm):
    """Start time of `arm`'s currently-running job, or None.

    THE POINT OF THIS. A per-type file existing does not mean it belongs to the run you are
    watching -- the previous attempt's output sits at the same path until the current one
    overwrites it, and a partial score read off those is a number from a different pipeline
    version entirely. That is the same class of mistake as trusting a "phase scan_pass:
    25665.4s" line from a run that died: the artifact is there, it just is not yours.
    """
    try:
        out = subprocess.run(["squeue", "-h", "-u", os.environ.get("USER", ""),
                              "-o", "%j|%S"], capture_output=True, text=True, timeout=15).stdout
    except Exception:
        return None
    for line in out.splitlines():
        name, _, start = line.partition("|")
        if name.strip() != arm:
            continue
        try:
            return datetime.datetime.strptime(start.strip(), "%Y-%m-%dT%H:%M:%S").timestamp()
        except ValueError:
            return None
    return None


def load_calls(path, family):
    out = []
    with open(path) as fin:
        for line in fin:
            f = line.rstrip("\n").split("\t")
            if len(f) < 4:
                continue
            try:
                pos = resolve_pos(int(f[2]), int(f[3]))
            except ValueError:
                continue
            out.append(CallRecord(f[0], pos, family, None,
                                  f[32] if len(f) > 32 else None, {}))
    return out


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("--arm", help="HapMap arm name (bcm/broad/nygc/uwsc/washu)")
    ap.add_argument("--out-dir", help="sample root directly (overrides --arm)")
    ap.add_argument("--truth-dir", default=DEFAULT_TRUTH_DIR)
    ap.add_argument("--window", type=int, default=200)
    ap.add_argument("--no-freshness", action="store_true",
                    help="skip the stale-file check (see _job_start_epoch)")
    a = ap.parse_args(argv)

    if a.out_dir:
        root = a.out_dir
    elif a.arm:
        root = f"{OUT_BASE}/mosaic_hapmap_wgs_{a.arm}/hapmap_{a.arm}_wgs"
    else:
        ap.error("need --arm or --out-dir")

    job_start = None if a.no_freshness or not a.arm else _job_start_epoch(a.arm)

    hg005, masked = load_exclusions(a.truth_dir)
    print("=" * 78)
    print(f"PARTIAL mosaic score -- per-type files under {root}")
    print(f"match window: +/-{a.window} bp   (families with no file yet are skipped)")
    print("=" * 78)

    tot_tp = tot_truth = tot_fp = 0
    for fam in FAMILIES:
        path = os.path.join(root, fam, MOSAIC_NAME)
        if not os.path.isfile(path):
            print(f"\n--- {fam} --- not finished yet ({MOSAIC_NAME} absent)")
            continue
        calls = load_calls(path, fam)
        truth = load_truth(a.truth_dir, fam, tier=1)
        bins, o, unm = score(calls, truth, hg005, masked, a.window)
        tot_tp += o["tp"]; tot_truth += o["truth"]; tot_fp += o["fp"]
        mtime = os.path.getmtime(path)
        stamp = datetime.datetime.fromtimestamp(mtime).strftime("%m-%d %H:%M")
        stale = job_start is not None and mtime < job_start
        print(f"\n--- {fam} ---" + (f"   [written {stamp}]" if not stale else
              f"   *** STALE: written {stamp}, BEFORE this arm's running job started -- "
              f"these numbers are from a PREVIOUS attempt ***"))
        print(f"  truth {o['truth']}  calls {o['calls']}  TP {o['tp']}  FN {o['fn']}  FP {o['fp']}")
        print(f"  sensitivity {o['sens']:.1f}%   precision {o['prec']:.1f}%")
        print(f"  unmatched: {unm['germline']} germline (HG005), "
              f"{unm['ambiguous']} ambiguous (masked), {unm['fp']} counted as FP")
        if not truth:
            continue
        print(f"  {'VAF':<12}{'truth':>7}{'TP':>7}{'sens':>9}")
        for b in bins:
            if b["n"]:
                print(f"  {b['label']:<12}{b['n']:>7}{b['tp']:>7}{b['sens']:>8.1f}%")
    if tot_truth:
        print(f"\n{'ALL SCORED FAMILIES':<24}TP {tot_tp}  FP {tot_fp}  "
              f"of {tot_truth} truth  ->  sens {100.0*tot_tp/tot_truth:.1f}%"
              + (f", prec {100.0*tot_tp/(tot_tp+tot_fp):.1f}%" if (tot_tp + tot_fp) else ""))
    return 0


if __name__ == "__main__":
    sys.exit(main())
