#!/usr/bin/env python3
"""Accumulate per-phase wall-clock timings across EVERY mosaic attempt, for cost estimation.

    python3 docs/handoff/runs/mosaic_timing_ledger.py            # print the ledger
    python3 docs/handoff/runs/mosaic_timing_ledger.py --tsv out.tsv
    python3 docs/handoff/runs/mosaic_timing_ledger.py --estimate 5982   # cost a 5982 GB arm

WHY IT READS THE .out LOGS AND NOT step_timings.tsv. run_xtea.init_timing opens
step_timings.tsv with mode "w", so every rerun of an arm TRUNCATES the timings of the
attempt before it -- and with --resume in play an arm can be run several times. The
per-job SLURM stdout (mosaic_wgs_<arm>_<jobid>.out) is the durable record: one file per
attempt, never rewritten, and it carries the same `[xtea2 + T] phase NAME: Ns` lines.

WHAT COUNTS AS A COMPLETED PHASE, and why the job's exit State is the wrong test. Before
2026-08-13 `phase()` logged its duration from a bare `finally`, so a phase that DIED printed
the same line as one that succeeded -- uwsc (50447196) logged "phase scan_pass: 25665.4s" on
the way out of an OOM having never written scan_clip_realign.tsv. But keying on the job's
State instead would throw away every good measurement we have: bcm's scan_pass is CANCELLED
at the job level and unquestionably finished, because cns_emit ran for two hours afterwards.

The rule used here is the one that actually follows from how the driver exits: a phase is
CONFIRMED when ANY later event (another phase, or any step) is logged after it. The driver
dies at the first failure, so anything logged afterwards proves the phase before it returned.
Only the LAST event in an attempt is unprovable from the log -- and that is exactly the uwsc
case -- so only that one is marked SUSPECT and held out of the throughput numbers.
"""
import argparse
import glob
import os
import re
import subprocess
import sys
from collections import defaultdict

OUT_BASE = "/n/scratch/users/c/cos689/xtea2_testdata"
# this file lives at <repo>/docs/handoff/runs/, so the repo root is FOUR levels up
REPO = os.path.dirname(os.path.dirname(os.path.dirname(
    os.path.dirname(os.path.abspath(__file__)))))
MANIFEST = os.path.join(REPO, "hapmapmix_mosaictruthtest.tsv")
ARMS = ["broad", "nygc", "uwsc", "bcm", "washu"]     # shallow -> deep

PHASE_OK = re.compile(r"^\[xtea2 \+\s*([\d.]+)s\]\s*phase (\S+): ([\d.]+)s"
                      r"(?:\s*\(children peak RSS ([\d.]+) GB\))?")
ANY_EVENT = re.compile(r"^\[xtea2 \+\s*([\d.]+)s\]\s*(?:phase |step |=== )")
PHASE_FAIL = re.compile(r"^\[xtea2 \+\s*([\d.]+)s\]\s*phase (\S+): FAILED after ([\d.]+)s")
STEP_OK = re.compile(r"^\[xtea2 \+\s*([\d.]+)s\]\s*step (\S+): ([\d.]+)s")
RESUMED = re.compile(r"^\[resume\].*?(scan_pass|cns_emit).*?(reusing|ADOPTED)")
COV = re.compile(r"cov\(raw\)=([\d.]+)")


def sh(cmd):
    return subprocess.run(cmd, shell=True, capture_output=True, text=True).stdout


def job_states():
    """jobid -> State, for every mosaic job sacct still remembers."""
    out = sh("sacct --format=JobID,JobName%12,State -P --noheader -S 2026-08-01 2>/dev/null")
    st = {}
    for line in out.splitlines():
        f = line.split("|")
        if len(f) < 3 or "." in f[0]:
            continue
        if f[1].strip() in ARMS:
            st[f[0].strip()] = f[2].strip()
    return st


def arm_meta():
    """arm -> (n_bams, input_gb) from the manifest."""
    meta = {}
    try:
        with open(MANIFEST) as fin:
            next(fin)
            for line in fin:
                c = line.rstrip("\n").split("\t")
                if len(c) < 6:
                    continue
                n, gb = meta.get(c[1], (0, 0.0))
                meta[c[1]] = (n + 1, gb + float(c[5]))
    except OSError:
        pass
    return meta


def harvest():
    """One record per (arm, jobid) attempt."""
    states = job_states()
    attempts = []
    for path in sorted(glob.glob(f"{OUT_BASE}/mosaic_wgs_*_*.out")):
        base = os.path.basename(path)[len("mosaic_wgs_"):-len(".out")]
        arm, _, jobid = base.rpartition("_")
        if arm not in ARMS:
            continue
        rec = {"arm": arm, "jobid": jobid, "state": states.get(jobid, "?"),
               "phases": {}, "at": {}, "failed": {}, "resumed": set(), "cov": None,
               "steps": {}, "last_event_t": 0.0, "mtime": os.path.getmtime(path)}
        try:
            with open(path, errors="replace") as fin:
                for line in fin:
                    line = line.rstrip("\n")
                    m = ANY_EVENT.match(line)
                    if m:
                        rec["last_event_t"] = max(rec["last_event_t"], float(m.group(1)))
                    m = PHASE_OK.match(line)
                    if m:
                        rec["phases"][m.group(2)] = (float(m.group(3)),
                                                     float(m.group(4)) if m.group(4) else None)
                        rec["at"][m.group(2)] = float(m.group(1))
                        continue
                    m = STEP_OK.match(line)
                    if m:
                        rec["steps"][m.group(2)] = float(m.group(3))
                        continue
                    m = PHASE_FAIL.match(line)
                    if m:
                        rec["failed"][m.group(2)] = float(m.group(3))
                        continue
                    m = RESUMED.match(line.strip())
                    if m:
                        rec["resumed"].add(m.group(1))
                        continue
                    if rec["cov"] is None:
                        m = COV.search(line)
                        if m:
                            rec["cov"] = float(m.group(1))
        except OSError:
            continue
        if rec["phases"] or rec["failed"] or rec["resumed"]:
            attempts.append(rec)
    return attempts


def phase_confirmed(rec, name):
    """A phase is proven complete when something was logged AFTER it: the driver exits at the
    first failure, so a later event can only exist if this phase returned. The last event in
    an attempt is the one that cannot be proven this way -- the uwsc trap."""
    t = rec["at"].get(name)
    if t is None:
        return False
    return rec["last_event_t"] > t + 1e-9


def fmt(s):
    if s is None:
        return "-"
    h = s / 3600.0
    return f"{h:.2f}h" if h >= 1 else f"{s/60:.1f}m"


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("--tsv", help="also write the raw ledger here")
    ap.add_argument("--estimate", type=float,
                    help="estimate scan_pass+cns_emit for an arm of this input size (GB)")
    a = ap.parse_args(argv)

    meta = arm_meta()
    attempts = harvest()
    attempts.sort(key=lambda r: (ARMS.index(r["arm"]), r["mtime"]))

    rows = ["arm\tjobid\tstate\tphase\twall_s\tpeak_rss_gb\tstatus\tn_bams\tinput_gb\tcov"]
    print("=" * 100)
    print("MOSAIC PHASE TIMINGS -- every attempt, newest last per arm")
    print("=" * 100)
    print(f"{'arm':6s}{'job':>10s}{'state':>14s}  {'scan_pass':>11s}{'cns_emit':>11s}"
          f"{'downstream':>12s}{'shared':>9s}   note")
    for r in attempts:
        nb, gb = meta.get(r["arm"], (0, 0.0))
        note = []
        if r["resumed"]:
            note.append("adopted:" + ",".join(sorted(r["resumed"])))
        if r["failed"]:
            note.append("FAILED:" + ",".join(sorted(r["failed"])))
        unproven = [p for p in r["phases"] if not phase_confirmed(r, p)]
        if unproven:
            note.append(f"SUSPECT: {','.join(sorted(unproven))} ({r['state'][:20]})")
        def g(n):
            return r["phases"].get(n, (None, None))[0]
        print(f"{r['arm']:6s}{r['jobid']:>10s}{r['state'][:13]:>14s}  "
              f"{fmt(g('scan_pass')):>11s}{fmt(g('cns_emit')):>11s}"
              f"{fmt(g('downstream_serial')):>12s}{fmt(g('shared_tail')):>9s}   "
              f"{'; '.join(note)}")
        for ph, (sec, rss) in sorted(r["phases"].items()):
            status = "confirmed" if phase_confirmed(r, ph) else "SUSPECT"
            rows.append(f"{r['arm']}\t{r['jobid']}\t{r['state']}\t{ph}\t{sec:.1f}\t"
                        f"{rss if rss is not None else ''}\t{status}\t{nb}\t{gb:.0f}\t"
                        f"{r['cov'] if r['cov'] else ''}")
        for ph, sec in sorted(r["failed"].items()):
            rows.append(f"{r['arm']}\t{r['jobid']}\t{r['state']}\t{ph}\t{sec:.1f}\t\tFAILED\t"
                        f"{nb}\t{gb:.0f}\t{r['cov'] if r['cov'] else ''}")

    # throughput, from COMPLETED attempts only
    print("\n" + "=" * 100)
    print("THROUGHPUT (completed scans only -- the basis for any extrapolation)")
    print("=" * 100)
    good = []
    for r in attempts:
        if "scan_pass" not in r["phases"] or r["resumed"]:
            continue
        if not phase_confirmed(r, "scan_pass"):
            continue
        nb, gb = meta.get(r["arm"], (0, 0.0))
        good.append((r["arm"], nb, gb, r["cov"], r["phases"]["scan_pass"][0],
                     r["phases"].get("cns_emit", (None,))[0]))
    if good:
        print(f"  {'arm':6s}{'bams':>5s}{'input_GB':>10s}{'cov':>8s}{'scan':>9s}"
              f"{'s/GB':>8s}{'cns_emit':>10s}")
        for arm, nb, gb, cov, scan, emit in good:
            print(f"  {arm:6s}{nb:>5d}{gb:>10.0f}{(f'{cov:.0f}x' if cov else '-'):>8s}"
                  f"{fmt(scan):>9s}{(scan/gb if gb else float('nan')):>8.2f}{fmt(emit):>10s}")
        sgb = [s / g for _, _, g, _, s, _ in good if g]
        if not sgb:
            print("\n  (no per-GB rate: manifest sizes unavailable)")
            sgb = None
        mean = (sum(sgb) / len(sgb)) if sgb else None
        if sgb:
            print(f"\n  scan_pass: mean {mean:.2f} s/GB of input "
                  f"(range {min(sgb):.2f}-{max(sgb):.2f}, n={len(sgb)})")
        if a.estimate and sgb:
            print(f"  -> a {a.estimate:.0f} GB arm: scan ~{fmt(mean*a.estimate)} "
                  f"(range {fmt(min(sgb)*a.estimate)}-{fmt(max(sgb)*a.estimate)})")
    else:
        print("  (no COMPLETED scan yet to base an estimate on)")

    if a.tsv:
        with open(a.tsv, "w") as fout:
            fout.write("\n".join(rows) + "\n")
        print(f"\nwrote {a.tsv} ({len(rows)-1} rows)")
    return 0


if __name__ == "__main__":
    sys.exit(main())
