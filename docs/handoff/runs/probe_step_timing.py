#!/usr/bin/env python3
"""Where did a two-arm wall-clock delta actually go, and can the flag under test have caused it.

Offline, seconds, no BAM, no job. Reads two finished runs' own logs --

    python docs/handoff/runs/probe_step_timing.py --off-arm <run>/<sample> --on-arm <run>/<sample>

WHAT THIS IS FOR, AND THE TRAP IT EXISTS TO CLOSE. Every pipeline gate in this project reports a
`downstream_parallel` wall time against a band from prior runs, and every one of them carries the
same disclaimer: *"one arm on one node has no control."* That disclaimer is true of the PHASE
number and it is what made the decoy-concatenation gate's `603.6s` against a `462-512s` band
un-actionable -- handoff 2026-08-05 §5 filed it as "one rerun settles it", i.e. as 24 cores and an
hour of priopark.

**It does not need a rerun, because the control is already inside the job.** Three facts about
how these runs are logged turn one arm into a controlled experiment:

  1. **`--parallel-downstream` runs the three families concurrently in ONE process tree on ONE
     node.** So L1, Alu and SVA share a machine, a filesystem and a minute-by-minute load.
  2. **A family whose bwa reference is byte-identical between the arms is doing identical work.**
     When a flag changes only SVA's transduction reference, L1 is a same-node, same-job,
     same-phase NULL ARM. If L1 moved too, the machine moved; if L1 is flat and SVA is not, the
     change is convicted. This is the whole instrument.
  3. **`x_transduction.td_phase` already writes a 13-line sub-decomposition of the transduction
     step to each family's `.steps.log`** (`[TD-PHASE]<TAB>pid<TAB>rep_type<TAB>label<TAB>secs`).
     So "the step got slower" resolves to WHICH of flank_build / flank_bwa_index / native_gather /
     clip_bwa_realign / clip_sam_parse / disc_bwa_realign / ... got slower, with no instrumenting.

This is the third time this project's intermediates have outlived a handoff's belief that a
measurement needed a run (plan §10: `probe_nearby_clip` found the disc-to-consensus SAM, gate [A]
found `temp_transduction_clip.sam`). **Check the logs before scoping a rerun.**

MEASURED on the decoy-concatenation pair -- job `49792012` (OFF) vs `49879677` (ON):

  * **`downstream_parallel` IS SVA's serial chain.** OFF 63.0+354.3+51.0+2.2+0.9 = 471.4s against
    a reported 471.5s; ON 66.8+433.8+58.7+30.0+14.2 = 603.5s against 603.6s. So the phase is a
    critical path, not an average, and only SVA's steps can move it.
  * **40% of the +132.1s delta is in steps the flag cannot reach** -- `SVA:postF_all` 2.2 -> 30.0s
    and `SVA:postF_hc` 0.9 -> 14.2s, which are downstream of the call set (it differs by ONE call
    between these arms) and which the timeline below shows **ran with nothing else overlapping in
    BOTH arms**. A 13x and 16x slowdown on an uncontended flag-independent step is the machine.
  * **L1 -- the null arm -- got FASTER**: `TOTAL_call_candidate_transduction_v3` 169.95 -> 165.60s
    (-4.35s), `clip_bwa_realign` 75.97 -> 72.49s. Same job, same node, same phase, byte-identical
    reference. **So the node is exonerated and the remaining SVA cost is real.**
  * **The real cost is +55.4s and it is the two bwa realignment passes**: SVA `clip_bwa_realign`
    125.95 -> 164.15s (+38.20) and `disc_bwa_realign` 101.33 -> 118.51s (+17.18), out of a
    +57.40s move in the whole timed region.

AND THE MECHANISM IS NOT REFERENCE SIZE, WHICH IS WHY THE GATE'S OWN GUESS MISSED IT. The gate
reasoned "45 records out of 5,778 cannot plausibly account for it" and that is correct as stated:
SVA's reference grows 17,331,894 -> 17,367,290 bases, **+0.20%**, and `flank_bwa_index` actually
got *faster* (6.58 -> 5.73s). The 45 records are REPEAT CONSENSUS sequences -- full-length L1,
HERV-K, six SVA and 35 Alu -- so they are seed magnets: a read that used to fail fast now finds
many seed hits to extend. The output volume confirms it is work-per-read and not more reads:
`temp_transduction_clip.sam` 1,979,331,242 -> 1,978,591,405 bytes (**-0.04%**) and
`temp_transduction_disc.sam` 772,545,001 -> 772,187,258 (**-0.05%**). **Same reads, same output,
+30% alignment time.**

So the honest headline is **~55s (~12% of the phase), not 132s (28%)** -- and note that the gate's
stated prediction ("if anything expect SVA's disc parse to get slightly cheaper") was RIGHT and
was swamped: `clip_sam_parse` -2.41s and `disc_sam_parse` -0.78s, real savings, 20x smaller than
the alignment cost nobody predicted.

WHAT THIS PROBE DOES NOT DO. It cannot separate machine drift from real cost when EVERY family's
reference changed -- there is no null arm then, and a rerun is the only answer. It reports the
absence of a control rather than guessing past it.
"""
import argparse
import collections
import os
import re
import sys

STEP_END = re.compile(
    r"\[xtea2 \+ *([0-9.]+)s\] step ([A-Za-z0-9_]+):([A-Za-z0-9_]+): ([0-9.]+)s")
PHASE_END = re.compile(r"\[xtea2 \+ *([0-9.]+)s\] phase ([A-Za-z0-9_]+): ([0-9.]+)s")
TD_PHASE = re.compile(r"\[TD-PHASE\]\t(\d+)\t(\S+)\t(\S+)\t([0-9.]+)")

# The transduction reference each family's bwa passes are built against. Byte-equality of this
# file between the arms is what makes a family a null arm -- see the module docstring. It is
# deliberately compared as BYTES and not as a record count: the concatenation gate learned that
# "same number of records" and "same reference" are different claims.
REF_REL = "{fam}/tmp/transduction/all_with_polymerphic_flanks.fa"


def _driver_log(run_dir):
    """The driver log sits beside the sample dir, not inside it (a gate got this wrong once)."""
    for cand in (os.path.join(os.path.dirname(run_dir.rstrip("/")), "driver.log"),
                 os.path.join(run_dir, "driver.log")):
        if os.path.isfile(cand):
            return cand
    return None


def read_steps(run_dir):
    """(start, end, family, step, dur) for every per-family step, phase-relative."""
    log = _driver_log(run_dir)
    if log is None:
        return [], {}
    ev, phases = [], {}
    with open(log, errors="replace") as fh:
        for line in fh:
            m = STEP_END.search(line)
            if m:
                end, fam, step, dur = (float(m.group(1)), m.group(2),
                                       m.group(3), float(m.group(4)))
                if fam != "shared":
                    ev.append((end - dur, end, fam, step, dur))
                continue
            m = PHASE_END.search(line)
            if m:
                phases[m.group(2)] = float(m.group(3))
    return sorted(ev), phases


def read_td_phase(run_dir, sample, fam):
    """The transduction step's own 13-line sub-decomposition, in emission order."""
    path = os.path.join(run_dir, "{0}_{1}.steps.log".format(sample, fam))
    out = collections.OrderedDict()
    if not os.path.isfile(path):
        return out
    with open(path, errors="replace") as fh:
        for line in fh:
            m = TD_PHASE.search(line)
            if m:
                # Later emissions win: a resumed step would append, and the last run is the
                # one whose wall time the driver reported.
                out[m.group(3)] = float(m.group(4))
    return out


def ref_state(off_dir, on_dir, fam):
    """'null arm' | 'CHANGED' | 'unknown' -- byte equality of the family's bwa reference."""
    a = os.path.join(off_dir, REF_REL.format(fam=fam))
    b = os.path.join(on_dir, REF_REL.format(fam=fam))
    if not (os.path.isfile(a) and os.path.isfile(b)):
        return "unknown", None, None
    sa, sb = os.path.getsize(a), os.path.getsize(b)
    if sa != sb:
        return "CHANGED", sa, sb
    with open(a, "rb") as fa, open(b, "rb") as fb:
        same = fa.read() == fb.read()
    return ("null arm" if same else "CHANGED"), sa, sb


def critical_path(ev, phase_secs):
    """Which family's serial chain equals the phase wall time.

    Under --parallel-downstream the phase is a MAX over families, not a sum over steps, so only
    the critical family's steps can move it. Naming that family is what stops a reader summing
    every step's delta and over-attributing the phase.
    """
    per_fam = collections.defaultdict(float)
    for _s, _e, fam, _st, dur in ev:
        per_fam[fam] += dur
    if not per_fam:
        return None, {}
    fam = max(per_fam, key=lambda f: per_fam[f])
    return fam, per_fam


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--off-arm", required=True, help="<out>/<sample> of the baseline arm")
    ap.add_argument("--on-arm", required=True, help="<out>/<sample> of the arm under test")
    ap.add_argument("--sample", default=None,
                    help="sample name for the .steps.log files (default: basename of --on-arm)")
    ap.add_argument("--families", default="L1,Alu,SVA")
    args = ap.parse_args()

    off, on = args.off_arm.rstrip("/"), args.on_arm.rstrip("/")
    sample = args.sample or os.path.basename(on)
    fams = [f for f in args.families.split(",") if f]

    for d in (off, on):
        if not os.path.isdir(d):
            sys.stderr.write("FATAL: not a directory: {0}\n".format(d))
            return 2

    ev_off, ph_off = read_steps(off)
    ev_on, ph_on = read_steps(on)
    if not ev_off or not ev_on:
        sys.stderr.write("FATAL: no per-step timings; is driver.log beside the sample dir?\n")
        return 2

    print("=== two-arm step timing ===")
    print("  OFF : {0}".format(off))
    print("  ON  : {0}".format(on))
    print()

    # ---- 0. the null arm, which is what makes this a measurement -------------------------
    print("--- 0. which families are NULL ARMS (bwa reference byte-identical) ---")
    nulls = []
    for fam in fams:
        state, sa, sb = ref_state(off, on, fam)
        extra = ""
        if sa is not None:
            extra = "  ({0} -> {1} bytes)".format(sa, sb)
        print("  {0:<5} {1}{2}".format(fam, state, extra))
        if state == "null arm":
            nulls.append(fam)
    if nulls:
        print("  => {0} share the node, the job and the phase with the arm under test and do".format(
            ",".join(nulls)))
        print("     IDENTICAL work. If they moved, the machine moved.")
    else:
        print("  => NO NULL ARM. Every family's reference changed, so machine drift and real")
        print("     cost are not separable from these logs. A rerun IS required here.")
    print()

    # ---- 1. phases ------------------------------------------------------------------------
    print("--- 1. phases ---")
    print("  {0:<24}{1:>10}{2:>10}{3:>10}{4:>9}".format("phase", "OFF", "ON", "delta", "pct"))
    for name in [k for k in ph_off if k in ph_on]:
        o, n = ph_off[name], ph_on[name]
        print("  {0:<24}{1:>10.1f}{2:>10.1f}{3:>+10.1f}{4:>8.1f}%".format(
            name, o, n, n - o, 100.0 * (n - o) / o if o else 0.0))
    print()

    # ---- 2. the critical path ---------------------------------------------------------------
    cf_off, _ = critical_path(ev_off, ph_off)
    cf_on, _ = critical_path(ev_on, ph_on)
    print("--- 2. the critical path (the phase is a MAX over families, not a sum) ---")
    for label, ev, ph, cf in (("OFF", ev_off, ph_off, cf_off), ("ON", ev_on, ph_on, cf_on)):
        chain = sum(d for _s, _e, f, _st, d in ev if f == cf)
        rep = ph.get("downstream_parallel")
        agree = ("== reported {0:.1f}s".format(rep)
                 if rep is not None and abs(chain - rep) < 2.0
                 else "reported {0}".format("n/a" if rep is None else "{0:.1f}s".format(rep)))
        print("  {0}: longest family = {1}, its serial chain = {2:.1f}s  {3}".format(
            label, cf, chain, agree))
    print("  => only the critical family's steps can move the phase.")
    print()

    # ---- 3. per-step, with contention ------------------------------------------------------
    dur_off = {(f, s): d for _a, _b, f, s, d in ev_off}
    dur_on = {(f, s): d for _a, _b, f, s, d in ev_on}

    def overlapped(ev, fam, step):
        for s, e, f, st, _d in ev:
            if f == fam and st == step:
                others = set()
                for s2, e2, f2, _st2, _d2 in ev:
                    if f2 != fam and s2 < e and e2 > s:
                        others.add(f2)
                return ",".join(sorted(others)) or "ALONE"
        return "?"

    print("--- 3. per-step deltas, largest first, with what overlapped each step ---")
    print("  {0:<20}{1:>9}{2:>9}{3:>9}   {4:<14}{5}".format(
        "step", "OFF", "ON", "delta", "OFF ovlp", "ON ovlp"))
    rows = []
    for key in dur_on:
        if key in dur_off:
            rows.append((dur_on[key] - dur_off[key], key))
    for delta, (fam, step) in sorted(rows, reverse=True):
        o_ov = overlapped(ev_off, fam, step)
        n_ov = overlapped(ev_on, fam, step)
        flag = ""
        # The sharpest single line this probe can print: a step that ran uncontended in BOTH
        # arms and still moved is the machine, whatever the flag did.
        if o_ov == "ALONE" and n_ov == "ALONE" and abs(delta) >= 5.0:
            flag = "  <-- UNCONTENDED IN BOTH ARMS: machine, not work"
        print("  {0:<20}{1:>9.1f}{2:>9.1f}{3:>+9.1f}   {4:<14}{5}{6}".format(
            fam + ":" + step, dur_off[(fam, step)], dur_on[(fam, step)], delta,
            o_ov, n_ov, flag))
    print()

    # ---- 4. inside the transduction step ----------------------------------------------------
    print("--- 4. inside the transduction step (x_transduction.td_phase, already logged) ---")
    for fam in fams:
        a = read_td_phase(off, sample, fam)
        b = read_td_phase(on, sample, fam)
        if not a or not b:
            print("  {0}: no [TD-PHASE] lines in {1}_{0}.steps.log".format(fam, sample))
            continue
        state, _sa, _sb = ref_state(off, on, fam)
        print("  === {0}  (reference: {1}) ===".format(fam, state))
        print("    {0:<38}{1:>9}{2:>9}{3:>9}".format("td-phase", "OFF", "ON", "delta"))
        bwa_delta = 0.0
        for lab in a:
            if lab not in b:
                continue
            d = b[lab] - a[lab]
            if "bwa_realign" in lab:
                bwa_delta += d
            print("    {0:<38}{1:>9.2f}{2:>9.2f}{3:>+9.2f}{4}".format(
                lab, a[lab], b[lab], d, "  <<<" if abs(d) >= 5.0 else ""))
        print("    ---- the bwa realignment passes together: {0:+.2f}s ----".format(bwa_delta))
    print()

    print("--- how to read this ---")
    if nulls:
        print("  A null arm that is flat or faster while the arm under test is slower means the")
        print("  cost is REAL and the node is exonerated. A null arm that moved with it means the")
        print("  machine moved. Attribute only what the flag can reach, and quote the step, not")
        print("  the phase: 40% of this pair's phase delta was in flag-independent steps.")
    else:
        print("  No null arm -- this run cannot separate drift from cost. Rerun both arms.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
