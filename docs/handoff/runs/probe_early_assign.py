#!/usr/bin/env python3
"""Stage 2c.3's prediction: what moves if family is decided at cns_filter instead
of at the postF_hc merge.

    python docs/handoff/runs/probe_early_assign.py --run <sample_root>

Runs offline against a completed run's kept intermediates -- no BAM, no compute --
so the prediction exists BEFORE the gate does, which is the only order in which
"calls may move, but only in a predicted direction" means anything.

WHY cns_filter IS THE ANSWER TO "how far up front", measured on BASELINE-0:

  * The rule needs `two_side_tprt_both`, which decides 13 of the merge's 22
    conflicts. That flag is written BY cns_filter (x_clip_disc_filter.py:1845).
    Nothing earlier has it, so anything earlier is a different rule, not an
    earlier one.
  * There is no family signal earlier to substitute. `scan_clip_realign.tsv`
    looks like per-family evidence -- three rows per position, one per family,
    each from its own bwa (realign.hpp) -- but the counts come out identical: of
    155,615 positions with >=5 clip support, only 504 differ across families at
    all. The three clip candidate lists are the same 577K positions (577293 /
    577312 / 577347).
  * Family-discriminating evidence first appears on the DISC side, where the mate
    is tested against each family's own annotation: 2244 / 9468 / 1183 candidates,
    union 12,223, only 672 rows redundant. Real, but only 5% overlap -- so an
    up-front decision there would buy 5% of cns_filter and cost the TPRT override.

WHAT THIS PRINTS, and why both halves are needed:

  (a) DISAGREEMENT -- the locus is contended at cns level and the cns-level rule
      names a different family than the postF_hc-level rule does.
  (b) LOSS -- the cns-level winner is later dropped by its own postF while the
      loser would have survived. Today the loser is called; with the decision
      moved up, nobody calls it. (a) alone would report these as mere flips and
      hide that a call disappears.

ONE APPROXIMATION, stated because it bounds what the gate may assume: this holds
each family's transduction/sibling output fixed. Those steps CREATE sites (L1 cns
300 -> cns2 315), so a real run may move loci this does not name. The gate
therefore reports every moved locus rather than asserting this exact list.
"""
import argparse
import os
import sys

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)),
                                "..", "..", "..", "xtea"))

from locus_table import LocusTable, load_rows, TYPE_PRIORITY  # noqa: E402

CNS = "candidate_disc_filtered_cns.txt"
HCPF = "candidate_disc_filtered_cns.txt.high_confident.post_filtering.txt"


def assign(rows_by_type, window):
    t = LocusTable(window=window)
    loci = t.assign_rows(rows_by_type)
    return {(l.chrm, l.pos): l for l in loci}, t


def nearest(key, by_chrm, window):
    chrm, pos = key
    best = None
    for p in by_chrm.get(chrm, ()):
        d = abs(p - pos)
        if d <= window and (best is None or d < best[0]):
            best = (d, p)
    return None if best is None else (chrm, best[1])


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--run", required=True,
                    help="sample_root of a completed --keep-intermediates run")
    ap.add_argument("--window", type=int, default=50)
    a = ap.parse_args()

    cns = {t: load_rows(os.path.join(a.run, t, CNS), t) for t in TYPE_PRIORITY}
    hcpf = {t: load_rows(os.path.join(a.run, t, HCPF), t) for t in TYPE_PRIORITY}

    print("input population")
    for t in TYPE_PRIORITY:
        print(f"  {t:4s} cns={len(cns[t]):6d}  postF_hc={len(hcpf[t]):6d}")

    late, late_tbl = assign(hcpf, a.window)     # what the pipeline does today
    early, early_tbl = assign(cns, a.window)    # what a decision at cns_filter does

    print(f"\nassignment populations: late(postF_hc)={len(late)}  early(cns)={len(early)}")
    for label, tbl in (("late (the 22 the plan gates on)", late_tbl), ("early", early_tbl)):
        print(f"conflict summary -- {label}")
        for k, v in sorted(tbl.conflict_summary().items(), key=lambda x: -x[1]):
            print(f"  {k[0]:4s} over {k[1]:4s} by {k[2]:14s} {v}")

    early_by_chrm = {}
    for (c, p) in early:
        early_by_chrm.setdefault(c, []).append(p)
    hcpf_pos = {t: {(c, p) for (c, p, _, _, _) in
                    ((r.chrm, r.pos, 0, 0, 0) for r in hcpf[t])} for t in TYPE_PRIORITY}

    agree = unmatched = 0
    moved, lost = [], []
    for key, locus in late.items():
        ek = nearest(key, early_by_chrm, a.window)
        if ek is None:
            unmatched += 1
            continue
        w = early[ek]
        if w.type == locus.type:
            agree += 1
            continue
        survived = any(c == w.chrm and abs(p - w.pos) <= a.window
                       for (c, p) in hcpf_pos[w.type])
        (moved if survived else lost).append((key, locus.type, w.type))

    print(f"\n=== (a) DISAGREEMENT over the {len(late)} called loci ===")
    print(f"  same family early vs late : {agree}")
    print(f"  DIFFERENT family          : {len(moved) + len(lost)}")
    print(f"  no early locus within {a.window}bp: {unmatched}"
          f"   (transduction/sibling created these after cns)")
    for (key, lt, et) in sorted(moved):
        print(f"    FLIP {key[0]}:{key[1]}  {lt} -> {et}")

    print(f"\n=== (b) LOSS: early winner does not survive its own postF ===")
    print(f"  loci called today that a cns-level decision would lose: {len(lost)}")
    for (key, lt, et) in sorted(lost):
        print(f"    LOST {key[0]}:{key[1]}  today={lt}; early gave it to {et}, "
              f"which postF then drops")


if __name__ == "__main__":
    main()
