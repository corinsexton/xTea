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

WHAT THIS PRINTS -- three populations, because the interesting one is the third:

  * conflicts BOTH stages see, and whether they agree. A disagreement here is a
    call that moves; zero disagreements means the relocation is output-neutral.
  * conflicts resolved ONLY at cns -- the barrier pre-empts them.
  * conflicts resolved ONLY at the merge -- the losing family has no cns row,
    because transduction/sibling create its record after the barrier has run. This
    contention is not reachable from cns at all, which is why a resolution point
    after sibling remains necessary and Stage 2c cannot end with one decision until
    site creation is unified too.

WHAT THIS SCRIPT GOT WRONG THE FIRST TIME, corrected here after job 49502115 ran
the real thing and moved ZERO calls against this script's predicted four:

  * **The two resolutions are not alternatives; both run.** The barrier resolves
    family over cns candidates and the postF_hc merge then resolves again over
    whatever survives. So the barrier changes the OUTPUT only where it removes a
    row that would have WON at the merge. On HG002 WGS it never does: all 7
    conflicts it pre-empts are decided identically at both ends.
  * **"late=L1, early=Alu" is usually not a flip.** The old comparison matched a
    late locus to the NEAREST early locus of ANY family. When L1 has no cns row at
    that position -- because transduction/sibling create the L1 record later -- a
    nearby Alu or SVA cns row got matched to it and reported as a family flip. It
    is not one: the barrier cannot take a locus away from a family that has not
    claimed it yet. Measured: the barrier left L1 at 300 of 300 rows.

So the comparison below is between the two CONFLICT POPULATIONS, not between two
assignment maps. That is the question that predicts output: of the loci both
stages see, do they decide the same way?
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

    # A conflict's identity is its winner position; both stages record it the same way.
    def by_pos(tbl):
        return {(c.chrm, c.winner_pos): (c.winner_type, c.loser_type, c.rule)
                for c in tbl.conflicts}

    early_c, late_c = by_pos(early_tbl), by_pos(late_tbl)
    both = sorted(set(early_c) & set(late_c))
    only_early = sorted(set(early_c) - set(late_c))
    only_late = sorted(set(late_c) - set(early_c))

    print(f"\n=== conflicts both stages see: {len(both)} ===")
    disagree = 0
    for k in both:
        same = early_c[k] == late_c[k]
        disagree += not same
        print(f"  {'SAME ' if same else 'DIFF '} {k[0]}:{k[1]}  "
              f"cns={'/'.join(early_c[k])}  merge={'/'.join(late_c[k])}")
    print(f"  --> disagreements: {disagree}"
          f"   (each one is a call that MOVES; zero means output is unchanged)")

    print(f"\n=== resolved only at cns: {len(only_early)} ===")
    print("  the barrier pre-empts these; the merge never sees them. Harmless unless the")
    print("  loser would have won at the merge -- check the losers against the merge list.")
    for k in only_early:
        print(f"    {k[0]}:{k[1]}  {'/'.join(early_c[k])}")

    print(f"\n=== resolved only at the merge: {len(only_late)} ===")
    print("  the losing family has no cns row here -- transduction/sibling create it after")
    print("  the barrier has already run, so this contention is NOT reachable from cns and")
    print("  a resolution point after sibling stays necessary.")
    for k in only_late:
        print(f"    {k[0]}:{k[1]}  {'/'.join(late_c[k])}")


if __name__ == "__main__":
    main()
