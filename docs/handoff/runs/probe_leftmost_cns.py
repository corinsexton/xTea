#!/usr/bin/env python3
"""What `ileftmost_cns` decides, swept through the SHIPPING predicate over a finished run.

    python docs/handoff/runs/probe_leftmost_cns.py --run-dir <out>/<sample>

WHAT THIS IS FOR. `HANDOFF_2026-08-04_sva_decoy_fix.md` §7 item 4 calls this "the ONLY one of
plan §5's four 'interrogate' thresholds still worth an item", and the reason given is a 10x
difference between the families: `family_params` carries `leftmost_cns` 500 for L1 and 5000 for
SVA. Its own docstring states the consequence:

    "`pick_full_length_insertion` keeps a candidate whose cns cluster start is
     `0 < start < leftmost_cns`. The longest record in `consensus/SVA.fa` is 3183 bp, so a 5000
     threshold admits EVERY candidate with a positive cluster start. The polymorphic-flank
     library would go from a filtered subset to essentially the whole candidate list, and every
     downstream bwa realignment with it. This moves calls in bulk; it is a threshold to calibrate
     against the truth set, not to flip blind."

**THAT WAS AN INFERENCE FROM A SEQUENCE LENGTH, NEVER A MEASUREMENT ON A RUN.** This probe is the
measurement. Verified independently: `consensus/SVA.fa`'s longest record IS 3183 bp, and
`consensus/LINE1.fa`'s is 6342 bp -- which is why the saturation is a fact about SVA and not
about the threshold, and why L1 at 5000 does NOT saturate.

TWO THINGS TO READ BEFORE READING THE NUMBERS, both of which reframe the agenda item.

1. **SVA's 5000 IS DORMANT, and that is documented, not discovered.** `x_TEA_main.py:311` resolves
   `x_family = family_params.from_sva_flag(options.sva)`, and `run_xtea.build_type_steps` passes
   `--sva` to **cns_filter alone** (matching the golden `run_xTEA_pipeline.sh`), so the
   transduction step resolves `GENERIC`, whose `leftmost_cns` is **500**.
   `family_params.py:18` says so outright: "`leftmost_cns` was therefore **always 500** -- SVA's
   5000 branch". So the live 10x difference is **zero**: today both families run 500. The item is
   not "why do they differ", it is "should SVA's dormant 5000 be woken", and waking it is
   documented as a one-word change (resolve from `--rtype` instead of `--sva`).

2. **A LARGER value is MORE PERMISSIVE.** `_reach_5_prime_end` keeps a candidate when a cluster
   start is `0 < start < i_leftmost`, so raising the bound admits more. Easy to read backwards.

IT DOES NOT REIMPLEMENT THE PREDICATE. It instantiates `x_transduction.XTransduction` and calls
`pick_full_length_insertion` itself, which is the same rule `probe_td_pool.py` and
`probe_td_source_classes.py` follow, and it exists because **a probe with its own copy of a rule
can agree with a scorer that disagrees with it.** An earlier version of this file DID carry its
own copy of `_reach_5_prime_end` and its own candidate parser while its docstring claimed
otherwise, and the two were not equivalent: the copy swallowed `ValueError` on a malformed
cluster field and skipped rows with fewer than 21 columns, where the shipping code raises and
`load_in_candidate_list_str_version` does not. A probe that is more forgiving than the step
under-reports exactly the rows worth looking at.

WHAT IT CANNOT SAY. Whether the extra candidates are GOOD. It counts what crosses the threshold;
only a truth-set-scored arm can say whether they should. See the handoff for what that arm is.
"""
import argparse
import os
import sys

REPO_ROOT = os.path.dirname(os.path.dirname(os.path.dirname(
    os.path.dirname(os.path.abspath(__file__)))))
sys.path.insert(0, os.path.join(REPO_ROOT, "xtea"))

import x_transduction  # noqa: E402
from x_intermediate_sites import XIntemediateSites  # noqa: E402

# The cns_filter step's output, which is `pick_full_length_insertion`'s `sf_te_candidates`.
CANDIDATE_FILE = "candidate_disc_filtered_cns.txt"
DEFAULT_SWEEP = (100, 500, 1000, 2000, 5000, 10000)
# family -> (canonical value, is it LIVE in the transduction step today)
CANONICAL = {"L1": (500, True), "SVA": (5000, False)}


def main():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--run-dir", required=True,
                   help="a finished run's sample directory, e.g. <out>/HG002.wgs")
    p.add_argument("--families", default="L1,SVA")
    p.add_argument("--sweep", default=",".join(str(v) for v in DEFAULT_SWEEP))
    a = p.parse_args()

    sweep = sorted(set(int(s) for s in a.sweep.split(",") if s.strip())
                   | set(v for v, _live in CANONICAL.values()))

    print("=" * 78)
    print("ileftmost_cns, swept through the SHIPPING pick_full_length_insertion")
    print("  run dir : {0}".format(a.run_dir))
    print("  L1      : canonical 500,  and it is LIVE")
    print("  SVA     : canonical 5000, and it is DORMANT -- the transduction step resolves")
    print("            GENERIC (500), see family_params.py:18 and x_TEA_main.py:305-311")
    print("  larger = MORE permissive")
    print("=" * 78)

    # One instance, reused: `pick_full_length_insertion` reads only its two arguments, and the
    # constructor's other work (the reference path, the decoy policy) is irrelevant here. Passing
    # a reference that does not exist is deliberate -- if a future version of the method starts
    # touching it, this probe should fail loudly rather than quietly measure something else.
    xtd = x_transduction.XTransduction(a.run_dir + "/", 1, "no_such_ref.fa")
    xintmdt = XIntemediateSites()

    problems = []
    for fam in [s.strip() for s in a.families.split(",") if s.strip()]:
        sf_cand = os.path.join(a.run_dir, fam, CANDIDATE_FILE)
        print()
        print("--- {0} ---".format(fam))
        if not os.path.isfile(sf_cand):
            print("  SKIP: no {0}".format(sf_cand))
            continue
        # The denominator comes from the SHIPPING loader too, so it cannot disagree with the
        # numerator about which rows are candidates at all.
        m_all = xintmdt.load_in_candidate_list_str_version(sf_cand)
        n_total = sum(len(v) for v in m_all.values())
        print("  candidates: {0}".format(n_total))
        print("  {0:>13} {1:>8} {2:>9}".format("ileftmost_cns", "picked", "% of all"))
        pct_at = {}
        for thresh in sweep:
            n_picked = len(xtd.pick_full_length_insertion(sf_cand, thresh))
            pct = 100.0 * n_picked / n_total if n_total else float("nan")
            pct_at[thresh] = pct
            note = ""
            canon, live = CANONICAL.get(fam, (None, None))
            if thresh == canon:
                note = "  <- canonical, {0}".format("LIVE" if live else "DORMANT")
            print("  {0:>13} {1:>8} {2:>8.1f}%{3}".format(thresh, n_picked, pct, note))

        # THE CHECK, and it is the docstring's claim rather than a round number: SVA's 5000 must
        # admit essentially everything, because 5000 exceeds SVA.fa's longest record (3183 bp).
        # L1 must NOT, because LINE1.fa's longest is 6342 bp -- so this asserts the ASYMMETRY,
        # which is the part that would break if the libraries or the column layout changed.
        if fam == "SVA" and 5000 in pct_at:
            if pct_at[5000] >= 99.9:
                print("  check: SVA at 5000 admits {0:.1f}% -- PASS, the docstring's claim holds"
                      .format(pct_at[5000]))
            else:
                problems.append("SVA at 5000 admits {0:.1f}%, not ~100%. family_params' claim "
                                "rests on 5000 > SVA.fa's longest record (3183 bp); if this is "
                                "below 100 that premise or the column layout has moved"
                                .format(pct_at[5000]))
        if fam == "L1" and 5000 in pct_at:
            if pct_at[5000] < 99.9:
                print("  check: L1 at 5000 admits {0:.1f}% -- PASS, it does NOT saturate, which "
                      "is the asymmetry (LINE1.fa's longest record is 6342 bp)"
                      .format(pct_at[5000]))
            else:
                problems.append("L1 at 5000 admits {0:.1f}%. It should not saturate -- LINE1.fa's "
                                "longest record is 6342 bp, above 5000. If it does, the "
                                "saturation is not about consensus length after all and the "
                                "reasoning in family_params' docstring needs re-deriving"
                                .format(pct_at[5000]))

    print()
    print("=" * 78)
    if problems:
        print("CHECK FAILED:")
        for s in problems:
            print("  {0}".format(s))
        return 1
    print("CHECK PASSED -- and note what this does NOT say: it counts what crosses the")
    print("threshold, not whether those candidates should. Only a truth-set-scored arm can")
    print("say that, and family_params' own docstring says so: 'not to flip blind'.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
