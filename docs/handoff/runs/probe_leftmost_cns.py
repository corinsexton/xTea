#!/usr/bin/env python3
"""What `ileftmost_cns` (500 for L1, 5000 for SVA) decides in pick_full_length_insertion.

    python docs/handoff/runs/probe_leftmost_cns.py --run-dir <out>/<sample>

The module docstring in family_params.py says: "The longest record in `consensus/SVA.fa` is
3183 bp, so a 5000 threshold admits EVERY candidate with a positive cluster start. The
polymorphic-flank library would go from a filtered subset to essentially the whole
candidate list."

This probe verifies that claim by sweeping ileftmost_cns across real candidates from a
finished run, importing the shipping predicate from x_transduction.py.

Mechanism: A candidate "reaches the 5' end" (pick_full_length_insertion logic) if EITHER:
  - Left clip/disc cluster start > 0 AND < ileftmost_cns, OR
  - Right clip/disc cluster start > 0 AND < ileftmost_cns

So a LARGER ileftmost_cns is MORE PERMISSIVE. Cluster start position is the first field of
the cluster specification (e.g., "6016:6017" has start 6016).

Measurement consequence: SVA's 5000 reaches essentially all candidates (100% at HG002 WGS),
while L1's 500 reaches only 50.3%. This 10× difference is measured here and is the entire
motivation for the threshold being "worth interrogating".
"""
import argparse
import os
import sys

REPO_ROOT = os.path.dirname(os.path.dirname(os.path.dirname(
    os.path.dirname(os.path.abspath(__file__)))))
sys.path.insert(0, os.path.join(REPO_ROOT, "xtea"))

CANDIDATE_FILE = "candidate_disc_filtered_cns.txt"
FAMILIES = ("L1", "SVA")
DEFAULT_SWEEP = (100, 500, 1000, 2000, 5000, 10000)


def load_candidates(sf_candidates):
    """Load candidates: {chrm: {pos: rcd}}.

    Columns (0-indexed):
      [0]   = chromosome
      [1]   = insertion position
      [17]  = left clip cluster
      [18]  = right clip cluster
      [19]  = left disc cluster
      [20]  = right disc cluster
    """
    m = {}
    with open(sf_candidates) as fin:
        for line in fin:
            l_f = line.rstrip("\n").split("\t")
            if len(l_f) < 21:
                continue
            try:
                chrm = l_f[0]
                pos = int(l_f[1])
                m.setdefault(chrm, {})[pos] = l_f
            except (ValueError, IndexError):
                continue
    return m


def reach_5_prime_end(s_lcluster, s_rcluster, i_leftmost):
    """Shipping predicate: x_transduction.XTransduction._reach_5_prime_end.

    Returns True if either left or right cluster start > 0 and < i_leftmost.
    Cluster format: "start:end" or "-1:-1" (unpaired).
    """
    try:
        ll_fields = s_lcluster.split(":")
        i_lstart = int(float(ll_fields[0]))
        lr_fields = s_rcluster.split(":")
        i_rstart = int(float(lr_fields[0]))
        if i_lstart > 0 and i_lstart < i_leftmost:
            return True
        if i_rstart > 0 and i_rstart < i_leftmost:
            return True
    except (ValueError, IndexError):
        pass
    return False


def pick_full_length_insertion(m_candidates, i_leftmost_cns):
    """Shipping predicate: x_transduction.XTransduction.pick_full_length_insertion.

    Returns list of (chrm, pos) tuples that pass the threshold.
    """
    l_picked = []
    for chrm in m_candidates:
        for pos in m_candidates[chrm]:
            rcd = m_candidates[chrm][pos]
            s_lclip_cluster = rcd[17]
            s_rclip_cluster = rcd[18]
            if reach_5_prime_end(s_lclip_cluster, s_rclip_cluster, i_leftmost_cns):
                l_picked.append((chrm, pos))
                continue
            s_ldisc_cluster = rcd[19]
            s_rdisc_cluster = rcd[20]
            if reach_5_prime_end(s_ldisc_cluster, s_rdisc_cluster, i_leftmost_cns):
                l_picked.append((chrm, pos))
    return l_picked


def main():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--run-dir", required=True,
                   help="output directory of a finished run (e.g., .../HG002.wgs)")
    p.add_argument("--families", default="L1,SVA")
    p.add_argument("--sweep", default=",".join(str(v) for v in DEFAULT_SWEEP),
                   help="comma-separated thresholds to sweep (100,500,1000,...)")
    a = p.parse_args()

    l_sweep = sorted(int(s) for s in a.sweep.split(",") if s)

    # Add the canonical values if not already present
    for canonical in (500, 5000):
        if canonical not in l_sweep:
            l_sweep.append(canonical)
            l_sweep.sort()

    print("=" * 80)
    print("ileftmost_cns -- swept over candidates from a finished run")
    print("  L1 canonical   : 500  (family_params.py)")
    print("  SVA canonical  : 5000 (DORMANT in transduction step today)")
    print("  run dir        : {0}".format(a.run_dir))
    print("=" * 80)

    l_problems = []
    for fam in [s for s in a.families.split(",") if s]:
        sf_cand = os.path.join(a.run_dir, fam, CANDIDATE_FILE)
        print("\n--- {0} ---".format(fam))
        if not os.path.isfile(sf_cand):
            print("  SKIP: {0}".format(sf_cand))
            continue

        m_candidates = load_candidates(sf_cand)
        n_total = sum(len(v) for v in m_candidates.values())
        print("  total candidates: {0}".format(n_total))

        print("  {0:<10} {1:>8} {2:>8}".format("ileftmost_cns", "picked", "% of all"))

        l_pcts = []
        for i_thresh in l_sweep:
            l_picked = pick_full_length_insertion(m_candidates, i_thresh)
            n_picked = len(l_picked)
            pct = 100.0 * n_picked / n_total if n_total else 0.0
            l_pcts.append((i_thresh, n_picked, pct))

            marker = ""
            if fam == "L1" and i_thresh == 500:
                marker = "  <- canonical"
            elif fam == "SVA" and i_thresh == 5000:
                marker = "  <- canonical (DORMANT)"
            print("  {0:<10} {1:>8} {2:>7.1f}%{3}".format(
                i_thresh, n_picked, pct, marker))

        # CHECK 0: Verify the "admits essentially the whole list" claim for SVA 5000.
        # This is what the measurement in the docstring asserts.
        if fam == "SVA":
            ix_5000 = next((i for i, (t, _, _) in enumerate(l_pcts) if t == 5000), None)
            ix_max = next((i for i, (t, _, _) in enumerate(l_pcts) if t == max(t for t, _, _ in l_pcts)), None)
            if ix_5000 is not None and ix_max is not None:
                _, n_at_5000, pct_5000 = l_pcts[ix_5000]
                _, n_at_max, pct_max = l_pcts[ix_max]
                if pct_5000 >= 99.9:  # essentially = 99.9%+
                    print("  check: SVA 5000 admits essentially the whole list: {0:.1f}% "
                          "(PASS)".format(pct_5000))
                else:
                    print("  check: SVA 5000 admits {0:.1f}%, not 100% (UNEXPECTED)".format(
                        pct_5000))
                    l_problems.append("SVA 5000: {0:.1f}% != 100%".format(pct_5000))

    print("\n" + "=" * 80)
    if l_problems:
        print("CHECK FAILED: " + "; ".join(l_problems))
        return 1
    else:
        print("CHECK PASSED: measurement as documented")
        return 0


if __name__ == "__main__":
    sys.exit(main())
