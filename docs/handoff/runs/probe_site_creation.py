#!/usr/bin/env python3
"""Where does a call's record first exist? — the measurement Stage 2c's gate said was missing.

WHY THIS EXISTS. The Stage 2c gate (job 49502115) found that the family-assignment barrier moves
zero calls, and the reason was structural: 15 of the 22 family conflicts the postF_hc merge
resolves are UNREACHABLE from `candidate_disc_filtered_cns.txt`, because the losing family has no
cns row there at all — transduction and sibling create that record *after* the barrier runs. So
"decide family once" cannot collapse to one decision point until **site creation** is unified, and
that dependency had never been quantified. This script quantifies it, offline, from BASELINE-0's
kept intermediates. No BAM, no pipeline run.

WHAT IT ANSWERS, per family:

  1. how many final calls (postF_hc) have a cns row at all — i.e. existed before the barrier;
  2. of the rest, which step invented them: transduction (`candidate_disc_filtered_cns2.txt`) or
     sibling (`candidate_sibling_transduction2.txt`);
  3. for the merge's conflict population, whether BOTH sides already exist at cns — which is
     exactly the condition for a pre-transduction decision point to be able to resolve it.

MEASURED on the Stage 2c gate's OFF arm (HG002 WGS, `stage2c_family_assign_49502115/OFF`):

    postF_hc calls whose record first exists in ...   cns    cns2 (transduction)   sibling
      L1   200                                        150     50                   0
      Alu 1230                                       1230      0                   0
      SVA   85                                         83      2                   0

  So **52 of 1,515 per-family final records — every one of them L1 or SVA, every one of them from
  the transduction step, none from sibling — do not exist when the barrier runs.** And of the 22
  family conflicts, **exactly 7 have both sides present at cns**; the other 15 have a winner or a
  loser that transduction had not yet invented.

  That 7 is the gate's own number, reached from the opposite direction: the gate counted which
  conflicts the barrier *actually* resolved and got 7 with 0 disagreements, this counts which ones
  it *could* have and gets the same 7. Two independent derivations, so the ceiling on any decision
  point placed before transduction is 7 of 22 — not a tuning question, a site-creation one.

USAGE
    python docs/handoff/runs/probe_site_creation.py <run_root>/<sample>
e.g.
    python docs/handoff/runs/probe_site_creation.py \
        /n/scratch/users/c/cos689/xtea2_testdata/baseline0_hg002_L1-Alu-SVA/HG002.wgs

READ THE "sibling" COLUMN AS "ONLY in sibling". The stages are tried in chain order — cns, then
cns2, then sibling — so a record sibling also carries, having read it out of cns2, is attributed to
cns2. That is the right attribution for the question being asked (when does this record first
exist?), and it is why sibling's column can be 0 while its file has 48 rows. It is NOT evidence
that the sibling step creates nothing; it is evidence that nothing reaches postF_hc *through*
sibling alone — consistent with the chain, whose postF_hc input is
`candidate_disc_filtered_cns2.txt.high_confident`, not sibling's own output.

MATCHING. Positions are matched with locus_table's own ±MERGE_WINDOW rule, not by equality: the
chain shifts a call's coordinate (postF re-anchors, transduction re-selects), so exact matching
would report invented sites that are merely moved. Nearest-within-window, one-to-one, greedy in
position order — the same discipline benchmark_truth_set.py uses for truth matching, and for the
same reason: without it, one cns row can absorb several late calls and the count is meaningless.

    MEASURED 2026-08-03: on THIS hop the window turns out to do no work, and that strengthens
    the result rather than weakening it. Every postF_hc position is present at distance 0 in
    `candidate_disc_filtered_cns2.txt` — the post-cns chain does not re-anchor at all — and the
    nearest cns row to each of the 52 unmatched records is >= 3,293 bp away, median 2.4 Mb. So
    the 52 are literally created, and no wider window recovers them. Their origin is now known
    too: all 52 are already rows in `<disc>.clip_sites_raw_disc.txt`, which the disc adapter
    writes BEFORE cns_filter runs. See xtea/site_universe.py and
    docs/handoff/runs/gate_site_universe.py.
"""
import os
import sys
from bisect import bisect_left

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)),
                                "..", "..", "..", "xtea"))
import locus_table  # noqa: E402

TYPES = ("L1", "Alu", "SVA")
WIN = locus_table.MERGE_WINDOW


def read_positions(path):
    """{chrom: sorted [pos]} from any of the chain's raw-result files (col 0, col 1)."""
    out = {}
    if not os.path.exists(path):
        return None
    with open(path) as fh:
        for line in fh:
            if not line.strip() or line.startswith("#"):
                continue
            f = line.split()
            if len(f) < 2:
                continue
            try:
                out.setdefault(f[0], []).append(int(f[1]))
            except ValueError:
                continue
    for v in out.values():
        v.sort()
    return out


def matched(late, early, window=WIN):
    """Split `late` into (has an early record within window, does not). One-to-one: an early
    position is consumed by the first late call that claims it, so two late calls 10 bp apart
    cannot both be explained by one early row."""
    hit, miss = [], []
    for chrom in sorted(late):
        cand = sorted(early.get(chrom, []))
        used = [False] * len(cand)
        for p in sorted(late[chrom]):
            i = bisect_left(cand, p - window)
            best, bd = -1, window + 1
            while i < len(cand) and cand[i] <= p + window:
                if not used[i] and abs(cand[i] - p) < bd:
                    best, bd = i, abs(cand[i] - p)
                i += 1
            if best >= 0:
                used[best] = True
                hit.append((chrom, p))
            else:
                miss.append((chrom, p))
    return hit, miss


def main():
    if len(sys.argv) != 2:
        sys.exit(__doc__)
    root = sys.argv[1]

    stages = {}   # type -> {stage: {chrom: [pos]}}
    for t in TYPES:
        d = os.path.join(root, t)
        if not os.path.isdir(d):
            continue
        stages[t] = {
            "cns":    read_positions(os.path.join(d, "candidate_disc_filtered_cns.txt")),
            "cns2":   read_positions(os.path.join(d, "candidate_disc_filtered_cns2.txt")),
            "sib":    read_positions(os.path.join(d, "candidate_sibling_transduction2.txt")),
            "hc_pf":  read_positions(os.path.join(
                d, "candidate_disc_filtered_cns.txt.high_confident.post_filtering.txt")),
        }
    if not stages:
        sys.exit(f"no per-type directories under {root}")

    def n(d):
        return 0 if d is None else sum(len(v) for v in d.values())

    print("=" * 78)
    print("SITE CREATION — how many of each family's final calls existed before the barrier")
    print("=" * 78)
    print(f"{'':6s} {'cns':>8s} {'cns2':>8s} {'sibling':>8s} {'postF_hc':>9s}")
    for t, s in stages.items():
        print(f"{t:6s} {n(s['cns']):8d} {n(s['cns2']):8d} {n(s['sib']):8d} {n(s['hc_pf']):9d}")

    print()
    print("Of each family's postF_hc calls, where the record FIRST exists "
          f"(+-{WIN} bp, one-to-one):")
    print(f"{'':6s} {'total':>7s} {'in cns':>8s} {'in cns2':>8s} {'in sib':>7s} {'nowhere':>8s}")
    invented = {}
    for t, s in stages.items():
        if s["hc_pf"] is None:
            continue
        hc = s["hc_pf"]
        in_cns, rest = matched(hc, s["cns"] or {})
        rest_by_chrom = {}
        for c, p in rest:
            rest_by_chrom.setdefault(c, []).append(p)
        in_cns2, rest2 = matched(rest_by_chrom, s["cns2"] or {})
        rest2_by_chrom = {}
        for c, p in rest2:
            rest2_by_chrom.setdefault(c, []).append(p)
        in_sib, nowhere = matched(rest2_by_chrom, s["sib"] or {})
        invented[t] = set(rest)
        print(f"{t:6s} {n(hc):7d} "
              f"{len(in_cns):8d} {len(in_cns2):8d} {len(in_sib):7d} {len(nowhere):8d}")

    # --- the claim the Stage 2c gate rested on -------------------------------------------------
    conf = os.path.join(root, "merged", "locus_family_conflicts.tsv")
    if not os.path.exists(conf):
        print(f"\n(no {conf}; skipping the conflict-population leg)")
        return
    print()
    print("=" * 78)
    print("THE CONFLICT POPULATION — could the barrier have seen BOTH sides at cns?")
    print("=" * 78)
    rows = []
    with open(conf) as fh:
        # locus_table writes the header commented ("#chrm\t..."), so strip the marker or every
        # lookup below silently misses and the leg reports a confident zero.
        head = fh.readline().lstrip("#").rstrip("\n").split("\t")
        for line in fh:
            rows.append(dict(zip(head, line.rstrip("\n").split("\t"))))
    if not rows:
        print("no conflicts recorded")
        return

    # A conflict is resolvable at the barrier only if BOTH families already have a record there.
    # The gate reported this from the other end -- which conflicts the barrier actually resolved --
    # and the two must agree: a conflict with both sides at cns is one the barrier saw.
    n_both, by_rule = 0, {}
    for r in rows:
        chrom = r["chrm"]
        win_t, lose_t = r["winner_type"], r["loser_type"]
        w_inv = (chrom, int(r["winner_pos"])) in invented.get(win_t, ())
        l_inv = (chrom, int(r["loser_pos"])) in invented.get(lose_t, ())
        where = ("both at cns" if not (w_inv or l_inv) else
                 "loser invented" if l_inv and not w_inv else
                 "winner invented" if w_inv and not l_inv else "both invented")
        n_both += not (w_inv or l_inv)
        k = (win_t, lose_t, r["rule"], where)
        by_rule[k] = by_rule.get(k, 0) + 1
    for k in sorted(by_rule):
        win, lose, rule, where = k
        print(f"  {win:4s} over {lose:4s} by {rule:14s} {where:16s} {by_rule[k]}")
    print(f"\n{n_both} of {len(rows)} conflicts have BOTH sides present at cns — those and only")
    print("those are the ones a decision point before transduction can resolve. The other")
    print(f"{len(rows) - n_both} need unified SITE CREATION, because one side does not exist yet.")


if __name__ == "__main__":
    main()
