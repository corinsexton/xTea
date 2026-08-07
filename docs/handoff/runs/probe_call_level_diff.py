#!/usr/bin/env python3
"""CALL-LEVEL diff: master vs xTea2, and each arm against its own replicates.

`probe_repro_comparison.py` answers "how many". This answers "which ones, and what did the caller
believe about them" -- the actual VCF records behind array 50062602's counts. Offline, seconds.

THREE POPULATIONS, KEPT APART, because conflating them is how a reproducibility claim goes wrong:

  A. WITHIN-ARM instability -- calls an arm makes in some replicates and not others. This is the
     arm's own noise, and it is the yardstick for B.
  B. BETWEEN-ARM difference over each arm's STABLE set (calls present in all replicates). Only a
     difference here can be attributed to the port rather than to scheduling.
  C. TRUTH-anchored differences -- of B, which ones the truth set adjudicates. A call one arm
     lacks is only a sensitivity regression if a truth record is behind it.

MATCHING IS POSITIONAL AT +-200 bp, the same window `benchmark_truth_set.py` scores at, because
the two arms do not place a breakpoint identically -- an exact-coordinate diff would report
hundreds of "differences" that are the same insertion called a few bases apart. Sections A and B
therefore report both: matched-but-moved (jitter) separately from unmatched (a real difference in
what was called).

USAGE
    python3 docs/handoff/runs/probe_call_level_diff.py
    python3 docs/handoff/runs/probe_call_level_diff.py --family SVA --limit 40
    python3 docs/handoff/runs/probe_call_level_diff.py --info-fields AF,LCLIP,RCLIP,LDISC,RDISC
"""
import argparse
import glob
import os
import sys
from collections import Counter, defaultdict

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.dirname(os.path.dirname(os.path.dirname(HERE)))
sys.path.insert(0, os.path.join(REPO, "benchmark"))
sys.path.insert(0, HERE)
import benchmark_truth_set as bts          # noqa: E402
import probe_repro_comparison as prc       # noqa: E402

FAMILIES = ("L1", "Alu", "SVA")
ARMS = ("master", "xtea2")
# The evidence a reader actually wants next to a disputed call: how much support, of what kind,
# and what the caller concluded. Deliberately short -- the full INFO is in the VCF.
DEFAULT_INFO = "SUBTYPE,TD_SRC,AF,LCLIP,RCLIP,LDISC,RDISC"


def load(arm, rep, family):
    """[(chrom, pos, info_dict)] for one replicate/family, or None if the tree is absent."""
    tree = os.path.join(prc.SCRATCH, "repro_%s_rep%s" % (arm, rep))
    vcf = prc.family_vcf(tree, arm, family)
    if not vcf:
        return None
    return [(c.chrom, c.pos, c.info) for c in bts.read_calls_vcf(vcf, [family])]


def key_set(calls):
    return {(c, p) for c, p, _ in calls}


def cluster_match(a_keys, b_keys, window):
    """Greedy nearest pairing between two key sets. Returns (pairs, a_only, b_only).

    Reuses the shipping matcher's discipline -- one-to-one, by increasing distance -- rather than
    "is anything nearby", so a single call cannot absorb several of the other arm's.
    """
    by_chrom = defaultdict(list)
    for c, p in b_keys:
        by_chrom[c].append(p)
    cands = []
    for c, p in a_keys:
        for q in by_chrom.get(c, ()):
            d = abs(p - q)
            if d <= window:
                cands.append((d, (c, p), (c, q)))
    cands.sort()
    ta, tb, pairs = set(), set(), []
    for d, ak, bk in cands:
        if ak in ta or bk in tb:
            continue
        ta.add(ak)
        tb.add(bk)
        pairs.append((ak, bk, d))
    return pairs, sorted(a_keys - ta), sorted(b_keys - tb)


def fmt(key, info, fields, truth_hit=None):
    bits = []
    for f in fields:
        v = info.get(f, "-")
        if f == "TD_SRC" and v not in ("-", "not_transduction") and len(v) > 34:
            v = v[:31] + "..."
        bits.append("%s=%s" % (f, v))
    tag = ""
    if truth_hit is True:
        tag = "  [ON-TRUTH]"
    elif truth_hit is False:
        tag = "  [off-truth]"
    return "%s:%-11d %s%s" % (key[0], key[1], " ".join(bits), tag)


def truth_lookup(family, window):
    """(chrom,pos) -> True if within `window` of a truth record. Built once per family."""
    truth = bts.load_truth_family(prc.TRUTH_DIR, family)
    by_chrom = defaultdict(list)
    for t in truth:
        by_chrom[t.chrom].append(t.pos)
    for c in by_chrom:
        by_chrom[c].sort()

    def on_truth(chrom, pos):
        return any(abs(pos - q) <= window for q in by_chrom.get(chrom, ()))
    return on_truth


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--family", default="L1,Alu,SVA")
    ap.add_argument("--reps", default="1,2,3")
    ap.add_argument("--window", type=int, default=200)
    ap.add_argument("--limit", type=int, default=25, help="max records printed per list")
    ap.add_argument("--info-fields", default=DEFAULT_INFO)
    a = ap.parse_args()
    fams = [f.strip() for f in a.family.split(",") if f.strip()]
    reps = [r.strip() for r in a.reps.split(",") if r.strip()]
    fields = [f.strip() for f in a.info_fields.split(",") if f.strip()]

    print("=" * 96)
    print("CALL-LEVEL DIFF  --  master vs xTea2, HG002 WGS, array 50062602, %d replicates each"
          % len(reps))
    print("  matching: positional, +-%d bp, one-to-one   info shown: %s"
          % (a.window, ",".join(fields)))
    print("=" * 96)

    pooled = {arm: dict(stable=0, truth_all=0, wins=0, unstable=0) for arm in ARMS}

    for family in fams:
        print("\n" + "#" * 96)
        print("# %s" % family)
        print("#" * 96)
        on_truth = truth_lookup(family, a.window)

        loaded = {}
        for arm in ARMS:
            for rep in reps:
                c = load(arm, rep, family)
                if c is None:
                    print("  (no %s VCF for %s rep %s)" % (family, arm, rep))
                else:
                    loaded[(arm, rep)] = c
        if not loaded:
            continue

        # ---------- A. within-arm instability ----------
        print("\n--- A. WITHIN-ARM: calls an arm makes in SOME replicates but not all ---")
        stable, info_of = {}, {}
        for arm in ARMS:
            sets = [key_set(loaded[(arm, r)]) for r in reps if (arm, r) in loaded]
            if len(sets) < 2:
                continue
            inter, union = set.intersection(*sets), set.union(*sets)
            stable[arm] = inter
            info_of[arm] = {}
            for r in reps:
                for c, p, i in loaded.get((arm, r), []):
                    info_of[arm].setdefault((c, p), i)
            wob = sorted(union - inter)
            print("\n  %s: %d call(s) in all %d reps, %d in the union -> %d UNSTABLE"
                  % (arm, len(inter), len(sets), len(union), len(wob)))
            if not wob:
                print("    (none -- every replicate produced the same call set)")
                continue
            # A call absent from one replicate but present within the window in that replicate is
            # JITTER, not an appearing/disappearing call. Split them, because they mean opposite
            # things about reproducibility.
            jitter, real = [], []
            for k in wob:
                present_in, near_in = [], []
                for r in reps:
                    if (arm, r) not in loaded:
                        continue
                    ks = key_set(loaded[(arm, r)])
                    if k in ks:
                        present_in.append(r)
                    elif any(c == k[0] and abs(p - k[1]) <= a.window for c, p in ks):
                        near_in.append(r)
                (jitter if near_in and len(present_in) + len(near_in) == len(sets)
                 else real).append((k, present_in, near_in))
            print("    %d are COORDINATE JITTER (the call exists in every replicate, moved <=%d bp)"
                  % (len(jitter), a.window))
            for k, pin, nin in jitter[:a.limit]:
                print("      %s   exact in rep %s, moved in rep %s"
                      % (fmt(k, info_of[arm][k], fields, on_truth(*k)),
                         ",".join(pin), ",".join(nin)))
            print("    %d APPEAR/DISAPPEAR between replicates (no nearby call in the others)"
                  % len(real))
            for k, pin, nin in real[:a.limit]:
                print("      %s   present in rep %s ONLY"
                      % (fmt(k, info_of[arm][k], fields, on_truth(*k)), ",".join(pin)))

        # ---------- B. between-arm, over the stable sets ----------
        if len(stable) < 2:
            continue
        print("\n--- B. BETWEEN-ARM, over each arm's STABLE set (present in all replicates) ---")
        pairs, m_only, x_only = cluster_match(stable["master"], stable["xtea2"], a.window)
        print("  master stable %d   xTea2 stable %d   matched %d   master-only %d   xTea2-only %d"
              % (len(stable["master"]), len(stable["xtea2"]), len(pairs), len(m_only), len(x_only)))

        off = Counter(d for _, _, d in pairs)
        exact = off.get(0, 0)
        print("\n  of the %d matched, the breakpoint is IDENTICAL in %d (%.1f%%); offsets: %s"
              % (len(pairs), exact, bts._rate(exact, len(pairs)),
                 ", ".join("%dbp:%d" % (d, n) for d, n in sorted(off.items())[:8])
                 + (" ..." if len(off) > 8 else "")))
        # Where the two arms agree a call exists, do they agree what it IS? This separates
        # "different calls" from "same call, different interpretation".
        for f in ("SUBTYPE", "TD_SRC"):
            diff = [(ak, bk) for ak, bk, _ in pairs
                    if info_of["master"][ak].get(f) != info_of["xtea2"][bk].get(f)]
            print("  matched calls whose %-7s DISAGREES: %d of %d (%.1f%%)"
                  % (f, len(diff), len(pairs), bts._rate(len(diff), len(pairs))))
            for ak, bk in diff[:min(a.limit, 6)]:
                print("      %s:%d  master %s=%s" % (ak[0], ak[1], f,
                                                     info_of["master"][ak].get(f, "-")))
                print("      %s  xTea2  %s=%s" % (" " * (len("%s:%d" % ak)), f,
                                                  info_of["xtea2"][bk].get(f, "-")))

        # "Only in arm X's stable set" is NOT the same as "arm Y never called it": arm Y may make
        # the call in some replicates and not others, in which case the difference is arm Y's
        # INSTABILITY rather than a disagreement about the locus. Measured, and it happens -- on L1,
        # chr17:6438313 and chr9:41047840 are called by master in 2 of 3 reps and by xTea2 in all 3,
        # so listing them as "xTea2-only" without this annotation would overstate the difference.
        for label, keys, arm, other in (("MASTER-ONLY", m_only, "master", "xtea2"),
                                        ("xTEA2-ONLY", x_only, "xtea2", "master")):
            other_sometimes = set()
            for r in reps:
                if (other, r) not in loaded:
                    continue
                ks = key_set(loaded[(other, r)])
                for k in keys:
                    if any(c == k[0] and abs(p - k[1]) <= a.window for c, p in ks):
                        other_sometimes.add(k)
            ont = [k for k in keys if on_truth(*k)]
            print("\n  %s stable calls: %d, of which %d are ON-TRUTH and %d off-truth"
                  % (label, len(keys), len(ont), len(keys) - len(ont)))
            print("    of those %d, the OTHER arm made %d in SOME (not all) of its replicates --"
                  % (len(keys), len(other_sometimes)))
            print("    those are the other arm's INSTABILITY, not a disagreement about the locus;")
            print("    the remaining %d it never called at all." % (len(keys) - len(other_sometimes)))
            print("    (an off-truth call the other arm lacks is a PRECISION difference; an")
            print("     on-truth one is a SENSITIVITY difference -- section C adjudicates those)")
            shown = ont + [k for k in keys if not on_truth(*k)]
            for k in shown[:a.limit]:
                mark = "  <- %s calls it in SOME reps" % other if k in other_sometimes else ""
                print("      %s%s" % (fmt(k, info_of[arm][k], fields, on_truth(*k)), mark))
            if len(shown) > a.limit:
                print("      ... %d more (raise --limit)" % (len(shown) - a.limit))

        # ---------- C. truth-anchored ----------
        print("\n--- C. TRUTH-ANCHORED: which truth records each arm recovers in ALL replicates ---")
        truth = bts.load_truth_family(prc.TRUTH_DIR, family)
        hits = {}
        for arm in ARMS:
            per = []
            for r in reps:
                if (arm, r) not in loaded:
                    continue
                calls = [bts.CallRecord(c, p, family, i.get("TD_SRC"), i.get("SUBTYPE"), i)
                         for c, p, i in loaded[(arm, r)]]
                m = bts.match_one_to_one(truth, calls, a.window)
                per.append({(truth[ti].chrom, truth[ti].pos) for ti in m})
            hits[arm] = set.intersection(*per) if per else set()
        m_win, x_win = sorted(hits["master"] - hits["xtea2"]), sorted(hits["xtea2"] - hits["master"])
        print("  truth records: %d total | master-all-reps %d | xTea2-all-reps %d"
              % (len(truth), len(hits["master"]), len(hits["xtea2"])))
        print("  RECOVERED BY MASTER IN EVERY REP AND BY xTea2 IN NONE: %d" % len(m_win))
        for k in m_win[:a.limit]:
            near = [(c, p) for c, p in stable["master"] if c == k[0] and abs(p - k[1]) <= a.window]
            src = info_of["master"][near[0]] if near else {}
            print("      truth %s:%-11d  master called it: %s"
                  % (k[0], k[1], " ".join("%s=%s" % (f, src.get(f, "-")) for f in fields)
                     if src else "(matched a non-stable call)"))
        pooled["master"]["truth_all"] += len(hits["master"])
        pooled["xtea2"]["truth_all"] += len(hits["xtea2"])
        pooled["master"]["wins"] += len(m_win)
        pooled["xtea2"]["wins"] += len(x_win)
        for arm in ARMS:
            pooled[arm]["stable"] += len(stable.get(arm, ()))
            sets_ = [key_set(loaded[(arm, r)]) for r in reps if (arm, r) in loaded]
            if len(sets_) > 1:
                pooled[arm]["unstable"] += len(set.union(*sets_) - set.intersection(*sets_))
        print("  RECOVERED BY xTea2 IN EVERY REP AND BY MASTER IN NONE: %d" % len(x_win))
        for k in x_win[:a.limit]:
            near = [(c, p) for c, p in stable["xtea2"] if c == k[0] and abs(p - k[1]) <= a.window]
            src = info_of["xtea2"][near[0]] if near else {}
            print("      truth %s:%-11d  xTea2 called it: %s"
                  % (k[0], k[1], " ".join("%s=%s" % (f, src.get(f, "-")) for f in fields)
                     if src else "(matched a non-stable call)"))
    if len(fams) > 1:
        print("\n" + "=" * 96)
        print("POOLED OVER %s -- THE REPRODUCIBLE CALL SETS ONLY" % "+".join(fams))
        print("=" * 96)
        for arm in ARMS:
            d = pooled[arm]
            print("  %-7s stable calls %5d   unstable %3d   truth records recovered in EVERY rep "
                  "%5d   loci the other arm never gets %3d"
                  % (arm, d["stable"], d["unstable"], d["truth_all"], d["wins"]))
        m, x = pooled["master"], pooled["xtea2"]
        if m["truth_all"] == x["truth_all"]:
            print("\n  THE TWO ARMS RECOVER THE SAME NUMBER OF TRUTH RECORDS REPRODUCIBLY: %d each."
                  % m["truth_all"])
            print("  They differ in WHICH: master gets %d xTea2 never does, xTea2 gets %d master"
                  % (m["wins"], x["wins"]))
            print("  never does. So on the calls each caller can actually repeat, SENSITIVITY IS")
            print("  EQUAL -- and master's apparent advantage in the per-replicate means comes")
            print("  entirely from calls it makes INCONSISTENTLY (%d unstable vs %d)."
                  % (m["unstable"], x["unstable"]))
            print("\n  STATE BOTH HALVES, because only one of them favours each arm. Master DOES")
            print("  recover more truth records PER RUN (probe_repro_comparison: mean 1324.7 vs")
            print("  1320.0) -- that part is real and is not an artifact. What it does not do is")
            print("  recover the SAME ones twice: only %d are common to all its replicates, and the"
                  % m["truth_all"])
            print("  extra few are a different few each time. Which framing matters depends on the")
            print("  use: a one-off call set gains from master's lottery tickets, a cohort or any")
            print("  comparison across samples is corrupted by them.")
        else:
            lead = "master" if m["truth_all"] > x["truth_all"] else "xTea2"
            print("\n  Reproducibly recovered truth records: master %d, xTea2 %d -- %s ahead by %d,"
                  % (m["truth_all"], x["truth_all"], lead, abs(m["truth_all"] - x["truth_all"])))
            print("  against master %d unstable calls and xTea2 %d."
                  % (m["unstable"], x["unstable"]))
    return 0


if __name__ == "__main__":
    sys.exit(main())
