#!/usr/bin/env python3
"""Where the transduction step's loci go, and which drops the ledger does NOT explain.

Offline, seconds, no BAM: reads one finished run's own intermediates plus its decision ledgers
(`xtea/transduction_evidence.py`) and reconciles them by SET DIFFERENCE, not subtraction.

    python docs/handoff/runs/probe_decision_funnel.py --run-dir <out>/<sample>

WHY THIS EXISTS. The ledger attributes every drop that one of the five gate groups made. It
cannot attribute a locus that never entered a gate's DOMAIN -- and after the WGS gate (job
`49634290`) the funnel showed those populations are the large ones:

    L1  .slct 34,197 -> clip-selected 2,802 -> sources 159 -> cluster-consistent 24 -> ...
    SVA .slct 39,452 -> clip-selected 3,233 -> sources 105 -> cluster-consistent 56 -> ...

with gate A recording 20,332 L1 drops and gate B 676. Those do not add up to the attrition, and
the gaps are the point of this probe. It names three populations the ledger is silent about:

  1. **`.slct` sites that never appear in the clip parse at all** — no clipped read realigned to
     the flank library, so the locus produces no evidence and no gate ever sees it.
  2. **clip-selected sites with no qualifying disc alignment** — gate B's domain is the disc
     SAM's own map, so a site whose disc reads all failed the mapping-quality or
     mapped-fraction tests is absent from it rather than dropped by it.
  3. **loci dropped by `check_seprt_disc_consistency`**, between gates B and C. This is the one
     that is a THRESHOLD rather than an absence of data: it is where plan §5's `NEARBY_CLIP=50`
     ("arbitrary proxy -- replace with a measured clip-to-insertion distance distribution")
     actually lives, and it is the step's largest threshold-driven rejection after gate A.

(1) and (2) are "no evidence was produced", which is a different statement from "a threshold
rejected it" and should never be folded together with it. (3) is un-instrumented and is the
obvious next extension of the ledger -- with the caveat that
`check_seprt_disc_consistency` is SHARED with `cns_filter`, so instrumenting it touches a path
this refactor deliberately did not.

The file names are the step's own hard-coded ones; a rename upstream makes this probe report
zeros, so it fails loud on a missing input instead.
"""
import argparse
import os
import sys
from collections import Counter

REPO_ROOT = os.path.dirname(os.path.dirname(os.path.dirname(
    os.path.dirname(os.path.abspath(__file__)))))
sys.path.insert(0, os.path.join(REPO_ROOT, "xtea"))
import transduction_evidence as te  # noqa: E402

# The transduction step's output stem, per family directory.
STEM = "candidate_disc_filtered_cns2.txt.all_non_sibling_td.txt"
SLCT = "candidate_list_from_disc.txt.clip_sites_raw_disc.txt.slct"
# Its intermediates, in the order the chain writes them.
CLIP_SELECTED = STEM + ".unique_trsdct_half_clip"
DISC_SOURCES = STEM + ".unique_trsdct_disc_only"
DISC_ORPHANS = STEM + ".unique_trsdct_disc_orphan"
CLUSTER_OK = STEM + ".unique_trsdct_disc_only_half_clip_half_disc"
CANDIDATES = STEM + ".unique_trsdct_disc_only_half_clip_half_disc_polyA"
AFTER_SANITY = STEM + ".unique_trsdct_disc_only_half_clip_half_disc_polyA_after_filter"
NEW_RECORDS = STEM + ".new_sites"
LEDGER = STEM + ".decision_ledger.tsv"


def read_positions(path, required=True):
    """{(chrom, pos)} from any of the step's two-column-or-more site files."""
    if not os.path.isfile(path):
        if required:
            raise SystemExit("missing input: {0}\n(a rename upstream would make this probe "
                             "report zeros; it fails instead)".format(path))
        return set()
    out = set()
    with open(path) as fin:
        for line in fin:
            if line.startswith("#"):
                continue
            f = line.split()
            if len(f) < 2:
                continue
            try:
                out.add((f[0], int(f[1])))
            except ValueError:
                continue
    return out


def read_ledger_positions(path):
    """({(chrom,pos)} per phase, Counter of reasons)."""
    by_phase, reasons = {}, Counter()
    if not os.path.isfile(path):
        raise SystemExit("missing ledger: {0}".format(path))
    with open(path) as fin:
        for line in fin:
            if line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 4:
                continue
            by_phase.setdefault(f[2], set()).add((f[0], int(f[1])))
            reasons[f[3]] += 1
    return by_phase, reasons


def main():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--run-dir", required=True,
                   help="a finished run's SAMPLE directory (holding L1/ Alu/ SVA/)")
    p.add_argument("--families", default="L1,SVA",
                   help="flanked families only -- the step does not run for Alu")
    a = p.parse_args()

    print("=" * 78)
    print("THE TRANSDUCTION DECISION'S FUNNEL, and what the ledger does not explain")
    print("  run: {0}".format(a.run_dir))
    print("=" * 78)

    for fam in [f for f in a.families.split(",") if f]:
        d = os.path.join(a.run_dir, fam)
        if not os.path.isdir(d):
            print("\n{0}: no such directory, skipped".format(fam))
            continue
        slct = read_positions(os.path.join(d, SLCT))
        clip_ok = read_positions(os.path.join(d, CLIP_SELECTED))
        sources = read_positions(os.path.join(d, DISC_SOURCES))
        orphans = read_positions(os.path.join(d, DISC_ORPHANS))
        cluster_ok = read_positions(os.path.join(d, CLUSTER_OK))
        cands = read_positions(os.path.join(d, CANDIDATES))
        sanity_ok = read_positions(os.path.join(d, AFTER_SANITY))
        records = read_positions(os.path.join(d, NEW_RECORDS), required=False)
        phases, reasons = read_ledger_positions(os.path.join(d, LEDGER))

        a_drops = phases.get("A_clip", set())
        b_drops = phases.get("B_source", set())
        c_drops = phases.get("C_candidate", set())
        d_drops = phases.get("D_sanity", set())
        e_drops = phases.get("E_genotype", set())

        print("\n" + "-" * 78)
        print("{0}".format(fam))
        print("  {0:<44}{1:>8}".format("sites in .slct (the step's domain)", len(slct)))
        print("  {0:<44}{1:>8}   gate A dropped {2}".format(
            "-> clip-selected", len(clip_ok), len(a_drops)))
        print("  {0:<44}{1:>8}   gate B dropped {2}".format(
            "-> disc sources picked", len(sources), len(b_drops)))
        print("  {0:<44}{1:>8}   NOT INSTRUMENTED".format(
            "-> cluster-consistent", len(cluster_ok)))
        print("  {0:<44}{1:>8}   gate C dropped {2}".format(
            "-> candidates", len(cands), len(c_drops)))
        print("  {0:<44}{1:>8}   gate D dropped {2}".format(
            "-> after sanity", len(sanity_ok), len(d_drops)))
        print("  {0:<44}{1:>8}   gate E dropped {2}".format(
            "-> new records written", len(records), len(e_drops)))
        print("  {0:<44}{1:>8}".format("(orphan-labelled loci)", len(orphans)))

        # --- the reconciliation, by set difference ---
        # NOT A NESTED FUNNEL, and this is the thing to know before reading the numbers.
        # Gate A screens on CLIP evidence and gate B on DISC evidence, and both domains are
        # derived from the SAME `.slct` list independently: the clip gather writes
        # `all_disc_pos` over the whole `.slct`, and `collect_disc_reads_of_given_list` reads
        # its site-list argument into a variable it then never uses (verified in
        # `x_clip_disc_filter.py`) -- the step's own comment says so too. So gate B's domain is
        # NOT a subset of gate A's survivors, and subtracting one from the other (the obvious
        # reading) produces a number that means nothing. Each absence is measured against
        # `.slct`, and the INTERSECTION is reported separately, because that intersection is
        # the real domain of everything downstream of gate B.
        a_domain = a_drops | clip_ok
        b_domain = b_drops | sources
        no_clip_evidence = slct - a_domain
        no_disc_evidence = slct - b_domain
        both = clip_ok & b_domain
        disc_not_clip = b_domain - clip_ok
        cluster_dropped = sources - cluster_ok
        print("\n  WHAT THE LEDGER DOES NOT EXPLAIN -- two INDEPENDENT screens over .slct")
        print("  {0:<44}{1:>8}  {2:>5.1f}% of .slct".format(
            "1. no clip read realigned at all", len(no_clip_evidence),
            100.0 * len(no_clip_evidence) / len(slct) if slct else 0.0))
        print("  {0:<44}{1:>8}  {2:>5.1f}% of .slct".format(
            "2. no qualifying disc alignment", len(no_disc_evidence),
            100.0 * len(no_disc_evidence) / len(slct) if slct else 0.0))
        print("  {0:<44}{1:>8}  {2:>5.1f}% of picked sources".format(
            "3. check_seprt_disc_consistency", len(cluster_dropped),
            100.0 * len(cluster_dropped) / len(sources) if sources else 0.0))
        print("     (3) is the only one of the three that is a THRESHOLD rather than an")
        print("     absence of evidence -- it is where plan sec 5's NEARBY_CLIP=50 lives.")
        print("\n  HOW THE TWO SCREENS INTERSECT (the later gates' real domain)")
        print("  {0:<44}{1:>8}".format("clip-selected (gate A survivors)", len(clip_ok)))
        print("  {0:<44}{1:>8}".format("in gate B's domain (had disc evidence)", len(b_domain)))
        print("  {0:<44}{1:>8}".format("BOTH -- clip-selected AND disc evidence", len(both)))
        print("  {0:<44}{1:>8}".format("disc evidence but NOT clip-selected", len(disc_not_clip)))
        print("     A locus in the second-to-last row can still reach gate C and be dropped")
        print("     there as `no_clip_evidence`, which is why that reason exists at all.")

        # Consistency checks: if any fails, a domain assumption above is wrong and none of the
        # numbers mean what the labels say.
        problems = []
        if not a_drops.isdisjoint(clip_ok):
            problems.append("gate A's drops overlap the clip-selected set")
        if not b_drops.isdisjoint(sources):
            problems.append("gate B's drops overlap the picked sources")
        if not cluster_ok <= sources:
            problems.append("cluster-consistent is NOT a subset of the picked sources -- the "
                            "containment call_candidates relies on is broken")
        if not a_domain <= slct:
            problems.append("gate A saw loci absent from .slct")
        if not b_domain <= slct:
            problems.append("gate B saw loci absent from .slct")
        print("\n  reconciliation checks: {0}".format(
            "OK" if not problems else "FAILED -- " + "; ".join(problems)))

        print("  ledger by reason: {0}".format(
            " ".join("{0}={1}".format(k, reasons[k]) for k in sorted(reasons))))
        unknown = [r for r in reasons if r not in te.PHASE_OF_REASON]
        if unknown:
            raise SystemExit("ledger has reasons the module does not know: {0}".format(unknown))


if __name__ == "__main__":
    main()
