#!/usr/bin/env python3
"""Why the decoy flip's new SVA calls are all off-truth -- and what the +3 accounting missed.

    python docs/handoff/runs/probe_sva_flip_offtruth.py \
        --before <BASELINE-0.vcf> --after <arm SVA.vcf> --truth-dir <truth> \
        [--verdicts <verify_on/verdicts.tsv>]

THE QUESTION. The flip (`--td-decoy-structural` on by default, `3d1437a`) bought 3 SVA calls and
moved position precision 93.9% -> 90.6%. TP and FN did not move, so all three are off-truth. The
2026-08-05 handoff §6 item 7 filed the honest options as "correct transductions the MEI-position
truth set cannot represent, or wrong" -- and this is what can be decided offline.

WHAT IT CHECKS, and each line answers a different explanation:

 1. **Strict superset.** Is the after-set the before-set plus N, or did positions also move? A
    shifted call and a new call are different events and precision cannot tell them apart.
 2. **Is it a WINDOWING near-miss?** The scorer matches at +-200 bp. If a new call sits 300 bp from
    a truth record, "off-truth" is the window's verdict and not the caller's. So: distance to the
    nearest truth record of EVERY family on the same contig, not just its own.
 3. **Can the position be on-truth AT ALL?** A contig with zero truth records for that family has
    no matchable record, so a call there is off-truth by construction of the denominator. Reported
    per family so a blanket contig exclusion can be told from a family-specific absence.
 4. **Is the site on the same chromosome as its source?** An intra-chromosomal call near its own
    source is the classic false-positive mode -- the source copy's own reads mis-assigned locally.
    Inter-chromosomal calls cannot be explained that way.
 5. **What did the READS say?** The verifier's per-entry verdict, which is a different instrument
    from the truth set and can disagree with it without either being wrong.
 6. **AND THE ONE THE +3 ACCOUNTING MISSES.** A position present in BOTH arms whose `TD_SRC` went
    from `not_transduction` to an assertion is a real effect of the flip that moves NO precision
    number, because the position was already called. Counting only new positions hides it, and it
    is why the verification pool grew by 4 while the call set grew by 3.

WHAT IT CANNOT SAY: whether an off-truth call is a real insertion. Nothing offline can. What it can
do is remove the explanations that are checkable, and say which of the remaining ones the available
instruments are even able to address.
"""
import argparse
import os
import sys

REPO_ROOT = os.path.dirname(os.path.dirname(os.path.dirname(
    os.path.dirname(os.path.abspath(__file__)))))
sys.path.insert(0, os.path.join(REPO_ROOT, "benchmark"))
import benchmark_truth_set as bts  # noqa: E402

NOT_TD = "not_transduction"
FAMILIES = ("L1", "Alu", "SVA")


def load_verdicts(path):
    """{(chrom,pos): (verdict, detail, td_src)} from a `verify_transduction` verdicts TSV."""
    m = {}
    with open(path) as fin:
        for line in fin:
            if line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 17:
                continue
            try:
                i_pos = int(f[2])
            except ValueError:
                continue        # the column-name row, which does NOT start with '#'
            m[(f[1], i_pos)] = (f[15], f[16], f[17] if len(f) > 17 else "")
    return m


def nearest_truth(truth, chrom, pos):
    """(distance, n_records_on_contig) -- None distance when the contig has no record."""
    l_d = [abs(t.pos - pos) for t in truth if t.chrom == chrom]
    return (min(l_d) if l_d else None), len(l_d)


def main():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--before", required=True, help="the BEFORE VCF (BASELINE-0)")
    p.add_argument("--after", required=True, help="the AFTER VCF (the flipped arm)")
    p.add_argument("--truth-dir", required=True)
    p.add_argument("--verdicts", default=None,
                   help="verify_transduction's verdicts.tsv for the AFTER arm")
    p.add_argument("--family", default="SVA", choices=list(FAMILIES))
    p.add_argument("--window", type=int, default=200, help="the scorer's position window")
    a = p.parse_args()

    truth = {fam: bts.load_truth_family(a.truth_dir, fam) for fam in FAMILIES}
    before = bts.read_calls_vcf(a.before, families=[a.family])
    after = bts.read_calls_vcf(a.after, families=[a.family])
    m_before = {(c.chrom, c.pos): c for c in before}
    m_after = {(c.chrom, c.pos): c for c in after}
    m_verdicts = load_verdicts(a.verdicts) if a.verdicts else {}

    print("=" * 78)
    print("WHY THE FLIP'S NEW {0} CALLS ARE OFF-TRUTH".format(a.family))
    print("  before : {0}  ({1} {2} calls)".format(a.before, len(before), a.family))
    print("  after  : {0}  ({1} {2} calls)".format(a.after, len(after), a.family))
    print("  truth  : {0}  ({1} {2} records)".format(
        a.truth_dir, len(truth[a.family]), a.family))
    print("=" * 78)

    # --- 1. strict superset ---------------------------------------------------------------
    gained = sorted(set(m_after) - set(m_before))
    lost = sorted(set(m_before) - set(m_after))
    print("\n1. IS IT ADDITIVE?  gained {0}, lost {1}  -> {2}".format(
        len(gained), len(lost),
        "STRICT SUPERSET" if not lost else "NOT a superset -- positions also moved"))
    if lost:
        print("   lost: {0}".format(lost))

    # --- 6. the effect the +N accounting misses -------------------------------------------
    converted = []
    for key, c in m_after.items():
        if key not in m_before:
            continue
        s_b = m_before[key].td_src or NOT_TD
        s_a = c.td_src or NOT_TD
        if s_b != s_a:
            converted.append((key, s_b, s_a))
    print("\n2. POSITIONS PRESENT IN BOTH ARMS WHOSE TD_SRC CHANGED: {0}".format(len(converted)))
    print("   These move NO precision number -- the position was already called -- so counting")
    print("   only new positions understates the flip. This is where a pool delta bigger than")
    print("   the call delta comes from.")
    for key, s_b, s_a in sorted(converted):
        s_note = ""
        if s_b == NOT_TD and s_a != NOT_TD:
            s_note = "   <- GAINED a transduction assertion"
        elif s_a == NOT_TD:
            s_note = "   <- LOST one"
        print("   {0}:{1}\n     before {2}\n     after  {3}{4}".format(
            key[0], key[1], s_b, s_a, s_note))
        if key in m_verdicts:
            v, d, _ = m_verdicts[key]
            print("     reads: {0} ({1})".format(v, d))

    # --- 2..5, per gained call ------------------------------------------------------------
    print("\n3. THE GAINED CALLS, one block each")
    for key in gained:
        c = m_after[key]
        chrom, pos = key
        print("\n   --- {0}:{1} ---".format(chrom, pos))
        print("     SUBTYPE  {0}".format(c.subtype))
        print("     TD_SRC   {0}   [shape={1} call_class={2} source_class={3}]".format(
            c.td_src, bts.td_src_shape(c.td_src), bts.td_src_call_class(c.td_src),
            bts.td_src_source_class(c.td_src)))
        print("     support  LCLIP={0} RCLIP={1} LDISC={2} RDISC={3} LPOLYA={4} RPOLYA={5} "
              "AF={6}".format(c.info.get("LCLIP"), c.info.get("RCLIP"), c.info.get("LDISC"),
                              c.info.get("RDISC"), c.info.get("LPOLYA"),
                              c.info.get("RPOLYA"), c.info.get("AF")))

        # 4. inter- or intra-chromosomal, and how far from its own source
        parsed = bts.parse_td_src(c.td_src)
        if parsed and parsed[0] == "region":
            src_chrom, src_start = parsed[1], parsed[2]
            b_same = (src_chrom == chrom)
            print("     source   {0} -- {1}{2}".format(
                src_chrom, "SAME contig as the site" if b_same else "DIFFERENT contig",
                "; distance {0} bp".format(abs(src_start - pos)) if b_same else
                ", so this cannot be the source copy's own reads mis-assigned locally"))

        # 2 + 3. the truth set's ability to match this position at all
        for fam in FAMILIES:
            d, n = nearest_truth(truth[fam], chrom, pos)
            if n == 0:
                print("     truth {0:<4} NO records on {1} at all -> a {0} call here CANNOT be "
                      "on-truth".format(fam, chrom))
            else:
                print("     truth {0:<4} nearest {1:,} bp ({2} records on {3}) -> {4}".format(
                    fam, d, n, chrom,
                    "WITHIN the {0} bp window".format(a.window) if d <= a.window
                    else "far outside the {0} bp window".format(a.window)))

        # 5. what the reads said
        if key in m_verdicts:
            v, d, _ = m_verdicts[key]
            print("     reads    {0} ({1})".format(v, d))
        elif m_verdicts:
            print("     reads    not in the verdicts file")

    # --- the summary the handoff has to quote ----------------------------------------------
    if m_verdicts and gained:
        m_tally = {}
        for key in gained:
            v = m_verdicts.get(key, ("NOT_IN_POOL", "", ""))[0]
            m_tally[v] = m_tally.get(v, 0) + 1
        print("\n4. THE GAINED CALLS' READ VERDICTS: {0}".format(
            " ".join("{0}={1}".format(k, v) for k, v in sorted(m_tally.items()))))
        print("   Quote this split, not the pool's -- the pool contains entries at positions that")
        print("   were ALREADY called, and those are not part of the precision cost.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
