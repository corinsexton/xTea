#!/usr/bin/env python3
"""PROBE, not a pipeline component: can truth transductions be derived by aligning the
truth ALTs to the GENOME instead of to the family consensus?

WHY THIS EXISTS. `benchmark/benchmark_truth_set.py derive-td` follows plan §7.1: align the
truth ALT to `consensus/LINE1.fa`, take the 3' residual, map it back, require adjacency to
a full-length reference L1. That yields a transduction denominator of **0** for L1, and the
reason is not a threshold — 11 of the 20 residuals ≥100 bp map inside a reference L1 and
every one is an older subfamily (L1PA6/7/8 at 4-7% divergence, L1MB/MD/ME at 16-25%).
`consensus/LINE1.fa` is 4 L1HS-centric records, so an older element's body partly fails to
align and surfaces as "residual". Loosening the criterion admits body, not transductions.

THE IDEA. Skip the consensus. Align the ALT straight to hg38 and read its segmentation:
an assembled MEI insertion decomposes into (a) segments landing on reference repeat copies
of its family -- the element body, whatever the subfamily, because the reference contains
every subfamily by construction -- and (b) at most one segment landing on unique sequence
somewhere else. A 3' transduction is exactly case (b) adjacent to a full-length reference
copy, and the subfamily-coverage problem cannot arise.

ANSWERED 2026-08-02, job 49482779 — NO, and that CONFIRMS the consensus route rather than
rescuing it. The genome explains the truth ALTs almost completely (median coverage 0.99;
161/196 L1 and 77/90 SVA above 0.95), and after that **195 of 197 L1 ALTs contain no unique
non-L1 segment of 100 bp or more at all**. One L1 candidate survives and it is not adjacent
to a full-length L1; SVA yields 10 candidates of which 3 are adjacent. So two independent
methods -- consensus residual and genome segmentation -- agree that the HG002 truth set
contains essentially no derivable 3' transductions: **L1 zero, SVA two or three.**

A plausible explanation, offered as a hypothesis and not as a measurement: the benchmark
was manually curated by verifying "a polyA tail and target-site duplication/deletion" per
insertion, and that selects for clean canonical TPRT products, which may be exactly the
insertions that do NOT carry a transduced 3' segment.

**Do not run a third variant of this.** The ALT route is exhausted on HG002; Stage 3's
transduction gate needs the paper's supplement or a different sample. Kept in the tree
because the negative result is the finding, and because the segmentation output is a
useful way to look at any assembly-derived truth set.

ONE BUG WORTH NOT REPEATING, found by the first run of this script: testing whether a
segment is "inside a repeat" at its START POSITION calls a 6 kb segment beginning 7 bp
before an annotated L1 unique, and the candidate list fills with 6 kb "transductions" that
are the element body. Overlap the whole reference span -- see repeat_overlap_frac.

    python docs/handoff/runs/probe_alt_vs_genome.py --family L1 \\
        --truth-dir .../test_data/HG002/truth --reference .../hg38_no_alt.fa \\
        --rep-lib-dir ~/xTea/rep_lib_annotation --workdir /n/scratch/.../probe

Submit it; it loads the hg38 bwa index. See docs/handoff/runs/derive_truth_td.sbatch for
the environment.
"""
import argparse
import os
import sys

# .../xTea/docs/handoff/runs/<this file> -- four levels up, not three.
REPO = os.path.dirname(os.path.dirname(os.path.dirname(
    os.path.dirname(os.path.abspath(__file__)))))
sys.path.insert(0, os.path.join(REPO, "benchmark"))
import benchmark_truth_set as bts  # noqa: E402

RMSK = {"L1": "LINE/{g}/{g}_L1.fa.out", "SVA": "SVA/{g}/{g}_SVA.out"}
FULL = {"L1": 5000, "SVA": 1000}


def repeat_overlap_frac(copies, chrom, beg, end):
    """What fraction of [beg, end) lies inside a same-family reference repeat.

    The whole span, not the start point. Element body overlaps ~fully; a genuine
    transduced block, which sits downstream of its source, overlaps ~not at all.
    """
    span = max(1, end - beg)
    hit = 0
    for (s, e) in copies.get(chrom, ()):
        if e < beg:
            continue
        if s >= end:
            break
        hit += min(e, end) - max(s, beg)
    return min(1.0, hit / float(span))


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--family", required=True, choices=["L1", "SVA"])
    p.add_argument("--truth-dir", required=True)
    p.add_argument("--reference", required=True)
    p.add_argument("--rep-lib-dir", required=True)
    p.add_argument("--genome", default="hg38")
    p.add_argument("--bwa", default="bwa")
    p.add_argument("--threads", type=int, default=8)
    p.add_argument("--workdir", required=True)
    p.add_argument("--min-segment", type=int, default=100,
                   help="ignore ALT segments shorter than this")
    p.add_argument("--adjacency", type=int, default=10000)
    p.add_argument("--max-repeat-overlap", type=float, default=0.2,
                   help="a candidate transduced segment may overlap same-family repeat by "
                        "at most this fraction of its reference span; above it the segment "
                        "is element body, not transduced sequence")
    a = p.parse_args()

    os.makedirs(a.workdir, exist_ok=True)
    truth = bts.load_truth_family(a.truth_dir, a.family)
    alt_fa = os.path.join(a.workdir, "alt.fa")
    alt_by_key = {}
    with open(alt_fa, "w") as fout:
        for r in truth:
            key = "{0}:{1}".format(r.chrom, r.pos)
            alt_by_key[key] = r.inserted_sequence()
            fout.write(">{0}\n{1}\n".format(key, alt_by_key[key]))
    print("[probe] {0} {1} truth ALTs -> genome".format(len(truth), a.family), flush=True)
    sam = bts.run_bwa(a.bwa, a.reference, alt_fa, a.threads,
                      os.path.join(a.workdir, "alt_vs_genome.sam"))
    aln = bts.all_alignments(sam)

    rmsk = os.path.join(a.rep_lib_dir, RMSK[a.family].format(g=a.genome))
    all_copies = bts.load_full_length_copies(rmsk, 0)
    full_copies = bts.load_full_length_copies(rmsk, FULL[a.family])
    print("[probe] {0} reference {1} intervals, {2} at >= {3} bp".format(
        sum(len(v) for v in all_copies.values()), a.family,
        sum(len(v) for v in full_copies.values()), FULL[a.family]), flush=True)

    n_unmapped = n_no_unique = 0
    candidates = []
    covfracs = []
    for r in truth:
        key = "{0}:{1}".format(r.chrom, r.pos)
        alt = alt_by_key[key]
        recs = aln.get(key)
        if not recs:
            n_unmapped += 1
            continue
        segs = []          # (alt_start, alt_end, chrom, pos, mapq, repeat_overlap_frac)
        covered = 0
        for flag, rn, ppos, mapq, cigar in recs:
            span = bts.query_span_in_alt(flag, cigar, len(alt))
            if span is None:
                continue
            covered += span[1] - span[0]
            if span[1] - span[0] < a.min_segment:
                continue
            # Overlap of the WHOLE reference span with same-family repeat, not containment
            # of its start point. A 6 kb segment whose start sits 7 bp before an annotated
            # L1 is the element BODY; a start-point test calls it unique and the candidate
            # list fills with 6 kb "transductions". Observed, on the first run of this probe.
            ref_len = sum(n for n, op in bts.parse_cigar(cigar) if op in "MD=XN")
            segs.append((span[0], span[1], rn, ppos, mapq,
                         repeat_overlap_frac(all_copies, rn, ppos, ppos + ref_len)))
        covfracs.append(covered / float(len(alt)) if alt else 0.0)
        # case (b): a long segment on UNIQUE sequence that is NOT inside a same-family
        # repeat, i.e. a candidate transduced block rather than element body.
        uniq = [s for s in segs if s[5] <= a.max_repeat_overlap and s[4] >= 30]
        if not uniq:
            n_no_unique += 1
            continue
        for s in uniq:
            hit = bts.adjacent_full_length_copy(full_copies, s[2], s[3], a.adjacency)
            candidates.append((key, s[1] - s[0], "{0}:{1}".format(s[2], s[3]), s[4],
                               "" if hit is None else "{0}-{1}".format(*hit), s[5]))

    covfracs.sort()
    print()
    print("[probe] ALT coverage by the genome: median {0:.2f}, >=0.95: {1}/{2}".format(
        covfracs[len(covfracs) // 2] if covfracs else 0,
        sum(1 for x in covfracs if x >= 0.95), len(covfracs)))
    print("[probe] {0} ALT(s) unmapped, {1} with no unique non-{2} segment >= {3} bp".format(
        n_unmapped, n_no_unique, a.family, a.min_segment))
    print("[probe] {0} candidate transduced segment(s); {1} adjacent to a full-length {2}"
          .format(len(candidates), sum(1 for c in candidates if c[4]), a.family))
    print()
    print("%-22s %-7s %-24s %-6s %-6s %s" % ("truth site", "seglen", "maps to", "mapq",
                                             "rpt%", "adjacent full-length copy"))
    for key, seglen, where, mapq, hit, ovl in sorted(candidates, key=lambda c: -c[1]):
        print("%-22s %-7d %-24s %-6d %-6.2f %s" % (key, seglen, where, mapq, ovl, hit or "-"))
    print()
    print("[probe] READ THIS AS: the last column is the transduction denominator this "
          "method would give. If it is still ~0, the ALT route is exhausted on HG002 and "
          "Stage 3 needs the paper's supplement or a different sample.")


if __name__ == "__main__":
    main()
