#!/usr/bin/env python
"""Score a mosaic call set against the HapMap-mixture truth set, STRATIFIED BY VAF.

    python benchmark/benchmark_hapmap_mosaic.py --calls run/sample.vcf [--region chr1]

WHY A CURVE AND NOT A NUMBER
----------------------------
This truth set is a contrived mix of six known donors at known proportions, so every truth
locus has a *predicted allele fraction*, and they span 0.25%-8.25%. A single sensitivity
number over that range answers no question anyone has: it blends loci a 466x run can see
easily with loci it mathematically cannot. What the mode has to report is where detection
falls off, so the output here is sensitivity as a function of VAF.

It replaces SMHT005-3G, whose 15 sMEI loci were measured on 2026-08-12 to carry NO short-read
signal above matched background (plan sec 7.2). This set is better on every axis that matters:
2,996 tier1 somatic loci instead of 15, all three families instead of L1 only, and a known VAF
per locus instead of none.

VAF IS DERIVED, AND THE DERIVATION IS CHECKED AGAINST THE SHIPPED FILE
----------------------------------------------------------------------
The truth set ships `tier1_L1_VAFs.tsv` -- L1 only. But the per-donor columns plus
`raw_mixture_percents.txt` determine the VAF outright:

    vaf = sum(mixture_fraction[d] for d in donors carrying the insertion) / 2

(the /2 because a donor carrying it is heterozygous). `derive_vaf` implements that, and
`--self-check` asserts it reproduces all 422 shipped L1 values EXACTLY -- which it does. That
is what licenses using it for Alu and SVA, which ship no VAF file and would otherwise have no
stratification at all. Run the self-check before trusting an Alu or SVA curve; if the mixture
proportions are ever revised, it is the thing that will notice.

WHAT COUNTS AS A FALSE POSITIVE, which is where a mosaic benchmark usually goes wrong
-------------------------------------------------------------------------------------
Not every unmatched call is an error, and treating it as one would make this set punish a
caller for being right. Two exclusions ship WITH the truth set and both are applied:

  * `hg005_blacklist.bed` -- every insertion carried by HG005, the 83.5% background donor,
    across all three families. Those are real, and they are GERMLINE at ~42% VAF. A mosaic
    caller finding them is working correctly; counting them as false positives would make
    precision a measure of how much germline signal the run saw.
  * `masked_regions.txt` -- tier3-but-not-tier1 loci (its own generator script:
    `grep -f l1_...tier1.txt -v l1_...tier3.txt`). Lower-confidence truth. A call there is
    neither credited nor penalised, because the truth set does not claim to know.

Unmatched calls are therefore reported in three buckets -- germline (blacklist), ambiguous
(masked), and genuine FP -- rather than collapsed into one, so precision is interpretable
instead of merely small.

Reuses `read_calls_vcf` and `match_one_to_one` from benchmark_truth_set rather than
reimplementing them (plan sec 10): the matching rule -- one-to-one, greedy by distance -- is
the thing that must not differ between the germline and mosaic scorers, or their numbers stop
being comparable.
"""
import argparse
import os
import sys
from collections import defaultdict

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from benchmark_truth_set import read_calls_vcf, match_one_to_one  # noqa: E402

DEFAULT_TRUTH_DIR = ("/n/data1/hms/dbmi/park/corinne/hapmap_mei_calls/"
                     "hapmap_mix_calls/truth_set_10_18_2024")

# family -> the tiered mixture file basename stem
FAMILY_FILE = {"L1": "l1", "Alu": "alu", "SVA": "sva"}

# From raw_mixture_percents.txt. Parsed at load time rather than hardcoded (see load_mixture);
# this is the fallback and the thing --self-check validates against.
MIXTURE = {"HG005": 0.835, "HG02622": 0.10, "HG02486": 0.02,
           "HG02257": 0.02, "HG002": 0.02, "HG00438": 0.005}

# Bin edges as FRACTIONS, placed strictly BETWEEN the ladder's rungs.
#
# The mixture is discrete, so predicted VAF takes exactly 15 values across the three
# families -- it is not a continuum, and binning it as one invents structure that is not
# there. The 15 rungs, and the donor arithmetic that produces each:
#
#   0.0025 (227)  HG00438 alone            0.05   (544)  HG02622
#   0.01   (1242) one 2% donor             0.0525 (23)   HG02622 + HG00438
#   0.0125 (102)  2% + HG00438             0.06   (274)  HG02622 + one 2%
#   0.02   (195)  two 2% donors            0.0625 (64)   + HG00438
#   0.0225 (47)   two 2% + HG00438         0.07   (132)  HG02622 + two 2%
#   0.03   (18)   three 2% donors          0.0725 (43)   + HG00438
#   0.0325 (20)   three 2% + HG00438       0.08   (37)   HG02622 + three 2%
#                                          0.0825 (28)   everything
#
# The bins below group rungs that share a donor structure, and every edge sits in a gap.
# An edge ON a rung is not merely inelegant: these values come out of floating-point
# arithmetic, so 0.03 may be 0.030000000000000002, and which side of a 0.03 edge those 18
# loci land on would depend on summation order. test_bin_edges_do_not_split_a_mixture_rung
# pins this -- it caught exactly that, with the previous 0.03 edge.
VAF_BINS = [(0.0,   0.005, "0.25%"),      # HG00438 only -- below detection at <500x
            (0.005, 0.015, "1-1.25%"),    # a single 2% donor
            (0.015, 0.040, "2-3.25%"),    # multiple 2% donors
            (0.040, 0.055, "5-5.25%"),    # HG02622, the 10% donor
            (0.055, 0.090, "6-8.25%")]    # HG02622 plus others


class MosaicTruth:
    __slots__ = ("chrom", "pos", "family", "donors", "vaf")

    def __init__(self, chrom, pos, family, donors, vaf):
        self.chrom, self.pos, self.family, self.donors, self.vaf = chrom, pos, family, donors, vaf


def load_mixture(truth_dir):
    """raw_mixture_percents.txt -> {donor: fraction}. Lines look like '2%:  HG002 (M, EUR)'."""
    path = os.path.join(truth_dir, "raw_mixture_percents.txt")
    if not os.path.exists(path):
        return dict(MIXTURE)
    out = {}
    for line in open(path):
        if ":" not in line:
            continue
        pct, rest = line.split(":", 1)
        pct = pct.strip().rstrip("%")
        donor = rest.strip().split()[0] if rest.strip() else ""
        try:
            out[donor] = float(pct) / 100.0
        except ValueError:
            continue
    return out or dict(MIXTURE)


def derive_vaf(donor_cells, donor_names, mixture):
    """Predicted mixture VAF: half the summed fraction of the donors carrying the insertion.

    Half because a carrying donor is heterozygous. Validated to reproduce all 422 shipped
    tier1_L1_VAFs.tsv values exactly (--self-check), which is what licenses using it for Alu
    and SVA -- neither ships a VAF file.
    """
    return sum(mixture.get(d, 0.0) for d, cell in zip(donor_names, donor_cells)
               if cell not in (".", "")) / 2.0


def load_truth(truth_dir, family, tier=1, somatic_only=True, region=None):
    stem = FAMILY_FILE[family]
    path = os.path.join(truth_dir, f"{stem}_hapmapmixture_tier{tier}.txt")
    mixture = load_mixture(truth_dir)
    out = []
    with open(path) as fin:
        hdr = fin.readline().rstrip("\n").split("\t")
        donor_names = hdr[2:8]
        for line in fin:
            c = line.rstrip("\n").split("\t")
            if len(c) < 9:
                continue
            if somatic_only and c[8] != "somatic":
                continue
            if region and c[0] != region:
                continue
            out.append(MosaicTruth(c[0], int(c[1]), family, tuple(donor_names[i]
                                   for i, v in enumerate(c[2:8]) if v not in (".", "")),
                                   derive_vaf(c[2:8], donor_names, mixture)))
    return out


def _load_positions(path, chrom_col=0, pos_col=1):
    out = defaultdict(list)
    if not os.path.exists(path):
        return out
    for line in open(path):
        f = line.rstrip("\n").split("\t")
        if len(f) <= max(chrom_col, pos_col):
            continue
        try:
            out[f[chrom_col]].append(int(f[pos_col]))
        except ValueError:
            continue                # header row (generate_blacklist.sh does not skip it)
    for c in out:
        out[c].sort()
    return out


def load_exclusions(truth_dir):
    """(hg005 germline positions, masked/ambiguous positions), both {chrom: sorted[pos]}."""
    return (_load_positions(os.path.join(truth_dir, "hg005_blacklist.bed")),
            _load_positions(os.path.join(truth_dir, "masked_regions.txt")))


def _near(index, chrom, pos, window):
    import bisect
    arr = index.get(chrom)
    if not arr:
        return False
    i = bisect.bisect_left(arr, pos - window)
    return i < len(arr) and arr[i] <= pos + window


def score(calls, truth, hg005, masked, window):
    """-> (per-bin rows, overall dict, unmatched-call breakdown)."""
    matched = match_one_to_one(truth, calls, window)
    bins = []
    for lo, hi, label in VAF_BINS:
        idx = [i for i, t in enumerate(truth) if lo <= t.vaf < hi]
        tp = sum(1 for i in idx if i in matched)
        bins.append({"label": label, "n": len(idx), "tp": tp,
                     "sens": (100.0 * tp / len(idx)) if idx else float("nan")})
    tp = len(matched)
    unmatched = [c for i, c in enumerate(calls) if i not in {ci for ci, _ in matched.values()}]
    germline = [c for c in unmatched if _near(hg005, c.chrom, c.pos, window)]
    ambiguous = [c for c in unmatched
                 if c not in germline and _near(masked, c.chrom, c.pos, window)]
    fp = [c for c in unmatched if c not in germline and c not in ambiguous]
    overall = {"truth": len(truth), "calls": len(calls), "tp": tp, "fn": len(truth) - tp,
               "fp": len(fp),
               "sens": (100.0 * tp / len(truth)) if truth else float("nan"),
               "prec": (100.0 * tp / (tp + len(fp))) if (tp + len(fp)) else float("nan")}
    return bins, overall, {"germline": len(germline), "ambiguous": len(ambiguous),
                           "fp": len(fp)}


def self_check(truth_dir):
    """The derivation must reproduce the shipped L1 VAFs EXACTLY, or Alu/SVA are guesses."""
    given = {}
    path = os.path.join(truth_dir, "tier1_L1_VAFs.tsv")
    with open(path) as fin:
        fin.readline()
        for line in fin:
            k, v = line.rstrip("\n").split("\t")
            given[k] = float(v)
    derived = {f"{t.chrom}:{t.pos}": t.vaf for t in load_truth(truth_dir, "L1")}
    bad = [k for k in given if k not in derived or abs(derived[k] - given[k]) > 1e-9]
    print(f"self-check: {len(given)} shipped L1 VAFs, {len(derived)} derived, "
          f"{len(given) - len(bad)} exact matches, {len(bad)} mismatches")
    if bad:
        for k in bad[:5]:
            print(f"  MISMATCH {k}: shipped {given[k]} derived {derived.get(k)}")
        return 1
    print("  the VAF derivation is sound; Alu/SVA stratification is licensed")
    return 0


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--calls", help="xTea VCF to score")
    ap.add_argument("--truth-dir", default=DEFAULT_TRUTH_DIR)
    ap.add_argument("--families", default="L1,Alu,SVA")
    ap.add_argument("--window", type=int, default=200,
                    help="match window in bp (default 200, same as the germline scorer)")
    ap.add_argument("--region", default=None,
                    help="restrict truth to one contig, e.g. chr1 for the chr1 slice. "
                         "Without it a chr1-only run scores against genome-wide truth and "
                         "its sensitivity is meaningless.")
    ap.add_argument("--tier", type=int, default=1)
    ap.add_argument("--self-check", action="store_true",
                    help="verify the VAF derivation against the shipped L1 file and exit")
    a = ap.parse_args(argv)

    if a.self_check:
        return self_check(a.truth_dir)
    if not a.calls:
        ap.error("--calls is required (or use --self-check)")

    families = [f.strip() for f in a.families.split(",") if f.strip()]
    all_calls = read_calls_vcf(a.calls, families)
    hg005, masked = load_exclusions(a.truth_dir)

    print("=" * 78)
    print(f"xTea2 mosaic vs HapMap mixture, tier{a.tier} somatic"
          + (f", {a.region} only" if a.region else ", genome-wide"))
    print(f"calls: {a.calls}")
    print(f"match window: +/-{a.window} bp")
    print("=" * 78)

    for fam in families:
        truth = load_truth(a.truth_dir, fam, tier=a.tier, region=a.region)
        calls = [c for c in all_calls if c.family == fam
                 and (a.region is None or c.chrom == a.region)]
        bins, overall, unm = score(calls, truth, hg005, masked, a.window)
        print(f"\n--- {fam} ---")
        print(f"  truth {overall['truth']}  calls {overall['calls']}  "
              f"TP {overall['tp']}  FN {overall['fn']}  FP {overall['fp']}")
        print(f"  sensitivity {overall['sens']:.1f}%   precision {overall['prec']:.1f}%")
        print(f"  unmatched calls: {unm['germline']} germline (HG005 blacklist), "
              f"{unm['ambiguous']} ambiguous (masked), {unm['fp']} counted as FP")
        if not truth:
            continue
        print(f"  {'VAF':<12}{'truth':>7}{'TP':>7}{'sens':>9}")
        for b in bins:
            if b["n"] == 0:
                continue
            print(f"  {b['label']:<12}{b['n']:>7}{b['tp']:>7}{b['sens']:>8.1f}%")
    return 0


if __name__ == "__main__":
    sys.exit(main())
