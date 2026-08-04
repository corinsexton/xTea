#!/usr/bin/env python3
"""Join a run's transduction decision ledgers to the truth set: which gate dropped each miss.

DOWNSTREAM_OVERHAUL_PLAN.md §5 on the `reason` field: "it is what makes the benchmark
interpretable, because every missed truth call gets attributed to the specific gate that
dropped it." `xtea/transduction_evidence.py` writes those reasons
(`<out>.decision_ledger.tsv`, one row per dropped locus); this joins them to
`benchmark_truth_set.py`'s false negatives so the sentence above becomes a table.

    python benchmark/ledger_attribution.py \
        --run-dir /path/to/out/HG002.wgs \
        --calls   /path/to/out/HG002.wgs/HG002.wgs.vcf \
        --truth-dir /path/to/HG002/truth

READ THIS BEFORE QUOTING ANY NUMBER IT PRINTS
---------------------------------------------
**The ledger covers the TRANSDUCTION step's decisions only, and that step is not where most
calls are made.** `cns_filter` produces the candidate list; the transduction step *adds*
records at sites it can support (the 52 of BASELINE-0's 1,515 that have no cns row) and
*updates* existing ones. So a false negative with **no ledger row** was not "let through by
these gates" -- it never reached them, having been dropped upstream or never been a
candidate at all. That population is expected to dominate, and it is reported as its own
line (`no_ledger_row`) rather than folded into an "unattributed" bucket, because the two
statements are completely different.

**A ledger row near a TRUE POSITIVE is not a contradiction either.** The step can reject its
own candidate at a position that `cns_filter` already called, and it can reject a different
locus that happens to sit inside the match window. Both are printed (`tp_with_ledger_row`)
precisely so the FN table is not read as causal on its own.

**Alu has no ledger, by construction** -- no transduction flank, so the step never runs for
it (`family_params.Alu.has_flank is False`). An Alu FN can never be attributed here, and
that is a property of the biology, not a gap in the tool.

MATCHING. The ledger's positions are candidate positions, which can sit up to the match
window away from a truth position, so rows are joined by proximity with the same default
window the scorer uses (200 bp). Where several rows fall in one window, ALL are reported and
the per-FN attribution takes the EARLIEST phase -- a locus dropped at gate A never had disc
evidence produced for it, so the earliest gate is the one that determined the outcome.
"""
import argparse
import os
import sys
from collections import Counter, defaultdict

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from benchmark_truth_set import (  # noqa: E402
    DEFAULT_POS_WINDOW, load_truth_family, match_one_to_one, read_calls_vcf)

sys.path.insert(0, os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))),
                                "xtea"))
import transduction_evidence as te  # noqa: E402

LEDGER_SUFFIX = ".decision_ledger.tsv"
# Phase order, for "the earliest gate wins" -- the module's own names, sorted, are already
# A_clip < B_source < C_candidate < D_sanity < E_genotype.
NO_ROW = "no_ledger_row"


def find_ledgers(run_dir):
    """{family_dir_name: [ledger paths]} under a run's sample directory."""
    out = defaultdict(list)
    for dirpath, _dirnames, filenames in os.walk(run_dir):
        for name in filenames:
            if name.endswith(LEDGER_SUFFIX):
                fam = os.path.basename(dirpath)
                out[fam].append(os.path.join(dirpath, name))
    return out


def read_ledger(path):
    """[(chrom, pos, phase, reason, detail)] -- and validate the reason vocabulary.

    A reason the module does not know means the file was written by a different version, so
    every count downstream would be silently mixing vocabularies. Raise instead.
    """
    rows = []
    with open(path) as fin:
        for line in fin:
            if line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 4:
                continue
            if f[3] not in te.PHASE_OF_REASON:
                raise ValueError(
                    "{0}: reason {1!r} is not in transduction_evidence's vocabulary -- the "
                    "ledger and the module disagree; do not mix their counts".format(
                        path, f[3]))
            rows.append((f[0], int(f[1]), f[2], f[3], f[4] if len(f) > 4 else ""))
    return rows


def index_by_chrom(rows):
    by_chrom = defaultdict(list)
    for r in rows:
        by_chrom[r[0]].append(r)
    for c in by_chrom:
        by_chrom[c].sort(key=lambda r: r[1])
    return by_chrom


def rows_near(by_chrom, chrom, pos, window):
    return [r for r in by_chrom.get(chrom, ()) if abs(r[1] - pos) <= window]


def _pct(num, den):
    return 0.0 if not den else 100.0 * num / den


def attribute(rows):
    """The reason that determined the outcome: the earliest phase among the rows."""
    if not rows:
        return NO_ROW, ""
    best = min(rows, key=lambda r: (r[2], r[3]))
    return best[3], best[4]


def main():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--run-dir", required=True,
                   help="a run's SAMPLE directory (the one holding L1/ Alu/ SVA/ merged/)")
    p.add_argument("--calls", required=True, help="the run's xTea VCF (merged or per-family)")
    p.add_argument("--truth-dir", required=True)
    p.add_argument("--families", default="L1,Alu,SVA")
    p.add_argument("--chrom", default="",
                   help="restrict truth AND calls to one chromosome. Needed to read a "
                        "single-contig run (the truth sets are whole-genome, so every "
                        "off-contig truth record would otherwise count as a miss)")
    p.add_argument("--window", type=int, default=DEFAULT_POS_WINDOW)
    p.add_argument("--details", type=int, default=25,
                   help="print up to N attributed false negatives per family (0 to suppress)")
    a = p.parse_args()

    families = [f for f in a.families.split(",") if f]
    ledgers = find_ledgers(a.run_dir)
    all_calls = read_calls_vcf(a.calls)
    if a.chrom:
        all_calls = [c for c in all_calls if c.chrom == a.chrom]

    print("=" * 78)
    print("LEDGER ATTRIBUTION -- which named gate dropped each missed truth call")
    print("  run    : {0}".format(a.run_dir))
    print("  calls  : {0}".format(a.calls))
    print("  truth  : {0}".format(a.truth_dir))
    print("  window : {0} bp".format(a.window))
    print("=" * 78)
    print()
    print("The ledger covers the TRANSDUCTION step only. `{0}` therefore means the locus "
          "never".format(NO_ROW))
    print("reached these gates -- not that they let it through. See the module docstring.")

    grand = Counter()
    for fam in families:
        paths = ledgers.get(fam, [])
        rows = []
        for path in sorted(paths):
            rows.extend(read_ledger(path))
        by_chrom = index_by_chrom(rows)

        truth = load_truth_family(a.truth_dir, fam)
        if a.chrom:
            truth = [t for t in truth if t.chrom == a.chrom]
        calls = [c for c in all_calls if c.family == fam]
        matched = match_one_to_one(truth, calls, a.window)
        fn_idx = [i for i in range(len(truth)) if i not in matched]

        print()
        print("-" * 78)
        print("{0}: truth {1}, calls {2}, TP {3}, FN {4}   ledger rows {5} from {6} file(s)"
              .format(fam, len(truth), len(calls), len(matched), len(fn_idx), len(rows),
                      len(paths)))
        if not paths:
            print("  no ledger for this family -- the transduction step does not run for it "
                  "(no flank),")
            print("  so no FN here can be attributed to a gate. Expected for Alu.")
        counts = Counter()
        shown = 0
        for i in fn_idx:
            t = truth[i]
            near = rows_near(by_chrom, t.chrom, t.pos, a.window)
            reason, detail = attribute(near)
            counts[reason] += 1
            grand[(fam, reason)] += 1
            if reason != NO_ROW and shown < a.details:
                shown += 1
                print("  FN {0}:{1}  <- {2} ({3}) {4}".format(
                    t.chrom, t.pos, reason, te.PHASE_OF_REASON.get(reason, "?"), detail))
                if len(near) > 1:
                    print("      {0} ledger rows in the window: {1}".format(
                        len(near), ", ".join("{0}@{1}".format(r[3], r[1]) for r in near)))
        print("  false negatives by determining gate:")
        for reason, n in sorted(counts.items(), key=lambda kv: (-kv[1], kv[0])):
            phase = NO_ROW if reason == NO_ROW else te.PHASE_OF_REASON.get(reason, "?")
            print("    {0:<24} {1:>6}   ({2})".format(reason, n, phase))

        # THE COUNTER-CHECK, and it is not decoration: a dropped locus lands within the
        # window of a TRUE POSITIVE at some base rate, so the same thing happening at a false
        # negative is only evidence to the extent that it happens MORE often. Without this
        # line the attribution table reads as causal when part of it is coincidence.
        tp_with_row = 0
        for i, (_ci, _d) in matched.items():
            if rows_near(by_chrom, truth[i].chrom, truth[i].pos, a.window):
                tp_with_row += 1
        n_fn_attributed = sum(n for r, n in counts.items() if r != NO_ROW)
        print("  coincidence control:")
        print("    true positives with a dropped locus in the window: {0}/{1} = {2:.1f}%"
              .format(tp_with_row, len(matched), _pct(tp_with_row, len(matched))))
        print("    false negatives with one (i.e. attributed):        {0}/{1} = {2:.1f}%"
              .format(n_fn_attributed, len(fn_idx), _pct(n_fn_attributed, len(fn_idx))))
        if matched and fn_idx:
            expected = _pct(tp_with_row, len(matched)) / 100.0 * len(fn_idx)
            print("    at the TP base rate you would expect {0:.1f} of the {1} FN attributed "
                  "by chance;".format(expected, len(fn_idx)))
            if expected:
                print("    observed {0} -> enrichment {1:.2f}x, EXCESS over chance ~{2:.0f} "
                      "call(s)".format(n_fn_attributed, n_fn_attributed / expected,
                                       n_fn_attributed - expected))
                print("    Quote the EXCESS, not the raw count -- and with this N, an "
                      "enrichment near 1.0 means nothing.")
            else:
                # No TP anywhere near a dropped locus: no base rate to compare against, so an
                # enrichment ratio would be a division by zero dressed up as a finding.
                print("    observed {0}; no TP has a nearby drop, so there is no base rate to "
                      "compare against".format(n_fn_attributed))
        if rows:
            print("  the ledger's own totals, by phase: {0}".format(
                " ".join("{0}={1}".format(k, v) for k, v in sorted(
                    Counter(r[2] for r in rows).items()))))
            print("  the ledger's own totals, by reason: {0}".format(
                " ".join("{0}={1}".format(k, v) for k, v in sorted(
                    Counter(r[3] for r in rows).items()))))

    print()
    print("=" * 78)
    print("ALL FAMILIES -- false negatives by determining gate")
    per_reason = Counter()
    for (_fam, reason), n in grand.items():
        per_reason[reason] += n
    total = sum(per_reason.values())
    for reason, n in sorted(per_reason.items(), key=lambda kv: (-kv[1], kv[0])):
        phase = NO_ROW if reason == NO_ROW else te.PHASE_OF_REASON.get(reason, "?")
        print("  {0:<24} {1:>6}  {2:>5.1f}%   ({3})".format(
            reason, n, 100.0 * n / total if total else 0.0, phase))
    print("  {0:<24} {1:>6}".format("TOTAL FN", total))
    print("=" * 78)


if __name__ == "__main__":
    main()
