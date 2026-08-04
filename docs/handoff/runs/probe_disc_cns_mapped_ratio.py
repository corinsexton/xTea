#!/usr/bin/env python3
"""What `MIN_DISC_MAPPED_RATIO = 0.7` costs: the knob that produces every `no_cns_disc_reads` locus.

    python docs/handoff/runs/probe_disc_cns_mapped_ratio.py --run-dir <out>/<sample>

OFFLINE, SECONDS, NO BAM AND NO RERUN -- and that last part is the point. Gate [B2]'s WGS split
(job 49703322) showed that **146 of its 184 drops are an ABSENCE of evidence, not a threshold
rejection**: 113 L1 and 33 SVA loci where the third bwa pass produced no usable disc-to-consensus
alignment at all. The gate for that is `parse_disc_algnmt_consensus`'s two filters -- unmapped, or
mapped fraction below `global_values.MIN_DISC_MAPPED_RATIO` -- which makes that constant the
largest single knob between gate [B] and gate [C], by a factor of four over everything [B2] itself
decides. It has never been examined and it is not in plan §5's interrogate list.

**THE HANDOFF CHAIN SAID THIS NEEDED AN INSTRUMENTED RUN. IT DOES NOT.** The claim was that the
disc-to-consensus SAM "is deleted at the end of the step"; it is deleted by
`run_xtea.cleanup_successful_run`, i.e. at the end of the RUN and only when `--keep-intermediates`
was not passed. **Every gate under `docs/handoff/runs/` passes it**, so
`<fam>/tmp/transduction/temp_transduction_disc_cns.sam` is sitting in every one of their output
trees -- 1.4 MB per family at WGS. So this probe reads what already exists.

WHAT IT MEASURES, per family:

  * the record-level funnel: total alignments, unmapped, and how many survive the mapped-fraction
    predicate at 0.7 and at each alternative cutoff;
  * **the locus-level funnel, which is the one that matters** -- how many loci have at least one
    surviving record. A locus with zero survivors is exactly a `no_cns_disc_reads` row, so the
    difference between two cutoffs' locus counts is what loosening the knob would RECOVER;
  * a **reconciliation against the ledger**: picked sources minus loci-with-a-survivor at 0.7 must
    equal the ledger's `no_cns_disc_reads` count. Two derivations from different files -- a SAM and
    a TSV -- and if they disagree, this probe's model of the filter is what is wrong.

The predicate is **imported, not reimplemented**: `is_clipped_part_qualified_algnmt` is called on a
real `XClipDiscFilter`, whose constructor only stores its arguments. A probe with its own copy of a
cutoff can agree with itself while the step disagrees.

WHAT IT CANNOT SAY. Recovering a locus at gate [B2] is not the same as calling it: a recovered
locus still has to cluster on a side and clear the two-sided support cutoff, then survive gates
[C], [D] and [E]. So read the locus deltas as an UPPER BOUND on what the knob can reach, and note
that plan §4's FN attribution measured gate [B2] as attributing **zero** false negatives -- so the
expected payoff here is precision-shaped or nothing, not sensitivity.
"""
import argparse
import os
import sys

REPO_ROOT = os.path.dirname(os.path.dirname(os.path.dirname(
    os.path.dirname(os.path.abspath(__file__)))))
sys.path.insert(0, os.path.join(REPO_ROOT, "xtea"))
import pysam            # noqa: E402
import global_values    # noqa: E402
import x_clip_disc_filter  # noqa: E402

DISC_CNS_SAM = os.path.join("tmp", "transduction", "temp_transduction_disc_cns.sam")
LEDGER = ("candidate_disc_filtered_cns2.txt.all_non_sibling_td.txt.decision_ledger.tsv")
SOURCES = ("candidate_disc_filtered_cns2.txt.all_non_sibling_td.txt.unique_trsdct_disc_only")
FAMILIES = ("L1", "SVA")
R_ABSENCE = "no_cns_disc_reads"
DEFAULT_CUTOFFS = (0.5, 0.6, 0.65, 0.7, 0.8, 0.9)


def _filter():
    """A real `XClipDiscFilter`. Its constructor only stores its arguments, and the predicate
    under test reads nothing else, so nonexistent paths are honest here rather than a stub."""
    return x_clip_disc_filter.XClipDiscFilter("no_such_bam_list.txt", "/tmp/", 1,
                                              "no_such_ref.fa")


def scan_sam(sf_sam, l_cutoffs):
    """(n_all, n_unmapped, {cutoff: (n_records, set_of_loci)}, all_loci_seen, per_locus).

    `all_loci_seen` counts a locus as seen if ANY record names it, mapped or not -- so
    "loci with a record but no surviving record" is separable from "loci with no record".
    `per_locus` is {locus: [n_records, n_unmapped]}, which is what separates the two ways a
    locus can end up with no evidence -- see the "CEILING" line in the report.
    """
    xf = _filter()
    n_all = n_unmapped = 0
    m_by_cutoff = dict((c, [0, set()]) for c in l_cutoffs)
    s_seen = set()
    m_per_locus = {}
    samfile = pysam.AlignmentFile(sf_sam, "r")
    for algnmt in samfile.fetch(until_eof=True):
        n_all += 1
        l_fields = algnmt.query_name.split(global_values.SEPERATOR)
        if len(l_fields) < 3:
            continue
        locus = (l_fields[-3], l_fields[-2])
        s_seen.add(locus)
        l_cnt = m_per_locus.setdefault(locus, [0, 0])
        l_cnt[0] += 1
        if algnmt.is_unmapped:
            n_unmapped += 1
            l_cnt[1] += 1
            continue
        l_cigar = algnmt.cigar
        for c in l_cutoffs:
            b_ok, _n = xf.is_clipped_part_qualified_algnmt(l_cigar, c)
            if b_ok:
                m_by_cutoff[c][0] += 1
                m_by_cutoff[c][1].add(locus)
    samfile.close()
    return n_all, n_unmapped, m_by_cutoff, s_seen, m_per_locus


def count_ledger_absences(sf_ledger):
    n = 0
    with open(sf_ledger) as fin:
        for line in fin:
            if line.startswith("#"):
                continue
            l_f = line.rstrip("\n").split("\t")
            if len(l_f) >= 4 and l_f[3] == R_ABSENCE:
                n += 1
    return n


def main():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--run-dir", required=True, help="a finished run's SAMPLE directory")
    p.add_argument("--cutoffs", default=",".join(str(c) for c in DEFAULT_CUTOFFS))
    a = p.parse_args()

    l_cutoffs = sorted(float(s) for s in a.cutoffs.split(",") if s)
    f_ship = global_values.MIN_DISC_MAPPED_RATIO
    if f_ship not in l_cutoffs:
        l_cutoffs.append(f_ship)
        l_cutoffs.sort()

    print("=" * 78)
    print("MIN_DISC_MAPPED_RATIO -- the filter behind every `no_cns_disc_reads` locus")
    print("  shipping value: {0} (global_values.MIN_DISC_MAPPED_RATIO)".format(f_ship))
    print("  run dir       : {0}".format(a.run_dir))
    print("=" * 78)

    l_problems = []
    for fam in FAMILIES:
        sf_sam = os.path.join(a.run_dir, fam, DISC_CNS_SAM)
        print("\n--- {0} ---".format(fam))
        if not os.path.isfile(sf_sam):
            print("  no disc-to-consensus SAM at {0}".format(sf_sam))
            print("  (the run needs --keep-intermediates; every gate under docs/handoff/runs "
                  "passes it)")
            continue
        n_all, n_unmapped, m_by_cutoff, s_seen, m_per_locus = scan_sam(sf_sam, l_cutoffs)
        print("  {0} ({1:.1f} MB)".format(sf_sam, os.path.getsize(sf_sam) / 1e6))
        print("  records {0}   unmapped {1} ({2:.1f}%)   loci named by any record {3}".format(
            n_all, n_unmapped, 100.0 * n_unmapped / n_all if n_all else 0.0, len(s_seen)))

        # The domain: gate [B2] sees exactly the picked sources (prep_representative_pos is
        # lossless, one entry per (chrm, pos)), so that file is the denominator.
        sf_src = os.path.join(a.run_dir, fam, SOURCES)
        n_picked = None
        if os.path.isfile(sf_src):
            with open(sf_src) as fin:
                n_picked = sum(1 for line in fin if line.strip())
            print("  picked sources (gate [B2]'s domain): {0}".format(n_picked))

        print("  {0:<8} {1:>10} {2:>10} {3:>12} {4:>14}".format(
            "cutoff", "records", "kept %", "loci kept", "vs shipping"))
        n_ship_loci = len(m_by_cutoff[f_ship][1])
        for c in l_cutoffs:
            n_rec, s_loci = m_by_cutoff[c]
            s_mark = " <- SHIPPING" if c == f_ship else ""
            print("  {0:<8} {1:>10} {2:>9.1f}% {3:>12} {4:>+14d}{5}".format(
                c, n_rec, 100.0 * n_rec / n_all if n_all else 0.0, len(s_loci),
                len(s_loci) - n_ship_loci, s_mark))

        # A locus can fail for two different reasons and they are not the same finding: no record
        # at all reached it (the re-selection or bwa produced nothing), or records reached it and
        # every one failed the ratio. Only the second is what this constant decides.
        s_absent = s_seen - m_by_cutoff[f_ship][1]
        print("  loci with a record but NO survivor at {0}: {1}".format(f_ship, len(s_absent)))
        if n_picked is not None:
            print("  loci in the domain with no record at all: {0}".format(
                n_picked - len(s_seen)))

        # THE CEILING, and it is the number that decides whether this constant is worth
        # interrogating at all. A locus with no surviving record got there in one of two ways:
        # every record for it came back UNMAPPED from bwa, in which case no value of this cutoff
        # can ever recover it; or at least one record mapped and failed the ratio, in which case
        # a looser cutoff might. Only the second population is this constant's to decide.
        n_all_unmapped = sum(1 for k in s_absent if m_per_locus[k][1] == m_per_locus[k][0])
        print("  CEILING: of those {0}, {1} have EVERY record unmapped -- unreachable at any "
              "cutoff --".format(len(s_absent), n_all_unmapped))
        print("           and {0} have >=1 mapped record that failed the ratio, which is the "
              "most".format(len(s_absent) - n_all_unmapped))
        print("           this constant can EVER reach. Compare it to the `loci kept` deltas "
              "above.")

        # THE RECONCILIATION. Two derivations, different files.
        sf_led = os.path.join(a.run_dir, fam, LEDGER)
        if os.path.isfile(sf_led) and n_picked is not None:
            n_absence = count_ledger_absences(sf_led)
            n_derived = n_picked - n_ship_loci
            s_verdict = "OK" if n_absence == n_derived else "MISMATCH"
            print("  ledger `{0}` rows: {1}   this probe's derivation: {2}   {3}".format(
                R_ABSENCE, n_absence, n_derived, s_verdict))
            if s_verdict == "MISMATCH":
                l_problems.append(
                    "{0}: ledger says {1} absences, the SAM says {2}".format(
                        fam, n_absence, n_derived))
        else:
            print("  (no ledger or no source file -- reconciliation SKIPPED, not passed)")

    print("\n  reconciliation checks: {0}".format(
        "OK" if not l_problems else "FAILED -- " + "; ".join(l_problems)))
    return 1 if l_problems else 0


if __name__ == "__main__":
    sys.exit(main())
