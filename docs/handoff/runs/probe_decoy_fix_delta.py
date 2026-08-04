#!/usr/bin/env python3
"""What the SVA decoy fix will do to gate [B], derived from a PRE-FIX run's own SAM.

    python docs/handoff/runs/probe_decoy_fix_delta.py --run-dir <out>/<sample> [--family SVA]

WHY THIS CAN PREDICT A RUN IT DID NOT PERFORM. `--td-decoy-structural` changes one thing: how
`parse_transduction_disc_algnmt` classifies a hit's reference name. It does not change the bwa
reference (`construct_novel_flanks` concatenates the family flank library with the polymorphic
flanks, and neither depends on the flag), it does not change the reads that reach that pass, and
it runs strictly after the alignment. **So a pre-fix run's `temp_transduction_disc.sam` is
byte-identical to the one the fixed run will produce**, and replaying the shipping parse over it
with `b_decoy_structural` set both ways is the fix's exact effect on gate [B] -- from a different
source than the run, which is what makes it a prediction rather than a restatement.

That SAM survives whenever `--keep-intermediates` was passed, which every gate under
`docs/handoff/runs/` does (772 MB for SVA at WGS, 393 MB for L1). It is not deleted at the end of
the step; `run_xtea.cleanup_successful_run` deletes it at the end of the RUN.

It drives the SHIPPING method -- `XTransduction.parse_transduction_disc_algnmt` -- at the arguments
the STEP ran with, which are NOT the CLI defaults (see the constants below: `n_disc_cutoff` 4, not
5, because `b_automatic` derives it from coverage). The OFF arm reproducing the run's own
`unique_trsdct_disc_only` and its gate-[B] ledger split is the check that they are right; treat a
mismatch there as this probe being misconfigured, not as a finding. It reports, per arm:

  * picked-source loci, split into REFERENCE copies and POLYMORPHIC ones. `chrom:start-end~rc`
    with start == end is polymorphic (`gnrt_flank_regions_of_polymerphic_insertions` writes
    `ins_pos` for both), so the split needs no annotation;
  * orphan loci;
  * the gate-[B] ledger's populations, so "the fix adds N loci to the decision" is separable from
    "the fix adds N picked sources".

WHAT IT CANNOT SAY: everything after gate [B]. A locus that gains a source still has to clear
[B2]'s cluster consistency, [C]'s conjunction, [D] and [E]. So the picked-source delta here is the
INPUT to the rest of the cascade, not a call count.
"""
import argparse
import os
import sys

REPO_ROOT = os.path.dirname(os.path.dirname(os.path.dirname(
    os.path.dirname(os.path.abspath(__file__)))))
sys.path.insert(0, os.path.join(REPO_ROOT, "xtea"))
import global_values          # noqa: E402
import transduction_evidence  # noqa: E402
import x_transduction         # noqa: E402

DISC_SAM = os.path.join("tmp", "transduction", "temp_transduction_disc.sam")
# THE PRODUCTION ARGUMENTS ARE NOT THE CLI DEFAULTS, and getting this wrong is how the first
# version of this probe reproduced 94 of the run's 105 picked sources. `x_TEA_main` sets
# `b_automatic = True` unless `--user`, which nothing passes, so it OVERRIDES `--nfclip/--nfdisc`
# with `Parameters.get_par_by_cov(coverage)`. HG002 WGS is 31.816x -> `m_cov_par[30] = (3, 4, 1)`,
# so the step runs at **n_clip_cutoff 3 and n_disc_cutoff 4**, i.e. `n_half_disc` **2.0**, not the
# 2.5 that `--nfdisc 5` implies. `run_xtea.py`'s locus-universe docstring says the same thing
# ("wrong by one at HG002 (4, not 5)"). The default here is the automatic value; --n-disc
# overrides it, and the OFF arm reproducing the run's own numbers is the check that it is right.
FLANK_LEN = 3000
N_DISC_DEFAULT = 4
F_MIN_UNIQ = 0.5
RTYPE = {"L1": 1, "SVA": 4}


def classify(m_picked):
    """(n_reference_copy, n_polymorphic) over {chrm: {pos: (source, hits, n)}}."""
    n_ref = n_poly = 0
    for chrm in m_picked:
        for pos in m_picked[chrm]:
            s_source = m_picked[chrm][pos][0]
            s_id = s_source.split(global_values.SEPERATOR)[0]
            try:
                _c, s_span = s_id.rsplit(":", 1)
                s_start, s_end = s_span.split("-")
                if int(s_start) == int(s_end):
                    n_poly += 1
                else:
                    n_ref += 1
            except ValueError:
                pass
    return n_ref, n_poly


def run_arm(sf_sam, fam, b_structural, s_workdir, n_disc):
    xtd = x_transduction.XTransduction(s_workdir, 1, "no_such_ref.fa")
    xtd.b_decoy_structural = b_structural
    xtd.f_min_uniq_ratio = F_MIN_UNIQ
    i_rep_type = RTYPE[fam]
    params = transduction_evidence.TransductionParams.for_disc_sources(
        n_disc / 2.0, F_MIN_UNIQ, i_rep_type)
    ledger = transduction_evidence.DecisionLedger()
    m_picked, m_orphan = xtd.parse_transduction_disc_algnmt(
        sf_sam, FLANK_LEN, n_disc / 2.0, i_rep_type, params=params, ledger=ledger)
    n_loci = sum(len(v) for v in m_picked.values())
    n_orphan = sum(len(v) for v in m_orphan.values())
    m_reasons = {}
    for row in ledger.rows:
        m_reasons[row[3]] = m_reasons.get(row[3], 0) + 1
    return n_loci, classify(m_picked), n_orphan, m_reasons


def main():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--run-dir", required=True, help="a finished run's SAMPLE directory")
    p.add_argument("--families", default="SVA,L1")
    p.add_argument("--workdir", default="/tmp/")
    p.add_argument("--n-disc", type=float, default=N_DISC_DEFAULT,
                   help="n_disc_cutoff as the STEP saw it -- the automatic, coverage-derived "
                        "value (4 at HG002's 31.8x), not --nfdisc")
    a = p.parse_args()

    print("=" * 78)
    print("THE SVA DECOY FIX's EFFECT ON GATE [B], replayed over a pre-fix run's own SAM")
    print("  run dir: {0}".format(a.run_dir))
    print("  args   : flklen {0}, n_disc_cutoff {1} -> n_half_disc {2}, f_min_uniq_ratio {3}"
          .format(FLANK_LEN, a.n_disc, a.n_disc / 2.0, F_MIN_UNIQ))
    print("  (n_disc_cutoff is the AUTOMATIC coverage-derived value, not --nfdisc; see the top"
          " of this file)")
    print("=" * 78)

    for fam in [s for s in a.families.split(",") if s]:
        sf_sam = os.path.join(a.run_dir, fam, DISC_SAM)
        print("\n--- {0} ---".format(fam))
        if not os.path.isfile(sf_sam):
            print("  no disc-to-flank SAM at {0} (needs --keep-intermediates)".format(sf_sam))
            continue
        print("  {0} ({1:.0f} MB)".format(sf_sam, os.path.getsize(sf_sam) / 1e6))
        m_arms = {}
        for (s_arm, b_structural) in (("OFF (shipping)", False), ("ON  (--td-decoy-structural)",
                                                                 True)):
            n_loci, (n_ref, n_poly), n_orphan, m_reasons = run_arm(
                sf_sam, fam, b_structural, a.workdir, a.n_disc)
            m_arms[b_structural] = (n_loci, n_ref, n_poly, n_orphan, m_reasons)
            print("  {0:<30} picked {1:>5}  (reference {2:>5}, polymorphic {3:>5})  "
                  "orphan {4:>4}".format(s_arm, n_loci, n_ref, n_poly, n_orphan))
            print("  {0:<30} gate [B] ledger: {1}".format("", " ".join(
                "{0}={1}".format(k, v) for k, v in sorted(m_reasons.items()))))
        off, on = m_arms[False], m_arms[True]
        print("  DELTA: picked {0:+d}   reference copies {1:+d}   polymorphic {2:+d}   "
              "orphan {3:+d}".format(on[0] - off[0], on[1] - off[1], on[2] - off[2],
                                     on[3] - off[3]))
        n_rows_off = sum(off[4].values())
        n_rows_on = sum(on[4].values())
        print("         gate [B] ledger rows {0} -> {1} ({2:+d}) -- loci that now REACH the"
              " decision".format(n_rows_off, n_rows_on, n_rows_on - n_rows_off))


if __name__ == "__main__":
    main()
