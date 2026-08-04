#!/usr/bin/env python3
"""What `NEARBY_CLIP = 50` decides in gate [B2], swept offline over a finished run's own SAM.

    python docs/handoff/runs/probe_nearby_clip.py --run-dir <out>/<sample>

Plan §5 lists `NEARBY_CLIP=50` as an "arbitrary proxy -- replace with a measured clip-to-insertion
distance distribution", and the handoff chain filed it as needing an instrumented run because the
disc-to-consensus SAM was believed deleted. **It is not** (see `probe_disc_cns_mapped_ratio.py`), and
the two quantities the constant compares are both in that SAM's read names:

    clip_pos    = `prep_representative_pos(m_disc_transd)[chrm][pos]`, which is the IDENTITY --
                  `m_represent[chrm][ins_pos] = ins_pos`, x_transduction.py:710. So clip_pos IS the
                  insertion position, despite the name.
    anchor_pos  = the read name's `[-4]` field, which `parse_disc_algnmt_consensus` stores verbatim.

So `abs(anchor_pos - clip_pos) < NEARBY_CLIP` -- the test that decides whether a disc read that is
ALSO clipped contributes its side from the clip orientation -- is computable from the SAM alone, and
this probe sweeps it by driving the **shipping** gate at each value:

  * `parse_disc_algnmt_consensus` builds `m_disc` (the real one, imported);
  * `m_represent_pos` is rebuilt from the run's own picked-source file, exactly as the step does;
  * `check_seprt_disc_consistency` runs at each `NEARBY_CLIP`, with a real `DecisionLedger`.

**CHECK 0 IS WHAT MAKES THE SWEEP MEANINGFUL**: at the shipping value the reason split must equal
the run's own ledger B2 rows. If it does not, the arguments below are wrong -- `i_concord_dist` and
`n_half_disc` are both derived, and `n_half_disc` in particular is NOT `--nfdisc`/2 (the step's
cutoffs come from `Parameters.get_par_by_cov`; 4 → 2.0 at HG002's 31.8x). Treat a mismatch as this
probe being misconfigured rather than as a finding.

WHAT IT CANNOT SAY: whether a recovered locus becomes a CALL. Gate [B2]'s survivors still face [C],
[D] and [E], and the FN attribution measured [B2] as attributing zero false negatives (job
`49703322`), so expect a precision-shaped effect or none.
"""
import argparse
import os
import sys

REPO_ROOT = os.path.dirname(os.path.dirname(os.path.dirname(
    os.path.dirname(os.path.abspath(__file__)))))
sys.path.insert(0, os.path.join(REPO_ROOT, "xtea"))
import global_values          # noqa: E402
import transduction_evidence  # noqa: E402
import x_clip_disc_filter     # noqa: E402

DISC_CNS_SAM = os.path.join("tmp", "transduction", "temp_transduction_disc_cns.sam")
STEM = "candidate_disc_filtered_cns2.txt.all_non_sibling_td.txt"
SOURCES = STEM + ".unique_trsdct_disc_only"
LEDGER = STEM + ".decision_ledger.tsv"
FAMILIES = ("L1", "SVA")
B2_PHASE = "B2_cns_cluster"
DEFAULT_SWEEP = (0, 10, 25, 50, 100, 200, 500, 10 ** 9)


def load_loci(sf_sources):
    """m_represent_pos, as `prep_representative_pos` builds it: {chrm: {pos: pos}}."""
    m = {}
    with open(sf_sources) as fin:
        for line in fin:
            l_f = line.split()
            if len(l_f) < 2:
                continue
            m.setdefault(l_f[0], {})[int(l_f[1])] = int(l_f[1])
    return m


def ledger_b2_split(sf_ledger):
    m = {}
    with open(sf_ledger) as fin:
        for line in fin:
            if line.startswith("#"):
                continue
            l_f = line.rstrip("\n").split("\t")
            if len(l_f) >= 4 and l_f[2] == B2_PHASE:
                m[l_f[3]] = m.get(l_f[3], 0) + 1
    return m


def main():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--run-dir", required=True)
    p.add_argument("--families", default="L1,SVA")
    p.add_argument("--sweep", default=",".join(str(v) for v in DEFAULT_SWEEP))
    # Both derived, both overridable, and check 0 is what says they are right. i_concord_dist is
    # `max(550, int(mean_is + 3*std))` -- 954 for HG002 WGS from its own basic_cov_is_rlth_info.
    p.add_argument("--concord-dist", type=int, default=954)
    p.add_argument("--concord-ratio", type=float, default=0.25)
    p.add_argument("--n-half-disc", type=float, default=2.0,
                   help="the step's n_disc_cutoff/2 at the AUTOMATIC coverage-derived cutoff "
                        "(4 at HG002's 31.8x), not --nfdisc/2")
    a = p.parse_args()

    l_sweep = sorted(int(s) for s in a.sweep.split(",") if s)
    n_ship = global_values.NEARBY_CLIP
    if n_ship not in l_sweep:
        l_sweep.append(n_ship)
        l_sweep.sort()

    print("=" * 78)
    print("NEARBY_CLIP -- swept through the SHIPPING gate [B2] on a finished run's own SAM")
    print("  shipping value: {0} (global_values.NEARBY_CLIP)".format(n_ship))
    print("  gate args     : concord_dist {0}, ratio {1}, n_half_disc {2}".format(
        a.concord_dist, a.concord_ratio, a.n_half_disc))
    print("  run dir       : {0}".format(a.run_dir))
    print("=" * 78)

    l_problems = []
    xf = x_clip_disc_filter.XClipDiscFilter("no_such_bam_list.txt", "/tmp/", 1, "no_such_ref.fa")
    for fam in [s for s in a.families.split(",") if s]:
        sf_sam = os.path.join(a.run_dir, fam, DISC_CNS_SAM)
        sf_src = os.path.join(a.run_dir, fam, SOURCES)
        print("\n--- {0} ---".format(fam))
        if not (os.path.isfile(sf_sam) and os.path.isfile(sf_src)):
            print("  missing the SAM or the source file (needs --keep-intermediates)")
            continue
        m_represent = load_loci(sf_src)
        n_domain = sum(len(v) for v in m_represent.values())
        m_disc, _m_sample, _m_polya = xf.parse_disc_algnmt_consensus(
            sf_sam, global_values.MIN_DISC_MAPPED_RATIO)
        print("  gate [B2] domain (picked sources): {0}".format(n_domain))

        # The distance distribution plan §5 asks for, over the records the test can reach at all:
        # a disc read that is clipped on exactly one side. The others never consult the constant.
        l_dist = []
        n_records = n_one_sided_clip = 0
        for chrm in m_disc:
            for pos in m_disc[chrm]:
                if chrm not in m_represent or pos not in m_represent[chrm]:
                    continue
                for (anchor_pos, b_l, b_r, _cns, _same, _rc) in m_disc[chrm][pos]:
                    n_records += 1
                    if b_l != b_r:                      # exactly one side clipped
                        n_one_sided_clip += 1
                        l_dist.append(abs(anchor_pos - m_represent[chrm][pos]))
        l_dist.sort()
        print("  records in the domain: {0}   of which clipped on exactly one side: {1} "
              "({2:.1f}%)".format(n_records, n_one_sided_clip,
                                  100.0 * n_one_sided_clip / n_records if n_records else 0.0))
        if l_dist:
            def _q(f):
                return l_dist[min(len(l_dist) - 1, int(f * len(l_dist)))]
            print("  |anchor_pos - ins_pos| over those: min {0} p25 {1} median {2} p75 {3} "
                  "p90 {4} max {5}".format(l_dist[0], _q(0.25), _q(0.5), _q(0.75), _q(0.9),
                                           l_dist[-1]))
            print("  under the shipping {0}: {1} of {2} ({3:.1f}%)".format(
                n_ship, sum(1 for d in l_dist if d < n_ship), len(l_dist),
                100.0 * sum(1 for d in l_dist if d < n_ship) / len(l_dist)))

        print("  {0:<12} {1:>10} {2:>10}   {3}".format("NEARBY_CLIP", "survivors", "vs ship",
                                                       "reasons"))
        n_ship_survivors = None
        m_ship_reasons = None
        for v in l_sweep:
            global_values.NEARBY_CLIP = v
            params = transduction_evidence.TransductionParams.for_cns_cluster(a.n_half_disc)
            ledger = transduction_evidence.DecisionLedger()
            m_kept = xf.check_seprt_disc_consistency(
                m_represent, m_disc, a.concord_dist, a.concord_ratio, a.n_half_disc, False,
                ledger=ledger)
            del params
            n_kept = sum(len(x) for x in m_kept.values())
            m_reasons = {}
            for row in ledger.rows:
                m_reasons[row[3]] = m_reasons.get(row[3], 0) + 1
            if v == n_ship:
                n_ship_survivors, m_ship_reasons = n_kept, m_reasons
            print("  {0:<12} {1:>10} {2:>+10}   {3}".format(
                "INF" if v >= 10 ** 9 else v, n_kept,
                0 if n_ship_survivors is None else n_kept - n_ship_survivors,
                " ".join("{0}={1}".format(k, n) for k, n in sorted(m_reasons.items()))))
        global_values.NEARBY_CLIP = n_ship

        # CHECK 0: the shipping value must reproduce the run's own ledger.
        sf_led = os.path.join(a.run_dir, fam, LEDGER)
        if os.path.isfile(sf_led) and m_ship_reasons is not None:
            m_run = ledger_b2_split(sf_led)
            b_ok = m_run == m_ship_reasons
            print("  check 0: the run's ledger B2 rows {0}".format(
                " ".join("{0}={1}".format(k, v) for k, v in sorted(m_run.items()))))
            print("           this probe at NEARBY_CLIP={0}: {1}".format(
                n_ship, "MATCH" if b_ok else "MISMATCH -- the gate arguments are wrong"))
            if not b_ok:
                l_problems.append("{0}: probe {1} vs ledger {2}".format(
                    fam, m_ship_reasons, m_run))
        else:
            print("  check 0: no ledger -- SKIPPED, not passed")

    print("\n  check 0 overall: {0}".format(
        "OK" if not l_problems else "FAILED -- " + "; ".join(l_problems)))
    return 1 if l_problems else 0


if __name__ == "__main__":
    sys.exit(main())
