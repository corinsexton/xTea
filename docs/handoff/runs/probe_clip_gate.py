#!/usr/bin/env python3
"""What gate [A] decides, swept offline through the SHIPPING gate on a finished run's own clip SAM.

    python docs/handoff/runs/probe_clip_gate.py --run-dir <out>/<sample>

Gate [A] does **97% of the rejecting** in the transduction step (plan §4 Stage 3's WGS table:
20,332 L1 and 27,590 SVA drops, against 1-2 each at gates D and E), and plan §5's semantic half
has been filed against it through four handoffs without a measurement. This is the measurement,
and it is offline: `tmp/transduction/temp_transduction_clip.sam` survives a finished run, so gate
[A]'s entire input -- the clip realignments AND the polyA map, which the same parse builds -- is
on disk.

WHAT IS DRIVEN, AND WHAT IS NOT REIMPLEMENTED. `XTransduction.parse_transduction_clip_algnmt` is
called as the step calls it (`x_transduction.py:281`) and `transduction_evidence.select_clip_sites`
is the gate itself (`:283`). Nothing here re-derives a predicate -- see plan §10's
`probe_leftmost_cns` entry for the 2.6x-wrong numbers a private copy of a shipping rule produced.

FOUR THINGS IT MEASURES

 1. **The clip-support sweep.** `n_clip_cutoff` is coverage-derived (`Parameters.get_par_by_cov`),
    so it is not a constant anyone chose per sample; the sweep says what the step's population
    would be at other values. `n_polya_cutoff = n_clip_cutoff/2 - 1` is derived from it inside
    `TransductionParams`, so one swept value moves both of gate [A]'s thresholds -- which is
    itself a finding about the parameterisation and is reported as such.

 2. **Whether `low_polya_clip` is REACHABLE at the shipping cutoff.** The polyA record only comes
    into existence when a polyA/polyT/ambiguous clip is seen, and the counter is incremented in
    the same block (`x_transduction.py:929-975`), so `sum(rcd[:6]) >= 1` for every entry that
    exists. At `n_clip_cutoff = 3` the derived `n_polya_cutoff` is 0.5, so `< 0.5` needs a sum of
    zero. The prediction is that the *threshold* decides nothing and gate [A]'s polyA decision is
    a **membership** test; the probe counts zero-sum entries and reports the first swept cutoff at
    which the reason appears at all.

 3. **The `n_clip` distribution of the `no_polya_clip` population**, which the ledger cannot
    supply: `_drop(R_NO_POLYA_CLIP)` records no detail, so a locus dropped for polyA-map absence
    is indistinguishable in the ledger from any other. That population is 6,056 (L1) and 11,304
    (SVA), and how much clip support it carries is the question of whether the membership test is
    discarding well-supported loci or noise.

 4. **How much of gate [A]'s evidence is hits to the DECOY CONSENSUS.** The clip path's decoy skip
    is commented out (`x_transduction.py:872-875`) -- only the *disc* path tallies cross-family
    decoy hits -- so a clipped read that realigns to the `>LINE1`/`>SVA`/`>ALU` consensus records
    counts as clip support. Such a record's `s_transduct` is written by the parse itself as
    `":-1--1~0~<pos>"` (`x_transduction.py:994` with `len(flank_id_fields) <= 1`), so this is read
    off the shipping writer's own output rather than re-classified here.

AND ONE SCOPE CORRECTION IT CAN MAKE. `NEARBY_CLIP` was measured as moving "exactly one locus in
the whole genome" (`probe_nearby_clip.py`) -- but that was measured **inside gate [B2]**, and the
constant ALSO filters the clip SAM parse at `x_transduction.py:865`, before the unmapped and
clip-qualified tests, which is gate [A]'s domain. So `--nearby 50,1000000000` re-parses at
unbounded distance and reports what gate [A]'s domain and survivor count do. That is a different
population from the one the earlier sweep bounded.

CHECK 0 IS TWO CHECKS AND BOTH ARE STRONGER THAN A COUNT. At the shipping `(NEARBY_CLIP,
n_clip_cutoff)` pair the probe must reproduce (a) the run's own ledger `A_clip` reason split and
(b) the run's own `.unique_trsdct_half_clip` file as a SET of loci, not merely its line count.
(b) is what makes the configuration -- including `i_flank_length`, which gate [A] cannot read
because it only sums list lengths -- verified rather than assumed. A mismatch means this probe is
misconfigured, not that a finding was made.

WHAT IT CANNOT SAY: whether a recovered locus becomes a CALL. Gate [A]'s survivors still face
[B], [B2], [C], [D] and [E], and `ledger_attribution.py` measured 306 of BASELINE-0's 322 false
negatives as never reaching this decision at all, with the 16 that do enriched by only ~8 calls
over chance. Expect a precision-shaped effect, and a truth-set-scored two-arm run is what would
settle any move.
"""
import argparse
import os
import re
import sys

REPO_ROOT = os.path.dirname(os.path.dirname(os.path.dirname(
    os.path.dirname(os.path.abspath(__file__)))))
sys.path.insert(0, os.path.join(REPO_ROOT, "xtea"))
import global_values          # noqa: E402
import transduction_evidence  # noqa: E402
import x_transduction         # noqa: E402

CLIP_SAM = os.path.join("tmp", "transduction", "temp_transduction_clip.sam")
STEM = "candidate_disc_filtered_cns2.txt.all_non_sibling_td.txt"
SURVIVORS = STEM + ".unique_trsdct_half_clip"
LEDGER = STEM + ".decision_ledger.tsv"
A_PHASE = "A_clip"
FAMILIES = ("L1", "SVA")
DEFAULT_SWEEP = (1, 2, 3, 4, 5, 6, 8, 10, 15, 20)
# The parse writes this for a hit to a record with no `~`-separated structure -- i.e. one of the
# decoy consensus records. x_transduction.py:994, with source_start/source_end left at -1.
NON_FLANK_PREFIX = ":-1--1~"
# `_drop(R_LOW_CLIP)`'s detail, which is where the run's own coverage-derived cutoff is written.
RE_LOW_CLIP_DETAIL = re.compile(r"^n_clip=\d+ < (\d+(?:\.\d+)?)$")


def quantiles(l_sorted):
    """min p25 median p75 p90 max, for a pre-sorted list. Empty list -> None."""
    if not l_sorted:
        return None

    def _q(f):
        return l_sorted[min(len(l_sorted) - 1, int(f * len(l_sorted)))]
    return (l_sorted[0], _q(0.25), _q(0.5), _q(0.75), _q(0.9), l_sorted[-1])


def ledger_a_split(sf_ledger):
    """{reason: count} over the ledger's A_clip rows, plus the cutoff its detail records."""
    m = {}
    s_cutoff = None
    with open(sf_ledger) as fin:
        for line in fin:
            if line.startswith("#"):
                continue
            l_f = line.rstrip("\n").split("\t")
            if len(l_f) < 4 or l_f[2] != A_PHASE:
                continue
            m[l_f[3]] = m.get(l_f[3], 0) + 1
            if s_cutoff is None and len(l_f) >= 5:
                mo = RE_LOW_CLIP_DETAIL.match(l_f[4])
                if mo is not None:
                    s_cutoff = mo.group(1)
    return m, s_cutoff


def load_survivor_set(sf_survivors):
    """{(chrm, pos)} from the run's own gate [A] output file."""
    s = set()
    with open(sf_survivors) as fin:
        for line in fin:
            l_f = line.split()
            if len(l_f) >= 2:
                s.add((l_f[0], int(l_f[1])))
    return s


def clip_counts(m_lclip, m_rclip):
    """{(chrm,pos): n_clip}, summed over both sides exactly as `select_clip_sites` sums it."""
    m = {}
    for m_side in (m_lclip, m_rclip):
        for chrm in m_side:
            for pos in m_side[chrm]:
                n = 0
                for ref_mpos in m_side[chrm][pos]:
                    n += len(m_side[chrm][pos][ref_mpos])
                m[(chrm, pos)] = m.get((chrm, pos), 0) + n
    return m


def non_flank_support(m_lclip, m_rclip):
    """(n_records, n_non_flank, {(chrm,pos): (n_records, n_non_flank)}).

    A non-flank record is one the parse wrote with no source coordinates -- see NON_FLANK_PREFIX.
    """
    n_all = n_nf = 0
    m_per_locus = {}
    for m_side in (m_lclip, m_rclip):
        for chrm in m_side:
            for pos in m_side[chrm]:
                for ref_mpos in m_side[chrm][pos]:
                    for s_td in m_side[chrm][pos][ref_mpos]:
                        b_nf = s_td.startswith(NON_FLANK_PREFIX)
                        n_all += 1
                        n_nf += 1 if b_nf else 0
                        a, b = m_per_locus.get((chrm, pos), (0, 0))
                        m_per_locus[(chrm, pos)] = (a + 1, b + (1 if b_nf else 0))
    return n_all, n_nf, m_per_locus


def sweep_family(xt, sf_sam, i_flank_length, l_sweep, n_ship_cutoff):
    """Parse once, then run the shipping gate at every cutoff. Returns a result dict."""
    m_lclip, m_rclip, m_polya = xt.parse_transduction_clip_algnmt(sf_sam, i_flank_length)

    m_counts = clip_counts(m_lclip, m_rclip)
    n_records, n_non_flank, m_nf_per_locus = non_flank_support(m_lclip, m_rclip)

    # Measurement 2: is `low_polya_clip` reachable? Every entry the parse creates is incremented
    # in the same block, so a zero sum should be impossible.
    n_polya_loci = 0
    n_polya_zero = 0
    for chrm in m_polya:
        for pos in m_polya[chrm]:
            n_polya_loci += 1
            if sum(m_polya[chrm][pos][:6]) == 0:
                n_polya_zero += 1

    l_rows = []
    m_ship = None
    for n_cut in l_sweep:
        # n_clip_cutoff ALONE, which is what gate [A] reads. Every disc-side threshold stays None
        # on purpose -- see TransductionParams' "a threshold a caller does not have is None, not
        # zero": if this sweep ever reached a gate it has no business in, the comparison raises
        # instead of silently admitting everything.
        params = transduction_evidence.TransductionParams(n_clip_cutoff=n_cut)
        ledger = transduction_evidence.DecisionLedger()
        m_slct = transduction_evidence.select_clip_sites(m_lclip, m_rclip, m_polya, params,
                                                        ledger=ledger)
        s_survivors = set()
        for chrm in m_slct:
            for pos in m_slct[chrm]:
                s_survivors.add((chrm, pos))
        m_reasons = {}
        for row in ledger.rows:
            m_reasons[row[3]] = m_reasons.get(row[3], 0) + 1
        row = {"cutoff": n_cut, "polya_cutoff": params.n_polya_cutoff,
               "n_survivors": len(s_survivors), "reasons": m_reasons,
               "survivors": s_survivors if n_cut == n_ship_cutoff else None}
        l_rows.append(row)
        if n_cut == n_ship_cutoff:
            m_ship = row

    # Measurement 3: the n_clip distribution of the no-polyA population, at the shipping cutoff.
    l_nopolya_nclip = []
    l_kept_nclip = []
    for (chrm, pos), n_clip in m_counts.items():
        if n_clip < n_ship_cutoff:
            continue
        if m_polya.get(chrm, {}).get(pos) is None:
            l_nopolya_nclip.append(n_clip)
        else:
            l_kept_nclip.append(n_clip)
    l_nopolya_nclip.sort()
    l_kept_nclip.sort()

    # Measurement 4, at the locus level: survivors whose clip support is ENTIRELY non-flank.
    n_all_nf_survivors = 0
    n_any_nf_survivors = 0
    if m_ship is not None:
        for key in m_ship["survivors"]:
            n_loc, n_loc_nf = m_nf_per_locus.get(key, (0, 0))
            if n_loc_nf > 0:
                n_any_nf_survivors += 1
                if n_loc_nf == n_loc:
                    n_all_nf_survivors += 1

    return {"n_loci": len(m_counts), "n_records": n_records, "n_non_flank": n_non_flank,
            "n_polya_loci": n_polya_loci, "n_polya_zero": n_polya_zero,
            "rows": l_rows, "ship": m_ship,
            "nopolya_nclip": l_nopolya_nclip, "kept_nclip": l_kept_nclip,
            "n_any_nf_survivors": n_any_nf_survivors,
            "n_all_nf_survivors": n_all_nf_survivors}


def main():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--run-dir", required=True,
                   help="<out>/<sample> -- the directory holding L1/ and SVA/")
    p.add_argument("--families", default=",".join(FAMILIES))
    p.add_argument("--sweep", default=",".join(str(v) for v in DEFAULT_SWEEP),
                   help="n_clip_cutoff values")
    p.add_argument("--nearby", default=str(global_values.NEARBY_CLIP),
                   help="NEARBY_CLIP values to re-parse at. The shipping value must be one of "
                        "them or check 0 cannot run. Each extra value costs a full SAM parse")
    p.add_argument("--n-clip-cutoff", type=int, default=None,
                   help="the run's own coverage-derived cutoff. Default: read from the ledger's "
                        "low_clip_support detail, which is where the run recorded it")
    p.add_argument("--flank-length", type=int, default=3000,
                   help="i_flank_length. Gate [A] cannot read it -- it sums list lengths -- and "
                        "check 0b is what proves that rather than assuming it")
    p.add_argument("--ref", default=None, help="reference FASTA; a SAM with a header needs none")
    a = p.parse_args()

    l_sweep = sorted(set(int(s) for s in a.sweep.split(",") if s))
    l_nearby = sorted(set(int(s) for s in a.nearby.split(",") if s))
    n_ship_nearby = global_values.NEARBY_CLIP

    print("=" * 78)
    print("GATE [A] -- swept through the SHIPPING gate on a finished run's own clip SAM")
    print("  run dir            : {0}".format(a.run_dir))
    print("  NEARBY_CLIP arms   : {0}   (shipping {1})".format(
        ", ".join("INF" if v >= 10 ** 9 else str(v) for v in l_nearby), n_ship_nearby))
    print("  n_clip_cutoff sweep: {0}".format(", ".join(str(v) for v in l_sweep)))
    print("=" * 78)

    l_problems = []
    for fam in [s for s in a.families.split(",") if s]:
        sf_sam = os.path.join(a.run_dir, fam, CLIP_SAM)
        sf_led = os.path.join(a.run_dir, fam, LEDGER)
        sf_srv = os.path.join(a.run_dir, fam, SURVIVORS)
        print("\n" + "=" * 78)
        print("--- {0} ---".format(fam))
        if not os.path.isfile(sf_sam):
            print("  no clip SAM at {0} (needs --keep-intermediates)".format(sf_sam))
            continue

        m_run_reasons, s_run_cutoff = (None, None)
        if os.path.isfile(sf_led):
            m_run_reasons, s_run_cutoff = ledger_a_split(sf_led)
        n_cut = a.n_clip_cutoff
        if n_cut is None:
            if s_run_cutoff is None:
                print("  no ledger low_clip_support detail to read the cutoff from; pass "
                      "--n-clip-cutoff")
                l_problems.append("{0}: cutoff unknown".format(fam))
                continue
            n_cut = int(float(s_run_cutoff))
            print("  n_clip_cutoff {0}, read from the run's own ledger detail".format(n_cut))
        else:
            print("  n_clip_cutoff {0}, given on the command line".format(n_cut))
        if n_cut not in l_sweep:
            l_sweep = sorted(set(l_sweep) | {n_cut})

        xt = x_transduction.XTransduction(os.path.join(a.run_dir, fam) + "/", 1, a.ref)

        for n_nearby in l_nearby:
            global_values.NEARBY_CLIP = n_nearby
            s_arm = "INF" if n_nearby >= 10 ** 9 else str(n_nearby)
            b_ship_arm = (n_nearby == n_ship_nearby)
            print("\n  NEARBY_CLIP = {0}{1}".format(
                s_arm, "   <- SHIPPING, check 0 runs on this arm" if b_ship_arm else ""))
            r = sweep_family(xt, sf_sam, a.flank_length, l_sweep, n_cut)

            print("    gate [A] domain: {0} loci from {1} qualified clip records".format(
                r["n_loci"], r["n_records"]))
            print("    of those records, hits to a DECOY CONSENSUS record (no source coords): "
                  "{0} ({1:.2f}%)".format(
                      r["n_non_flank"],
                      100.0 * r["n_non_flank"] / r["n_records"] if r["n_records"] else 0.0))
            print("    polyA-map loci: {0}, of which sum(rcd[:6]) == 0: {1}   "
                  "[0 means `low_polya_clip` needs n_polya_cutoff > 1, i.e. cutoff >= 5]".format(
                      r["n_polya_loci"], r["n_polya_zero"]))

            print("    {0:>6} {1:>7} {2:>10} {3:>9}   {4}".format(
                "cutoff", "polyA", "survivors", "vs ship", "reasons"))
            n_ship_srv = r["ship"]["n_survivors"] if r["ship"] else None
            for row in r["rows"]:
                print("    {0:>6} {1:>7} {2:>10} {3:>+9}   {4}{5}".format(
                    row["cutoff"], row["polya_cutoff"], row["n_survivors"],
                    row["n_survivors"] - n_ship_srv if n_ship_srv is not None else 0,
                    " ".join("{0}={1}".format(k, v) for k, v in sorted(row["reasons"].items())),
                    "   <- SHIPPING" if row["cutoff"] == n_cut else ""))

            l_first_low = [row["cutoff"] for row in r["rows"]
                           if row["reasons"].get("low_polya_clip", 0) > 0]
            print("    `low_polya_clip` first fires at cutoff: {0}".format(
                l_first_low[0] if l_first_low else "NEVER over this sweep"))

            q_no = quantiles(r["nopolya_nclip"])
            q_kept = quantiles(r["kept_nclip"])
            print("    at the shipping cutoff, loci with n_clip >= {0}:".format(n_cut))
            print("      dropped for polyA-map ABSENCE: {0}   n_clip min/p25/med/p75/p90/max "
                  "{1}".format(len(r["nopolya_nclip"]), q_no))
            print("      kept (in the polyA map)      : {0}   n_clip min/p25/med/p75/p90/max "
                  "{1}".format(len(r["kept_nclip"]), q_kept))
            n_strong = sum(1 for n in r["nopolya_nclip"] if n >= 10)
            print("      of the dropped, n_clip >= 10 : {0} ({1:.1f}%)".format(
                n_strong,
                100.0 * n_strong / len(r["nopolya_nclip"]) if r["nopolya_nclip"] else 0.0))
            print("    survivors with ANY decoy-consensus clip support: {0}; with ONLY that: "
                  "{1}".format(r["n_any_nf_survivors"], r["n_all_nf_survivors"]))

            if not b_ship_arm:
                continue

            # CHECK 0a: the run's own ledger A_clip rows.
            if m_run_reasons is None:
                print("    check 0a: no ledger -- SKIPPED, not passed")
                l_problems.append("{0}: no ledger".format(fam))
            else:
                b_ok = m_run_reasons == r["ship"]["reasons"]
                print("    check 0a: ledger A_clip {0}".format(
                    " ".join("{0}={1}".format(k, v) for k, v in sorted(m_run_reasons.items()))))
                print("              this probe  {0}   -> {1}".format(
                    " ".join("{0}={1}".format(k, v)
                             for k, v in sorted(r["ship"]["reasons"].items())),
                    "MATCH" if b_ok else "MISMATCH -- this probe is misconfigured"))
                if not b_ok:
                    l_problems.append("{0}: reasons {1} vs ledger {2}".format(
                        fam, r["ship"]["reasons"], m_run_reasons))

            # CHECK 0b: the run's own survivor file, as a SET.
            if not os.path.isfile(sf_srv):
                print("    check 0b: no survivor file -- SKIPPED, not passed")
                l_problems.append("{0}: no survivor file".format(fam))
            else:
                s_run = load_survivor_set(sf_srv)
                s_probe = r["ship"]["survivors"]
                b_ok = s_run == s_probe
                print("    check 0b: the run's .unique_trsdct_half_clip {0} loci; this probe "
                      "{1}; identical as SETS: {2}".format(
                          len(s_run), len(s_probe), "YES" if b_ok else "NO"))
                if not b_ok:
                    print("              probe-only {0}, run-only {1}".format(
                        len(s_probe - s_run), len(s_run - s_probe)))
                    l_problems.append("{0}: survivor sets differ".format(fam))

        global_values.NEARBY_CLIP = n_ship_nearby

    print("\n" + "=" * 78)
    print("check 0 overall: {0}".format(
        "OK" if not l_problems else "FAILED -- " + "; ".join(l_problems)))
    return 1 if l_problems else 0


if __name__ == "__main__":
    sys.exit(main())
