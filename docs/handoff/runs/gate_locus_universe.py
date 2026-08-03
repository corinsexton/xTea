#!/usr/bin/env python3
"""WGS-scale gate for the DRIVER-SIDE locus universe. Offline: no BAM, no bwa, no sbatch.

WHAT THIS GATES THAT `gate_site_universe.py` DOES NOT. That gate proves the lifted rule
reproduces a run's `.slct` when handed **the arguments that run's own step log records**. It
deliberately parses the log so it cannot drift into testing the driver instead of the run.

`xtea/downstream_unified.py` needs the opposite check. Its job is to build those arguments
**from the driver's own state** -- REP_TYPES, the rep-lib layout, the coverage-derived cutoff,
the step's boundary extend -- because the per-locus loop has no step log to read. Every one of
those is a place to be quietly wrong:

  - `annotation1` (`LINE/{g}/{g}_L1.fa.out`) is what build_type_steps hands transduction.
    `annotation` (`{g}_L1_larger_500_with_all_L1HS.out`) is cns_filter's, exists, and is a
    plausible-looking L1 annotation. Picking it changes the rule's third rejection.
  - the cutoff is NOT `--nd`: `b_automatic` is True unless `--user_specific`, so the step uses
    `Parameters.get_par_by_cov(purity-scaled coverage)[1]` -- 4 at HG002's 31.816x, not 5.
  - the rmsk boundary extend is XTransduction's own 450 (the transduction branch never calls
    `set_boundary_extend`; the sibling branch does, which is one reason `.slct2` differs).
  - `I_MIN_COPY_LEN` is 225.

So this rebuilds each argument the driver's way, **prints it next to what the run's log says
the step was actually given, and fails on any mismatch**, then runs the enumeration and byte-
compares. Both halves matter: agreeing bytes with a disagreeing annotation path would mean the
rule is insensitive to it on this data, not that the driver is right.

It also exercises `downstream_unified`'s union code at real scale rather than the ad-hoc union
`gate_site_universe.py` computes inline -- 41,604 loci over ~74K per-family claims, where the
unit tests see a handful.

USAGE
    python docs/handoff/runs/gate_locus_universe.py <run_root>/<sample> [rep_lib_dir] [genome]
e.g.
    python docs/handoff/runs/gate_locus_universe.py \\
        /n/scratch/users/c/cos689/xtea2_testdata/stage2c_family_assign_49502115/OFF/HG002.wgs

Exit 0 only if every flanked family reproduced byte-for-byte, every derived argument matched
the run's, AND the number of comparisons equals the number of flanked families present. That
last counter is job 49516776's lesson: a gate whose only check runs at comparison time cannot
tell "the code is wrong" from "the harness compared nothing".
"""
import os
import sys
import tempfile

REPO_ROOT = os.path.dirname(os.path.dirname(os.path.dirname(
    os.path.dirname(os.path.abspath(__file__)))))   # runs/ -> handoff/ -> docs/ -> repo root
sys.path.insert(0, os.path.join(REPO_ROOT, "xtea"))
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import downstream_unified                  # noqa: E402
from run_xtea import REP_TYPES             # noqa: E402
from x_parameter import Parameters         # noqa: E402
from gate_site_universe import (           # noqa: E402  (one copy of the log parsing)
    find_steps_log, parse_step_argv, automatic_disc_cutoff)

DEFAULT_REP_LIB = "/home/cos689/xTea/rep_lib_annotation"
DEFAULT_GENOME = "hg38"


def driver_annotation(rep_lib_dir, genome, family):
    """`run_xtea.resolve_annot`'s result for the annotation transduction is given.

    Kept as the same two-step (template from REP_TYPES, `{g}` -> genome) rather than a
    hard-coded path, so a change to the layout table is caught here instead of read past.
    """
    return os.path.join(rep_lib_dir, REP_TYPES[family]["annotation1"].format(g=genome))


def main():
    if len(sys.argv) < 2:
        print(__doc__)
        return 2
    sample_dir = sys.argv[1].rstrip("/")
    rep_lib_dir = sys.argv[2] if len(sys.argv) > 2 else DEFAULT_REP_LIB
    genome = sys.argv[3] if len(sys.argv) > 3 else DEFAULT_GENOME

    n_disc_cutoff, f_cov = automatic_disc_cutoff(sample_dir)
    print("driver-derived disc cutoff: Parameters.get_par_by_cov({0:.3f})[1] = {1}".format(
        f_cov, n_disc_cutoff))
    assert n_disc_cutoff == Parameters().get_par_by_cov(f_cov)[1]

    tmpdir = tempfile.mkdtemp(prefix="locus_universe_gate_")
    expected = 0          # flanked families present in this run
    compared = 0
    identical = 0
    mismatched_args = 0
    slct_by_family = {}

    for family in REP_TYPES:
        sf_log = find_steps_log(sample_dir, family)
        if sf_log is None:
            print("{0}: no steps log -- not in this run, skipping".format(family))
            continue
        argv = parse_step_argv(sf_log, "--transduction")
        if argv is None:
            print("{0}: steps log has no --transduction command; skipping".format(family))
            continue

        b_flank = bool(REP_TYPES[family]["flank"])
        # --- the driver's arguments, built the driver's way -------------------------------
        type_dir = os.path.join(sample_dir, family)
        d_raw_disc = os.path.join(type_dir,
                                  "candidate_list_from_disc.txt.clip_sites_raw_disc.txt")
        d_cns = os.path.join(type_dir, "candidate_disc_filtered_cns.txt")
        d_annot = driver_annotation(rep_lib_dir, genome, family)
        d_rtype = REP_TYPES[family]["rtype"]
        # --- and what the run's step was actually given ----------------------------------
        for (label, got, want) in (("raw_disc", d_raw_disc, argv["--input2"]),
                                   ("cns list", d_cns, argv["-i"]),
                                   ("annotation", d_annot, argv["-a"]),
                                   ("rtype", str(d_rtype), argv["--rtype"])):
            if os.path.realpath(str(got)) != os.path.realpath(str(want)) and got != want:
                print("  {0}: MISMATCH {1}\n      driver: {2}\n      run   : {3}".format(
                    family, label, got, want))
                mismatched_args += 1
            else:
                print("  {0}: {1:10s} matches the run ({2})".format(family, label, got))

        sf_out = downstream_unified.family_slct_path(tmpdir, family)
        n_cand, n_new = downstream_unified.enumerate_family_sites(
            d_raw_disc, d_cns, d_annot, d_rtype, n_disc_cutoff,
            bool(argv.get("--tumor", False)), sf_out, b_flank=b_flank,
            working_folder=tmpdir, log=lambda *a, **k: None)
        slct_by_family[family] = sf_out
        print("  {0}: {1} cns candidates + {2} created = {3} sites".format(
            family, n_cand, n_new, n_cand + n_new))

        sf_want = d_raw_disc + ".slct"
        if not b_flank:
            # No transduction flank -> the step copies its input through and writes no
            # `.slct` at all, so there is nothing to compare. A pass, not a skip: its
            # universe contribution is the cns list, which probe_site_creation.py confirms
            # from the other side (all 1,230 Alu records first exist at cns).
            print("  {0}: flankless -- step writes no .slct; nothing to compare".format(family))
            if os.path.isfile(sf_want):
                print("  {0}: FATAL: a .slct exists for a family that should have none".format(
                    family))
                mismatched_args += 1
            continue
        expected += 1
        if not os.path.isfile(sf_want):
            print("  {0}: FATAL: flanked family has no .slct at {1} -- was the run made "
                  "with --keep-intermediates?".format(family, sf_want))
            continue
        compared += 1
        with open(sf_want, "rb") as fa, open(sf_out, "rb") as fb:
            a, b = fa.read(), fb.read()
        same = (a == b)
        identical += 1 if same else 0
        print("  {0}: vs the run's own .slct: {1} ({2} vs {3} bytes)".format(
            family, "IDENTICAL" if same else "DIFFERS", len(a), len(b)))
        if not same:
            la, lb = a.decode().splitlines(), b.decode().splitlines()
            print("      rows {0} vs {1}".format(len(la), len(lb)))
            for i in range(min(len(la), len(lb))):
                if la[i] != lb[i]:
                    print("      first differing row {0}:\n        run: {1}\n        drv: {2}"
                          .format(i, la[i], lb[i]))
                    break

    # --- the union, through the module the loop will use ----------------------------------
    universe = downstream_unified.build_locus_universe(slct_by_family)
    per_fam = universe.per_family()
    n_claims = sum(per_fam.values())
    print("\nlocus universe: {0} loci from {1} per-family claims ({2:.2f}x)".format(
        universe.n_loci, n_claims,
        (n_claims / universe.n_loci) if universe.n_loci else 0.0))
    for tag in sorted(per_fam):
        print("    {0:4s} {1:8d}".format(tag, per_fam[tag]))
    by_n = universe.by_n_families()
    print("    claimed by " + ", ".join(
        "{0} famil{1}: {2}".format(n, "y" if n == 1 else "ies", by_n[n])
        for n in sorted(by_n)))
    print("    reachable only by creation (no family has a cns row there): {0}".format(
        universe.n_created_only()))
    sf_table = os.path.join(tmpdir, "locus_universe.tsv")
    universe.write(sf_table)
    n_rows = sum(1 for line in open(sf_table) if not line.startswith("#"))
    print("    table {0}: {1} rows".format(sf_table, n_rows))
    if n_rows != universe.n_loci:
        print("    FATAL: the written table has {0} rows for {1} loci".format(
            n_rows, universe.n_loci))
        mismatched_args += 1

    print("\nproducts compared: {0} (expected {1})   identical: {2}   "
          "argument mismatches: {3}".format(compared, expected, identical, mismatched_args))
    ok = (compared == expected and identical == expected and expected > 0
          and mismatched_args == 0)
    print("GATE: {0}".format("PASS" if ok else "FAIL"))
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
