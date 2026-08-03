#!/usr/bin/env python3
"""WGS-scale gate for the lifted site-creation rule -- and the union locus universe it makes
measurable. Offline: reads a finished run's kept intermediates. No BAM, no bwa, no sbatch.

WHAT IT GATES. `xtea/site_universe.write_slct_file` is the rule that used to live inside
`XTransduction.re_slct_with_clip_raw_disc_sites`, and it is what CREATES the 52 of
BASELINE-0's 1,515 final per-family records that have no `candidate_disc_filtered_cns.txt`
row (50 L1, 2 SVA -- `probe_site_creation.py`). `tests/test_site_universe.py` pins it against
a verbatim copy of the replaced loop, but on synthetic rows: a handful of chromosomes, a few
dozen sites, a fake annotation. Everything that could still differ is scale- or data-shaped:

  - the real rmsk annotation is 978,710 L1 intervals loaded at extnd 450 / min copy 225,
    and the rule's third rejection reads `get_div_subfamily_with_min_div` off it;
  - the real raw_disc is ~39,600 rows against 300 cns candidates spread over every contig,
    including decoys and the proximity window's edges as the data happens to place them;
  - the real cutoff is not `--nd`: `b_automatic` is True by default, so it comes from the
    run's own coverage through `Parameters.get_par_by_cov` (4 at HG002's 31.8x, not 5).

So this reruns the rule with the arguments the run itself used -- parsed out of the run's
`*.steps.log`, not retyped -- and byte-compares against the `.slct` file that run produced.

WHAT IT ALSO REPORTS, because it is free once the selections are in memory: the size of the
**union locus universe** across families, which is the per-locus loop's iteration domain
(plan sec 3, "one pass over loci"). Nobody had that number; the 53,945-of-54K figure in plan
sec 4 2b is about the td gather's fetch WINDOWS, which is a different thing.

Alu is expected to have no `.slct` at all -- no transduction flank, so x_TEA_main copies its
input through and the rule never runs. A missing Alu file is a pass, not a skip.

USAGE
    python docs/handoff/runs/gate_site_universe.py <run_root>/<sample>
e.g.
    python docs/handoff/runs/gate_site_universe.py \
        /n/scratch/users/c/cos689/xtea2_testdata/stage2c_family_assign_49502115/OFF/HG002.wgs

Exit status is 0 only if every family that has a `.slct` reproduced byte-for-byte AND the
number of comparisons matches the number of families that should have one. That second
counter is the lesson job 49516776 cost: a gate whose only completeness check runs at
comparison time cannot tell "the code is wrong" from "the harness compared nothing".
"""
import os
import shlex
import sys
import tempfile

REPO_ROOT = os.path.dirname(os.path.dirname(os.path.dirname(
    os.path.dirname(os.path.abspath(__file__)))))   # runs/ -> handoff/ -> docs/ -> repo root
sys.path.insert(0, os.path.join(REPO_ROOT, "xtea"))
import site_universe                       # noqa: E402
from x_parameter import Parameters         # noqa: E402
from x_transduction import XTransduction   # noqa: E402

FAMILIES = ("L1", "Alu", "SVA")
I_MIN_COPY_LEN = 225          # x_TEA_main.py:519, the transduction branch's own value
TD_NON_SIBLING_SUFFIX = ".all_non_sibling_td.txt"   # global_values.TD_NON_SIBLING_SUFFIX


def find_steps_log(sample_dir, family):
    """The run's per-family step log, whatever the sample is called."""
    for name in os.listdir(sample_dir):
        if name.endswith("_{0}.steps.log".format(family)):
            return os.path.join(sample_dir, name)
    return None


def parse_step_argv(sf_log, s_step_flag):
    """The command line this run actually executed for one step, as a dict of flags.

    Taken from the log rather than rebuilt from build_type_steps: the point of this gate is
    to reproduce THAT run, and a rebuilt argv would silently track any later change to the
    driver instead of failing.
    """
    with open(sf_log) as fin:
        for line in fin:
            if s_step_flag not in line or "x_TEA_main.py" not in line:
                continue
            toks = shlex.split(line.strip())
            argv = {}
            i = 0
            while i < len(toks):
                t = toks[i]
                if t.startswith("-") and i + 1 < len(toks) and not toks[i + 1].startswith("-"):
                    argv[t] = toks[i + 1]
                    i += 2
                else:
                    argv[t] = True
                    i += 1
            return argv
    return None


def read_mean_insert_size(sample_dir):
    """The run's mean insert size, which the SIBLING caller passes as the rmsk boundary
    extend (`xorphan.set_boundary_extend(mean_is)`, x_TEA_main.py:615) instead of
    XTransduction's default 450. Same cached basic-info file as the cutoff."""
    sf_basic = os.path.join(sample_dir, "tmp_basic", "basic_cov_is_rlth_info.txt")
    with open(sf_basic) as fin:
        fields = fin.readline().split()
    return float(fields[3])


def automatic_disc_cutoff(sample_dir):
    """The disc cutoff the run used. `b_automatic` is True unless --user_specific, so --nd is
    NOT it: x_TEA_main takes rcd[1] from Parameters.get_par_by_cov(coverage), and the
    coverage is in the basic-info file the run cached."""
    sf_basic = os.path.join(sample_dir, "tmp_basic", "basic_cov_is_rlth_info.txt")
    with open(sf_basic) as fin:
        fields = fin.readline().split()
    f_cov = float(fields[1])
    par = Parameters().get_par_by_cov(f_cov)
    return par[1], f_cov


def main():
    if len(sys.argv) != 2:
        print(__doc__)
        return 2
    sample_dir = sys.argv[1].rstrip("/")
    n_disc_cutoff, f_cov = automatic_disc_cutoff(sample_dir)
    print("run coverage {0:.3f} -> automatic disc cutoff {1} "
          "(the --nd on the command line is NOT what the step used)".format(f_cov, n_disc_cutoff))

    mean_is = read_mean_insert_size(sample_dir)
    print("run mean insert size {0} -- the sibling caller's rmsk boundary extend "
          "(the transduction caller's is XTransduction's default 450)".format(mean_is))

    expected = 0          # per-caller selections that exist and must be reproduced
    compared = 0
    identical = 0
    universes = {}
    tmpdir = tempfile.mkdtemp(prefix="site_universe_gate_")

    for family in FAMILIES:
        sf_log = find_steps_log(sample_dir, family)
        if sf_log is None:
            print("{0}: no steps log -- family not in this run, skipping".format(family))
            continue
        argv = parse_step_argv(sf_log, "--transduction")
        if argv is None:
            print("{0}: steps log has no --transduction command; skipping".format(family))
            continue
        sf_raw_disc = argv["--input2"]
        sf_cns = argv["-i"]
        sf_rmsk = argv["-a"]
        i_rep_type = int(argv["--rtype"])
        b_tumor = bool(argv.get("--tumor", False))

        if not os.path.isfile(sf_raw_disc + ".slct"):
            # Alu: no flank -> x_TEA_main copies its input through and never re-selects, so
            # its locus universe is exactly its cns candidate list. Reported, not gated.
            print("{0}: no .slct (rtype {1}) -- expected for a family with no transduction "
                  "flank; nothing to reproduce. Its universe is the cns list alone".format(
                      family, i_rep_type))
            if os.path.isfile(sf_cns):
                universes[family] = site_universe.read_positions(sf_cns)
            continue

        # BOTH production callers, because they differ in every argument that could hide a
        # bug: the sibling one passes b_with_candidate=False, a different candidate file
        # (the transduction step's non-sibling output) and a different rmsk boundary extend
        # (mean insert size, a float, instead of 450 -- x_TEA_main.py:615).
        #   (label, candidates, b_with_candidate, boundary_extnd, produced file)
        sf_sib_argv = parse_step_argv(sf_log, "--sibling")
        callers = [(".slct", sf_cns, True, None, sf_raw_disc + ".slct")]
        if sf_sib_argv is not None and os.path.isfile(sf_raw_disc + ".slct2"):
            callers.append((".slct2",
                            sf_sib_argv["-i"] + TD_NON_SIBLING_SUFFIX, False, mean_is,
                            sf_raw_disc + ".slct2"))

        for (label, sf_candidates, b_with_candidate, extnd, sf_want) in callers:
            expected += 1
            xtd = XTransduction(tmpdir, 1, None)
            if extnd is not None:
                xtd.set_boundary_extend(extnd)
            xannotation = xtd.prep_annotation_interval_tree(sf_rmsk, I_MIN_COPY_LEN)
            sf_out = os.path.join(tmpdir, "{0}{1}".format(family, label))
            n_cand, n_new = site_universe.write_slct_file(
                sf_raw_disc, sf_candidates, n_disc_cutoff, xannotation, i_rep_type, b_tumor,
                sf_out, b_with_candidate=b_with_candidate, log=lambda *a, **k: None)
            compared += 1

            with open(sf_want, "rb") as fa, open(sf_out, "rb") as fb:
                a, b = fa.read(), fb.read()
            same = (a == b)
            identical += 1 if same else 0
            print("{0}{1}: rtype {2}, tumor {3} -- {4} candidates + {5} new sites = {6} rows; "
                  "vs the run's own file: {7} ({8} vs {9} bytes)".format(
                      family, label, i_rep_type, b_tumor, n_cand, n_new, n_cand + n_new,
                      "IDENTICAL" if same else "DIFFERS", len(a), len(b)))
            if not same:
                la = a.decode().splitlines()
                lb = b.decode().splitlines()
                print("    rows {0} vs {1}".format(len(la), len(lb)))
                for i in range(min(len(la), len(lb))):
                    if la[i] != lb[i]:
                        print("    first differing row {0}:\n      run: {1}\n      lib: {2}".format(
                            i, la[i], lb[i]))
                        break
            if label == ".slct":
                universes[family] = site_universe.read_positions(sf_out)

    # --- the union locus universe: the per-locus loop's iteration domain -------------------
    if universes:
        print("\nunion locus universe across {0}".format("+".join(sorted(universes))))
        all_pos = {}
        for family, m in universes.items():
            n = sum(len(v) for v in m.values())
            print("    {0:4s} {1:8d} sites".format(family, n))
            for chrm in m:
                for p in m[chrm]:
                    all_pos.setdefault((chrm, p), set()).add(family)
        by_count = {}
        for fams in all_pos.values():
            by_count[len(fams)] = by_count.get(len(fams), 0) + 1
        n_per_family = sum(sum(len(v) for v in m.values()) for m in universes.values())
        print("    union (exact position) {0} sites: {1}".format(
            len(all_pos), ", ".join("{0} claimed by {1} family/ies".format(v, k)
                                    for k, v in sorted(by_count.items()))))
        print("    so a loop over the union visits {0} loci where the per-family chains visit "
              "{1} -- {2:.2f}x".format(len(all_pos), n_per_family,
                                       n_per_family / float(len(all_pos))))

    print("\nproducts compared: {0} (expected {1})".format(compared, expected))
    ok = (expected > 0) and (compared == expected) and (identical == expected)
    print("GATE: {0} -- {1}/{2} byte-identical".format("PASS" if ok else "FAIL",
                                                       identical, expected))
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
