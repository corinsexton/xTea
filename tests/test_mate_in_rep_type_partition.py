"""`--types` must partition cleanly: a family's call set may not depend on which OTHER
families were requested.

WHY THIS GATE EXISTS
--------------------
`mate_in_rep_clip_cnt` is the number of clipped reads at a site whose mate lands inside any
selected repeat family's annotation -- a repeat-anchored support signal, a union over the
families named by --mate-in-rep-families, emitted identically for every type at a site, and
consumed only by the SVA cns path (spec §4 DC1). It is xtea-scan's 6th realign-out column,
which is the only reason it used to be called `col5`: a column position, not a concept.

run_xtea.py used to build those annotation BEDs by iterating the *requested* types and
skipping any outside the mate-in-rep family set, i.e. the BED set was
`requested_types & mate_in_rep_families`:

    --types L1,Alu,SVA  (default families alu,sva)  ->  ALU + SVA  ->  mate-in-(Alu | SVA)
    --types SVA                                    ->  SVA only   ->  mate-in-SVA only
    --types L1                                     ->  no BEDs    ->  0 at every site

Because union >= subset always, and because this value gates the SVA cns consistency check,
`--types SVA` could emit FEWER SVA calls than the same sample run as `--types L1,Alu,SVA`.
That is a live correctness bug, not a redesign artifact: asking for less must not change what
you get for what you did ask for.

WHAT THIS PINS
--------------
The derivation, not a timing or a call count -- the family set is known from
--mate-in-rep-families before any read is touched, so a unit test is the right altitude and a
WGS run is not.

1. `mate_in_rep_bed_plan` does not accept the requested types AT ALL. This is the structural
   guarantee: a function that cannot see `types` cannot branch on it. If someone threads
   `types` back in, this test fails at the signature, before any behavioural assertion.
2. The plan depends only on --mate-in-rep-families, over every possible requested type set.
3. The default plan is unchanged (ALU + SVA, SVA extended by 200), so three-type output does
   not move -- the other half of the fix's signature.
4. `--col5-families` still works as a deprecated alias, so existing run scripts (several
   under docs/handoff/runs/ pass it) do not break.

DOWNSTREAM_OVERHAUL_PLAN.md §6.1/§6.2, Stage 0.
"""
import inspect
import itertools
import os
import subprocess
import sys

REPO_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(REPO_ROOT, "xtea"))
import run_xtea  # noqa: E402  (needs the path insert above)
import scan_adapter  # noqa: E402

ALL_TYPES = ["L1", "Alu", "SVA"]


def _all_requested_type_sets():
    """Every non-empty subset of the callable types, i.e. every legal --types value."""
    for n in range(1, len(ALL_TYPES) + 1):
        for combo in itertools.combinations(ALL_TYPES, n):
            yield list(combo)


def test_bed_plan_cannot_see_requested_types():
    """The structural guarantee: `types` is not a parameter, so it cannot be branched on.

    This is the assertion that actually prevents regression. The behavioural tests below
    would all still pass if someone added an unused `types` argument -- and then the next
    edit could start using it.
    """
    params = list(inspect.signature(run_xtea.mate_in_rep_bed_plan).parameters)
    assert params == ["mate_in_rep_fams"], (
        f"mate_in_rep_bed_plan must take ONLY the family set, got {params}. The count is "
        "defined by --mate-in-rep-families alone; giving this function access to the "
        "requested --types is how `--types SVA` came to emit fewer SVA calls than "
        "`--types L1,Alu,SVA`."
    )


def test_plan_is_identical_across_every_requested_type_set():
    """The behaviour: one --mate-in-rep-families value, one plan, whatever --types says.

    Iterating requested type sets here is deliberately vacuous -- mate_in_rep_bed_plan is
    not given them. That is the point: the test states the property the caller must
    preserve, so it keeps meaning if the call site ever regains a type-dependent path.
    """
    for fams in ({"alu", "sva"}, {"sva"}, {"l1"}, {"l1", "alu", "sva"}, {"l1", "sva"}):
        expected = run_xtea.mate_in_rep_bed_plan(fams)
        for types in _all_requested_type_sets():
            got = run_xtea.mate_in_rep_bed_plan(fams)
            assert got == expected, (
                f"BED plan for --mate-in-rep-families={sorted(fams)} changed when "
                f"--types={types}: {got} != {expected}"
            )
            # The contributing families must be exactly those requested as mate-in-rep
            # families -- never narrowed to the intersection with `types`.
            assert {t for t, _fam, _e in got} == {
                run_xtea.MATE_IN_REP_KEY_TO_TYPE[k] for k in fams}


def test_single_type_run_still_gets_the_full_union():
    """The specific regression: SVA-only must see the same Alu|SVA union as a three-type run.

    This is the case that was silently wrong. Under the old code the SVA-only plan lost the
    ALU BED, and since mate-in-(Alu|SVA) >= mate-in-SVA, the SVA cns consistency check saw a
    weaker signal and could drop calls a three-type run kept.
    """
    plan = run_xtea.mate_in_rep_bed_plan({"alu", "sva"})
    families = [fam for _t, fam, _e in plan]
    assert families == ["ALU", "SVA"], families
    assert "ALU" in families, (
        "an SVA-only run must still count mates landing in Alu: SVA is a composite element "
        "with an Alu-derived domain, so mate-in-Alu is genuine SVA evidence (spec §4 DC1)"
    )


def test_l1_only_run_is_not_silently_zeroed():
    """`--types L1` used to generate no BEDs at all, making the count 0 at every site."""
    plan = run_xtea.mate_in_rep_bed_plan({"alu", "sva"})
    assert plan, ("--mate-in-rep-families=alu,sva must yield BEDs regardless of --types "
                  "being L1 only")


def test_default_plan_is_unchanged_so_three_type_output_does_not_move():
    """The other half of the fix's signature: three-type calls must be byte-identical.

    Pins the exact tuples, including SVA's +200 annotation extend -- the one asymmetry that
    would change SVA's BED content if it were lost -- and the documented default family set.
    """
    assert run_xtea.MATE_IN_REP_FAMILIES_DEFAULT == "alu,sva"
    default_fams = set(run_xtea.MATE_IN_REP_FAMILIES_DEFAULT.split(","))
    assert run_xtea.mate_in_rep_bed_plan(default_fams) == [
        ("Alu", "ALU", 0),
        ("SVA", "SVA", scan_adapter.SVA_ANNOTATION_EXTND),
    ]
    assert scan_adapter.SVA_ANNOTATION_EXTND == 200


def test_every_family_has_a_resolvable_annotation():
    """The precondition that makes the fix possible at all.

    The old call site read the annotation out of `type_paths`, which is built only for
    REQUESTED types -- that dict is the mechanism by which the bug was expressible. The fix
    resolves from REP_TYPES instead, which requires that every family nameable in
    --mate-in-rep-families carries a real annotation entry (not the `None` -> "null"
    sentinel that Alu's transduction flank uses). If a future family were added to
    MATE_IN_REP_KEY_TO_TYPE without one, the non-requested case would resolve to "null" and
    silently contribute an empty BED -- the same class of silent zeroing this fix removed.
    """
    for key, t in run_xtea.MATE_IN_REP_KEY_TO_TYPE.items():
        rel = run_xtea.REP_TYPES[t]["annotation"]
        assert rel is not None, (
            f"--mate-in-rep-families accepts {key!r} but {t} has no annotation")
        path = run_xtea.resolve_annot("/rep_lib", "hg38", rel)
        assert path != "null"
        assert path.startswith("/rep_lib/") and "hg38" in path


def test_deprecated_col5_families_alias_still_exists():
    """`--col5-families` must keep working: run scripts under docs/handoff/runs/ pass it.

    Checked through `--help` rather than by importing the parser, because the parser is built
    inside main(); this asserts the user-visible contract, which is the thing that would
    break a run script. Renaming a CLI flag out from under committed scripts is a loud
    failure at best -- argparse rejects the old flag and the job dies at startup -- and this
    keeps that from happening at all.
    """
    out = subprocess.run(
        [sys.executable, os.path.join(REPO_ROOT, "xtea", "run_xtea.py"), "--help"],
        capture_output=True, text=True, timeout=180)
    assert out.returncode == 0, out.stderr
    text = out.stdout
    assert "--mate-in-rep-families" in text, "the new flag must be documented"
    assert "--col5-families" in text, (
        "the deprecated alias must stay: docs/handoff/runs/*.sbatch pass it")
    # The alias must be marked, so nobody reads it as the current name. rindex, not index:
    # argparse prints the flag twice, first in the usage synopsis (which carries no help
    # text) and again in the options section -- the synopsis hit would always fail.
    col5_help = text[text.rindex("--col5-families"):]
    assert "DEPRECATED" in col5_help[:400], (
        "--col5-families must be labelled DEPRECATED in its help text")
