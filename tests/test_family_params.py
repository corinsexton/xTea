"""The per-family parameter bundle that replaced `global_values.IS_CALL_SVA`.

WHY THIS GATE EXISTS
--------------------
`IS_CALL_SVA` was a process-wide mutable flag, set from `x_TEA_main --sva` and read in
three modules. DOWNSTREAM_OVERHAUL_PLAN.md §4 Stage 2 retires it: after Stage 2c one
process handles all three families, so "which family am I" cannot be ambient.

The replacement is agreement-gated -- Stage 2 permits calls to move ONLY as
family-assignment flips, and a parameter-plumbing change is not one. So the interesting
assertions here are not "SVA gets SVA's numbers"; they are the two ways this refactor
could silently change behaviour:

  1. **The dormancy.** `--sva` is passed to the cns_filter step ONLY (the golden
     `run_xTEA_pipeline.sh` puts it on -C/-D/-N, not on transduction/sibling/postF), so
     `IS_CALL_SVA` was False inside the transduction subprocess even for SVA, and
     `leftmost_cns` was 500 there. Resolving the bundle from `--rtype` instead -- the
     obvious-looking cleanup -- would wake SVA's 5000 and move calls in bulk: the longest
     record in consensus/SVA.fa is 3183 bp, so a 5000 threshold makes
     `_reach_5_prime_end` true for every candidate with a positive cluster start, and the
     polymorphic-flank library stops being a filtered subset. Test 3 pins the resolution
     rule so that change has to be deliberate.

  2. **The default.** Every object that used to read the global must now default to the
     GENERIC bundle, whose numbers are exactly what `IS_CALL_SVA == False` meant. A
     constructor that defaulted to the L1 *or* SVA bundle instead would be a silent
     behaviour change at the five call sites that pass nothing.
"""
import inspect
import os
import sys

REPO_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(REPO_ROOT, "xtea"))
import family_params  # noqa: E402  (needs the path insert above)
import global_values  # noqa: E402
import run_xtea  # noqa: E402


# --- 1. the values ---------------------------------------------------------------------

def test_family_values():
    """The three bundles carry the numbers the plan names, per family."""
    l1, alu, sva = (family_params.for_name(t) for t in ("L1", "Alu", "SVA"))

    # flank: None for Alu -- a non-autonomous SINE cannot 3'-read-through.
    assert (l1.has_flank, alu.has_flank, sva.has_flank) == (True, False, True)
    # annotation_extend: 0 / 0 / 200.
    assert (l1.annotation_extend, alu.annotation_extend, sva.annotation_extend) == (0, 0, 200)
    # leftmost_cns: 500 / 500 / 5000.
    assert (l1.leftmost_cns, alu.leftmost_cns, sva.leftmost_cns) == (500, 500, 5000)
    # The two cns-path relaxations are SVA's alone.
    assert [f.relaxed_cns_consistency for f in (l1, alu, sva)] == [False, False, True]
    assert [f.one_side_flank_high_confidence for f in (l1, alu, sva)] == [False, False, True]
    # rtype bitmask agrees with x_TEA_main's --rtype and run_xtea's REP_TYPES.
    assert (l1.rtype, alu.rtype, sva.rtype) == (1, 2, 4)
    for t in ("L1", "Alu", "SVA"):
        assert family_params.for_name(t).rtype == run_xtea.REP_TYPES[t]["rtype"]
        assert family_params.for_name(t).vcf_alt == run_xtea.REP_TYPES[t]["vcf_alt"]
        assert family_params.for_name(t).has_flank == bool(run_xtea.REP_TYPES[t]["flank"])


def test_generic_is_what_the_absent_flag_meant():
    """GENERIC must be numerically identical to `IS_CALL_SVA == False`.

    Five call sites construct these objects without naming a family and rely on this.
    """
    g = family_params.GENERIC
    assert g.name is None                      # distinguishable from L1, deliberately
    assert g.leftmost_cns == 500
    assert g.annotation_extend == 0
    assert g.relaxed_cns_consistency is False
    assert g.one_side_flank_high_confidence is False


def test_the_global_is_gone():
    """Not merely unused -- absent. A leftover attribute would let a stale
    `global_values.IS_CALL_SVA = True` assignment (docs/handoff/runs had one) keep
    parsing while silently doing nothing.
    """
    assert not hasattr(global_values, "IS_CALL_SVA")
    assert not hasattr(global_values, "turn_on_sva")


# --- 2. the defaults at every construction site ------------------------------------------

def test_constructors_default_to_generic():
    """XClipDiscFilter / XTransduction / XOrphanTransduction take `family` and default it
    to GENERIC, so an un-updated caller behaves exactly as it did under the global."""
    from x_clip_disc_filter import XClipDiscFilter
    from x_orphan_transduction import XOrphanTransduction
    from x_transduction import XTransduction

    for cls in (XClipDiscFilter, XTransduction, XOrphanTransduction):
        params = inspect.signature(cls.__init__).parameters
        assert "family" in params, "{0} lost its family parameter".format(cls.__name__)
        assert params["family"].default is None

    assert XTransduction("/tmp", 1, None).family is family_params.GENERIC
    assert XOrphanTransduction("/tmp", 1, None).family is family_params.GENERIC
    assert XClipDiscFilter(None, "/tmp", 1, None).family is family_params.GENERIC
    sva = family_params.for_name("SVA")
    assert XTransduction("/tmp", 1, None, family=sva).family is sva


# --- 3. the resolution rule, i.e. the dormancy ------------------------------------------

def test_resolution_from_the_sva_flag_only():
    """`from_sva_flag` is the faithful stand-in for the old `if options.sva:` -- and it is
    what x_TEA_main uses for EVERY step, including transduction, which never receives the
    flag. So transduction resolves GENERIC even when its --rtype says SVA."""
    assert family_params.from_sva_flag(True) is family_params.for_name("SVA")
    assert family_params.from_sva_flag(False) is family_params.GENERIC


def test_transduction_step_is_not_given_the_sva_flag():
    """The structural fact behind the dormancy, read off the argv the driver actually
    builds rather than asserted in prose. If someone adds --sva to the transduction or
    sibling step, SVA's leftmost_cns=5000 wakes up and calls move -- this test is where
    that decision has to be made consciously."""
    steps = _build_sva_steps()
    for name in ("transduction", "sibling", "postF_all", "postF_hc"):
        assert "--sva" not in steps[name], (
            "{0} now receives --sva. That wakes SVA's leftmost_cns=5000 in "
            "pick_full_length_insertion and is a Stage 3 change, not a plumbing one -- "
            "see xtea/family_params.py.".format(name))
    # ...and the one step that DOES get it, so this test cannot pass vacuously.
    assert "--sva" in steps["cns_filter"]


def test_sva_leftmost_cns_would_admit_every_candidate():
    """Why waking 5000 is not a small change: it is larger than any SVA consensus record,
    so `_reach_5_prime_end` becomes true for every positive cluster start.

    Pinned as an inequality against the parsed consensus when it is available, and as the
    bare numeric claim when it is not, so the test is meaningful on a bare checkout.
    """
    sva = family_params.for_name("SVA")
    assert sva.leftmost_cns == 5000
    sf = os.environ.get("XTEA2_REP_LIB_DIR", os.path.expanduser("~/xTea/rep_lib_annotation"))
    sf = os.path.join(sf, "consensus", "SVA.fa")
    if not os.path.isfile(sf):
        return  # bare checkout: the numeric assertion above is what remains
    longest = cur = 0
    with open(sf) as fin:
        for line in fin:
            if line.startswith(">"):
                longest = max(longest, cur)
                cur = 0
            else:
                cur += len(line.strip())
    longest = max(longest, cur)
    assert 0 < longest < sva.leftmost_cns, (
        "consensus/SVA.fa's longest record is {0} bp; leftmost_cns=5000 no longer "
        "trivially admits every candidate, so the reasoning in family_params' docstring "
        "needs rechecking".format(longest))


def _build_sva_steps():
    """The SVA per-type step chain's argv, keyed by step name."""
    opt = _Opt()
    cfg = {"reference": "/ref.fa", "blacklist": "/bl.bed", "gene_gff3": "/g.gff3",
           "model": "/m.onnx"}
    paths = {"rtype": 4, "flank": "/flank.fa", "cns": "/cns.fa",
             "annotation": "/a.out", "annotation1": "/a1.out"}
    dirs = {"type": "/out/SVA", "bam_list": "/bl.txt", "tmp_cns": "/tmp_cns",
            "tmp_td": "/tmp_td", "control_bam_list": "/cbl.txt"}
    return {name: argv
            for name, argv in run_xtea.build_type_steps("SVA", cfg, paths, dirs, opt)}


class _Opt:
    cores = 1
    nfclip = 3
    nfdisc = 5
    flklen = 3000
    tumor = False
    purity = 0.45
    case_control = False
    sample = "s"
