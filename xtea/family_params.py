"""One explicit per-family parameter bundle, replacing `global_values.IS_CALL_SVA`.

DOWNSTREAM_OVERHAUL_PLAN.md §4 Stage 2: "Retire `global_values.IS_CALL_SVA` global
state for an explicit per-family parameter bundle carrying `flank` (None for Alu),
`leftmost_cns` (500 / 5000 for SVA), `annotation_extend` (0 / 200 for SVA)."

Why a bundle rather than a boolean. `IS_CALL_SVA` was a process-wide mutable flag set
from a CLI flag (`x_TEA_main.py --sva`) and read four places, three modules apart. It
made the family an ambient property of the *process*, which is exactly the shape the
unified downstream has to get rid of: after Stage 2c one process handles all three
families, so "which family am I" cannot be a global. It also made the flag's reach
invisible -- which is how the dormancy documented below went unnoticed.

**THE ONE THING TO KNOW BEFORE CHANGING ANY RESOLUTION SITE.** `--sva` is passed to
the **cns_filter step only** (`run_xtea.py:build_type_steps`, matching the golden
`run_xTEA_pipeline.sh`, which puts it on -C/-D/-N and not on transduction/sibling/postF).
So in the transduction subprocess `IS_CALL_SVA` was **always False, including for SVA**,
and `leftmost_cns` was therefore **always 500** -- SVA's 5000 branch
(`x_transduction.py`) has never executed in this pipeline, in xTea2 or in master.

That dormancy is preserved here, deliberately and visibly: every call site resolves its
bundle from the same input the old global had, so this module is byte-identical in
behaviour. `SVA.leftmost_cns` carries the *intended* 5000 so the value lives with its
family rather than in a dead branch, and `GENERIC` is what a step that was never told
its family resolves to. Waking the SVA value is a one-word change (resolve from
`--rtype` instead of `--sva` in `x_TEA_main.py`'s transduction branch) and it is a
**Stage 3** change, not a Stage 2 one, because it is not a family-assignment flip:

    `pick_full_length_insertion` keeps a candidate whose cns cluster start is
    `0 < start < leftmost_cns`. The longest record in `consensus/SVA.fa` is 3183 bp,
    so a 5000 threshold admits EVERY candidate with a positive cluster start. The
    polymorphic-flank library would go from a filtered subset to essentially the whole
    candidate list, and every downstream bwa realignment with it. This moves calls in
    bulk; it is a threshold to calibrate against the truth set, not to flip blind.

`tests/test_family_params.py` pins both the values and the resolution rule.
"""


class FamilyParams:
    """The family-varying behaviour of one MEI family. Immutable by convention."""

    __slots__ = ("name", "rtype", "vcf_alt", "has_flank", "annotation_extend",
                 "leftmost_cns", "relaxed_cns_consistency",
                 "one_side_flank_high_confidence")

    def __init__(self, name, rtype, vcf_alt, has_flank, annotation_extend,
                 leftmost_cns, relaxed_cns_consistency,
                 one_side_flank_high_confidence):
        self.name = name
        self.rtype = rtype
        self.vcf_alt = vcf_alt
        # Alu is a non-autonomous SINE with no internal Pol II promoter, so it cannot
        # 3'-read-through: no transduction flank exists. Biologically correct, not a gap.
        self.has_flank = has_flank
        # Boundary widening applied to this family's repeat intervals for the
        # mate-in-repeat test (disc step and mate_in_rep_clip_cnt must agree -- an
        # inconsistent disc=200/cns=0 is what silently dropped SVA before).
        self.annotation_extend = annotation_extend
        # "at least hit the first N bp of the consensus" in pick_full_length_insertion.
        # See the module docstring: SVA's 5000 is DORMANT in today's pipeline.
        self.leftmost_cns = leftmost_cns
        # SVA's VNTR expansion makes copy length vary, so a site with neither side
        # clip/disc distance-consistent is NOT filtered out for SVA.
        self.relaxed_cns_consistency = relaxed_cns_consistency
        # A one-side-flanking transduction call is promoted to high-confidence for SVA.
        self.one_side_flank_high_confidence = one_side_flank_high_confidence

    def __repr__(self):
        return "FamilyParams({0})".format(self.name)


# The bundle a step resolves to when it was never told which family it is running.
# Numerically identical to L1's, which is what `IS_CALL_SVA == False` meant.
GENERIC = FamilyParams(
    name=None, rtype=0, vcf_alt=None, has_flank=True, annotation_extend=0,
    leftmost_cns=500, relaxed_cns_consistency=False,
    one_side_flank_high_confidence=False)

FAMILIES = {
    "L1": FamilyParams(
        name="L1", rtype=1, vcf_alt="LINE1", has_flank=True, annotation_extend=0,
        leftmost_cns=500, relaxed_cns_consistency=False,
        one_side_flank_high_confidence=False),
    "Alu": FamilyParams(
        name="Alu", rtype=2, vcf_alt="ALU", has_flank=False, annotation_extend=0,
        leftmost_cns=500, relaxed_cns_consistency=False,
        one_side_flank_high_confidence=False),
    "SVA": FamilyParams(
        name="SVA", rtype=4, vcf_alt="SVA", has_flank=True, annotation_extend=200,
        leftmost_cns=5000, relaxed_cns_consistency=True,
        one_side_flank_high_confidence=True),
}

# Lowest set bit wins, mirroring x_rep_type.RepType.get_rep_type. A merged bitmask
# (e.g. 7 = L1|Alu|SVA, which the unified tail uses) therefore resolves to L1 -- which
# is why the unified tail must NOT resolve a family this way; it carries a per-locus
# type map instead.
_RTYPE_ORDER = ("L1", "Alu", "SVA")


def for_name(name):
    """Bundle for a family key ("L1"/"Alu"/"SVA"). Raises on an unknown key rather
    than silently handing back GENERIC -- a typo'd family must not read as "no family"."""
    try:
        return FAMILIES[name]
    except KeyError:
        raise KeyError("unknown MEI family {0!r}; known: {1}".format(
            name, ", ".join(sorted(FAMILIES))))


def for_rtype(i_rep_type):
    """Bundle for x_TEA_main's `--rtype` bitmask (L1=1, Alu=2, SVA=4).

    Returns GENERIC for 0 and for masks naming only families this module does not
    model (HERV/MIT/MSTA/pseudogene), because those run with L1's numbers today.
    """
    for key in _RTYPE_ORDER:
        fam = FAMILIES[key]
        if i_rep_type & fam.rtype:
            return fam
    return GENERIC


def from_sva_flag(b_sva):
    """Bundle a step resolves from `x_TEA_main --sva`, which is the ONLY family signal
    the pre-Stage-2 steps had.

    This is the faithful replacement for `if options.sva: global_values.turn_on_sva()`.
    It returns GENERIC -- not the L1 bundle -- when the flag is absent, because absent
    means "this step was not told", and the two are worth telling apart even though
    their numbers coincide today.
    """
    return FAMILIES["SVA"] if b_sva else GENERIC
