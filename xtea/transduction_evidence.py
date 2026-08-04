"""The transduction decision: one per-locus evidence record, and one named gate per threshold.

DOWNSTREAM_OVERHAUL_PLAN.md sec 5 asks for exactly this -- "replace the six-phase veto cascade
with one explicit per-site evidence record feeding one decision function, so every threshold is
a named parameter rather than a buried constant" -- and sec 8 names this module, `SiteEvidence`,
`decide()` and `TransductionParams`. This is that module, built **answer-preserving**: every
predicate here is a faithful copy of the one it replaces, including two defects (marked
DEFECT-PRESERVED below). The thresholds move into `TransductionParams` and the *reasons* become
data; no call moves. Changing any of these numbers is the Stage 3 work this unblocks, and it is
truth-set gated -- see sec 5's "interrogate" list, which names four of the gates below.

WHY IT IS FIVE GATE GROUPS AND NOT ONE FUNCTION, and this is a property of the data flow
-----------------------------------------------------------------------------------------
`decide(evidence, params)` as sec 5 draws it assumes the evidence for a locus is all available
at once. It is not, and no refactor makes it so: **two of the five groups' inputs are produced
by batched BAM passes over sets that earlier groups define.**

    clip SAM parse  -> [A] clip gate ------------.
    disc SAM parse  -> [B] source gate ---> [C] candidate gate
                                                 |
                            (ReadDepth over the surviving candidates ONLY)
                                                 v
                                            [D] sanity gate
                                                 |
                            (the genotyper over the surviving sites ONLY)
                                                 v
                                            [E] genotype gate

`filter_transduction`'s `calc_coverage_of_two_regions` and `_collect_gntp_features`' genotyper
each open the BAM for the set handed to them. Hoisting them ahead of [C] to make one flat
`decide()` would compute local coverage and genotype features for the whole `.slct` list --
34-40K loci per family instead of a few thousand -- which is not a refactor either. So the shape
is five ordered gate groups with two batched fetches between them, and the honest reading is
that "one decision function" means *one place per threshold*, not one call.

WHAT THE GATES ARE FOR (the reason ledger)
------------------------------------------
Every gate returns a `Verdict` carrying a stable `reason` string, and `DecisionLedger` writes
them to a TSV beside the step's output. sec 5: "`reason` is not decoration -- it is what makes
the benchmark interpretable, because every missed truth call gets attributed to the specific
gate that dropped it." Today those reasons exist only as `print()` calls in the step's stdout
("[Transduction filter]:Filtered out: ... fall in repetitive region"), unparseable and
unjoined to the truth set; BASELINE-0's 322 false negatives have no attribution at all. The
ledger is the same information as a joinable file. It is **diagnostic**: nothing reads it, and
the reason strings are the only part of this module that is a public contract.

ORDER IS PART OF THE ANSWER. Every driver here iterates the maps it is handed and builds its
output in that same order, because the step's output files are written in map iteration order
and two runs that agree on every locus still differ in BYTES if they built a map in a different
order -- job 49596361 failed on exactly that (see `downstream_unified.LocusUniverse`). The
differential test in `tests/test_transduction_evidence.py` asserts key ORDER as well as content
for that reason.
"""
from collections import namedtuple

import global_values


# --------------------------------------------------------------------------------------
# The reason vocabulary. These strings are written to the ledger, so they are a contract:
# rename one and every benchmark join against an older run's ledger silently misses.
# --------------------------------------------------------------------------------------
# [A] clip evidence, from the clip-to-flank SAM parse
R_LOW_CLIP = "low_clip_support"            # fewer clipped reads realigning than n_clip_cutoff
R_NO_POLYA_CLIP = "no_polya_clip"          # the locus has no polyA/polyT clip record at all
R_LOW_POLYA_CLIP = "low_polya_clip"        # polyA clip reads below n_polya_cutoff
# [B] disc source, from the disc-to-flank SAM parse
R_DECOY_CROSS_FAMILY = "decoy_cross_family"  # enough disc reads hit ANOTHER family's decoy
R_LOW_DISC = "low_disc_support"             # total disc support over all sources too low
R_NO_DOMINANT_SOURCE = "no_dominant_source"  # no single source dominates, or it is not unique
# [C] the candidate conjunction
R_NO_CLIP_EVIDENCE = "no_clip_evidence"     # not in the clip-selected set
R_NO_DISC_POLYA = "no_disc_polya"           # no polyA among the disc reads realigned to cns
R_NO_DISC_CLUSTER = "no_disc_cluster"       # the dominant side's hits do not form a cluster
# [D] site sanity, needs local coverage
R_SAME_FAMILY_REPEAT = "same_family_repeat"  # sits in a low-divergence copy of its own family
R_COVERAGE_ISLAND = "coverage_island"        # local coverage above i_max_cov on any window
R_TWO_SIDED_POLYA = "two_sided_polya"        # polyA clips on BOTH sides
R_TWO_SIDED_POLYT = "two_sided_polyt"        # polyT clips on BOTH sides
# [E] genotype
R_LARGE_INDEL = "large_indel_excess"         # too many large-indel disc pairs for the coverage

# Which gate group each reason belongs to, for the ledger's `phase` column. A phase is
# meaningful on its own: a locus dropped at A never had disc evidence produced for it, so
# "dropped at A" and "dropped at D" are not comparable strengths of evidence against.
PHASE_OF_REASON = {
    R_LOW_CLIP: "A_clip", R_NO_POLYA_CLIP: "A_clip", R_LOW_POLYA_CLIP: "A_clip",
    R_DECOY_CROSS_FAMILY: "B_source", R_LOW_DISC: "B_source",
    R_NO_DOMINANT_SOURCE: "B_source",
    R_NO_CLIP_EVIDENCE: "C_candidate", R_NO_DISC_POLYA: "C_candidate",
    R_NO_DISC_CLUSTER: "C_candidate",
    R_SAME_FAMILY_REPEAT: "D_sanity", R_COVERAGE_ISLAND: "D_sanity",
    R_TWO_SIDED_POLYA: "D_sanity", R_TWO_SIDED_POLYT: "D_sanity",
    R_LARGE_INDEL: "E_genotype",
}

# `keep` is the answer; `reason` and `detail` are why. detail carries the numbers that decided
# it, so a ledger row is readable without re-running anything.
Verdict = namedtuple("Verdict", "keep reason detail")
KEPT = Verdict(True, "", "")


def _drop(reason, detail=""):
    return Verdict(False, reason, detail)


class TransductionParams:
    """Every threshold the transduction decision uses, named, in one object.

    The derived cutoffs are computed here with the step's own expressions, floats and all:
    `n_polya_cutoff = n_clip_cutoff / 2 - 1` is 1.0 for `--nclip 4`, and `n_half_disc =
    n_disc_cutoff / 2` is 2.5 for `--nd 5`. They are NOT rounded, because `>= 2.5` and `>= 3`
    admit different loci and the step has always used the float. `n_half_disc`'s `<= 1` fallback
    is the step's too (`n_disc_cutoff` itself when halving would reach 1 or less).

    WHAT IS NOT HERE. `i_concord_dist` / `f_concord_ratio` are handed to the cluster checker
    rather than compared against, so they stay arguments of the driver that calls it, and the
    family bundle (`family_params`) stays separate because it is not a threshold set.

    **A THRESHOLD A CALLER DOES NOT HAVE IS `None`, NOT ZERO.** Two of the three constructors
    are partial -- `filter_transduction` and `add_new_rcd_with_cutoff` are shared with the
    sibling step (`x_orphan_transduction`), which reaches those gates from its own per-site
    parse and never runs gates A-C -- and a partial object leaves every threshold it does not
    own as None. A gate reached with a None threshold raises `TypeError` on the comparison
    instead of silently admitting or rejecting everything, which is what a 0 default would do.
    """

    __slots__ = ("n_clip_cutoff", "n_polya_cutoff", "n_disc_cutoff", "n_half_disc",
                 "n_decoy_cutoff", "n_half_clip", "f_multi_source_min_ratio",
                 "f_min_uniq_ratio", "i_max_cov", "ave_cov", "max_indel_ratio",
                 "i_rep_type", "i_orphan_anchor_slack", "i_rep_div_cutoff")

    def __init__(self, n_clip_cutoff=None, n_disc_cutoff=None, i_max_cov=None, ave_cov=None,
                 max_indel_ratio=None, f_min_uniq_ratio=None, i_rep_type=-1,
                 f_multi_source_min_ratio=None, i_rep_div_cutoff=None, n_half_clip=None):
        self.n_clip_cutoff = n_clip_cutoff
        self.n_disc_cutoff = n_disc_cutoff
        # x_transduction.call_candidate_transduction_v3's own two lines.
        if n_clip_cutoff is None:
            self.n_polya_cutoff = None
        else:
            n_polya_cutoff = n_clip_cutoff / 2 - 1
            self.n_polya_cutoff = 0 if n_polya_cutoff < 0 else n_polya_cutoff
        if n_disc_cutoff is None:
            self.n_half_disc = None
            self.n_decoy_cutoff = None
        else:
            n_half_disc = n_disc_cutoff / 2
            self.n_half_disc = n_disc_cutoff if n_half_disc <= 1 else n_half_disc
            # `_chk_decoy_cns_cnt` is called with `n_half_cutoff*2` -- which is n_disc_cutoff
            # again for the normal case and NOT for the `<= 1` fallback, so it is derived
            # from the fallback-adjusted value rather than from n_disc_cutoff.
            self.n_decoy_cutoff = self.n_half_disc * 2
        # The two-sided polyA/polyT threshold: `filter_transduction` is called with
        # `n_clip_cutoff / 2`, a different quantity from n_polya_cutoff above. It is
        # overridable because the sibling step passes its own (0).
        if n_half_clip is not None:
            self.n_half_clip = n_half_clip
        else:
            self.n_half_clip = None if n_clip_cutoff is None else n_clip_cutoff / 2
        self.f_multi_source_min_ratio = (global_values.TRANSDCT_MULTI_SOURCE_MIN_RATIO
                                         if f_multi_source_min_ratio is None
                                         else f_multi_source_min_ratio)
        self.f_min_uniq_ratio = f_min_uniq_ratio
        self.i_max_cov = i_max_cov
        self.ave_cov = ave_cov
        self.max_indel_ratio = max_indel_ratio
        self.i_rep_type = i_rep_type
        # `_slct_orphan_transduction`'s `(anchor_pos + 5) < ins_pos`. Unexplained in the
        # original and left as it is; sec 5's "interrogate" list is where it belongs.
        self.i_orphan_anchor_slack = 5
        self.i_rep_div_cutoff = (global_values.TD_REP_DIVERGENT_CUTOFF
                                 if i_rep_div_cutoff is None else i_rep_div_cutoff)

    @classmethod
    def for_site_sanity(cls, i_max_cov, n_half_clip, i_rep_type):
        """Gate D only -- what `filter_transduction`'s two callers have.

        The sibling step passes `n_half_clip = 0` AND an empty polyA map, so its two-sided
        polyA/polyT gates cannot fire either way; passing its 0 through keeps that visible
        rather than implicit.
        """
        return cls(i_max_cov=i_max_cov, n_half_clip=n_half_clip, i_rep_type=i_rep_type)

    @classmethod
    def for_genotype_gate(cls, ave_cov, max_indel_ratio):
        """Gate E only -- what `add_new_rcd_with_cutoff`'s two callers have."""
        return cls(ave_cov=ave_cov, max_indel_ratio=max_indel_ratio)

    @classmethod
    def for_disc_sources(cls, n_half_disc, f_min_uniq_ratio, i_rep_type):
        """Gate B only, from the ALREADY-HALVED cutoff its caller holds.

        `parse_transduction_disc_algnmt` is handed `n_half_cutoff`, not `n_disc_cutoff`, and
        derives the decoy cutoff as `n_half_cutoff * 2`. Taking the halved value directly
        reproduces that exactly instead of un-halving and re-halving it, which is not a
        round trip: with `--nd 1` the `<= 1` fallback fires on the way back and turns 1 into 2.
        """
        p = cls(f_min_uniq_ratio=f_min_uniq_ratio, i_rep_type=i_rep_type)
        p.n_half_disc = n_half_disc
        p.n_decoy_cutoff = n_half_disc * 2
        return p


class SiteEvidence:
    """One locus's evidence, as the gates need it. Groups are sec 5's, filled per phase.

    A field left None means "not produced yet at this phase" and any gate that reads it says
    which phase produces it. That is deliberate rather than sloppy: the two batched BAM fetches
    (see the module docstring) make partial evidence the real state of a locus for most of the
    decision, and a record that pretended otherwise would have to be built twice.

    NO pysam, NO annotation index, NO BAM in this module. `b_in_same_family_repeat` arrives as
    a bool and the four coverages as ints, so every gate is unit-testable with plain numbers --
    which is what lets the differential test drive thousands of randomized loci in a second.
    """

    __slots__ = ("chrm", "pos",
                 # clip: from parse_transduction_clip_algnmt. Its polyA record is six
                 # counters -- [0,1] left/right polyA, [2,3] left/right polyT, [4,5]
                 # left/right ambiguous (both) -- of which gate A sums all six and gate D
                 # reads only the first four.
                 "n_clip", "n_polya_clip", "n_lpolya", "n_rpolya", "n_lpolyt", "n_rpolyt",
                 # disc/source: from parse_transduction_disc_algnmt
                 "sources", "n_non_unique_by_source", "n_decoy_alu", "n_decoy_l1",
                 "n_decoy_sva", "b_disc_polya",
                 # cluster
                 "b_l_transduct", "b_r_transduct", "n_disc_transduct", "s_source",
                 # sanity
                 "i_lcov1", "i_rcov1", "i_lcov2", "i_rcov2", "b_in_same_family_repeat",
                 # genotype
                 "n_disc_large_indel")

    def __init__(self, chrm=None, pos=None, n_clip=None, n_polya_clip=None,
                 n_lpolya=0, n_rpolya=0, n_lpolyt=0, n_rpolyt=0,
                 sources=None, n_non_unique_by_source=None, n_decoy_alu=0, n_decoy_l1=0,
                 n_decoy_sva=0, b_disc_polya=False, b_l_transduct=False,
                 b_r_transduct=False, n_disc_transduct=0, s_source=None,
                 i_lcov1=0, i_rcov1=0, i_lcov2=0, i_rcov2=0,
                 b_in_same_family_repeat=False, n_disc_large_indel=0):
        self.chrm = chrm
        self.pos = pos
        self.n_clip = n_clip
        self.n_polya_clip = n_polya_clip
        self.n_lpolya = n_lpolya
        self.n_rpolya = n_rpolya
        self.n_lpolyt = n_lpolyt
        self.n_rpolyt = n_rpolyt
        self.sources = sources
        self.n_non_unique_by_source = n_non_unique_by_source
        self.n_decoy_alu = n_decoy_alu
        self.n_decoy_l1 = n_decoy_l1
        self.n_decoy_sva = n_decoy_sva
        self.b_disc_polya = b_disc_polya
        self.b_l_transduct = b_l_transduct
        self.b_r_transduct = b_r_transduct
        self.n_disc_transduct = n_disc_transduct
        self.s_source = s_source
        self.i_lcov1 = i_lcov1
        self.i_rcov1 = i_rcov1
        self.i_lcov2 = i_lcov2
        self.i_rcov2 = i_rcov2
        self.b_in_same_family_repeat = b_in_same_family_repeat
        self.n_disc_large_indel = n_disc_large_indel

    def __repr__(self):
        return "SiteEvidence({0}:{1})".format(self.chrm, self.pos)


# --------------------------------------------------------------------------------------
# [A] The clip gate -- was XTransduction.filter_by_clip_algnmt's tail
# --------------------------------------------------------------------------------------

def gate_clip_evidence(ev, params):
    """Enough clipped reads realigning to the flank library, with a polyA signature.

    `n_polya_clip` is the sum of all SIX counters the clip parse keeps -- left/right polyA,
    left/right polyT, and left/right AMBIGUOUS (a clip that reads as both polyA and polyT,
    which `parse_transduction_clip_algnmt` tallies separately in fields 4 and 5) -- not just
    the two polyA ones. That breadth is what makes this a *signature* test rather than a
    strand-specific one. `n_polya_clip is None` means the locus is absent from the polyA map
    entirely, which the original treats as a drop rather than as zero, and the two are worth
    telling apart in the ledger.
    """
    if ev.n_clip < params.n_clip_cutoff:
        return _drop(R_LOW_CLIP, "n_clip={0} < {1}".format(ev.n_clip, params.n_clip_cutoff))
    if ev.n_polya_clip is None:
        return _drop(R_NO_POLYA_CLIP)
    if ev.n_polya_clip < params.n_polya_cutoff:
        return _drop(R_LOW_POLYA_CLIP, "n_polya={0} < {1}".format(
            ev.n_polya_clip, params.n_polya_cutoff))
    return KEPT


def select_clip_sites(m_lclip, m_rclip, m_polya, params, ledger=None):
    """{chrm: {pos: (n_clip, n_polya)}} -- the clip-selected loci, in the original's order.

    The two-pass accumulation over m_lclip then m_rclip is kept exactly: it is what decides
    the output's key order (left-map chromosomes first, then chromosomes only the right map
    has), and the counts are summed across both sides before the cutoff, not tested per side.
    """
    m_clip_cnt = {}
    for m_side in (m_lclip, m_rclip):
        for ins_chrm in m_side:
            if ins_chrm not in m_clip_cnt:
                m_clip_cnt[ins_chrm] = {}
            for ins_pos in m_side[ins_chrm]:
                n_side = 0
                for ref_mpos in m_side[ins_chrm][ins_pos]:
                    n_side += len(m_side[ins_chrm][ins_pos][ref_mpos])
                if ins_pos not in m_clip_cnt[ins_chrm]:
                    m_clip_cnt[ins_chrm][ins_pos] = 0
                m_clip_cnt[ins_chrm][ins_pos] += n_side

    m_slct = {}
    for ins_chrm in m_clip_cnt:
        for ins_pos in m_clip_cnt[ins_chrm]:
            rcd = m_polya.get(ins_chrm, {}).get(ins_pos)
            ev = SiteEvidence(chrm=ins_chrm, pos=ins_pos,
                              n_clip=m_clip_cnt[ins_chrm][ins_pos],
                              n_polya_clip=None if rcd is None else sum(rcd[:6]))
            verdict = gate_clip_evidence(ev, params)
            if ledger is not None:
                ledger.record(ev, verdict)
            if not verdict.keep:
                continue
            if ins_chrm not in m_slct:
                m_slct[ins_chrm] = {}
            m_slct[ins_chrm][ins_pos] = (ev.n_clip, ev.n_polya_clip)
    return m_slct


# --------------------------------------------------------------------------------------
# [B] The source gate -- was _chk_decoy_cns_cnt + _pick_trsdct_rescued
# --------------------------------------------------------------------------------------

def gate_decoy_cross_family(ev, params):
    """Reject a locus whose disc reads prefer ANOTHER family's decoy sequence.

    Only L1 (bit 1) and SVA (bit 4) have arms, and each rejects on the other two families --
    there is no Alu arm because Alu has no transduction flank and this step never runs for it
    (`is_matched_rep_type` has the same asymmetry; `tests/test_site_universe.py` pins it).
    `i_rep_type < 0` means "family unknown" and disables the gate entirely, which is how the
    original's early return behaves.

    sec 5 calls this "the three-way decoy tally" and lists it as a likely implementation
    artifact to interrogate. It is also the evidence that SVA/Alu mapping competition is real
    (SVA is SINE-VNTR-*Alu*), so it is not a candidate for silent deletion.
    """
    if params.i_rep_type < 0:
        return KEPT
    n = params.n_decoy_cutoff
    if params.i_rep_type & 1 != 0:                      # L1
        if ev.n_decoy_alu >= n:
            return _drop(R_DECOY_CROSS_FAMILY, "ALU={0} >= {1}".format(ev.n_decoy_alu, n))
        if ev.n_decoy_sva >= n:
            return _drop(R_DECOY_CROSS_FAMILY, "SVA={0} >= {1}".format(ev.n_decoy_sva, n))
    elif params.i_rep_type & 4 != 0:                    # SVA
        if ev.n_decoy_alu >= n:
            return _drop(R_DECOY_CROSS_FAMILY, "ALU={0} >= {1}".format(ev.n_decoy_alu, n))
        if ev.n_decoy_l1 >= n:
            return _drop(R_DECOY_CROSS_FAMILY, "L1={0} >= {1}".format(ev.n_decoy_l1, n))
    return KEPT


def select_by_decoy(m_transduction, m_alu, m_l1, m_sva, params, ledger=None):
    """The decoy tally, per locus. Returns m_transduction itself when the family is unknown.

    Returning the same object rather than a copy is the original's behaviour and matters: the
    caller then hands one map to two consumers.
    """
    if params.i_rep_type < 0:
        return m_transduction
    m_slct = {}
    for ins_chrm in m_transduction:
        for ins_pos in m_transduction[ins_chrm]:
            ev = SiteEvidence(chrm=ins_chrm, pos=ins_pos,
                              n_decoy_alu=m_alu.get(ins_chrm, {}).get(ins_pos, 0),
                              n_decoy_l1=m_l1.get(ins_chrm, {}).get(ins_pos, 0),
                              n_decoy_sva=m_sva.get(ins_chrm, {}).get(ins_pos, 0))
            verdict = gate_decoy_cross_family(ev, params)
            if ledger is not None:
                ledger.record(ev, verdict)
            if not verdict.keep:
                continue
            if ins_chrm not in m_slct:
                m_slct[ins_chrm] = {}
            m_slct[ins_chrm][ins_pos] = m_transduction[ins_chrm][ins_pos]
    return m_slct


def dominant_source(sources):
    """(max_source, max_cnt, n_all) over {source: [hit records]}.

    Strictly-greater comparison, so the FIRST source wins a tie -- which makes the answer
    depend on `sources`' insertion order, i.e. on SAM record order. Preserved as-is; the
    alternative (sort by source name) is a call-moving change.
    """
    max_cnt = 0
    max_source = ""
    n_all = 0
    for tmp_source in sources:
        tmp_cnt = len(sources[tmp_source])
        n_all += tmp_cnt
        if tmp_cnt > max_cnt:
            max_cnt = tmp_cnt
            max_source = tmp_source
    return max_source, max_cnt, n_all


def gate_disc_source(ev, params):
    """One source must dominate the locus's disc support, and be uniquely mapped enough.

    Three named numbers, all of which sec 5 keeps as load-bearing: the total support cutoff
    (`n_half_disc`), the source dominance ratio (0.45,
    `global_values.TRANSDCT_MULTI_SOURCE_MIN_RATIO`), and the unique/non-unique ratio (0.5,
    `f_min_uniq_ratio`, strict `>`). The non-unique count is the reads that named the SAME
    dominant source at MAPQ <= TRANSDCT_UNIQ_MAPQ; the disc parse keeps them aside precisely
    so this ratio can be formed.

    The original prints "no enough disc support" and "ratio is low" here; those two prints are
    this gate's two reasons.
    """
    max_source, max_cnt, n_all = dominant_source(ev.sources)
    if n_all < params.n_half_disc:
        return _drop(R_LOW_DISC, "n_disc={0} < {1}".format(n_all, params.n_half_disc))
    f_ratio = float(max_cnt) / float(n_all)
    n_non_unique = 0
    if ev.n_non_unique_by_source:
        n_non_unique = ev.n_non_unique_by_source.get(max_source, 0)
    if (max_cnt >= params.n_half_disc
            and f_ratio >= params.f_multi_source_min_ratio
            and float(max_cnt) / float(max_cnt + n_non_unique) > params.f_min_uniq_ratio):
        return KEPT
    return _drop(R_NO_DOMINANT_SOURCE,
                 "max={0}/{1} ratio={2:.3f} non_unique={3}".format(
                     max_cnt, n_all, f_ratio, n_non_unique))


def pick_dominant_sources(m_transduction, m_raw_td_src, params, ledger=None):
    """{chrm: {pos: (source, [hit records], n_support)}} for the loci with one dominant source.

    The original also built a second map of the *rejected* multi-source loci and its caller
    discards it (`m_rscued_no_use`); it is not rebuilt here. That map was also wrong -- its
    accumulation sat inside `if ins_chrm not in m_n_source_picked`, so only the first locus per
    chromosome was ever recorded -- which is worth knowing if anyone revives the idea. The
    ledger's `no_dominant_source` rows are the same population, correctly.
    """
    m_picked = {}
    for ins_chrm in m_transduction:
        for ins_pos in m_transduction[ins_chrm]:
            ev = SiteEvidence(chrm=ins_chrm, pos=ins_pos,
                              sources=m_transduction[ins_chrm][ins_pos],
                              n_non_unique_by_source=m_raw_td_src.get(
                                  ins_chrm, {}).get(ins_pos, {}))
            verdict = gate_disc_source(ev, params)
            if ledger is not None:
                ledger.record(ev, verdict)
            if not verdict.keep:
                continue
            max_source, max_cnt, _n_all = dominant_source(ev.sources)
            if ins_chrm not in m_picked:
                m_picked[ins_chrm] = {}
            l_pos = [rcd for rcd in ev.sources[max_source]]
            m_picked[ins_chrm][ins_pos] = (max_source, l_pos, max_cnt)
    return m_picked


def select_orphans(m_transduction, m_picked, params):
    """Loci whose dominant source is supported from BOTH sides -- an orphan transduction.

    Orphan = transduced sequence landed WITHOUT the element body, so the anchors bracket the
    insertion point instead of sitting on one side of it. Not a gate: a locus that fails this
    is still a normal transduction candidate, so nothing is dropped and nothing is recorded in
    the ledger. `(anchor_pos + slack) < ins_pos` puts an anchor on the left; everything else
    counts right, including anchors exactly at the position.

    The original's own comment -- "here miss those orphan promoted SV case!!!" -- is about this
    two-sided requirement. sec 5 folds orphan/sibling/regular into one call class with a typed
    source; this function is the piece that would become `source == inferred_region`.
    """
    m_orphan = {}
    for ins_chrm in m_transduction:
        for ins_pos in m_transduction[ins_chrm]:
            if (ins_chrm not in m_picked) or (ins_pos not in m_picked[ins_chrm]):
                continue
            n_left = 0
            n_right = 0
            tmp_source = m_picked[ins_chrm][ins_pos][0]
            for rcd in m_transduction[ins_chrm][ins_pos][tmp_source]:
                anchor_pos = int(rcd[1])
                if (anchor_pos + params.i_orphan_anchor_slack) < ins_pos:
                    n_left += 1
                else:
                    n_right += 1
            if (n_left >= params.n_half_disc) and (n_right >= params.n_half_disc):
                if ins_chrm not in m_orphan:
                    m_orphan[ins_chrm] = {}
                m_orphan[ins_chrm][ins_pos] = (tmp_source, n_left, n_right)
    return m_orphan


# --------------------------------------------------------------------------------------
# [C] The candidate gate -- was call_out_candidate_transduct
# --------------------------------------------------------------------------------------

def dominant_side_cluster(l_trsdct_pos, clip_pos, f_is_cluster):
    """(b_l, b_r, trsdct_l, trsdct_r, n_disc) for the side most of the hits are on.

    **DEFECT-PRESERVED, verbatim.** The cluster test sits INSIDE the per-read loop in the
    original, so:

      * it runs once per read instead of once per locus (O(n^2) sorting on a hot path), and
      * only the LAST read's evaluation survives for the side it picked, while the other
        side's flag keeps whatever an earlier iteration left there. A locus can therefore
        come out with `b_l_trsdct` from read 3 and `trsdct_l/trsdct_r` from read 9, i.e.
        boundaries that belong to the other side.

    This is copied exactly rather than lifted out of the loop because lifting it MOVES CALLS:
    the "both sides true" outcome, which `collect_ins_info` reads to decide which side's cns
    interval to blank, is reachable only through the stale flag. It is a Stage 3 correction
    with a truth-set measurement behind it, not a cleanup. `f_is_cluster(l_pos)` is the
    cluster checker with its distance/ratio already bound, and it SORTS ITS ARGUMENT IN PLACE
    -- also part of the reproduced behaviour, since the list is appended to between calls.
    """
    b_l = False
    b_r = False
    trsdct_l = -1
    trsdct_r = -1
    n_disc = 0
    l_l_pos = []
    l_r_pos = []
    for (hit_ref_pos, anchor_pos, _dist_from_rep) in l_trsdct_pos:
        if anchor_pos < clip_pos:
            l_l_pos.append(hit_ref_pos)
        else:
            l_r_pos.append(hit_ref_pos)
        if len(l_l_pos) > len(l_r_pos):
            b_l, trsdct_l, trsdct_r = f_is_cluster(l_l_pos)
            n_disc = len(l_l_pos)
        else:
            b_r, trsdct_l, trsdct_r = f_is_cluster(l_r_pos)
            n_disc = len(l_r_pos)
    return b_l, b_r, trsdct_l, trsdct_r, n_disc


def gate_candidate(ev, params):
    """The conjunction that turns evidence into a candidate: clip AND disc-polyA AND cluster.

    Nothing about the numbers here -- it is three memberships, and the reasons are what make
    the strata visible. `no_disc_polya` in particular is the polyA among the disc reads
    realigned to the CONSENSUS, a different observation from gate A's clip polyA, and it drops
    the largest share of loci at this phase.
    """
    if ev.n_clip is None:
        return _drop(R_NO_CLIP_EVIDENCE)
    if not ev.b_disc_polya:
        return _drop(R_NO_DISC_POLYA)
    if not (ev.b_l_transduct or ev.b_r_transduct):
        return _drop(R_NO_DISC_CLUSTER, "n_disc={0}".format(ev.n_disc_transduct))
    return KEPT


def call_candidates(m_repsnt_pos, m_disc_transduct, m_clip_slct, m_disc_polya, m_orphan,
                    f_is_cluster, params, ledger=None):
    """{chrm: {pos: (b_l, b_r, n_disc, "source~l-r")}} -- the candidate transductions.

    Two passes, as the original: the clustered candidates over `m_repsnt_pos`, then the
    orphans over `m_clip_slct`. An orphan is written unconditionally as `(False, True, ...)`
    with a `~orphan` source tag and can OVERWRITE a clustered record for the same locus --
    preserved, and it is why `collect_ins_info` reads a right-side record for every orphan.

    A locus in `m_repsnt_pos` that has no `m_disc_transduct` entry produces no record and no
    ledger row: it was already dropped at gate B, which recorded it there.
    """
    m_transduction = {}
    for ins_chrm in m_repsnt_pos:
        for ins_pos in m_repsnt_pos[ins_chrm]:
            clip_rcd = m_clip_slct.get(ins_chrm, {}).get(ins_pos)
            ev = SiteEvidence(chrm=ins_chrm, pos=ins_pos,
                              n_clip=None if clip_rcd is None else clip_rcd[0],
                              b_disc_polya=(ins_pos in m_disc_polya.get(ins_chrm, {})))
            if clip_rcd is None or not ev.b_disc_polya:
                # The original tests these two before touching the disc map, so a locus
                # failing them is never cluster-evaluated. Same order here.
                verdict = gate_candidate(ev, params)
                if ledger is not None:
                    ledger.record(ev, verdict)
                continue
            if (ins_chrm not in m_disc_transduct) or (ins_pos not in m_disc_transduct[ins_chrm]):
                # UNREACHABLE from the production caller, so it records nothing rather than
                # inventing a reason for a population that does not exist. `m_repsnt_pos` is
                # derived FROM `m_disc_transduct`: the step builds it as
                # `prep_representative_pos(m_disc_filtered)`, and `m_disc_filtered` comes from
                # `check_seprt_disc_consistency(prep_representative_pos(m_disc_transduct), …)`,
                # which only ever filters its first argument. So `m_repsnt_pos`'s keys are a
                # subset of this map's -- measured on the WGS gate (job 49634290): L1
                # 24-of-159 and SVA 56-of-105 loci survive to `m_disc_filtered` and **0** of
                # either are absent from `m_disc_transduct`.
                #
                # Kept as a guard because the containment is a property of the CALLER, not of
                # this function, and a future caller could hand in an independent site list.
                # `test_a_locus_absent_from_the_source_map_is_skipped_not_recorded` pins the
                # behaviour so it stays a deliberate skip rather than becoming a KeyError.
                continue
            trsdct_source, l_trsdct_pos, _n_unique = m_disc_transduct[ins_chrm][ins_pos]
            clip_pos = m_repsnt_pos[ins_chrm][ins_pos]
            (ev.b_l_transduct, ev.b_r_transduct, trsdct_l, trsdct_r,
             ev.n_disc_transduct) = dominant_side_cluster(l_trsdct_pos, clip_pos, f_is_cluster)
            verdict = gate_candidate(ev, params)
            if ledger is not None:
                ledger.record(ev, verdict)
            if not verdict.keep:
                continue
            if ins_chrm not in m_transduction:
                m_transduction[ins_chrm] = {}
            m_transduction[ins_chrm][ins_pos] = (
                ev.b_l_transduct, ev.b_r_transduct, ev.n_disc_transduct,
                "{0}~{1}-{2}".format(trsdct_source, trsdct_l, trsdct_r))

    for ins_chrm in m_clip_slct:
        for ins_pos in m_clip_slct[ins_chrm]:
            if (ins_chrm in m_orphan) and (ins_pos in m_orphan[ins_chrm]):
                if ins_chrm not in m_transduction:
                    m_transduction[ins_chrm] = {}
                s_source = m_orphan[ins_chrm][ins_pos][0]
                n_disc = m_orphan[ins_chrm][ins_pos][1] + m_orphan[ins_chrm][ins_pos][2]
                m_transduction[ins_chrm][ins_pos] = (
                    False, True, n_disc, "{0}~orphan".format(s_source))
    return m_transduction


# --------------------------------------------------------------------------------------
# [D] The sanity gate -- was filter_transduction's body
# --------------------------------------------------------------------------------------

def gate_site_sanity(ev, params):
    """Repeat context, coverage island, and the two-sided polyA/polyT veto.

    All four numbers are sec 5 items. The first two it keeps ("repeat-divergence guard",
    "coverage and large-indel sanity"); the polyA/polyT veto it flags as "a blunt tandem-dup
    filter -- make conditional on source strength", which is a Stage 3 change.

    The coverage test is `>` against `i_max_cov` on all FOUR windows, and the two window sizes
    are not what `filter_transduction`'s comments say: it passes `focal_win = LOCAL_COV_WIN`
    (**200**) commented "search for coverage island, by default 900" and `focal_win2 =
    COV_ISD_CHK_WIN` (**900**) commented "calculate the local coverage, by default 200". The
    names and the comments disagree; the VALUES are 200 for cov1 and 900 for cov2, and cov2 is
    what `collect_ins_info` reports as the site's local coverage. Fixing the comments there
    is safe; swapping the constants would move calls.
    """
    if ev.b_in_same_family_repeat:
        return _drop(R_SAME_FAMILY_REPEAT)
    if (ev.i_lcov1 > params.i_max_cov or ev.i_rcov1 > params.i_max_cov
            or ev.i_lcov2 > params.i_max_cov or ev.i_rcov2 > params.i_max_cov):
        return _drop(R_COVERAGE_ISLAND, "cov={0},{1},{2},{3} > {4}".format(
            ev.i_lcov1, ev.i_rcov1, ev.i_lcov2, ev.i_rcov2, params.i_max_cov))
    # Strict `>` on both sides, and only when the locus HAS a clip-polyA record: absent means
    # the gate does not apply, not that the counts are zero.
    if ev.n_lpolya is not None:
        if ev.n_lpolya > params.n_half_clip and ev.n_rpolya > params.n_half_clip:
            return _drop(R_TWO_SIDED_POLYA, "l={0} r={1} > {2}".format(
                ev.n_lpolya, ev.n_rpolya, params.n_half_clip))
        if ev.n_lpolyt > params.n_half_clip and ev.n_rpolyt > params.n_half_clip:
            return _drop(R_TWO_SIDED_POLYT, "l={0} r={1} > {2}".format(
                ev.n_lpolyt, ev.n_rpolyt, params.n_half_clip))
    return KEPT


def apply_site_filters(m_ori_ins, m_read_depth, m_clip_polya, f_in_same_family_repeat,
                       params, ledger=None):
    """Filter a candidate map by site sanity, carrying each locus's payload through unchanged.

    Payload-agnostic on purpose: the transduction step's values are
    `(b_l, b_r, n_disc, source)` tuples and the sibling step's are
    `(b_l, b_r, src_chrm, src_pos, nspot)`, and neither is inspected here.

    `m_read_depth` must cover every locus in `m_ori_ins` that is NOT in a same-family repeat --
    it is the batched `calc_coverage_of_two_regions` result, and a missing locus is a bug rather
    than a zero, so it raises the same KeyError the original would. The repeat test is
    evaluated FIRST and the depth read skipped when it fires, which is the original's
    short-circuit: a repeat-dropped locus is never looked up. Its four coverages are then left
    at 0 in the record, and nothing reads them -- the gate has already returned.
    """
    m_new_ins = {}
    for ins_chrm in m_ori_ins:
        for ins_pos in m_ori_ins[ins_chrm]:
            polya = m_clip_polya.get(ins_chrm, {}).get(ins_pos)
            b_in_rep = f_in_same_family_repeat(ins_chrm, ins_pos)
            depth = (0, 0, 0, 0) if b_in_rep else m_read_depth[ins_chrm][ins_pos]
            ev = SiteEvidence(
                chrm=ins_chrm, pos=ins_pos,
                i_lcov1=depth[0], i_rcov1=depth[1], i_lcov2=depth[2], i_rcov2=depth[3],
                n_lpolya=None if polya is None else polya[0],
                n_rpolya=None if polya is None else polya[1],
                n_lpolyt=None if polya is None else polya[2],
                n_rpolyt=None if polya is None else polya[3],
                b_in_same_family_repeat=b_in_rep)
            verdict = gate_site_sanity(ev, params)
            if ledger is not None:
                ledger.record(ev, verdict)
            if not verdict.keep:
                continue
            if ins_chrm not in m_new_ins:
                m_new_ins[ins_chrm] = {}
            m_new_ins[ins_chrm][ins_pos] = m_ori_ins[ins_chrm][ins_pos]
    return m_new_ins


# --------------------------------------------------------------------------------------
# [E] The genotype gate -- was add_new_rcd_with_cutoff's one `continue`
# --------------------------------------------------------------------------------------

def gate_genotype(ev, params):
    """Too many large-indel discordant pairs for this coverage means a structural artifact.

    `ave_cov * max_indel_ratio` (0.2) is sec 5's "large-indel sanity", kept. Note the
    threshold is the SAMPLE's average coverage, not the site's local coverage -- the two
    differ by an order of magnitude at a coverage island, and using the local one would make
    the gate self-defeating there.
    """
    if ev.n_disc_large_indel > (params.ave_cov * params.max_indel_ratio):
        return _drop(R_LARGE_INDEL, "n={0} > {1:.2f}".format(
            ev.n_disc_large_indel, params.ave_cov * params.max_indel_ratio))
    return KEPT


# --------------------------------------------------------------------------------------
# The ledger
# --------------------------------------------------------------------------------------

class DecisionLedger:
    """Every dropped locus and the named gate that dropped it, as a TSV.

    Drops only, deliberately: a locus that survives every gate is in the step's own output
    already, so recording it here would double the file to say nothing new. A locus can appear
    more than once -- gates A and D both see loci that gate C never reaches -- and the `phase`
    column is what makes those rows comparable.

    Diagnostic. Nothing in the pipeline reads it, and it is written next to the step's output
    rather than into `merged/` because the drops are per family and a locus's fate differs
    between families.
    """

    __slots__ = ("rows", "counts")

    def __init__(self):
        self.rows = []
        self.counts = {}

    def record(self, ev, verdict):
        if verdict.keep:
            return
        self.counts[verdict.reason] = self.counts.get(verdict.reason, 0) + 1
        self.rows.append((ev.chrm, ev.pos, PHASE_OF_REASON.get(verdict.reason, "?"),
                          verdict.reason, verdict.detail))

    def write(self, sf_out):
        with open(sf_out, "w") as fout:
            fout.write("#chrom\tpos\tphase\treason\tdetail\n")
            for (chrm, pos, phase, reason, detail) in self.rows:
                fout.write("{0}\t{1}\t{2}\t{3}\t{4}\n".format(chrm, pos, phase, reason, detail))
        return sf_out

    def summary(self):
        """"reason=n" per reason, phase order then reason name -- one log line."""
        return " ".join("{0}={1}".format(r, self.counts[r]) for r in sorted(
            self.counts, key=lambda r: (PHASE_OF_REASON.get(r, "?"), r)))
