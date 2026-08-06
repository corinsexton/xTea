"""The Stage 3 transduction gate's STRATIFICATION and its NO-REGRESSION rule.

WHY THESE AND NOT OTHERS
------------------------
`benchmark/verify_transduction.py` is mostly an alignment pipeline -- pysam fetches, two bwa
passes -- and those parts have a harness, not unit tests. What is testable is what decides
pass/fail, and there are exactly three places this gate can be quietly wrong:

1. **The strata.** `HANDOFF_2026-08-04_td_gate_built.md` sec 10.2 makes
   `family x call class x shape x sidedness` primary because the source-class axis is defined for
   7 of 115 entries. Both new facts are derived, so a derivation that silently disagrees with
   `benchmark_truth_set.td_src_shape` would restratify the gate without saying so.

2. **Reading the FROZEN baseline.** It was written before the shape axis existed, so it has no
   `shape` or `sidedness` column. `read_baseline_verdicts` must re-derive both from `td_src_raw`
   and `src_start`. If it ever read the columns instead, the frozen baseline would silently
   compare as all-empty strata and the gate would PASS on anything.

3. **The pass/fail rule.** sec 10.4 gates on ABSOLUTE VERIFIED counts per stratum, and the
   failure mode that matters is a gain in one stratum paying for a loss in another -- which is
   exactly what a total-count or percentage comparison would allow. `SVA regular reference_copy`
   holds ONE entry.

The operating-point refusal is tested through `format_op_point` / `_OP_POINT_RE` rather than
through `main`, because that round-trip IS the property: a header this program writes has to be
a header it can read back, or every baseline comparison is a comparison across thresholds.
"""
import os
import sys

REPO_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(REPO_ROOT, "benchmark"))
import benchmark_truth_set as bts  # noqa: E402
import verify_transduction as vt  # noqa: E402


# --- 1. the strata -----------------------------------------------------------------------

def _row(raw, src_start, family="L1", call_class="regular", source_class="disc_span"):
    """A PoolRow from the pool TSV's own column names, so the test exercises the real reader."""
    return vt.PoolRow(0, {
        "family": family, "ins_chrom": "chr1", "ins_pos": "1000",
        "call_class": call_class, "source_class": source_class,
        "src_chrom": "chr4", "src_start": str(src_start), "src_end": "138066924",
        "arms": "master", "td_src_raw": raw,
    })


def test_shape_comes_from_the_truth_set_harness_and_not_from_a_second_copy():
    """sec 1's rule: one implementation of the rule, or a probe can agree with a scorer that
    disagrees with it. Pinned by equality against `bts`, not by restating the shape names."""
    for raw in ("chr4:-1-138066924", "chr8:51817625-51818695",
                "chr6:32546845-32546845~0~orphan", "chrY:9591830~sibling"):
        assert _row(raw, 0).shape == bts.td_src_shape(raw)


def test_sidedness_is_src_start_and_it_is_the_axis_that_discriminates():
    """`src_start = -1` means one of the two disc-pair midpoints was never established. It
    verifies ~30 points worse (55.0% vs 83.9%), which is why it is a stratum and not a footnote."""
    assert _row("chr4:-1-138066924", -1).sidedness == "one_sided"
    assert _row("chr8:51817625-51818695", 51817625).sidedness == "two_sided"


def test_a_sibling_point_has_no_sidedness_because_its_bounds_were_invented():
    """The pool builder gives a sibling xTea's own [POS-100, POS+100]
    (`x_orphan_transduction.py:1190-1191`), so its `src_start` is a window bound and not a
    midpoint the caller established. Reading sidedness off it would report noise as an axis --
    and it would read `two_sided`, i.e. the CONFIDENT value, which is the dangerous direction."""
    r = _row("chrY:9591830~sibling", 9591730, call_class="sibling", source_class="-")
    assert r.shape == bts.SHAPE_SIBLING_PT
    assert r.sidedness == "-"
    assert r.shape_stratum == ("L1", "sibling", "sibling_point", "-")


def test_the_two_strata_are_independent_axes():
    """A bare one-sided entry and a bare two-sided one share a source class and differ in the
    primary stratum -- which is the whole reason sec 10.2 promoted shape x sidedness."""
    a = _row("chr4:-1-138066924", -1)
    b = _row("chr8:51817625-51818695", 51817625)
    assert a.source_stratum == b.source_stratum == ("L1", "regular", "disc_span")
    assert a.shape_stratum != b.shape_stratum


# --- 2. reading the frozen baseline -----------------------------------------------------

# The FROZEN format: no `shape` column, no `sidedness` column. Written verbatim in the column
# order `td_gate_baseline_hg002/verdicts.tsv` uses, because the point of this fixture is that a
# file predating the axis still stratifies.
FROZEN_HEADER = ("family\tins_chrom\tins_pos\tcall_class\tsource_class\tsrc_chrom\tsrc_start\t"
                 "src_end\tarms\tdepth\trep_min_div\tn_clip_support\tn_disc_support\tverdict\t"
                 "why\ttd_src_raw")
FROZEN_ROWS = (
    # bare + two-sided, VERIFIED
    "L1\tchr1\t100\tregular\tdisc_span\tchr4\t138066900\t138066924\tmaster\t50\tnone\t0\t3\t"
    "VERIFIED\tsupport=3\tchr4:138066900-138066924",
    # bare + ONE-sided, VERIFIED -- same source class, different primary stratum
    "L1\tchr1\t200\tregular\tdisc_span\tchr4\t-1\t138066924\tmaster\t50\tnone\t0\t3\t"
    "VERIFIED\tsupport=3\tchr4:-1-138066924",
    # REFUTED: must not reach either count
    "L1\tchr1\t300\tregular\tdisc_span\tchr4\t-1\t138066924\tmaster\t50\tnone\t0\t0\t"
    "REFUTED\tsupport=0\tchr4:-1-138066924",
    # tailed, the only shape whose source class means anything
    "SVA\tchr17\t45593489\torphan\treference_copy\tchr17\t46252374\t46253788\tdecoy\t40\t8.0\t0\t9"
    "\tVERIFIED\tsupport=9\tchr17:46252374-46253788~0~orphan",
)


def _frozen(tmp_path, point="N=3 f_cns=0.8 f_src=0.65 cov>=5 rep_div<5.0"):
    p = tmp_path / "verdicts.tsv"
    p.write_text("# per-entry verdicts at {0}.\n# THE SWEEP is the result.\n{1}\n{2}\n".format(
        point, FROZEN_HEADER, "\n".join(FROZEN_ROWS)))
    return str(p)


def test_the_frozen_baseline_stratifies_without_a_shape_column(tmp_path):
    """THE LOAD-BEARING ONE. The frozen file has no `shape`/`sidedness` column, so both are
    re-derived from `td_src_raw` and `src_start`. A reader that took them from columns would see
    every stratum as empty here and the gate would PASS on anything."""
    point, shape_v, source_v, n = vt.read_baseline_verdicts(_frozen(tmp_path))
    assert n == 4
    assert point == "N=3 f_cns=0.8 f_src=0.65 cov>=5 rep_div<5.0"
    assert shape_v == {
        ("L1", "regular", "bare", "two_sided"): 1,
        ("L1", "regular", "bare", "one_sided"): 1,
        ("SVA", "orphan", "tailed", "two_sided"): 1,
    }
    # the REFUTED row is in `n` and in NEITHER count: decision 5's denominator rule
    assert sum(shape_v.values()) == 3
    assert source_v == {("L1", "regular", "disc_span"): 2,
                        ("SVA", "orphan", "reference_copy"): 1}


def test_both_axes_partition_the_same_verified_entries(tmp_path):
    """A cheap invariant that catches a whole class of stratification bug: two partitions of one
    set have the same total. It is the check that confirmed the real frozen baseline's 50."""
    _p, shape_v, source_v, _n = vt.read_baseline_verdicts(_frozen(tmp_path))
    assert sum(shape_v.values()) == sum(source_v.values())


def test_the_operating_point_round_trips(tmp_path):
    """A header this program WRITES has to be one it can READ. If it cannot, `--baseline` either
    refuses every valid baseline or -- worse -- compares across thresholds, where the deltas are
    artifacts of N rather than facts about the caller."""
    mid = (3, 0.8, 0.65, 5, 5.0)
    written = vt.format_op_point(tuple(str(x) for x in mid))
    point, _s, _so, _n = vt.read_baseline_verdicts(_frozen(tmp_path, point=written))
    assert point == written


def test_a_baseline_with_no_stated_operating_point_reads_as_None(tmp_path):
    """Which is what lets `main` REFUSE it (exit 3) instead of comparing counts that are not
    comparable. Silence about the thresholds is not permission to ignore them."""
    p = tmp_path / "v.tsv"
    p.write_text("# some other tool's file\n{0}\n{1}\n".format(FROZEN_HEADER, FROZEN_ROWS[0]))
    point, _s, _so, n = vt.read_baseline_verdicts(str(p))
    assert point is None
    assert n == 1


# --- 3. the pass/fail rule --------------------------------------------------------------

def _strata(counts):
    """{stratum: n_verified} -> the {stratum: Counter} shape `tally_strata` returns."""
    out = {}
    for k, n in counts.items():
        c = vt.Counter()
        c[vt.VERIFIED] = n
        out[k] = c
    return out


A = ("L1", "regular", "bare", "two_sided")
B = ("SVA", "regular", "tailed", "two_sided")


def test_no_regression_passes_when_nothing_is_lost():
    assert vt.compare_no_regression("t", {A: 3, B: 1}, _strata({A: 3, B: 1})) is True


def test_growth_passes_and_is_not_required():
    assert vt.compare_no_regression("t", {A: 3}, _strata({A: 5})) is True
    assert vt.compare_no_regression("t", {}, _strata({B: 2})) is True


def test_losing_one_entry_in_one_stratum_fails():
    """sec 10.4 gates on ABSOLUTE counts, not a tolerance, because the strata that matter hold
    1-4 entries -- `SVA regular reference_copy` holds ONE, where a percentage is noise and one
    reclassification swings it by 100 points."""
    assert vt.compare_no_regression("t", {A: 3, B: 1}, _strata({A: 3, B: 0})) is False


def test_a_gain_elsewhere_cannot_pay_for_a_loss():
    """THE FAILURE MODE A TOTAL WOULD ALLOW. Totals are equal at 4 and the gate must still FAIL:
    the two strata are different hypotheses about different biology, so entries in one do not
    substitute for entries in the other."""
    base, cur = {A: 3, B: 1}, {A: 1, B: 3}
    assert sum(base.values()) == sum(cur.values())
    assert vt.compare_no_regression("t", base, _strata(cur)) is False


def test_a_stratum_absent_from_the_baseline_is_not_a_loss():
    """A new stratum has baseline 0, and 0 -> n is growth. This is what makes the gate usable on
    an arm that finds something no baseline arm did -- which is the decoy fix's whole shape."""
    assert vt.compare_no_regression("t", {A: 3}, _strata({A: 3, B: 2})) is True


# --- 6. the self-containing source window (decision 6) --------------------------------------
#
# `PoolRow.self_window` disqualifies an entry whose source window reaches over the insertion it
# is supposed to explain: reads are gathered AT the insertion and counted for aligning to the
# window, so a window covering that point verifies itself off ordinary local coverage.
#
# MUTATIONS RUN, all of which bite:
#   | mutation | caught by |
#   |---|---|
#   | the same-chromosome test dropped | test_a_window_on_another_contig_cannot_self_contain |
#   | `ws <= pos <= we` weakened to `ws < pos < we` | the boundary test |
#   | only the D window checked | test_the_upstream_mirror_is_checked_too |
#   | the guard placed after the repeat-window test | test_self_window_outranks_a_repeat_window |
#   | the guard ignoring --allow-self-window | test_the_escape_hatch_restores_the_old_verdict |

def _srow(ins_pos, src_start, src_end, ins_chrom="chr7", src_chrom="chr7"):
    return vt.PoolRow(0, {
        "family": "SVA", "ins_chrom": ins_chrom, "ins_pos": str(ins_pos),
        "call_class": "regular", "source_class": "polymorphic",
        "src_chrom": src_chrom, "src_start": str(src_start), "src_end": str(src_end),
        "arms": "master", "td_src_raw": "{0}:{1}-{2}".format(src_chrom, src_start, src_end),
    })


def test_the_downstream_window_containing_the_insertion_is_detected():
    """The `chr17:45585847` geometry: src_end below the insertion, window reaching over it."""
    r = _srow(ins_pos=45585847, src_start=45583245, src_end=45583245)
    assert r.self_window(3000, False) == ("D", 2602)


def test_the_upstream_mirror_is_checked_too():
    """The `chr6:31329020` geometry -- caught by the U window, which `--one-sided` would skip."""
    r = _srow(ins_pos=31329020, src_start=31330396, src_end=31330396)
    assert r.self_window(3000, False) == ("U", 1624)
    # ...and under --one-sided only the downstream window exists, so it is not reachable.
    assert r.self_window(3000, True) is None


def test_a_distant_source_on_the_same_chromosome_is_not_flagged():
    """Same contig is NOT sufficient -- this is geometry, not proximity."""
    r = _srow(ins_pos=1000000, src_start=1200000, src_end=1200000)
    assert r.self_window(3000, False) is None


def test_the_1000bp_caller_constant_would_have_missed_both_real_cases():
    """The finding, pinned: mirroring MIN_POLYMORPHIC_SOURCE_DIST would catch neither.

    If someone later "simplifies" the guard into a distance threshold, this says what breaks.
    """
    for ins, src in ((45585847, 45583245), (31329020, 31330396)):
        assert abs(ins - src) > 1000, "both real cases are past the caller's guard"
        assert _srow(ins, src, src).self_window(3000, False) is not None


def test_a_window_on_another_contig_cannot_self_contain():
    """Ordinary inter-chromosomal transduction, which is what this gate is for."""
    r = _srow(ins_pos=45585847, src_start=45583245, src_end=45583245,
              ins_chrom="chr17", src_chrom="chr5")
    assert r.self_window(3000, False) is None


def test_the_window_boundary_is_inclusive():
    """A point exactly ON the edge is inside the fetched window, so it must count."""
    assert _srow(ins_pos=5000, src_start=5000, src_end=5000).self_window(3000, False) is not None
    assert _srow(ins_pos=8000, src_start=5000, src_end=5000).self_window(3000, False) is not None
    # one past the downstream edge, and the upstream mirror does not reach it either
    assert _srow(ins_pos=8001, src_start=5000, src_end=5000).self_window(3000, False) is None


def test_the_window_width_changes_what_is_flagged():
    """No new constant: the guard is a function of the window decision 1 already builds."""
    r = _srow(ins_pos=7000, src_start=5000, src_end=5000)
    assert r.self_window(3000, False) is not None    # window reaches 8000
    assert r.self_window(1000, False) is None        # window reaches 6000


def _result(self_window=None, rep_div=None, n=143):
    res = vt.RowResult(_srow(45585847, 45583245, 45583245))
    res.depth = 100
    res.self_window = self_window
    res.rep_div = rep_div
    res.rep_name = "SVA_D"
    res.support = [("clip", "D", 0.0, 1.0)] * n
    return res


POINT = (3, 0.8, 0.65, 5, 5.0)


def test_a_self_containing_window_is_undecidable_however_much_support():
    state, why = _result(self_window=("D", 2602)).verdict(*POINT)
    assert state == vt.UNDECIDABLE
    assert why == "self_window(D@+2602)"


def test_the_escape_hatch_restores_the_old_verdict():
    """`--allow-self-window` leaves `self_window` unset, so the entry verifies as it used to.

    This is what reproduces a baseline frozen before the guard, and the 143 is the real value.
    """
    state, why = _result(self_window=None).verdict(*POINT)
    assert state == vt.VERIFIED
    assert why == "support=143"


def test_self_window_outranks_a_repeat_window():
    """Circularity disqualifies even in a pristine unique region, so it is checked first."""
    state, why = _result(self_window=("D", 2602), rep_div=1.0).verdict(*POINT)
    assert state == vt.UNDECIDABLE
    assert why.startswith("self_window")


def test_the_guard_does_not_disturb_ordinary_entries():
    """The 113 of 115 the guard must leave alone."""
    assert _result(self_window=None, n=3).verdict(*POINT)[0] == vt.VERIFIED
    assert _result(self_window=None, n=0).verdict(*POINT)[0] == vt.REFUTED
    assert _result(self_window=None, n=1).verdict(*POINT)[0] == vt.SUB_THRESHOLD
