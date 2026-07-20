"""Unit tests for the concordance comparator (synthetic call-sets, no external inputs).

These pin the matching semantics of the reproducibility gate: ±window, type-aware, one-to-one.
They run everywhere (pure Python) and must stay green.
"""
from tests.utils.vcf_io import MEICall, ALU, LINE1, SVA
from tests.utils.concordance import compare_callsets, dedupe_conflict_loci


def _c(chrom, pos, t):
    return MEICall(chrom=chrom, pos=pos, mei_type=t)


def test_identical_callsets_are_fully_concordant():
    calls = [_c("chr1", 100, ALU), _c("chr1", 5000, SVA), _c("chr2", 200, LINE1)]
    r = compare_callsets(calls, list(calls))
    assert r.sensitivity == 1.0 and r.ppv == 1.0
    assert not r.missed and not r.extra


def test_within_window_matches():
    truth = [_c("chr1", 1000, ALU)]
    query = [_c("chr1", 1000 + 25, ALU)]  # exactly at the window edge
    r = compare_callsets(truth, query, window=25)
    assert r.n_matched == 1 and r.sensitivity == 1.0


def test_outside_window_misses():
    truth = [_c("chr1", 1000, ALU)]
    query = [_c("chr1", 1000 + 26, ALU)]
    r = compare_callsets(truth, query, window=25)
    assert r.n_matched == 0
    assert r.missed and r.extra
    assert r.sensitivity == 0.0 and r.ppv == 0.0


def test_type_mismatch_does_not_match():
    truth = [_c("chr1", 1000, ALU)]
    query = [_c("chr1", 1000, SVA)]
    r = compare_callsets(truth, query)
    assert r.n_matched == 0


def test_type_mismatch_ignored_when_match_type_false():
    truth = [_c("chr1", 1000, ALU)]
    query = [_c("chr1", 1000, SVA)]
    r = compare_callsets(truth, query, match_type=False)
    assert r.n_matched == 1


def test_one_to_one_matching_no_double_count():
    # One query call must not satisfy two nearby truth calls.
    truth = [_c("chr1", 1000, ALU), _c("chr1", 1010, ALU)]
    query = [_c("chr1", 1005, ALU)]
    r = compare_callsets(truth, query, window=25)
    assert r.n_matched == 1
    assert len(r.missed) == 1  # the second truth call is unmatched


def test_nearest_wins():
    truth = [_c("chr1", 1000, ALU)]
    query = [_c("chr1", 1020, ALU), _c("chr1", 1002, ALU)]
    r = compare_callsets(truth, query, window=25)
    assert r.n_matched == 1
    # the far query (1020) is the leftover extra; nearest (1002) matched
    assert len(r.extra) == 1 and r.extra[0].pos == 1020


def test_extra_query_lowers_ppv():
    truth = [_c("chr1", 1000, ALU)]
    query = [_c("chr1", 1000, ALU), _c("chr1", 9000, ALU)]
    r = compare_callsets(truth, query)
    assert r.sensitivity == 1.0
    assert r.ppv == 0.5


def test_empty_truth_and_query_edge_cases():
    assert compare_callsets([], []).sensitivity == 1.0
    assert compare_callsets([], []).ppv == 1.0
    # L1==0 golden: empty truth, empty query → trivially concordant (documents the coverage gap)
    assert compare_callsets([], [_c("chr1", 1, ALU)]).ppv == 0.0


def _c_subtype(chrom, pos, t, subtype=None):
    info = {"SUBTYPE": subtype} if subtype else {}
    return MEICall(chrom=chrom, pos=pos, mei_type=t, info=info)


def test_dedupe_keeps_non_conflicting_calls_untouched():
    calls = [_c("chr1", 100, ALU), _c("chr2", 200, SVA), _c("chr3", 300, LINE1)]
    assert dedupe_conflict_loci(calls) == calls


def test_dedupe_collapses_same_locus_reported_under_two_types_by_priority():
    # L1 beats Alu at the same locus (load-order priority) when neither carries TPRT-both.
    l1_call = _c("chr1", 1000, LINE1)
    alu_call = _c("chr1", 1010, ALU)  # within the +-50bp merge window
    survivors = dedupe_conflict_loci([l1_call, alu_call])
    assert survivors == [l1_call]


def test_dedupe_tprt_both_overrides_load_order_priority():
    # Alu loads after L1 in priority order, but carries two_side_tprt_both -> Alu wins anyway.
    l1_call = _c("chr1", 1000, LINE1)
    alu_tprt_call = _c_subtype("chr1", 1010, ALU, subtype="two_side_tprt_both")
    survivors = dedupe_conflict_loci([l1_call, alu_tprt_call])
    assert survivors == [alu_tprt_call]


def test_dedupe_beyond_window_are_independent_loci():
    l1_call = _c("chr1", 1000, LINE1)
    alu_call_far = _c("chr1", 1000 + 51, ALU)  # just outside the +-50bp merge window
    survivors = dedupe_conflict_loci([l1_call, alu_call_far])
    assert survivors == [l1_call, alu_call_far]


def test_dedupe_fixes_the_real_wgs_chr1_241980550_three_way_pileup():
    # docs/handoff/MOVE_A_PLAN.md's documented 3-way conflict: L1 + SVA + Alu all fire near
    # chr1:~241980550; two_side_tprt_both on the SVA call makes SVA the correct winner.
    l1_call = _c("chr1", 241980558, LINE1)
    sva_call = _c_subtype("chr1", 241980544, SVA, subtype="two_side_tprt_both")
    alu_call = _c("chr1", 241980550, ALU)
    survivors = dedupe_conflict_loci([l1_call, sva_call, alu_call])
    assert survivors == [sva_call]


def test_xtea2_vs_golden_style_comparison_no_longer_penalizes_known_conflict_loci():
    # Reproduces the exact false-negative pattern seen in the WGS move-A.2 gate run: master's
    # undeduplicated golden lists one physical insertion twice (SVA + LINE1, 14bp apart); xtea2's
    # output correctly has it once, as SVA. Without dedup this scores as a missed truth call even
    # though xtea2 called the locus.
    golden_sva = _c_subtype("chr1", 241980544, SVA, subtype="two_side_tprt_both")
    golden_line1 = _c("chr1", 241980558, LINE1)
    xtea2_query = [_c("chr1", 241980544, SVA)]

    r_undeduped = compare_callsets([golden_sva, golden_line1], xtea2_query, window=25)
    assert r_undeduped.sensitivity < 1.0  # the bug this test documents

    truth = dedupe_conflict_loci([golden_sva, golden_line1])
    r_deduped = compare_callsets(truth, xtea2_query, window=25)
    assert r_deduped.sensitivity == 1.0
    assert r_deduped.ppv == 1.0
