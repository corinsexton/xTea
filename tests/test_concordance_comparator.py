"""Unit tests for the concordance comparator (synthetic call-sets, no external inputs).

These pin the matching semantics of the reproducibility gate: ±window, type-aware, one-to-one.
They run everywhere (pure Python) and must stay green.
"""
from tests.utils.vcf_io import MEICall, ALU, LINE1, SVA
from tests.utils.concordance import compare_callsets


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
