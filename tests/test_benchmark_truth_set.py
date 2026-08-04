"""The truth-set benchmark's parsers and matcher.

WHY THESE AND NOT OTHERS
------------------------
`benchmark/benchmark_truth_set.py` is the Stage 3 gate. A gate that is quietly wrong is
worse than no gate, and the three places it can be quietly wrong are all pure functions:

1. **`TD_SRC` parsing.** The field has three shapes and one booby trap -- the region
   delimiter is `-` and an unknown start is literally `-1`, so `chr10:-1-108181069` splits
   naively into ("", "1", "108181069"). Get it wrong and ~40 of 1493 real calls silently
   acquire an absurd source locus, which shows up as a plausible-looking
   TP_POS_WRONG_SOURCE rate rather than as an error.

2. **The 3' residual.** The whole truth-side transduction derivation rests on it, and the
   obvious implementation -- "take the primary's trailing soft clip" -- is wrong twice
   over. The ALT is written reference-forward, so a minus-strand insertion's 3' end is its
   LEFT end (half the records); and a supplementary alignment to the consensus is explained
   sequence, not residual (38 of 189 aligned L1 truth ALTs have one, usually a 5'
   inversion). The residual is the ALT minus the union of every alignment's query span,
   taken at the end the PRIMARY's strand designates.

3. **One-to-one matching.** At a 200 bp window over clustered MEI loci, "is any call near
   this truth record" lets one call satisfy several truth records. That inflates
   sensitivity and breaks TP + FN = |truth|.

The bwa-driven parts of `derive-td` are not unit-tested: they are an alignment pipeline,
and docs/handoff/runs/derive_truth_td.sbatch is their harness.
"""
import os
import sys

import pytest

REPO_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(REPO_ROOT, "benchmark"))
import benchmark_truth_set as bts  # noqa: E402


# --- 1. TD_SRC ---------------------------------------------------------------------------

def test_td_src_the_three_real_shapes():
    """Every shape measured on real output (BASELINE-0: 1403 / 40 / 50)."""
    assert bts.parse_td_src("not_transduction") is None
    assert bts.parse_td_src("") is None
    assert bts.parse_td_src(None) is None
    assert bts.parse_td_src("chr8:51817625-51818695") == ("region", "chr8", 51817625, 51818695)
    assert bts.parse_td_src("chrY:9591830~sibling") == ("sibling", "chrY", 9591830)


def test_td_src_unknown_start_is_minus_one_not_a_delimiter():
    """THE trap. `chrX:-1-45277389` is a region with an unknown start, not a malformed one."""
    assert bts.parse_td_src("chrX:-1-45277389") == ("region", "chrX", -1, 45277389)
    assert bts.parse_td_src("chr10:-1-108181069") == ("region", "chr10", -1, 108181069)
    # ...and the comparable point falls back to the END, the only coordinate that is known.
    assert bts.td_src_point(bts.parse_td_src("chrX:-1-45277389")) == ("chrX", 45277389)
    assert bts.td_src_point(bts.parse_td_src("chr8:51817625-51818695")) == ("chr8", 51817625)
    assert bts.td_src_point(bts.parse_td_src("chrY:9591830~sibling")) == ("chrY", 9591830)
    assert bts.td_src_point(None) is None


def test_td_src_negative_end_is_refused_not_guessed():
    """Splitting on the LAST `-` is only correct while the end is non-negative. If that
    ever changes the parser must fail, not return a nonsense region."""
    with pytest.raises(ValueError):
        bts.parse_td_src("chr1:100--1")


# --- 2. the 3' residual -------------------------------------------------------------------

def _rec(flag, cigar, rname="cns", pos=1, mapq=60):
    return (flag, rname, pos, mapq, cigar)


ALT = "AAAACCCCGGGGTTTT"      # 16 bp, all four ends distinguishable


def test_three_prime_residual_forward_is_the_right_end_of_the_alt():
    assert bts.three_prime_residual([_rec(0, "12M4S")], ALT) == "TTTT"
    assert bts.three_prime_residual([_rec(0, "3S10M3S")], ALT) == "TTT"
    # a leading clip alone is the element's 5' end -- not a 3' residual
    assert bts.three_prime_residual([_rec(0, "4S12M")], ALT) == ""
    assert bts.three_prime_residual([], ALT) == ""


def test_three_prime_residual_reverse_is_the_LEFT_end_revcomped():
    """The ALT is written reference-forward, so a minus-strand insertion's 3' end is the
    ALT's left end -- and the residual must come back in the element's own orientation.

    Taking the primary's trailing soft clip instead (the obvious reading, since bwa
    revcomps SEQ) silently returns the element's 5' end on every minus-strand insertion,
    which is about half of them.
    """
    # record covers ALT[4:16]; the unexplained part is ALT[0:4] = "AAAA" -> revcomp "TTTT"
    assert bts.three_prime_residual([_rec(16, "12M4S")], ALT) == "TTTT"
    # record covers ALT[0:12]; the element's 3' end (left) is fully explained
    assert bts.three_prime_residual([_rec(16, "4S12M")], ALT) == ""


def test_supplementary_coverage_counts_as_explained():
    """A 5' inversion makes bwa split the ALT instead of clipping it. The second hit is
    still consensus sequence, so it is not residual -- 38 of 189 L1 truth ALTs have one.
    """
    # primary aligns ALT[0:8], supplementary aligns ALT[8:16] (hard-clipped, as bwa emits)
    recs = [_rec(0, "8M8S"), _rec(0x800, "8H8M")]
    assert bts.three_prime_residual(recs, ALT) == ""
    # ...and with only the primary, the same input reads as an 8 bp residual
    assert bts.three_prime_residual([recs[0]], ALT) == "GGGGTTTT"


def test_query_span_in_alt():
    assert bts.query_span_in_alt(0, "4S12M", 16) == (4, 16)
    assert bts.query_span_in_alt(0, "12M4S", 16) == (0, 12)
    assert bts.query_span_in_alt(16, "12M4S", 16) == (4, 16)     # mirrored
    assert bts.query_span_in_alt(0x800, "8H8M", 16) == (8, 16)   # H counts like S
    # insertions consume query, deletions do not
    assert bts.query_span_in_alt(0, "6M2I8M", 16) == (0, 16)
    assert bts.query_span_in_alt(0, "8M2D8M", 16) == (0, 16)
    # a cigar that does not account for the whole query belongs to another sequence
    assert bts.query_span_in_alt(0, "8M", 16) is None


def test_trim_polya_strips_the_tail_that_is_not_source_sequence():
    assert bts.trim_polya("ACGTACGTAAAAAAAAAA") == "ACGTACGT"
    # Note the last real base goes with the run: a genuine sequence ending in T is
    # indistinguishable from the head of a polyT tail. Inherent, and harmless -- the
    # residual is being length-thresholded at 100 bp and mapped, not compared exactly.
    assert bts.trim_polya("ACGTACGTTTTTTTTTTT") == "ACGTACG"
    assert bts.trim_polya("ACGTACGT") == "ACGTACG"
    assert bts.trim_polya("ACGTACGA") == "ACGTACG"
    assert bts.trim_polya("ACGTACGC") == "ACGTACGC"
    assert bts.trim_polya("AAAAAAAA") == ""
    # only ONE of the two runs is trimmed: an internal A-run followed by T's is a real
    # sequence, and eating both would shorten genuine transduced sequence.
    assert bts.trim_polya("ACGTAAAATTTT") == "ACGTAAAA"


def test_parse_cigar():
    assert bts.parse_cigar("12M4S") == [(12, "M"), (4, "S")]
    assert bts.parse_cigar("3S10M2I5D3S") == [(3, "S"), (10, "M"), (2, "I"), (5, "D"), (3, "S")]


# --- 3. matching --------------------------------------------------------------------------

class _T:
    def __init__(self, chrom, pos):
        self.chrom, self.pos = chrom, pos


class _C:
    def __init__(self, chrom, pos):
        self.chrom, self.pos = chrom, pos


def test_matching_is_one_to_one():
    """Two truth records 100 bp apart and ONE call between them: exactly one match, the
    nearer. x_cmp.py's per-truth "any call within the window" test would score both."""
    truth = [_T("chr1", 1000), _T("chr1", 1100)]
    calls = [_C("chr1", 1040)]
    m = bts.match_one_to_one(truth, calls, 200)
    assert len(m) == 1 and 0 in m and m[0][0] == 0


def test_matching_prefers_the_nearest_pair_globally():
    truth = [_T("chr1", 1000), _T("chr1", 1100)]
    calls = [_C("chr1", 1099), _C("chr1", 1005)]
    m = bts.match_one_to_one(truth, calls, 200)
    assert len(m) == 2
    assert m[0][0] == 1 and m[1][0] == 0        # each took its nearest


def test_matching_respects_the_window_and_the_contig():
    assert bts.match_one_to_one([_T("chr1", 1000)], [_C("chr1", 1201)], 200) == {}
    assert len(bts.match_one_to_one([_T("chr1", 1000)], [_C("chr1", 1200)], 200)) == 1
    assert bts.match_one_to_one([_T("chr1", 1000)], [_C("chr2", 1000)], 200) == {}


def test_matching_makes_the_two_totals_add_up():
    """TP + FN = |truth| and TP + FP = |calls| -- the property that makes sensitivity and
    precision comparable, and the one many-to-one matching breaks."""
    truth = [_T("chr1", p) for p in (1000, 1100, 5000, 9000)]
    calls = [_C("chr1", p) for p in (1010, 1090, 1120, 20000)]
    m = bts.match_one_to_one(truth, calls, 200)
    tp = len(m)
    assert tp + (len(truth) - tp) == len(truth)
    assert tp + (len(calls) - tp) == len(calls)
    assert tp == 2      # 5000 and 9000 unmatched; only two of the three clustered calls used


# --- the truth VCF itself -------------------------------------------------------------------

def test_headerless_truth_vcf_and_the_inserted_sequence(tmp_path):
    """The truth VCFs have no `##` and no `#CHROM`. And the insertion is ALT minus the
    anchoring REF base -- stripped only when it is actually the same base."""
    p = tmp_path / "t.vcf"
    p.write_text(
        "chr12\t70660189\t.\tT\tTAAACCC\t60\t.\tQNAME=NA24385#S10-H1-000000F;QSTART=6945822\tGT\t1/1\tH1\n"
        "chr1\t100\t.\tG\tGACGT\t60\t.\tQNAME=x\tGT\t1/1\tH1,H2\n")
    recs = bts.read_truth_vcf(str(p))
    assert [(r.chrom, r.pos) for r in recs] == [("chr12", 70660189), ("chr1", 100)]
    assert recs[0].inserted_sequence() == "AAACCC"
    assert recs[1].inserted_sequence() == "ACGT"
    assert bts.parse_info(recs[0].info)["QSTART"] == "6945822"
    # A non-conventional record keeps its whole ALT rather than losing a real base.
    recs[0].ref, recs[0].alt = "T", "GAAACCC"
    assert recs[0].inserted_sequence() == "GAAACCC"


def test_tprt_must_stay_a_subset_of_all(tmp_path):
    """The `_tprt` flag is applied to `_all` records rather than loaded as its own set,
    which is only sound while `_tprt` nests. If it stops nesting, fail loudly."""
    d = tmp_path
    (d / "L1_all.vcf").write_text("chr1\t100\t.\tG\tGACGT\t60\t.\tQNAME=x\n")
    (d / "L1_tprt.vcf").write_text("chr1\t100\t.\tG\tGACGT\t60\t.\tQNAME=x\n")
    recs = bts.load_truth_family(str(d), "L1")
    assert len(recs) == 1 and recs[0].tprt is True

    (d / "L1_tprt.vcf").write_text(
        "chr1\t100\t.\tG\tGACGT\t60\t.\tQNAME=x\n"
        "chr2\t999\t.\tG\tGACGT\t60\t.\tQNAME=y\n")
    with pytest.raises(ValueError, match="absent from _all"):
        bts.load_truth_family(str(d), "L1")


def test_sva_ships_one_combined_file():
    """SVA has no separate `_all`; SVA_tprt_all.vcf is both the denominator and the
    stratum, so every SVA record must come back flagged tprt."""
    assert bts.TRUTH_FILES["SVA"]["all"] == bts.TRUTH_FILES["SVA"]["tprt"] == "SVA_tprt_all.vcf"


def test_full_length_copy_adjacency(tmp_path):
    """RepeatMasker .out parsing (fields 4/5/6, two header lines) plus the adjacency test
    on both sides -- the reference copy may be on either strand."""
    p = tmp_path / "r.out"
    p.write_text(
        "   SW   perc perc perc  query\n"
        "score   div. del. ins.  sequence\n"
        "\n"
        "  100  1.0  0.0  0.0  chr1  1000  7000  x  +  L1HS  LINE/L1  1 6000 0  1\n"
        "  100  1.0  0.0  0.0  chr1  20000 20500 x  +  L1PA2 LINE/L1  1  500 0  2\n")
    copies = bts.load_full_length_copies(str(p), 5000)
    assert copies == {"chr1": [(1000, 7000)]}      # the 500 bp copy is not full length
    assert bts.adjacent_full_length_copy(copies, "chr1", 7500, 10000) == (1000, 7000)
    assert bts.adjacent_full_length_copy(copies, "chr1", 500, 10000) == (1000, 7000)
    assert bts.adjacent_full_length_copy(copies, "chr1", 30000, 10000) is None
    assert bts.adjacent_full_length_copy(copies, "chr2", 7500, 10000) is None


def test_contains_is_inside_not_merely_near(tmp_path):
    """The divergent-body test. A 3' residual landing INSIDE a reference copy of its own
    family is far more likely to be element body the consensus does not cover than
    transduced sequence -- measured: 11 of 20 L1 residuals >=100 bp do, and every one is
    an older subfamily (L1PA6/7/8 at 4-7% divergence, L1MB/MD/ME at 16-25%). It must not
    be confused with the adjacency test, which is deliberately loose on both sides.
    """
    copies = {"chr1": [(1000, 2000), (5000, 6000)]}
    assert bts.contains(copies, "chr1", 1500) is True
    assert bts.contains(copies, "chr1", 1000) is True
    assert bts.contains(copies, "chr1", 2000) is True
    assert bts.contains(copies, "chr1", 2001) is False
    assert bts.contains(copies, "chr1", 999) is False
    assert bts.contains(copies, "chr1", 5500) is True     # past the early-exit boundary
    assert bts.contains(copies, "chr2", 1500) is False
    # ...and 2001 is NOT inside, while the adjacency test does accept it
    assert bts.adjacent_full_length_copy(copies, "chr1", 2001, 100) == (1000, 2000)


# --- 4. the FN attribution's phase order -------------------------------------------------

def test_the_earliest_gate_wins_and_phase_order_is_not_lexical():
    """`ledger_attribution.attribute` picks the gate that determined the outcome, and "earliest"
    is the CASCADE's order, not the phase strings' alphabetical order.

    Gate B2 (`check_seprt_disc_consistency`) is the reason this needs a test: its phase name
    `B2_cns_cluster` sorts BEFORE `B_source` in ASCII, because '2' < '_'. The tool used to take
    `min(rows)` on the phase string, which was correct only as long as A < B < C < D < E held by
    accident. With B2 present that would attribute a false negative to B2 when gate B rejected
    it first -- a wrong attribution that looks entirely plausible in the output table.
    """
    sys.path.insert(0, os.path.join(REPO_ROOT, "benchmark"))
    import ledger_attribution as la

    b_row = ("chr1", 100, "B_source", "low_disc_support", "n=1")
    b2_row = ("chr1", 101, "B2_cns_cluster", "no_cns_disc_cluster", "n_l=0 n_r=0")
    a_row = ("chr1", 102, "A_clip", "low_clip_support", "n_clip=1 < 4")
    assert la.attribute([b2_row, b_row])[0] == "low_disc_support"
    assert la.attribute([b_row, b2_row])[0] == "low_disc_support"
    assert la.attribute([b2_row, a_row])[0] == "low_clip_support"
    assert la.attribute([b2_row])[0] == "no_cns_disc_cluster"
    assert la.attribute([])[0] == la.NO_ROW
    # The hazard itself, so this test's premise cannot go stale silently.
    assert min(r[2] for r in (b_row, b2_row)) == "B2_cns_cluster"
