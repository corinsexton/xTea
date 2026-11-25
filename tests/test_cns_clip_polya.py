"""Poly-A counting in parse_clip_realignment_consensus (the cns step's clip realignment).

The counter used to accept only an A run in the realigned clip. The SAM stores the clip in
the consensus record's orientation, and a family's .fa has many records (SVA has 88), some
with internal T runs. So a poly-A clip that ties onto one of those reverse-strand reads
TTTTT and was dropped: HG002 chr5:56152167 SVA lost its only poly-A read that way and
failed post-filtering's poly-A requirement.

The obvious "if reverse, look for T" fix is also wrong: a poly-T clip that reverse-aligns
to a consensus poly-A tail is stored as AAAAA. These four cases pin both directions, plus
the A/T tally ([2]/[3]) that turns into the sign of nlpolyA/nrpolyA downstream.

Z2 (Zofia, SMaHT P25): accepting ANY T run also counted ordinary L1/SVA body sequence as
poly-A. A T-only slice now counts only when it is a homopolymer (>75% T), which the chr5
rescue clip is and internal sequence is not.
"""
import os
import sys

import pysam
import pytest

REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(REPO, "xtea"))
import global_values  # noqa: E402
import x_clip_disc_filter as xcdf  # noqa: E402

INS_POS = 1000
REVERSE = 16


def _name(i, side):
    # chrm~map_pos~L/R~is_rc~insertion_pos~rid~sample_id
    return global_values.SEPERATOR.join(
        ["chr5", str(INS_POS), side, "0", str(INS_POS), str(i), "0"])


def _count(tmp_path, records):
    """records: (seq, flag, side). Returns m_polyA['chr5'][INS_POS] or None."""
    sam = str(tmp_path / "temp_clip.sam")
    header = {"HD": {"VN": "1.6"}, "SQ": [{"SN": "SVA_s75", "LN": 2000}]}
    with pysam.AlignmentFile(sam, "w", header=header) as out:
        for i, (seq, flag, side) in enumerate(records):
            a = pysam.AlignedSegment(out.header)
            a.query_name = _name(i, side)
            a.flag = flag
            a.reference_id = 0
            a.reference_start = 500
            a.mapping_quality = 0
            a.cigartuples = [(0, len(seq))]
            a.query_sequence = seq
            a.query_qualities = pysam.qualitystring_to_array("I" * len(seq))
            out.write(a)
    obj = xcdf.XClipDiscFilter.__new__(xcdf.XClipDiscFilter)
    obj.sf_reference = None
    m_polyA = obj.parse_clip_realignment_consensus(sam, global_values.MIN_CLIP_MAPPED_RATIO)[5]
    return m_polyA.get("chr5", {}).get(INS_POS)


@pytest.mark.parametrize("seq, flag, genome_polyA", [
    ("GAAAAAAAAA", 0, True),         # forward, A in SAM: the case that always worked
    ("TTTTTTTTTC", REVERSE, True),   # chr5: poly-A clip reverse-tied onto an internal T run
    ("AAAAAAAAAG", REVERSE, False),  # poly-T clip reverse-aligned to the consensus A tail
    ("CTTTTTTTTT", 0, False),        # poly-T clip forward onto an internal T run
])
def test_polya_counted_on_either_strand(tmp_path, seq, flag, genome_polyA):
    got = _count(tmp_path, [(seq, flag, "L")])
    assert got is not None and got[0] == 1, "clip with an A or T run must be counted"
    assert got[2] == (1 if genome_polyA else 0), "A/T tally must follow the genome orientation"


def test_no_run_not_counted(tmp_path):
    assert _count(tmp_path, [("GACGTACGTA", 0, "L"), ("GACGTACGTA", REVERSE, "R")]) is None


def test_right_clip_tally(tmp_path):
    got = _count(tmp_path, [("TTTTTTTTTC", REVERSE, "R"), ("AAAAAAAAAG", REVERSE, "R")])
    assert got == [0, 2, 0, 1]


@pytest.mark.parametrize("seq", [
    "ATCCAGGAGCTGGTTTTTTG",   # HG002 L1 chr14:39405890 right clips: L1 ORF2 at cns 3071
    "GCAGTTTTTGTTTTGGCAAG",   # SVA body T run (Zofia chr3:177272377, chr13:53965793)
    "GGACTTATTGCCAAGACTGA",   # TTATT in ordinary L1 sequence (Zofia chr2:105546701)
])
@pytest.mark.parametrize("flag", [0, REVERSE])
def test_internal_t_run_not_counted(tmp_path, seq, flag):
    assert _count(tmp_path, [(seq, flag, "L"), (seq, flag, "R")]) is None


def test_a_run_rule_unchanged_for_mixed_slices(tmp_path):
    # master's rule: any A run counts, homopolymer or not
    got = _count(tmp_path, [("GGCAAAAAGCTCGGACGTCA", 0, "L")])
    assert got == [1, 0, 1, 0]
