"""Unit tests for xtea/locus_table.py (Phase C move A.1/A.2).

Pin the type-resolution rule (load order L1 -> SVA -> Alu, +-window bp cluster,
two_side_tprt_both overrides load-order priority) in isolation from master's
SampleVCFMerger, whose real implementation has a position-keyed-dict bug (B5,
docs/clipdisc_contract_spec.md sec 4) this module deliberately does not
reproduce: a losing conflict at a different position than the winner would be
written anyway (duplicate), and one at the exact same position silently
clobbers the winner regardless of priority/TPRT. These tests assert the
*correct* (documented, intended) outcome in both situations.

Move A.2 retargeted the resolution point from the three types' final VCFs to
their raw postF_hc candidate-list files (native xTea txt format), so the rows
here are synthetic postF_hc lines, not VCF lines.
"""
import os
import sys

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "xtea"))

from locus_table import LocusTable, TYPE_TO_RTYPE, assign_and_write, load_type_map  # noqa: E402


def _hc_pf(tmp_path, name, rows):
    """rows: list of (chrm, lclip_pos, rclip_pos, tprt) -> writes a minimal
    whitespace-separated postF_hc-format file (fields[0]=chrm, [2]=lclip_pos,
    [3]=rclip_pos, TPRT flag as a literal trailing token), returns its path."""
    sf = tmp_path / name
    with open(sf, "w") as fout:
        for chrm, lclip_pos, rclip_pos, tprt in rows:
            tprt_tok = "two_side_tprt_both" if tprt else "one_end_consistent"
            fout.write(f"{chrm}\t0\t{lclip_pos}\t{rclip_pos}\tfiller\t{tprt_tok}\n")
    return str(sf)


def test_no_conflict_keeps_everything(tmp_path):
    l1 = _hc_pf(tmp_path, "l1.txt", [("chr1", 100, 110, False)])
    alu = _hc_pf(tmp_path, "alu.txt", [("chr1", 5000, 5010, False)])
    loci = LocusTable().assign_type({"L1": l1, "Alu": alu})
    assert {(l.chrm, l.pos, l.type) for l in loci} == {("chr1", 100, "L1"), ("chr1", 5000, "Alu")}


def test_load_order_priority_wins_without_tprt(tmp_path):
    # L1 loads before Alu; neither is TPRT-both -> L1 (already placed) wins, Alu is dropped.
    l1 = _hc_pf(tmp_path, "l1.txt", [("chr1", 1000, 1010, False)])
    alu = _hc_pf(tmp_path, "alu.txt", [("chr1", 1010, 1020, False)])  # 10bp away, within window
    loci = LocusTable().assign_type({"L1": l1, "Alu": alu})
    assert len(loci) == 1
    assert (loci[0].chrm, loci[0].pos, loci[0].type) == ("chr1", 1000, "L1")


def test_tprt_both_overrides_load_order(tmp_path):
    # Alu loads after L1 but carries two_side_tprt_both while L1 doesn't -> Alu wins.
    l1 = _hc_pf(tmp_path, "l1.txt", [("chr1", 1000, 1010, False)])
    alu = _hc_pf(tmp_path, "alu.txt", [("chr1", 1010, 1020, True)])
    loci = LocusTable().assign_type({"L1": l1, "Alu": alu})
    assert len(loci) == 1
    assert (loci[0].chrm, loci[0].pos, loci[0].type) == ("chr1", 1010, "Alu")


def test_exact_position_tie_keeps_priority_winner(tmp_path):
    # Regression for B5: master's real merger silently lets the LAST-loaded type
    # (Alu) win an exact-position tie regardless of priority/TPRT. The fixed rule
    # must keep the higher-priority (or TPRT-both) call instead.
    sva = _hc_pf(tmp_path, "sva.txt", [("chr5", 142000088, 142000098, True)])
    alu = _hc_pf(tmp_path, "alu.txt", [("chr5", 142000088, 142000098, True)])
    loci = LocusTable().assign_type({"SVA": sva, "Alu": alu})
    assert len(loci) == 1
    assert loci[0].type == "SVA"


def test_out_of_window_is_not_a_conflict(tmp_path):
    l1 = _hc_pf(tmp_path, "l1.txt", [("chr1", 1000, 1010, False)])
    alu = _hc_pf(tmp_path, "alu.txt", [("chr1", 1051, 1061, False)])  # just outside 50bp window
    loci = LocusTable(window=50).assign_type({"L1": l1, "Alu": alu})
    assert len(loci) == 2


def test_assign_and_write_sorts_and_writes_type_map(tmp_path):
    l1 = _hc_pf(tmp_path, "l1.txt", [("chr2", 500, 510, False)])
    alu = _hc_pf(tmp_path, "alu.txt", [("chr1", 100, 110, False)])
    sf_merged = str(tmp_path / "merged.txt")
    sf_type_map = str(tmp_path / "type_map.tsv")
    loci = assign_and_write({"L1": l1, "Alu": alu}, sf_merged, sf_type_map)
    assert len(loci) == 2

    with open(sf_merged) as fin:
        lines = [ln.rstrip("\n") for ln in fin]
    assert len(lines) == 2
    assert lines[0].startswith("chr1\t")  # chr1 before chr2 -- sort -k1,1V -k2,2n
    assert lines[1].startswith("chr2\t")

    type_map = load_type_map(sf_type_map)
    assert type_map == {("chr1", 100): TYPE_TO_RTYPE["Alu"], ("chr2", 500): TYPE_TO_RTYPE["L1"]}
