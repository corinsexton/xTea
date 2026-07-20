"""Unit tests for xtea/locus_table.py (Phase C move A.1).

Pin the type-resolution rule (load order L1 -> SVA -> Alu, +-window bp cluster,
two_side_tprt_both overrides load-order priority) in isolation from master's
SampleVCFMerger, whose real implementation has a position-keyed-dict bug (B5,
docs/clipdisc_contract_spec.md sec 4) this module deliberately does not
reproduce: a losing conflict at a different position than the winner would be
written anyway (duplicate), and one at the exact same position silently
clobbers the winner regardless of priority/TPRT. These tests assert the
*correct* (documented, intended) outcome in both situations.
"""
import os
import sys

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "xtea"))

from locus_table import LocusTable, write_merged_vcf  # noqa: E402


def _vcf(tmp_path, name, rows):
    """rows: list of (chrm, pos, alt, tprt) -> writes a minimal VCF, returns its path."""
    sf = tmp_path / name
    with open(sf, "w") as fout:
        fout.write("##fileformat=VCFv4.2\n")
        fout.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n")
        for chrm, pos, alt, tprt in rows:
            info = "SVTYPE=INS" + (";two_side_tprt_both" if tprt else "")
            fout.write(f"{chrm}\t{pos}\t.\tN\t<INS:ME:{alt}>\t.\tPASS\t{info}\n")
    return str(sf)


def test_no_conflict_keeps_everything(tmp_path):
    l1 = _vcf(tmp_path, "l1.vcf", [("chr1", 100, "LINE1", False)])
    alu = _vcf(tmp_path, "alu.vcf", [("chr1", 5000, "ALU", False)])
    loci = LocusTable().assign_type({"L1": l1, "Alu": alu})
    assert {(l.chrm, l.pos, l.type) for l in loci} == {("chr1", 100, "L1"), ("chr1", 5000, "Alu")}


def test_load_order_priority_wins_without_tprt(tmp_path):
    # L1 loads before Alu; neither is TPRT-both -> L1 (already placed) wins, Alu is dropped.
    l1 = _vcf(tmp_path, "l1.vcf", [("chr1", 1000, "LINE1", False)])
    alu = _vcf(tmp_path, "alu.vcf", [("chr1", 1010, "ALU", False)])  # 10bp away, within window
    loci = LocusTable().assign_type({"L1": l1, "Alu": alu})
    assert len(loci) == 1
    assert (loci[0].chrm, loci[0].pos, loci[0].type) == ("chr1", 1000, "L1")


def test_tprt_both_overrides_load_order(tmp_path):
    # Alu loads after L1 but carries two_side_tprt_both while L1 doesn't -> Alu wins.
    l1 = _vcf(tmp_path, "l1.vcf", [("chr1", 1000, "LINE1", False)])
    alu = _vcf(tmp_path, "alu.vcf", [("chr1", 1010, "ALU", True)])
    loci = LocusTable().assign_type({"L1": l1, "Alu": alu})
    assert len(loci) == 1
    assert (loci[0].chrm, loci[0].pos, loci[0].type) == ("chr1", 1010, "Alu")


def test_exact_position_tie_keeps_priority_winner(tmp_path):
    # Regression for B5: master's real merger silently lets the LAST-loaded type
    # (Alu) win an exact-position tie regardless of priority/TPRT. The fixed rule
    # must keep the higher-priority (or TPRT-both) call instead.
    sva = _vcf(tmp_path, "sva.vcf", [("chr5", 142000088, "SVA", True)])
    alu = _vcf(tmp_path, "alu.vcf", [("chr5", 142000088, "ALU", True)])
    loci = LocusTable().assign_type({"SVA": sva, "Alu": alu})
    assert len(loci) == 1
    assert loci[0].type == "SVA"


def test_out_of_window_is_not_a_conflict(tmp_path):
    l1 = _vcf(tmp_path, "l1.vcf", [("chr1", 1000, "LINE1", False)])
    alu = _vcf(tmp_path, "alu.vcf", [("chr1", 1051, "ALU", False)])  # just outside 50bp window
    loci = LocusTable(window=50).assign_type({"L1": l1, "Alu": alu})
    assert len(loci) == 2


def test_write_merged_vcf_header_and_sort(tmp_path):
    l1 = _vcf(tmp_path, "l1.vcf", [("chr2", 500, "LINE1", False)])
    alu = _vcf(tmp_path, "alu.vcf", [("chr1", 100, "ALU", False)])
    out = str(tmp_path / "merged.vcf")
    loci = write_merged_vcf({"L1": l1, "Alu": alu}, out)
    assert len(loci) == 2
    with open(out) as fin:
        lines = [ln.rstrip("\n") for ln in fin]
    assert lines[0] == "##fileformat=VCFv4.2"
    assert lines[1] == "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO"
    data_lines = [ln for ln in lines if ln and not ln.startswith("#")]
    assert len(data_lines) == 2
    assert data_lines[0].startswith("chr1\t100\t")  # chr1 before chr2 -- sort -k1,1V -k2,2n
    assert data_lines[1].startswith("chr2\t500\t")
