"""Unit tests for xtea/x_gvcf.py's cvt_raw_rslt_to_gvcf_multi (Phase C move A.2).

gene annotation and genotyping are already type-agnostic (see MOVE_A_PLAN.md); the
one step that needs a per-record type instead of one global --rtype is gVCF, since
it stamps type into every record's ALT/SVTYPE. These tests check: (1) the new
fan-out function is byte-identical to the existing single-type function for an
equivalent single-type input (regression safety net for the untouched function),
(2) it correctly dispatches mixed-type records to their own per-type file plus one
merged file, and (3) it fails loud rather than silently mis-typing a record if the
type map and the record set ever disagree.
"""
import os
import sys

import pysam
import pytest

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "xtea"))

from x_gvcf import gVCF  # noqa: E402

# A real postF_hc row (49 whitespace-separated fields), used as a template: only
# chrm/lclip_pos/rclip_pos/TPRT-token vary between synthetic test rows, every other
# column is a fixed placeholder so every int()/float() cast in _cvt_to_vcf_rcd
# still succeeds.
_TEMPLATE = ["{chrm}", "0", "{lclip}", "{rclip}", "10", "6", "4", "14", "10", "6", "0",
             "35.4", "35.065", "6", "4", "14", "10", "6", "0", "5961:5974", "4619:4621",
             "4643.0:4916.0", "5846.0:6019.0", "not_transduction", "0", "0", "0", "14",
             "0", "10", "0", "Not-5prime-inversion", "{tprt}", "one_end_consistent",
             "hit_end_of_consensus", "10", "4", "7", "10", "34", "50", "0",
             "6:22:23:35:49:56:73:73:74:75", "0", "14", "10", "0", "1355",
             "Left_is_LINE1_copy_21.1"]
assert len(_TEMPLATE) == 49


def _row(chrm, lclip, rclip, tprt=False):
    tok = "two_side_tprt_both" if tprt else "one_end_consistent"
    vals = {"chrm": chrm, "lclip": lclip, "rclip": rclip, "tprt": tok}
    return "\t".join(f.format(**vals) if "{" in f else f for f in _TEMPLATE)


@pytest.fixture
def ref_fasta(tmp_path):
    sf = str(tmp_path / "ref.fa")
    with open(sf, "w") as fout:
        fout.write(">chr1\n")
        seq = "ACGT" * 1000  # 4000bp, comfortably covers every test position
        for i in range(0, len(seq), 60):
            fout.write(seq[i:i + 60] + "\n")
    pysam.faidx(sf)
    return sf


@pytest.fixture
def bam_file(tmp_path):
    sf = str(tmp_path / "reads.bam")
    header = {"HD": {"VN": "1.6"}, "SQ": [{"SN": "chr1", "LN": 4000}]}
    with pysam.AlignmentFile(sf, "wb", header=header):
        pass  # header only, no reads needed -- gnrt_vcf_head only reads m_header['SQ']
    return sf


def _data_lines(sf):
    with open(sf) as fin:
        return [ln for ln in fin if ln and not ln.startswith("#")]


def test_multi_matches_single_type_output(tmp_path, ref_fasta, bam_file):
    sf_raw = tmp_path / "hc_pf.txt"
    sf_raw.write_text(_row("chr1", 100, 110) + "\n")

    gvcf = gVCF()
    sf_single = str(tmp_path / "single.vcf")
    gvcf.cvt_raw_rslt_to_gvcf("SAMPLE", bam_file, str(sf_raw), 1, ref_fasta, sf_single)  # 1 = L1

    sf_type_map = tmp_path / "type_map.tsv"
    sf_type_map.write_text("chr1\t100\t1\n")
    sf_multi_l1 = str(tmp_path / "multi_L1.vcf")
    sf_merged = str(tmp_path / "merged.vcf")
    gvcf.cvt_raw_rslt_to_gvcf_multi("SAMPLE", bam_file, str(sf_raw), str(sf_type_map), ref_fasta,
                                     {1: sf_multi_l1}, sf_merged)

    single_lines = _data_lines(sf_single)
    assert len(single_lines) == 1
    assert single_lines == _data_lines(sf_multi_l1) == _data_lines(sf_merged)


def test_multi_dispatches_by_per_record_type(tmp_path, ref_fasta, bam_file):
    rows = [
        _row("chr1", 100, 110),    # -> L1  (rtype 1)
        _row("chr1", 1000, 1010),  # -> Alu (rtype 2)
        _row("chr1", 2000, 2010),  # -> SVA (rtype 4)
    ]
    sf_raw = tmp_path / "hc_pf.txt"
    sf_raw.write_text("\n".join(rows) + "\n")
    sf_type_map = tmp_path / "type_map.tsv"
    sf_type_map.write_text("chr1\t100\t1\nchr1\t1000\t2\nchr1\t2000\t4\n")

    gvcf = gVCF()
    m_sf_vcf = {1: str(tmp_path / "L1.vcf"), 2: str(tmp_path / "Alu.vcf"), 4: str(tmp_path / "SVA.vcf")}
    sf_merged = str(tmp_path / "merged.vcf")
    gvcf.cvt_raw_rslt_to_gvcf_multi("SAMPLE", bam_file, str(sf_raw), str(sf_type_map), ref_fasta,
                                     m_sf_vcf, sf_merged)

    def _alt(sf):
        lines = _data_lines(sf)
        assert len(lines) == 1
        return lines[0].split("\t")[4]

    assert _alt(m_sf_vcf[1]) == "<INS:ME:LINE1>"
    assert _alt(m_sf_vcf[2]) == "<INS:ME:ALU>"
    assert _alt(m_sf_vcf[4]) == "<INS:ME:SVA>"
    assert len(_data_lines(sf_merged)) == 3


def test_multi_missing_type_map_entry_raises(tmp_path, ref_fasta, bam_file):
    sf_raw = tmp_path / "hc_pf.txt"
    sf_raw.write_text(_row("chr1", 100, 110) + "\n")
    sf_type_map = tmp_path / "type_map.tsv"
    sf_type_map.write_text("")  # no entry for chr1:100

    gvcf = gVCF()
    with pytest.raises(RuntimeError):
        gvcf.cvt_raw_rslt_to_gvcf_multi("SAMPLE", bam_file, str(sf_raw), str(sf_type_map), ref_fasta,
                                         {1: str(tmp_path / "L1.vcf")}, str(tmp_path / "merged.vcf"))


def _info(line):
    return dict(kv.split("=", 1) for kv in line.split("\t")[7].split(";") if "=" in kv)


def test_other_families_tag_from_the_type_map(tmp_path, ref_fasta, bam_file):
    """A locus that beat other calls at the merge says so in INFO; one that beat
    nothing (4th column '.') gets no tag. The ##INFO header declares it."""
    sf_raw = tmp_path / "hc_pf.txt"
    sf_raw.write_text(_row("chr1", 100, 110) + "\n" + _row("chr1", 2000, 2010) + "\n")
    sf_type_map = tmp_path / "type_map.tsv"
    sf_type_map.write_text("chr1\t100\t4\tL1:100,Alu:120\nchr1\t2000\t4\t.\n")

    gvcf = gVCF()
    sf_sva, sf_merged = str(tmp_path / "SVA.vcf"), str(tmp_path / "merged.vcf")
    gvcf.cvt_raw_rslt_to_gvcf_multi("SAMPLE", bam_file, str(sf_raw), str(sf_type_map), ref_fasta,
                                     {4: sf_sva}, sf_merged)
    recs = _data_lines(sf_merged)
    assert _info(recs[0])["OTHER_FAMILIES"] == "L1:100,Alu:120"
    assert "OTHER_FAMILIES" not in _info(recs[1])
    assert _info(recs[0])["GENE_INFO"]          # the tag did not eat the next field
    assert any(ln.startswith("##INFO=<ID=OTHER_FAMILIES,") for ln in open(sf_merged))


def test_three_column_type_map_writes_no_tag(tmp_path, ref_fasta, bam_file):
    sf_raw = tmp_path / "hc_pf.txt"
    sf_raw.write_text(_row("chr1", 100, 110) + "\n")
    sf_type_map = tmp_path / "type_map.tsv"
    sf_type_map.write_text("chr1\t100\t1\n")
    gvcf = gVCF()
    sf_merged = str(tmp_path / "merged.vcf")
    gvcf.cvt_raw_rslt_to_gvcf_multi("SAMPLE", bam_file, str(sf_raw), str(sf_type_map), ref_fasta,
                                     {1: str(tmp_path / "L1.vcf")}, sf_merged)
    assert "OTHER_FAMILIES" not in _info(_data_lines(sf_merged)[0])
