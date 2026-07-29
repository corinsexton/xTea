"""`extract_mate_reads_of_region` must emit each disc mate ONCE per insertion site.

The disc FASTA it writes is the read-level evidence for a candidate: one record per
(mate read, insertion site). Emitting the same mate twice counts one read as two and
inflates that site's disc support, which moves calls across the cns filter's cutoffs.

Two things used to double-emit (both found on the NA12878 WGS run, where 28,173 of
583,489 Alu disc records -- 4.8% -- were duplicates):

  * a mate with a soft-clipped SUPPLEMENTARY alignment in the same bin matched twice
    (is_read1 is set on the supplementary too, and a soft-clipped supplementary carries
    the full query_sequence, so both records were byte-identical);
  * fetch() returns reads OVERLAPPING a bin while _parse_read_name_pos_of_region admits
    mates up to 100 bp beyond either edge, so a mate near a bin boundary was emitted by
    both adjacent bins.

These tests build a tiny synthetic BAM, so they need no external data and run in the fast
suite. The reference is a single 10 kb contig; sequences are arbitrary.
"""
import os
import sys

import pytest

pysam = pytest.importorskip("pysam")

# xtea/ modules use implicit relative imports (`from x_intermediate_sites import *`), so the
# package dir itself has to be importable -- same shim as test_downstream_from_cns.py.
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "xtea"))
from x_alignments import BamInfo  # noqa: E402

CHROM = "chr1"
CHROM_LEN = 10000
QNAME = "READ_A"
MATE_SEQ = "ACGT" * 25          # 100 bp -- the read1 mate whose seq lands in the FASTA
ANCHOR_SEQ = "TTGC" * 25
INS_POSITIONS = ("3000", "3500")  # two candidate sites this one mate supports


def _write_bam(path, records):
    """records: list of dicts with keys name/flag/pos/cigar/seq/mate_pos."""
    hdr = {"HD": {"VN": "1.6", "SO": "coordinate"},
           "SQ": [{"SN": CHROM, "LN": CHROM_LEN}]}
    with pysam.AlignmentFile(path, "wb", header=hdr) as out:
        for r in sorted(records, key=lambda x: x["pos"]):
            a = pysam.AlignedSegment(out.header)
            a.query_name = r["name"]
            a.flag = r["flag"]
            a.reference_id = 0
            a.reference_start = r["pos"]
            a.mapping_quality = 60
            a.cigarstring = r["cigar"]
            a.query_sequence = r["seq"]
            a.query_qualities = pysam.qualitystring_to_array("I" * len(r["seq"]))
            a.next_reference_id = 0
            a.next_reference_start = r["mate_pos"]
            out.write(a)
    pysam.index(path)


def _write_disc_pos(path, mate_pos, anchor_pos=3200):
    """One disc-pos line per insertion site, as collect_clipped_disc_reads_by_region writes them.

    Columns: ins_chrm, anchor_map_pos, anchor_rc, mate_rc, mate_chrm, mate_pos, qname,
             s_mate_first, ins_pos, b_lclip, b_rclip.
    s_mate_first=1 means "the mate is read1", which is what the extractor matches on.
    """
    with open(path, "w") as fh:
        for ins in INS_POSITIONS:
            fh.write("\t".join([CHROM, str(anchor_pos), "0", "1", CHROM, str(mate_pos),
                                QNAME, "1", ins, "0", "0"]) + "\n")


def _run_bins(tmpdir, bam, disc_pos, bins):
    """Run the extractor over `bins` and return the merged FASTA header lines."""
    info = BamInfo(bam, None)
    headers = []
    for bin_start, bin_end in bins:
        info.extract_mate_reads_of_region((CHROM, bin_start, bin_end, disc_pos, str(tmpdir) + "/"))
        fa = os.path.join(str(tmpdir), "{0}_{1}_{2}.tmp.disc.fa".format(CHROM, bin_start, bin_end))
        with open(fa) as fh:
            headers += [l.rstrip("\n") for l in fh if l.startswith(">")]
    return headers


def test_one_record_per_insertion_site(tmp_path):
    """Baseline: a lone primary mate yields exactly one record per insertion site."""
    bam = str(tmp_path / "plain.bam")
    _write_bam(bam, [
        {"name": QNAME, "flag": 83, "pos": 4000, "cigar": "100M", "seq": MATE_SEQ, "mate_pos": 3200},
        {"name": QNAME, "flag": 163, "pos": 3200, "cigar": "100M", "seq": ANCHOR_SEQ, "mate_pos": 4000},
    ])
    dp = str(tmp_path / "disc.pos")
    _write_disc_pos(dp, mate_pos=4000)

    headers = _run_bins(tmp_path, bam, dp, [(0, CHROM_LEN)])
    assert len(headers) == len(INS_POSITIONS) == len(set(headers))


def test_supplementary_mate_alignment_does_not_duplicate(tmp_path):
    """A soft-clipped supplementary alignment of the mate must not add a second record.

    This is the WGS case: flag 2113 = PAIRED|READ1|SUPPLEMENTARY, soft-clipped, so it
    reports is_read1 and carries the full query_sequence -- indistinguishable from the
    primary at the point the FASTA header is built.
    """
    bam = str(tmp_path / "supp.bam")
    _write_bam(bam, [
        {"name": QNAME, "flag": 83, "pos": 4000, "cigar": "100M", "seq": MATE_SEQ, "mate_pos": 3200},
        # read1 supplementary elsewhere in the same bin, soft-clipped -> full SEQ
        {"name": QNAME, "flag": 2113, "pos": 6000, "cigar": "60S40M", "seq": MATE_SEQ, "mate_pos": 3200},
        {"name": QNAME, "flag": 163, "pos": 3200, "cigar": "100M", "seq": ANCHOR_SEQ, "mate_pos": 4000},
    ])
    dp = str(tmp_path / "disc.pos")
    _write_disc_pos(dp, mate_pos=4000)

    headers = _run_bins(tmp_path, bam, dp, [(0, CHROM_LEN)])
    assert headers == sorted(set(headers), key=headers.index), \
        "supplementary alignment of the mate duplicated its FASTA record: %r" % (headers,)
    assert len(headers) == len(INS_POSITIONS)


def test_secondary_mate_alignment_does_not_duplicate(tmp_path):
    """Same for a secondary alignment (flag 0x100), which also carries the full SEQ."""
    bam = str(tmp_path / "sec.bam")
    _write_bam(bam, [
        {"name": QNAME, "flag": 83, "pos": 4000, "cigar": "100M", "seq": MATE_SEQ, "mate_pos": 3200},
        {"name": QNAME, "flag": 337, "pos": 6000, "cigar": "100M", "seq": MATE_SEQ, "mate_pos": 3200},
        {"name": QNAME, "flag": 163, "pos": 3200, "cigar": "100M", "seq": ANCHOR_SEQ, "mate_pos": 4000},
    ])
    dp = str(tmp_path / "disc.pos")
    _write_disc_pos(dp, mate_pos=4000)

    headers = _run_bins(tmp_path, bam, dp, [(0, CHROM_LEN)])
    assert len(headers) == len(INS_POSITIONS) == len(set(headers))


def test_mate_near_bin_boundary_is_emitted_exactly_once(tmp_path):
    """Adjacent bins must partition the mates, not overlap them.

    The mate sits 40 bp below the boundary: fetch() returns it for the lower bin (it starts
    there) and for the upper bin (it extends past the edge), and the +/-100 bp slack in
    _parse_read_name_pos_of_region puts it in both bins' name maps. Only the bin containing
    its start may emit it.
    """
    boundary = 5000
    mate_pos = boundary - 40
    bam = str(tmp_path / "edge.bam")
    _write_bam(bam, [
        {"name": QNAME, "flag": 83, "pos": mate_pos, "cigar": "100M", "seq": MATE_SEQ, "mate_pos": 3200},
        {"name": QNAME, "flag": 163, "pos": 3200, "cigar": "100M", "seq": ANCHOR_SEQ, "mate_pos": mate_pos},
    ])
    dp = str(tmp_path / "disc.pos")
    _write_disc_pos(dp, mate_pos=mate_pos)

    headers = _run_bins(tmp_path, bam, dp, [(0, boundary), (boundary, CHROM_LEN)])
    assert len(headers) == len(INS_POSITIONS) == len(set(headers)), \
        "mate straddling a bin boundary emitted by both bins: %r" % (headers,)


def test_mate_is_not_lost_when_it_starts_in_the_upper_bin(tmp_path):
    """Owning by start must not drop anything: a mate just above a boundary still appears."""
    boundary = 5000
    mate_pos = boundary + 10
    bam = str(tmp_path / "edge2.bam")
    _write_bam(bam, [
        {"name": QNAME, "flag": 83, "pos": mate_pos, "cigar": "100M", "seq": MATE_SEQ, "mate_pos": 3200},
        {"name": QNAME, "flag": 163, "pos": 3200, "cigar": "100M", "seq": ANCHOR_SEQ, "mate_pos": mate_pos},
    ])
    dp = str(tmp_path / "disc.pos")
    _write_disc_pos(dp, mate_pos=mate_pos)

    headers = _run_bins(tmp_path, bam, dp, [(0, boundary), (boundary, CHROM_LEN)])
    assert len(headers) == len(INS_POSITIONS) == len(set(headers))
    assert all(h.endswith(tuple(INS_POSITIONS)) for h in headers)
