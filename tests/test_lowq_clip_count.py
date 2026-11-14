"""Which reads the cns step counts as low-MAPQ clips (low_mapq_clip_cnt.txt).

A candidate is dropped when MAPQ<5 clipped reads within 35 bp are more than 65% of its clips
(MAX_LOWQ_CLIP_RATIO). Master counted these reads before its supplementary/secondary skip, so
a supplementary alignment counted too. Whether that matters depends on the aligner flags:

  * bwa mem without -Y (HG002, the mosaic HapMap BAMs) hard-clips supplementaries, and a
    hard-clipped read has no soft clip at either end, so it was never counted.
  * bwa mem -Y (the 1000 Genomes CRAMs) soft-clips them. A read whose primary alignment lies
    in the inserted element elsewhere then leaves a MAPQ-0 supplementary with its clip exactly
    at the junction -- evidence FOR the insertion, counted as background against it. On the
    1KG CRAMs 40-44% of the counted (read, candidate) pairs were such supplementaries, and the
    filter dropped true sites at about the rate it dropped false ones.

The count now takes primary alignments only (2026-11-13), which makes it the same with or
without -Y. Both code paths are pinned: the native scanner (cns_emit.hpp, single-file and
manifest modes) and the Python port in x_clip_disc_filter.py.
"""
import os
import subprocess
import sys

import pysam
import pytest

REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(REPO, "xtea"))
import global_values  # noqa: E402
import x_clip_disc_filter as xcdf  # noqa: E402

CHROM = "chr1"
CAND = 2000
RLEN = 150
SUPP, SEC = 2048, 256

# (name, flag, start, cigar, counted) -- every read MAPQ 0, clip anchored at CAND.
READS = [
    ("prim_left", 0, CAND, [(4, 50), (0, 100)], True),
    ("prim_right", 0, CAND - 100, [(0, 100), (4, 50)], True),
    ("supp_left", SUPP, CAND, [(4, 50), (0, 100)], False),
    ("supp_right", SUPP | 16, CAND - 100, [(0, 100), (4, 50)], False),
    ("sec_left", SEC, CAND, [(4, 50), (0, 100)], False),
]


def _segment(header, name, flag, start, cigar, mapq=0):
    a = pysam.AlignedSegment(header)
    a.query_name = name
    a.flag = flag
    a.reference_id = 0
    a.reference_start = start
    a.mapping_quality = mapq
    a.cigartuples = cigar
    n = sum(l for op, l in cigar if op in (0, 1, 4))
    a.query_sequence = "ACGT" * (n // 4) + "ACGT"[: n % 4]
    a.query_qualities = pysam.qualitystring_to_array("I" * n)
    a.next_reference_id = -1
    a.next_reference_start = -1
    return a


HEADER = {"HD": {"VN": "1.6", "SO": "coordinate"}, "SQ": [{"SN": CHROM, "LN": 10000}]}


# ---- Python port (XClipDisc) --------------------------------------------------------------

@pytest.mark.parametrize("name, flag, start, cigar, counted", READS)
def test_python_counts_primary_only(name, flag, start, cigar, counted):
    h = pysam.AlignmentHeader.from_dict(HEADER)
    a = _segment(h, name, flag, start, cigar)
    old = global_values.READ_LENGTH
    global_values.READ_LENGTH = RLEN
    try:
        pos = xcdf.low_mapq_clip_anchor(a)
    finally:
        global_values.READ_LENGTH = old
    if counted:
        assert pos is not None and abs(pos - CAND) < global_values.NEARBY_LOW_Q_CLIP
    else:
        assert pos is None, "supplementary/secondary low-MAPQ clips must not count"


def test_python_unclipped_not_counted():
    h = pysam.AlignmentHeader.from_dict(HEADER)
    assert xcdf.low_mapq_clip_anchor(_segment(h, "m", 0, CAND, [(0, 150)])) is None
    # a hard clip is not a soft clip: never counted (unchanged from master)
    assert xcdf.low_mapq_clip_anchor(_segment(h, "h", 0, CAND, [(5, 50), (0, 100)])) is None


# ---- native scanner (cns_emit.hpp) --------------------------------------------------------

def _write_inputs(tmp_path):
    bam = str(tmp_path / "t.bam")
    with pysam.AlignmentFile(bam, "wb", header=HEADER) as out:
        segs = [_segment(out.header, *r[:4]) for r in READS]
        for a in sorted(segs, key=lambda s: s.reference_start):
            out.write(a)
    pysam.index(bam)
    cand = tmp_path / "cand.txt"
    cand.write_text("%s %d\n" % (CHROM, CAND))
    cns = tmp_path / "cns.fa"
    cns.write_text(">c\n" + "ACGT" * 50 + "\n")
    return bam, str(cand), str(cns)


def _expected():
    return sum(1 for r in READS if r[4])


def _read_lowq(path):
    m = {}
    for ln in open(path):
        c, p, n = ln.split()
        m[(c, int(p))] = int(n)
    return m


def test_scanner_single_counts_primary_only(tmp_path, xtea2_scan_bin):
    bam, cand, cns = _write_inputs(tmp_path)
    lowq = str(tmp_path / "lowq.txt")
    # bwa is irrelevant here (no read passes the clip MAPQ gate); /bin/true stands in for it.
    subprocess.run([xtea2_scan_bin, "-b", bam, "--cns-candidates", cand, "--cns-ref", cns,
                    "--cns-clip-sam", str(tmp_path / "clip.sam"), "--cns-lowq", lowq,
                    "--cns-read-length", str(RLEN), "--cns-tmp-prefix", str(tmp_path / "tmp"),
                    "--bwa", "/bin/true"], check=True, capture_output=True)
    assert _read_lowq(lowq) == {(CHROM, CAND): _expected()}


def test_scanner_manifest_counts_primary_only(tmp_path, xtea2_scan_bin):
    bam, cand, cns = _write_inputs(tmp_path)
    lowq = str(tmp_path / "lowq.txt")
    man = tmp_path / "manifest.txt"
    man.write_text("Alu %s %s %s %s %s\n" % (cand, cns, tmp_path / "clip.sam",
                                             tmp_path / "disc.sam", lowq))
    subprocess.run([xtea2_scan_bin, "-b", bam, "--cns-emit-manifest", str(man),
                    "--cns-read-length", str(RLEN), "--cns-tmp-prefix", str(tmp_path / "tmp"),
                    "--bwa", "/bin/true"], check=True, capture_output=True)
    assert _read_lowq(lowq) == {(CHROM, CAND): _expected()}
