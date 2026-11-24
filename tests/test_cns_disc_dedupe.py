"""The cns step's disc FASTA holds one record per (read, mate, insertion site).

Z5 (Zofia, SMaHT P25, 2026-11): when a disc anchor read also has a hard-clipped supplementary
alignment near the same candidate, the native scanner (cns_emit.hpp) emitted the mate once per
alignment, so one template counted as two disc reads. Only the soft-clip blocks skip
supplementary/secondary alignments; a hard-clipped supplementary has no soft clip and reached the
disc block. The extra copies landed on the "missing" side of an Alu at chr10:133221839 and made it
look two-sided (3 FP across donors). Master writes one record per (qname~mate, site) because its
_parse_read_name_pos_of_region keys a dict on them, keeping the first occurrence's position and the
last occurrence's values.

bwa is replaced by a script that echoes the FASTA it is given, so disc.sam holds the FASTA records
the scanner wrote and the test needs no external data.
"""
import os
import stat
import subprocess
from collections import Counter

import pysam

CHROM = "chr1"
CAND = 2000
MATE_POS = 8000            # same chrom, far beyond the disc distance threshold
RLEN = 150
R1, R2, REV, SUPP = 64, 128, 16, 2048
HEADER = {"HD": {"VN": "1.6", "SO": "coordinate"}, "SQ": [{"SN": CHROM, "LN": 12000}]}


def _seg(header, name, flag, start, cigar, mate_start):
    a = pysam.AlignedSegment(header)
    a.query_name = name
    a.flag = flag | 1           # paired
    a.reference_id = 0
    a.reference_start = start
    a.mapping_quality = 60
    a.cigartuples = cigar
    n = sum(l for op, l in cigar if op in (0, 1, 4))
    a.query_sequence = ("ACGTTGCA" * 40)[:n]
    a.query_qualities = pysam.qualitystring_to_array("I" * n)
    a.next_reference_id = 0
    a.next_reference_start = mate_start
    return a


def _inputs(tmp_path, with_supp):
    bam = str(tmp_path / "t.bam")
    with pysam.AlignmentFile(bam, "wb", header=HEADER) as out:
        h = out.header
        segs = [
            _seg(h, "T1", R1, CAND - 60, [(0, RLEN)], MATE_POS),                 # anchor primary
            _seg(h, "T1", R2 | REV, MATE_POS, [(0, RLEN)], CAND - 60),           # its mate
            _seg(h, "T2", R1, CAND - 40, [(0, RLEN)], MATE_POS + 20),            # an unrelated template
            _seg(h, "T2", R2 | REV, MATE_POS + 20, [(0, RLEN)], CAND - 40),
        ]
        if with_supp:   # hard-clipped supplementary of T1's read1, still within the candidate window
            segs.append(_seg(h, "T1", R1 | SUPP, CAND + 10, [(5, 50), (0, 100)], MATE_POS))
        for a in sorted(segs, key=lambda s: s.reference_start):
            out.write(a)
    pysam.index(bam)
    cand = tmp_path / "cand.txt"
    cand.write_text("%s %d\n" % (CHROM, CAND))
    cns = tmp_path / "cns.fa"
    cns.write_text(">c\n" + "ACGT" * 50 + "\n")
    fake_bwa = tmp_path / "fake_bwa.sh"
    fake_bwa.write_text('#!/bin/sh\nfor a; do last="$a"; done\ncat "$last"\n')
    fake_bwa.chmod(fake_bwa.stat().st_mode | stat.S_IEXEC)
    return bam, str(cand), str(cns), str(fake_bwa)


def _records(path):
    """(qname, mate, ins_pos) of every disc FASTA header the scanner wrote."""
    out = []
    for ln in open(path):
        if ln.startswith(">"):
            f = ln[1:].split("~")
            out.append((f[0], f[1], int(f[8])))
    return out


def _run_single(tmp_path, with_supp, xtea2_scan_bin):
    bam, cand, cns, bwa = _inputs(tmp_path, with_supp)
    disc = str(tmp_path / "disc.sam")
    subprocess.run([xtea2_scan_bin, "-b", bam, "--cns-candidates", cand, "--cns-ref", cns,
                    "--cns-disc-sam", disc, "--cns-read-length", str(RLEN),
                    "--cns-tmp-prefix", str(tmp_path / "tmp"), "--bwa", bwa],
                   check=True, capture_output=True)
    return _records(disc)


def _run_manifest(tmp_path, with_supp, xtea2_scan_bin):
    bam, cand, cns, bwa = _inputs(tmp_path, with_supp)
    disc = tmp_path / "disc.sam"
    man = tmp_path / "manifest.txt"
    man.write_text("Alu %s %s %s %s %s\n" % (cand, cns, tmp_path / "clip.sam", disc,
                                             tmp_path / "lowq.txt"))
    subprocess.run([xtea2_scan_bin, "-b", bam, "--cns-emit-manifest", str(man),
                    "--cns-read-length", str(RLEN), "--cns-tmp-prefix", str(tmp_path / "tmp"),
                    "--bwa", bwa], check=True, capture_output=True)
    return _records(str(disc))


def _check(recs):
    c = Counter(recs)
    assert c == {("T1", "0", CAND): 1, ("T2", "0", CAND): 1}, c


def test_single_file_one_record_per_template_and_site(tmp_path, xtea2_scan_bin):
    _check(_run_single(tmp_path, True, xtea2_scan_bin))


def test_manifest_one_record_per_template_and_site(tmp_path, xtea2_scan_bin):
    _check(_run_manifest(tmp_path, True, xtea2_scan_bin))


def test_without_supplementary_unchanged(tmp_path, xtea2_scan_bin):
    # control: the fixture itself emits exactly one record per template
    (tmp_path / "s").mkdir(); (tmp_path / "m").mkdir()
    _check(_run_single(tmp_path / "s", False, xtea2_scan_bin))
    _check(_run_manifest(tmp_path / "m", False, xtea2_scan_bin))
