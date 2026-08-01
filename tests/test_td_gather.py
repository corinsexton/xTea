"""The native transduction gather must reproduce the Python one BYTE for byte.

WHY BYTES AND NOT RECORDS
-------------------------
The transduction gather (73-87% of the transduction step, itself ~93% of downstream wall time) is
ported to C++ in src/csrc/td_gather.hpp. Its two products, the clip FASTQ and the disc-mate FASTA,
are bwa inputs -- and bwa's tie-break among EQUAL-SCORING placements depends on a read's index in
its input (measured, cns_emit.hpp:798-817). A reordered FASTQ therefore silently moves which
repeat copy a multi-mapping clip lands on, and with it borderline calls. So "same reads" is not the
property that keeps Stage 1 agreement-gated; "same file" is.

Three ordering rules carry that, and each has a test here:
  * the clip FASTQ is candidate-major in the SITE FILE's order -- chromosomes in first-appearance
    order, sites in file order -- which is what the legacy per-chromosome Pool plus per-site file
    concatenation produces. The site file is NOT position-sorted (cns candidates are written
    first, then raw-disc lines), so `test_sites_file_order_is_the_emission_order` deliberately
    feeds an unsorted list;
  * n_cnt (read-name field 5) keeps master's GAPS, because master builds the right-clip name and
    increments the counter before applying that side's gates;
  * the disc FASTA is BAM-header-contig-major then coordinate-ordered, per mate read carrying its
    insertion sites in disc-pos order.

Everything skips cleanly when XTEA2_SCAN_BIN or the test BAM is absent.
"""
import os
import subprocess
import sys

import pysam
import pytest

REPO_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(REPO_ROOT, "xtea"))

EXTND = 550          # iextnd for a real run is max(400, mean_is + 3*std); 550 is representative
SITE_WINDOW = ("chr1", 23_000_000, 23_600_000)


def _clip_positions(bam, reference, chrom, start, end, want=40, min_gap=1500):
    """Positions carrying real soft-clip signal, so the gather has something to collect."""
    out = []
    with pysam.AlignmentFile(bam, "rb", reference_filename=reference) as bf:
        for a in bf.fetch(chrom, start, end):
            if a.is_duplicate or a.is_unmapped or not a.cigartuples:
                continue
            if a.mapping_quality < 12:
                continue
            if a.cigartuples[0][0] == 4 and a.cigartuples[0][1] >= 15:
                p = a.reference_start
            elif a.cigartuples[-1][0] == 4 and a.cigartuples[-1][1] >= 15:
                p = a.reference_end
            else:
                continue
            if out and p - out[-1] < min_gap:
                continue
            out.append(p)
            if len(out) >= want:
                break
    return out


@pytest.fixture(scope="module")
def site_positions(test_bam, reference):
    pos = _clip_positions(test_bam, reference, *SITE_WINDOW)
    if len(pos) < 10:
        pytest.skip(f"too few clip-bearing positions in {SITE_WINDOW} to exercise the gather")
    return pos


def _write_sites(path, chrom, positions):
    with open(path, "w") as f:
        for p in positions:
            f.write(f"{chrom}\t{p}\n")
    return str(path)


def _run_native(scan_bin, bam, reference, sites, outdir, threads=2, disc_fa=True):
    os.makedirs(outdir, exist_ok=True)
    clip_fq = os.path.join(outdir, "clip.fq")
    disc_pos = os.path.join(outdir, "all_disc_pos.disc_pos")
    cmd = [scan_bin, "--bam", bam, "--reference", reference, "--td-sites", sites,
           "--td-extnd", str(EXTND), "--td-clip-fastq", clip_fq,
           "--td-disc-pos", disc_pos, "--threads", str(threads)]
    fa = os.path.join(outdir, "disc.fa")
    if disc_fa:
        cmd += ["--td-disc-fa", fa]
    subprocess.run(cmd, check=True)
    return clip_fq, disc_pos, (fa if disc_fa else None)


def _run_python_clip(bam, reference, sites, outdir, threads=2):
    """The legacy clip gather. It also writes all_disc_pos.disc_pos as a side effect."""
    from x_clip_disc_filter import XClipDiscFilter
    import global_values
    wf = str(outdir) + "/"
    os.makedirs(wf, exist_ok=True)
    bam_list = os.path.join(wf, "bam_list.txt")
    with open(bam_list, "w") as f:
        f.write(bam + "\n")
    xf = XClipDiscFilter(bam_list, wf, threads, reference)
    clip_fq = wf + "clip.fq"
    xf.collect_clipped_reads(sites, EXTND, clip_fq)
    return clip_fq, wf + "all_disc_pos" + global_values.DISC_POS_SUFFIX, xf


def _assert_same(label, a, b):
    ba, bb = open(a, "rb").read(), open(b, "rb").read()
    if ba == bb:
        return
    la, lb = ba.splitlines(), bb.splitlines()
    first = next((i for i, (x, y) in enumerate(zip(la, lb)) if x != y), min(len(la), len(lb)))
    pytest.fail(f"{label} differs: python {len(la)} lines / {len(ba)} bytes, "
                f"native {len(lb)} lines / {len(bb)} bytes; first difference at line {first}:\n"
                f"  py: {la[first] if first < len(la) else '<eof>'!r}\n"
                f"  cpp: {lb[first] if first < len(lb) else '<eof>'!r}")


def test_clip_fastq_and_disc_pos_are_byte_identical(test_bam, reference, xtea2_scan_bin,
                                                    site_positions, tmp_path):
    sites = _write_sites(tmp_path / "sites.txt", SITE_WINDOW[0], site_positions)
    py_fq, py_dp, _ = _run_python_clip(test_bam, reference, sites, tmp_path / "py")
    cpp_fq, cpp_dp, _ = _run_native(xtea2_scan_bin, test_bam, reference, sites,
                                    str(tmp_path / "cpp"), disc_fa=False)
    assert os.path.getsize(py_fq) > 0, "the Python gather collected nothing — bad test window"
    _assert_same("clip FASTQ", py_fq, cpp_fq)
    _assert_same("all_disc_pos", py_dp, cpp_dp)


def test_sites_file_order_is_the_emission_order(test_bam, reference, xtea2_scan_bin,
                                                site_positions, tmp_path):
    """An UNSORTED site list, as sf_raw_slct really is: cns candidates first, raw-disc after.

    Sorting the list inside the gather would pass the previous test and fail this one, which is
    the failure mode worth pinning -- cns_emit's fan sorts, and copying that here would be wrong.
    """
    n = len(site_positions)
    shuffled = site_positions[n // 2:] + site_positions[:n // 2]   # deliberately non-monotonic
    sites = _write_sites(tmp_path / "sites.txt", SITE_WINDOW[0], shuffled)
    py_fq, py_dp, _ = _run_python_clip(test_bam, reference, sites, tmp_path / "py")
    cpp_fq, cpp_dp, _ = _run_native(xtea2_scan_bin, test_bam, reference, sites,
                                    str(tmp_path / "cpp"), disc_fa=False)
    _assert_same("clip FASTQ (unsorted sites)", py_fq, cpp_fq)
    _assert_same("all_disc_pos (unsorted sites)", py_dp, cpp_dp)
    # And the emission really is in the file's order, not the sorted one.
    seen = []
    with open(cpp_fq) as f:
        for i, line in enumerate(f):
            if i % 4:
                continue
            ins_pos = int(line.split("~")[4])
            if not seen or seen[-1] != ins_pos:
                seen.append(ins_pos)
    assert seen == [p for p in shuffled if p in set(seen)], \
        "clip FASTQ is not grouped by candidate in site-file order"


def test_native_gather_is_thread_count_invariant(test_bam, reference, xtea2_scan_bin,
                                                 site_positions, tmp_path):
    """Region stripes must not leak into the output; 1 thread and 4 must agree byte for byte."""
    sites = _write_sites(tmp_path / "sites.txt", SITE_WINDOW[0], site_positions)
    fq1, dp1, _ = _run_native(xtea2_scan_bin, test_bam, reference, sites,
                              str(tmp_path / "t1"), threads=1, disc_fa=False)
    fq4, dp4, _ = _run_native(xtea2_scan_bin, test_bam, reference, sites,
                              str(tmp_path / "t4"), threads=4, disc_fa=False)
    _assert_same("clip FASTQ (threads 1 vs 4)", fq1, fq4)
    _assert_same("all_disc_pos (threads 1 vs 4)", dp1, dp4)


def test_clip_fastq_n_cnt_keeps_masters_gaps(test_bam, reference, xtea2_scan_bin,
                                            site_positions, tmp_path):
    """Field 5 must be per-candidate and may skip values; it must never be a dense 0..k-1 run.

    The gaps come from master incrementing the counter before the right-clip gates. Nothing reads
    field 5 (parse_transduction_clip_algnmt uses 0,1,2,4), so this pins the mechanism rather than
    a consumed value -- but without it the byte-identity above would only hold by luck.
    """
    sites = _write_sites(tmp_path / "sites.txt", SITE_WINDOW[0], site_positions)
    cpp_fq, _, _ = _run_native(xtea2_scan_bin, test_bam, reference, sites,
                               str(tmp_path / "cpp"), disc_fa=False)
    per_site = {}
    with open(cpp_fq) as f:
        for i, line in enumerate(f):
            if i % 4:
                continue
            fields = line.rstrip().split("~")
            per_site.setdefault(int(fields[4]), []).append(int(fields[5]))
    assert per_site, "no clip records emitted"
    for site, idx in per_site.items():
        assert idx == sorted(idx), f"n_cnt is not ascending within candidate {site}"
        assert len(set(idx)) == len(idx), f"n_cnt repeats within candidate {site}"
    assert any(idx[-1] >= len(idx) for idx in per_site.values()), (
        "no candidate shows an n_cnt gap — the right-clip index burn is not being reproduced "
        "(or this window has no gated-out right clips, in which case widen SITE_WINDOW)")


@pytest.mark.slow
def test_disc_mate_fasta_is_byte_identical(test_bam, reference, xtea2_scan_bin,
                                           site_positions, tmp_path):
    """Marked slow: the legacy disc gather scans the WHOLE genome in 50 Mb bins (~2 min here,
    ~280s on WGS) regardless of how few sites it was given. That fixed cost is most of what the
    port removes."""
    sites = _write_sites(tmp_path / "sites.txt", SITE_WINDOW[0], site_positions)
    py_fq, py_dp, xf = _run_python_clip(test_bam, reference, sites, tmp_path / "py")
    py_fa = str(tmp_path / "py") + "/disc.fa"
    xf.collect_disc_reads(sites, EXTND, 50_000_000, py_fa)
    _, _, cpp_fa = _run_native(xtea2_scan_bin, test_bam, reference, sites, str(tmp_path / "cpp"))
    assert os.path.getsize(py_fa) > 0, "the Python disc gather collected nothing — bad test window"
    _assert_same("disc mate FASTA", py_fa, cpp_fa)
