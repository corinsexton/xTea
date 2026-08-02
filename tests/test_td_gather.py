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


# ---------------------------------------------------------------------------------------------
# Hermetic multi-contig case.  The real-data tests above all run on a chr1-only BAM, so they
# cannot see the two contig-ordering rules -- and the two are DIFFERENT orders, which is exactly
# the kind of thing that looks right until a whole-genome run:
#   * the clip FASTQ follows the SITE FILE's contig order (first appearance), while
#   * the disc FASTA follows the BAM HEADER's contig order, decoy/contig names skipped.
# The synthetic header therefore lists chr2 before chr1 while the site file names chr1 first, and
# every mate is placed on its anchor's OWN contig, so the two files must come out in opposite
# contig orders -- using either rule for the other file fails. Also covers: a site window clamped
# at position 1, overlapping windows (one read fanned to two candidates), a mate on a decoy contig
# (never a disc record), and a left clip too short to pass (which must NOT burn an n_cnt index).
# ---------------------------------------------------------------------------------------------
SYN_HEADER = {"HD": {"VN": "1.6", "SO": "coordinate"},
              "SQ": [{"SN": "chr2", "LN": 20000},        # header order, deliberately not chr1 first
                     {"SN": "chr1", "LN": 20000},
                     {"SN": "chrUn_decoy", "LN": 5000}]}
SYN_SITES = [("chr1", 500), ("chr2", 3000), ("chr1", 5000),   # interleaved + non-monotonic
             ("chr1", 12000), ("chr1", 12300), ("chr2", 9000)]  # 12000/12300 windows overlap
# Where the gathered mates live: same contig as their anchor (>2000 bp away is discordant on its
# own, DISC_THRESHOLD), outside every site window, and staggered so no two share a position.
SYN_MATE_POS = 18000


def _syn_read(header, name, tid, pos, cigar, *, is_read1=True, mtid=None, mpos=None,
              is_reverse=False, mate_reverse=False):
    a = pysam.AlignedSegment(header)
    a.query_name = name
    a.reference_id = tid
    a.reference_start = pos
    a.mapping_quality = 60
    a.cigarstring = cigar
    n = sum(l for op, l in a.cigartuples if op in (0, 1, 4))
    a.query_sequence = "ACGT" * (n // 4) + "ACGT"[:n % 4]
    a.query_qualities = pysam.qualitystring_to_array("I" * n)   # phred 40, passes the clip gate
    a.is_paired = True
    a.is_read1, a.is_read2 = is_read1, not is_read1
    a.is_reverse = is_reverse
    a.mate_is_reverse = mate_reverse
    if mtid is None:
        a.mate_is_unmapped = True
        a.next_reference_id, a.next_reference_start = -1, -1
    else:
        a.next_reference_id, a.next_reference_start = mtid, mpos
    return a


@pytest.fixture(scope="module")
def synthetic_bam(tmp_path_factory):
    hdr = pysam.AlignmentHeader.from_dict(SYN_HEADER)
    tid = {sq["SN"]: i for i, sq in enumerate(SYN_HEADER["SQ"])}
    recs = []
    for i, (chrom, S) in enumerate(SYN_SITES):
        mate_pos = SYN_MATE_POS + 10 * i
        tag = f"{chrom}{S}"
        recs += [
            _syn_read(hdr, tag + "Lclip", tid[chrom], S - 30, "20S80M"),
            _syn_read(hdr, tag + "Rclip", tid[chrom], S - 80, "80M20S"),
            _syn_read(hdr, tag + "Rburn", tid[chrom], S - 50, "97M3S"),   # burns an n_cnt index
            _syn_read(hdr, tag + "Lshort", tid[chrom], S - 40, "3S97M"),  # dropped, no burn
            # disc anchor + its mate, which the mate gather has to find by position
            _syn_read(hdr, tag + "Disc", tid[chrom], S - 10, "100M",
                      is_read1=True, mtid=tid[chrom], mpos=mate_pos),
            _syn_read(hdr, tag + "Disc", tid[chrom], mate_pos, "100M",
                      is_read1=False, mtid=tid[chrom], mpos=S - 10),
            # mate on a decoy contig: collected as a clip, never as a disc record
            _syn_read(hdr, tag + "Dec", tid[chrom], S + 10, "100M",
                      is_read1=True, mtid=tid["chrUn_decoy"], mpos=1000),
        ]
    recs.sort(key=lambda a: (a.reference_id, a.reference_start))
    d = tmp_path_factory.mktemp("td_syn")
    path = str(d / "syn.bam")
    with pysam.AlignmentFile(path, "wb", header=hdr) as out:
        for a in recs:
            out.write(a)
    pysam.index(path)
    return path


def test_multi_contig_gather_matches_python(synthetic_bam, xtea2_scan_bin, tmp_path):
    sites = str(tmp_path / "sites.txt")
    with open(sites, "w") as f:
        for chrom, p in SYN_SITES:
            f.write(f"{chrom}\t{p}\n")
    py_dir = tmp_path / "py"
    py_fq, py_dp, xf = _run_python_clip(synthetic_bam, "", sites, py_dir)
    py_fa = str(py_dir) + "/disc.fa"
    xf.collect_disc_reads(sites, EXTND, 50_000_000, py_fa)
    cpp_fq, cpp_dp, cpp_fa = _run_native(xtea2_scan_bin, synthetic_bam, "", sites,
                                         str(tmp_path / "cpp"))
    _assert_same("clip FASTQ (multi-contig)", py_fq, cpp_fq)
    _assert_same("all_disc_pos (multi-contig)", py_dp, cpp_dp)
    _assert_same("disc mate FASTA (multi-contig)", py_fa, cpp_fa)

    # The construction is only meaningful if it actually produced every signal class it claims
    # to. Clip records are named from the SITE, not the read, so the read-level classes have to
    # be identified by clip position and side.
    heads = [ln for i, ln in enumerate(open(cpp_fq).read().splitlines()) if i % 4 == 0]
    assert any("~L~" in h for h in heads) and any("~R~" in h for h in heads), \
        "synthetic clips were not collected on both sides"
    first_site = [h for h in heads if h.startswith("@chr1~") and h.split("~")[4] == "500"]
    assert [h.split("~")[5] for h in first_site] == ["0", "2"], (
        "expected n_cnt 0 then 2 at chr1:500 -- 0 for the right clip, 1 burned by the 3 bp right "
        f"clip that then failed the length gate, 2 for the left clip; got {first_site}")
    assert not any(h.split("~")[1] == "460" and h.split("~")[2] == "L" for h in heads), \
        "a 3 bp left clip should not survive the length gate"
    dp = open(cpp_dp).read()
    assert "Disc" in dp, "synthetic disc anchors were not collected"
    assert "Dec" not in dp, "a mate on chrUn_decoy must not become a disc record"
    # The two contig orders, which this fixture makes opposite: the site file names chr1 first,
    # the header lists chr2 first, and each mate sits on its anchor's own contig. Note the FASTQ
    # is GROUPED by contig (the legacy gather loops over an insertion-ordered {contig: [pos]}
    # dict), so site-file order means order of first appearance, not the literal line order.
    fq_contigs = [h.split("~")[0][1:] for h in heads]
    assert _dedup_runs(fq_contigs) == ["chr1", "chr2"], \
        f"clip FASTQ is not in site-file contig order: {_dedup_runs(fq_contigs)}"
    fa_ins_contigs = [n.split("~")[-3] for n in open(cpp_fa).read().splitlines()
                      if n.startswith(">")]
    assert _dedup_runs(fa_ins_contigs) == ["chr2", "chr1"], \
        f"disc FASTA is not in BAM-header contig order: {_dedup_runs(fa_ins_contigs)}"


def _dedup_runs(seq):
    out = []
    for x in seq:
        if not out or out[-1] != x:
            out.append(x)
    return out


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


# --------------------------------------------------------------------------------------------
# Stage 2c / plan §4 2b: several families out of ONE decode (--td-manifest).
#
# The saving is CPU, not wall clock -- L1's and SVA's gathers are already concurrent processes --
# so the property worth testing is not speed but that sharing the decode is INVISIBLE. Every test
# below is the same shape: run N entries together, run each one alone, require `cmp` on all three
# products. That is the whole safety argument for widening the fetch plan to the union of the
# entries' windows, and it is checked rather than reasoned about because the two things it rests
# on -- "holding is permissive, emission is exact" and the region ownership rule surviving any
# merge boundary -- are exactly the kind of invariant that reads as obviously true and isn't.
# --------------------------------------------------------------------------------------------
def _run_manifest(scan_bin, bam, reference, entries, outdir, threads=2):
    """entries: [(tag, sites_path, extnd, want_fa)] -> {tag: (clip_fq, disc_pos, disc_fa|None)}."""
    os.makedirs(outdir, exist_ok=True)
    manifest = os.path.join(outdir, "td_manifest.txt")
    out = {}
    with open(manifest, "w") as mf:
        mf.write("# tag sites extnd clip_fastq disc_fa disc_pos\n")
        for tag, sites, extnd, want_fa in entries:
            fq = os.path.join(outdir, f"{tag}.clip.fq")
            fa = os.path.join(outdir, f"{tag}.disc.fa") if want_fa else None
            dp = os.path.join(outdir, f"{tag}.all_disc_pos.disc_pos")
            mf.write(f"{tag}\t{sites}\t{extnd}\t{fq}\t{fa or '-'}\t{dp}\n")
            out[tag] = (fq, dp, fa)
    subprocess.run([scan_bin, "--bam", bam, "--reference", reference,
                    "--td-manifest", manifest, "--threads", str(threads)], check=True)
    return out


def _run_solo(scan_bin, bam, reference, tag, sites, extnd, want_fa, outdir, threads=2):
    """One entry, through the same manifest code path but on its own -- the comparison arm."""
    return _run_manifest(scan_bin, bam, reference, [(tag, sites, extnd, want_fa)],
                         outdir, threads)[tag]


def _assert_entry_matches_solo(label, shared, solo):
    for kind, a, b in zip(("clip FASTQ", "all_disc_pos", "disc FASTA"), solo, shared):
        if a is None or b is None:
            continue
        _assert_same(f"{label} {kind}", a, b)


def test_single_entry_manifest_equals_td_sites(test_bam, reference, xtea2_scan_bin,
                                               site_positions, tmp_path):
    """The manifest must not be a second implementation. --td-sites is a one-entry vector into
    the same function, so this pins that they really are the same path."""
    sites = _write_sites(tmp_path / "sites.txt", SITE_WINDOW[0], site_positions)
    flag_fq, flag_dp, flag_fa = _run_native(xtea2_scan_bin, test_bam, reference, sites,
                                            str(tmp_path / "flags"))
    man = _run_solo(xtea2_scan_bin, test_bam, reference, "L1", sites, EXTND, True,
                    str(tmp_path / "man"))
    assert os.path.getsize(flag_fq) > 0, "the gather collected nothing — bad test window"
    _assert_entry_matches_solo("--td-sites vs one-entry manifest",
                               man, (flag_fq, flag_dp, flag_fa))


def test_two_entries_share_one_decode_and_move_nothing(test_bam, reference, xtea2_scan_bin,
                                                       site_positions, tmp_path):
    """Overlapping-but-not-equal site sets, which is the real shape: on WGS 53,945 of ~54K
    transduction sites are shared by all three families, and the rest are not."""
    a_pos = site_positions[:len(site_positions) * 3 // 4]
    b_pos = site_positions[len(site_positions) // 4:]
    assert set(a_pos) & set(b_pos) and set(a_pos) ^ set(b_pos), \
        "this test is only meaningful with partially overlapping site sets"
    # Reversed on purpose for B: emission order is the ENTRY's own file order, and a shared
    # decode must not leak the union's order (or A's) into it.
    a_sites = _write_sites(tmp_path / "a.txt", SITE_WINDOW[0], a_pos)
    b_sites = _write_sites(tmp_path / "b.txt", SITE_WINDOW[0], list(reversed(b_pos)))

    shared = _run_manifest(xtea2_scan_bin, test_bam, reference,
                           [("L1", a_sites, EXTND, True), ("SVA", b_sites, EXTND, True)],
                           str(tmp_path / "shared"))
    solo_a = _run_solo(xtea2_scan_bin, test_bam, reference, "L1", a_sites, EXTND, True,
                       str(tmp_path / "solo_a"))
    solo_b = _run_solo(xtea2_scan_bin, test_bam, reference, "SVA", b_sites, EXTND, True,
                       str(tmp_path / "solo_b"))
    # Every product non-empty, or `cmp` on two empty files would pass and prove nothing. The
    # FASTA is small here (the chr1 slice holds few of the mates) but it must not be zero.
    for arm, name in ((solo_a, "L1"), (solo_b, "SVA")):
        for path, kind in zip(arm, ("clip FASTQ", "all_disc_pos", "disc FASTA")):
            assert os.path.getsize(path) > 0, \
                f"{name} {kind} is empty — bad test window, the comparison would be vacuous"
    _assert_entry_matches_solo("entry L1", shared["L1"], solo_a)
    _assert_entry_matches_solo("entry SVA", shared["SVA"], solo_b)
    # ...and the two entries are genuinely different outputs, or the test above proves nothing.
    assert open(shared["L1"][0], "rb").read() != open(shared["SVA"][0], "rb").read()


def test_entries_may_have_different_extnd(test_bam, reference, xtea2_scan_bin,
                                          site_positions, tmp_path):
    """The union fetch plan is built at the LARGEST extnd, so the narrower entry is fanned out of
    a collection holding reads it must not emit. If the hold test ever stopped being permissive
    -- or the fan stopped re-filtering exactly -- this is the test that fails."""
    sites = _write_sites(tmp_path / "sites.txt", SITE_WINDOW[0], site_positions)
    narrow, wide = 120, EXTND * 3
    shared = _run_manifest(xtea2_scan_bin, test_bam, reference,
                           [("narrow", sites, narrow, True), ("wide", sites, wide, True)],
                           str(tmp_path / "shared"))
    solo_n = _run_solo(xtea2_scan_bin, test_bam, reference, "narrow", sites, narrow, True,
                       str(tmp_path / "solo_n"))
    solo_w = _run_solo(xtea2_scan_bin, test_bam, reference, "wide", sites, wide, True,
                       str(tmp_path / "solo_w"))
    _assert_entry_matches_solo("narrow entry", shared["narrow"], solo_n)
    _assert_entry_matches_solo("wide entry", shared["wide"], solo_w)
    assert os.path.getsize(shared["wide"][0]) > os.path.getsize(shared["narrow"][0]), \
        "the wide window should collect strictly more clips — the extnd is not reaching the fan"


def test_multi_entry_is_thread_count_invariant(test_bam, reference, xtea2_scan_bin,
                                               site_positions, tmp_path):
    """Both fan-outs stripe over the fetch plan and concatenate part files in stripe order; with
    N entries the mate gather writes N part files per stripe, which is the new thing to get
    wrong."""
    a_sites = _write_sites(tmp_path / "a.txt", SITE_WINDOW[0], site_positions[::2])
    b_sites = _write_sites(tmp_path / "b.txt", SITE_WINDOW[0], site_positions[1::2])
    ents = [("L1", a_sites, EXTND, True), ("SVA", b_sites, EXTND, True)]
    one = _run_manifest(xtea2_scan_bin, test_bam, reference, ents, str(tmp_path / "t1"), threads=1)
    four = _run_manifest(xtea2_scan_bin, test_bam, reference, ents, str(tmp_path / "t4"), threads=4)
    for tag in ("L1", "SVA"):
        _assert_entry_matches_solo(f"{tag} threads 1 vs 4", four[tag], one[tag])


def test_multi_entry_multi_contig_matches_solo(synthetic_bam, xtea2_scan_bin, tmp_path):
    """The synthetic BAM's two contigs are in OPPOSITE order in the site file and the BAM header,
    which is what separates 'emission follows the entry's site file' from 'emission follows the
    header'. Here each entry gets one contig, so a shared decode that leaked the union's contig
    order would show up immediately."""
    both = str(tmp_path / "both.txt")
    with open(both, "w") as f:
        for chrom, p in SYN_SITES:
            f.write(f"{chrom}\t{p}\n")
    per = {}
    for chrom in ("chr1", "chr2"):
        per[chrom] = str(tmp_path / f"{chrom}.txt")
        with open(per[chrom], "w") as f:
            for c, p in SYN_SITES:
                if c == chrom:
                    f.write(f"{c}\t{p}\n")

    shared = _run_manifest(xtea2_scan_bin, synthetic_bam, "",
                           [("all", both, EXTND, True),
                            ("c1", per["chr1"], EXTND, True),
                            ("c2", per["chr2"], EXTND, True)],
                           str(tmp_path / "shared"))
    for tag, sites in (("all", both), ("c1", per["chr1"]), ("c2", per["chr2"])):
        solo = _run_solo(xtea2_scan_bin, synthetic_bam, "", tag, sites, EXTND, True,
                         str(tmp_path / f"solo_{tag}"))
        _assert_entry_matches_solo(f"entry {tag}", shared[tag], solo)
    # The single-contig entries partition the combined one; if the union decode had cross-fed
    # them, this is where an extra record would show up.
    n_all = sum(1 for ln in open(shared["all"][0]) if ln.startswith("@"))
    n_c1 = sum(1 for ln in open(shared["c1"][0]) if ln.startswith("@"))
    n_c2 = sum(1 for ln in open(shared["c2"][0]) if ln.startswith("@"))
    assert n_all == n_c1 + n_c2, f"{n_all} != {n_c1} + {n_c2}"


def test_manifest_rejects_a_shared_output_path(synthetic_bam, xtea2_scan_bin, tmp_path):
    """Two entries writing one file would interleave into it and neither would be readable.
    Caught up front, because the resulting file is undebuggable."""
    sites = str(tmp_path / "sites.txt")
    with open(sites, "w") as f:
        for chrom, p in SYN_SITES:
            f.write(f"{chrom}\t{p}\n")
    manifest = str(tmp_path / "m.txt")
    fq = str(tmp_path / "shared.fq")
    with open(manifest, "w") as mf:
        mf.write(f"a\t{sites}\t{EXTND}\t{fq}\t-\t-\n")
        mf.write(f"b\t{sites}\t{EXTND}\t{fq}\t-\t-\n")
    r = subprocess.run([xtea2_scan_bin, "--bam", synthetic_bam, "--td-manifest", manifest],
                       capture_output=True, text=True)
    assert r.returncode != 0, "two entries sharing an output path must not be accepted"
    assert "same output path" in r.stderr


def test_manifest_rejects_a_short_line(synthetic_bam, xtea2_scan_bin, tmp_path):
    """The cns manifest silently skips a malformed line. Here a dropped field means a family
    quietly does not run, which on a WGS gate reads as 'that family found nothing'."""
    manifest = str(tmp_path / "m.txt")
    with open(manifest, "w") as mf:
        mf.write("onlytag\tsites.txt\n")
    r = subprocess.run([xtea2_scan_bin, "--bam", synthetic_bam, "--td-manifest", manifest],
                       capture_output=True, text=True)
    assert r.returncode != 0 and "line 1" in r.stderr
