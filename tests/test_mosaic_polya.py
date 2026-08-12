"""M.0: the polyA additions to xtea-scan must not move a single germline number.

WHY THIS FILE EXISTS AT ALL
---------------------------
`--mosaic-polya` adds two counters to the realign pass that germline runs have been producing,
frozen and compared against for weeks. The feature is worth very little; NOT disturbing those
counts is worth a great deal, because the pipeline's strongest property is that three unpinned WGS
replicates produced byte-identical call sets (array 50062602). So most of this file is a
regression guard on the OLD behaviour, and only the last group is about the new feature.

THE FAILURE MODE, AND WHY IT TOOK TWO ATTEMPTS TO CATCH
-------------------------------------------------------
The polyA window is not the whole clip: master checks only the ~20 bp END that abuts the insertion
point, and WHICH end that is flips with both the clip side and the realignment strand
(clip_read.py:942-966, ported in realign.hpp:clip_has_polya). Get that backwards and the tally
reads the genomic flank instead of the polyA tract -- no crash, no error, just a quieter number
that looks like "this sample has few polyA clips".

The first version of this file tried to catch that with a real-data assertion that the tally is
non-zero. **It did not work.** Inverting the rule was MEASURED to leave every real-data test here
green: the aggregate moved 2370 -> 2287 polyA clips, 3.5%, because `clip_has_polya` short-circuits
to the whole sequence whenever the clip is <= CK_POLYA_SEQ_MAX (20 bp) and most real clips are.
A real-data test simply cannot see this bug.

What kills it is `test_window_rule_picks_the_breakpoint_end`: four synthetic 40 bp clips -- long
enough that the two ends differ -- with a homopolymer planted at one known physical end, both
strands. All four fail under the inverted rule. Recorded because the obvious test was the wrong
one, and the next person to touch clip_has_polya needs to know that.

The real-data tests still earn their place: they guard the germline counts, which is the property
that actually matters, and they exercise the shipping bwa-driven path rather than a fixture.
"""
import os
import subprocess

import pytest

REGION = "chr1:23000000-24000000"   # inside MEI_RICH_REGION, small enough for a fast bwa pass
CNS = ("LINE1.fa", "ALU.fa", "SVA.fa")


def _cns_paths(rep_lib_dir):
    out = []
    for name in CNS:
        p = os.path.join(rep_lib_dir, "consensus", name)
        if not os.path.exists(p) or not os.path.exists(p + ".bwt"):
            pytest.skip(f"consensus or its bwa index missing: {p}")
        out.append(p)
    return out


def _run_scan(binary, bam, reference, rep_lib_dir, tmp_path, tag, extra=()):
    l1, alu, sva = _cns_paths(rep_lib_dir)
    bwa = os.environ.get("XTEA2_BWA", "/home/cos689/software/bin/bwa")
    if not os.path.exists(bwa):
        pytest.skip("bwa not found; set XTEA2_BWA")
    out = tmp_path / f"realign.{tag}.txt"
    cmd = [binary, "--bam", bam, "--reference", reference, "--region", REGION,
           "--clip-fastq", str(tmp_path / f"clip.{tag}.fq"),
           "--realign-out", str(out), "--l1-cns", l1, "--alu-cns", alu, "--sva-cns", sva,
           "--bwa", bwa, "--bwa-threads", "2", "--threads", "2", *extra]
    r = subprocess.run(cmd, capture_output=True, text=True)
    assert r.returncode == 0, f"xtea_scan failed:\n{r.stderr[-3000:]}"
    assert out.exists(), "no realign-out produced"
    return [ln.split("\t") for ln in out.read_text().splitlines() if ln]


@pytest.fixture(scope="module")
def plain(xtea2_scan_bin, test_bam, reference, rep_lib_dir, tmp_path_factory):
    return _run_scan(xtea2_scan_bin, test_bam, reference, rep_lib_dir,
                     tmp_path_factory.mktemp("plain"), "plain")


@pytest.fixture(scope="module")
def withpa(xtea2_scan_bin, test_bam, reference, rep_lib_dir, tmp_path_factory):
    return _run_scan(xtea2_scan_bin, test_bam, reference, rep_lib_dir,
                     tmp_path_factory.mktemp("withpa"), "withpa", ("--mosaic-polya",))


def test_plain_output_still_has_exactly_six_columns(plain):
    """The default must stay the 6-column contract scan_adapter.py parses."""
    assert plain, "region produced no realign-out rows; the test region has no clip signal"
    assert {len(r) for r in plain} == {6}


def test_polya_columns_are_appended_not_interleaved(withpa):
    assert withpa, "region produced no realign-out rows"
    assert {len(r) for r in withpa} == {8}


def test_polya_flag_does_not_perturb_the_germline_counts(plain, withpa):
    """THE regression guard: chrom/pos/family/nL/nR/mate_in_rep identical with the flag on.

    Compared as full rows rather than counts, so a reordering or a dropped site fails here too.
    """
    assert [r[:6] for r in withpa] == plain


def test_polya_counts_never_exceed_their_side_counts(withpa):
    """npL/npR count a SUBSET of the same realigned clips nL/nR count."""
    bad = [r for r in withpa if int(r[6]) > int(r[3]) or int(r[7]) > int(r[4])]
    assert not bad, f"polyA count exceeds its side count on {len(bad)} rows, e.g. {bad[:3]}"


def test_polya_tally_is_not_dead(withpa):
    """Cheap smoke check: the feature produces SOMETHING on real MEI-rich sequence.

    NOT sufficient on its own. Flipping the window rule to the wrong end was measured to leave
    this assertion green (2370 -> 2287 polyA clips, a 3.5% move), because `clip_has_polya` uses
    the whole sequence whenever the clip is <= CK_POLYA_SEQ_MAX and most clips are. The window
    rule is pinned by test_window_rule_* below, which is where that mutant dies.
    """
    hits = sum(1 for r in withpa if int(r[6]) or int(r[7]))
    assert hits > 0, "no site has any polyA-bearing clip -- clip_has_polya is reading the wrong end"


def test_clip_polya_filter_only_ever_removes_clips(xtea2_scan_bin, test_bam, reference,
                                                   rep_lib_dir, tmp_path, plain):
    """--clip-polya-filter is a gate: it can drop clips from the realign FASTQ, never add them.

    Checked per (site, family) rather than in aggregate, so a filter that removed clips at one
    site while inventing them at another would still fail.
    """
    filtered = _run_scan(xtea2_scan_bin, test_bam, reference, rep_lib_dir, tmp_path,
                         "filtered", ("--clip-polya-filter",))
    base = {(r[0], r[1], r[2]): (int(r[3]), int(r[4])) for r in plain}
    grew = []
    for r in filtered:
        key = (r[0], r[1], r[2])
        bl, br = base.get(key, (0, 0))
        if int(r[3]) > bl or int(r[4]) > br:
            grew.append((key, (bl, br), (int(r[3]), int(r[4]))))
    assert not grew, f"the filter ADDED clips at {len(grew)} site/family keys, e.g. {grew[:3]}"


# ---------------------------------------------------------------------------------------------
# The window rule, pinned deterministically.
#
# Everything above runs on real data, where the window rule barely shows: `clip_has_polya` uses
# the WHOLE sequence when the clip is <= CK_POLYA_SEQ_MAX (20 bp) and most real clips are, so
# flipping first/last moved the aggregate by only 3.5% and left every assertion above green.
# These tests use clips of 40 bp -- long enough that the two ends genuinely differ -- with a
# homopolymer at one known physical end, so "which end" is the only thing under test.
#
# The physical end that matters is the one abutting the breakpoint: a LEFT clip's LAST bases, a
# RIGHT clip's FIRST bases. The b_rev half of the rule exists because a reverse realignment makes
# bwa store SEQ reverse-complemented, moving that end to the other side of the string. So the
# reverse case is built by reverse-complementing the consensus-matching half, which forces bwa to
# align reverse -- it is not simulated by setting a flag.
# ---------------------------------------------------------------------------------------------
_COMP = str.maketrans("ACGT", "TGCA")


def _rc(s):
    return s.translate(_COMP)[::-1]


def _synthetic_seq(n, seed):
    """Deterministic ACGT with no run of 5 -- so the ONLY polyA signal is the one we plant."""
    import random
    rng = random.Random(seed)
    out = []
    for _ in range(n):
        c = rng.choice("ACGT")
        while len(out) >= 4 and all(x == c for x in out[-4:]):
            c = rng.choice("ACGT")
        out.append(c)
    return "".join(out)


@pytest.fixture(scope="module")
def window_case(xtea2_scan_bin, tmp_path_factory):
    """One BAM whose only clips are four 40 bp cases with a planted 20 bp homopolymer."""
    pysam = pytest.importorskip("pysam")
    bwa = os.environ.get("XTEA2_BWA", "/home/cos689/software/bin/bwa")
    if not os.path.exists(bwa):
        pytest.skip("bwa not found; set XTEA2_BWA")
    d = tmp_path_factory.mktemp("window")

    cons = _synthetic_seq(200, seed=11)
    cons_fa = d / "L1.fa"
    cons_fa.write_text(">L1\n" + cons + "\n")
    r = subprocess.run([bwa, "index", str(cons_fa)], capture_output=True, text=True)
    assert r.returncode == 0, r.stderr[-2000:]

    ref = _synthetic_seq(4000, seed=99)
    ref_fa = d / "ref.fa"
    ref_fa.write_text(">chr1\n" + "\n".join(ref[i:i + 60] for i in range(0, len(ref), 60)) + "\n")
    subprocess.run(["samtools", "faidx", str(ref_fa)], capture_output=True)

    anchor, pa = cons[:20], "A" * 20
    # (name, position, clip as it appears IN THE READ, expected polyA under the correct rule)
    cases = [
        # forward realign (anchor matches consensus forward) -> SEQ as-is -> rule takes LAST 20
        ("fwd_polya_at_breakpoint", 600, anchor + pa, True),
        ("fwd_polya_at_far_end",    900, pa + anchor, False),
        # reverse realign (anchor reverse-complemented) -> bwa stores SEQ revcomp'd -> takes FIRST 20
        ("rev_polya_at_breakpoint", 1200, _rc(anchor) + pa, True),
        ("rev_polya_at_far_end",    1500, pa + _rc(anchor), False),
    ]

    hdr = {"HD": {"VN": "1.6", "SO": "coordinate"},
           "SQ": [{"SN": "chr1", "LN": len(ref)}]}
    bam = d / "clips.bam"
    with pysam.AlignmentFile(str(bam), "wb", header=hdr) as bf:
        for name, pos, clip, _ in cases:
            for i in range(4):          # 4 identical reads per case: distinct names, one site
                a = pysam.AlignedSegment()
                a.query_name = f"{name}_{i}"
                a.reference_id = 0
                a.reference_start = pos
                a.mapping_quality = 60
                a.cigartuples = [(4, len(clip)), (0, 60)]   # left soft clip + 60M
                a.query_sequence = clip + ref[pos:pos + 60]
                a.query_qualities = pysam.qualitystring_to_array("I" * (len(clip) + 60))
                a.flag = 0
                bf.write(a)
    subprocess.run(["samtools", "index", str(bam)], capture_output=True)

    out = d / "realign.txt"
    cmd = [xtea2_scan_bin, "--bam", str(bam), "--reference", str(ref_fa),
           "--region", f"chr1:0-{len(ref)}", "--clip-fastq", str(d / "clip.fq"),
           "--realign-out", str(out), "--l1-cns", str(cons_fa), "--bwa", bwa,
           "--bwa-threads", "1", "--threads", "1", "--mosaic-polya"]
    rr = subprocess.run(cmd, capture_output=True, text=True)
    assert rr.returncode == 0, f"scan failed:\n{rr.stderr[-3000:]}"
    rows = {}
    for ln in out.read_text().splitlines():
        f = ln.split("\t")
        if len(f) == 8:
            rows[int(f[1])] = {"nL": int(f[3]), "npL": int(f[6])}
    return cases, rows


@pytest.mark.parametrize("idx", range(4))
def test_window_rule_picks_the_breakpoint_end(window_case, idx):
    """A 40 bp clip is polyA-positive iff the homopolymer sits at the breakpoint-adjacent end.

    This is the assertion that fails when `take_first = (b_left == b_rev)` is inverted; the
    real-data tests above do not.
    """
    cases, rows = window_case
    name, pos, _clip, expect_polya = cases[idx]
    assert pos in rows, f"{name}: site {pos} produced no realigned clip (rows: {sorted(rows)})"
    got = rows[pos]
    assert got["nL"] > 0, f"{name}: clip did not realign at all, so the tally is untestable"
    if expect_polya:
        assert got["npL"] == got["nL"], f"{name}: expected every clip polyA-positive, got {got}"
    else:
        assert got["npL"] == 0, f"{name}: expected no polyA (it is at the far end), got {got}"
