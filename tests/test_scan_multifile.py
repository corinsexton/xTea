"""Multi-file input for the C++ xtea-scan engine: pooling exactness + the @SQ precondition.

Two invariants from the approved multi-BAM plan (`docs/handoff/HANDOFF_2026-07-21.md`,
decisions 1-2), neither previously covered:

  * **Pooling is exact.** Splitting one file's reads across two files and scanning them together
    must reproduce scanning the unsplit original. The split is by **read name**, not by region --
    that is what makes the comparison legitimate rather than circular. The engine's per-file
    read-collection state (disc-needs / anchors / mate-seq) is self-contained *because* mates
    never span physical files; a region split would violate that precondition and the test would
    be measuring its own construction error instead of the pooling logic.

  * **The identical-@SQ precondition fails loud.** Every `(tid,pos)` aggregation key in the
    engine is used unremapped across files, so a differing reference dictionary is silent
    corruption of the output rather than a degraded result. All three ways it can differ --
    sequence count, name/order, sequence length -- must exit non-zero and name the offending
    file, since an unattributed "header mismatch" leaves the user nothing to act on.

The pooling tests need the chr1 test BAM (skip cleanly without it); the @SQ tests are hermetic,
built from ~250-byte pysam BAMs. All skip when XTEA2_SCAN_BIN is unset.
"""
import os
import subprocess
import zlib

import pytest

from tests.utils.clipdisc_io import (
    run_xtea_scan, run_xtea_scan_clip_sites, run_xtea_scan_raw,
)

pysam = pytest.importorskip("pysam")

# A 2 Mb sub-window of the MEI-rich region. Pooling exactness is a property of the aggregation,
# not of window size, so this trades the full 7 Mb for a test that runs in seconds -- each test
# scans the same reads two or three times over.
_WINDOW_BP = 2_000_000


def _pair_bucket(qname, n_buckets=2):
    """Deterministic qname -> bucket, keeping both mates of a pair in the same file.

    crc32, not the builtin hash(): PYTHONHASHSEED randomizes str hashing per process, which would
    make the split -- and so any failure -- unreproducible between runs.
    """
    return zlib.crc32(qname.encode()) % n_buckets


@pytest.fixture(scope="module")
def split_bams(test_bam, mei_rich_region, tmp_path_factory):
    """Partition one window's reads by read name into 2 files, plus the unsplit combined file.

    Returns (combined_path, [partA_path, partB_path], stats).
    """
    chrom, start, _ = mei_rich_region
    end = start + _WINDOW_BP
    d = tmp_path_factory.mktemp("scan_multifile")
    combined = str(d / "combined.bam")
    parts = [str(d / "partA.bam"), str(d / "partB.bam")]

    src = pysam.AlignmentFile(test_bam)
    try:
        writers = [pysam.AlignmentFile(p, "wb", header=src.header) for p in parts]
        comb = pysam.AlignmentFile(combined, "wb", header=src.header)
        try:
            n_total = 0
            qnames_per_part = [set(), set()]
            for rec in src.fetch(chrom, start, end):
                b = _pair_bucket(rec.query_name)
                writers[b].write(rec)
                qnames_per_part[b].add(rec.query_name)
                comb.write(rec)
                n_total += 1
        finally:
            for w in writers:
                w.close()
            comb.close()
    finally:
        src.close()

    for p in [combined] + parts:
        pysam.index(p)

    stats = {"n_total": n_total, "qnames_per_part": qnames_per_part}
    return combined, parts, stats


def test_split_construction_is_sound(split_bams):
    """Guard the test's own fixture: both parts non-empty, mates never split across files.

    If this fails the pooling assertions below are meaningless, so it is asserted separately
    rather than folded into them.
    """
    combined, parts, stats = split_bams
    a, b = stats["qnames_per_part"]
    assert stats["n_total"] > 1000, f"window too sparse to be a real test: {stats['n_total']} reads"
    assert a and b, "one partition came out empty; the split proves nothing"
    assert not (a & b), (
        f"{len(a & b)} read names appear in BOTH files -- mates were split across files, "
        f"violating the precondition the engine's per-file state relies on"
    )
    n_split = sum(pysam.AlignmentFile(p).count(until_eof=True) for p in parts)
    assert n_split == stats["n_total"], f"split lost reads: {n_split} != {stats['n_total']}"


def test_pooled_clip_sites_match_unsplit(split_bams, xtea2_scan_bin, mei_rich_region):
    """Aggregated --clip-sites over "partA partB" must equal the single-file table exactly."""
    chrom, start, _ = mei_rich_region
    end = start + _WINDOW_BP
    combined, parts, _ = split_bams
    one = run_xtea_scan_clip_sites(xtea2_scan_bin, combined, chrom, start, end)
    pooled = run_xtea_scan_clip_sites(xtea2_scan_bin, " ".join(parts), chrom, start, end)
    assert pooled == one, (
        f"pooled clip-sites != unsplit: pooled={len(pooled)} single={len(one)}; "
        f"only_pooled={sorted(set(pooled) - set(one))[:5]} "
        f"only_single={sorted(set(one) - set(pooled))[:5]}; "
        f"count_diffs={[(k, pooled[k], one[k]) for k in sorted(set(pooled) & set(one)) if pooled[k] != one[k]][:5]}"
    )
    assert pooled, "expected a non-empty clip-site table"


def test_pooled_per_read_streams_match_unsplit(split_bams, xtea2_scan_bin, mei_rich_region):
    """The per-read CLIP/DISC streams must be the same *multiset* of records either way.

    Compared as sorted lines rather than as sets, deliberately: a set comparison would call a
    record emitted once per input file equal to one emitted once, so it would pass straight
    through the duplicate-emission bug that pooling is most likely to introduce. Sorting is what
    the comparison tolerates instead -- the interleaving of equal-position records across two
    input files is not part of the contract, but the record count is.
    """
    chrom, start, _ = mei_rich_region
    end = start + _WINDOW_BP
    combined, parts, _ = split_bams
    one = sorted(run_xtea_scan_raw(xtea2_scan_bin, combined, chrom, start, end).splitlines())
    pooled = sorted(run_xtea_scan_raw(xtea2_scan_bin, " ".join(parts), chrom, start, end).splitlines())
    assert len(pooled) == len(one), (
        f"pooled emitted {len(pooled)} records, unsplit {len(one)} "
        f"({'duplicates across files?' if len(pooled) > len(one) else 'records dropped?'})"
    )
    if pooled != one:
        diffs = [(p, o) for p, o in zip(pooled, one) if p != o][:5]
        pytest.fail(f"pooled per-read stream != unsplit; first mismatches (pooled, single): {diffs}")
    assert one, "expected non-empty per-read stream"


def test_pooled_qnames_carry_no_file_index_prefix(split_bams, xtea2_scan_bin, mei_rich_region):
    """Read names must not leak the internal `<file_idx>:<qname>` spill tag.

    Anchor qnames are tagged with their file index while spilled (so the pooled spill can tell
    two files' mates apart) and stripped at emit. A leaked prefix would not change record counts,
    so nothing else here would catch it -- but it would silently break every downstream consumer
    that joins on read name.
    """
    chrom, start, _ = mei_rich_region
    end = start + _WINDOW_BP
    _, parts, _ = split_bams
    clip_keys, disc_keys = run_xtea_scan(xtea2_scan_bin, " ".join(parts), chrom, start, end)
    qnames = [k[0] for k in clip_keys] + [k[0] for k in disc_keys]
    assert qnames, "expected a non-empty stream to check"
    tagged = [q for q in qnames if len(q) > 1 and q[0].isdigit() and q[1] == ":"]
    assert not tagged, f"{len(tagged)} read names still carry the file-index tag, e.g. {tagged[:5]}"


def test_pooling_is_threadsafe_across_files(split_bams, xtea2_scan_bin, mei_rich_region):
    """Multi-file pooling stays deterministic under threading, as the single-file path is."""
    chrom, start, _ = mei_rich_region
    end = start + _WINDOW_BP
    _, parts, _ = split_bams
    joined = " ".join(parts)
    t1 = run_xtea_scan_clip_sites(xtea2_scan_bin, joined, chrom, start, end,
                                  threads=1, tile_size=1_000_000)
    t4 = run_xtea_scan_clip_sites(xtea2_scan_bin, joined, chrom, start, end,
                                  threads=4, tile_size=250_000)
    assert t1 == t4, "multi-file clip-site aggregation is not thread/tile invariant"


# ============================================================================================
# The identical-@SQ precondition (hermetic)
# ============================================================================================
_SQ_BASE = [{"SN": "chrA", "LN": 1000}, {"SN": "chrB", "LN": 2000}]
_SQ_CASES = {
    "base":       _SQ_BASE,
    "base2":      _SQ_BASE,                                                   # identical twin
    "fewer_seqs": [{"SN": "chrA", "LN": 1000}],                               # count differs
    "renamed":    [{"SN": "chrA", "LN": 1000}, {"SN": "chrZ", "LN": 2000}],   # name/order differs
    "relength":   [{"SN": "chrA", "LN": 1000}, {"SN": "chrB", "LN": 9999}],   # length differs
}


@pytest.fixture(scope="module")
def sq_bams(tmp_path_factory):
    """{name: path} for tiny indexed BAMs differing only in their @SQ dictionary."""
    d = tmp_path_factory.mktemp("sq_bams")
    out = {}
    for name, sq in _SQ_CASES.items():
        path = str(d / f"{name}.bam")
        header = {"HD": {"VN": "1.6", "SO": "coordinate"}, "SQ": sq}
        with pysam.AlignmentFile(path, "wb", header=header) as bf:
            for i in range(2):
                a = pysam.AlignedSegment(bf.header)
                a.query_name = f"{name}_read{i}"
                a.query_sequence = "A" * 50
                a.query_qualities = pysam.qualitystring_to_array("I" * 50)
                a.flag = 0
                a.reference_id = 0
                a.reference_start = 100 + i * 100
                a.mapping_quality = 60
                a.cigarstring = "50M"
                bf.write(a)
        pysam.index(path)
        out[name] = path
    return out


def _run_scan(scan_bin, paths, region="chrA:0-1000"):
    """Invoke the binary over N paths in one --bam value; return the CompletedProcess."""
    return subprocess.run(
        [scan_bin, "--bam", " ".join(paths), "--region", region, "--out", "-"],
        capture_output=True, text=True,
    )


def test_identical_headers_are_accepted(sq_bams, xtea2_scan_bin):
    """Positive control: matching @SQ dictionaries must scan without complaint."""
    r = _run_scan(xtea2_scan_bin, [sq_bams["base"], sq_bams["base2"]])
    assert r.returncode == 0, f"identical headers rejected: rc={r.returncode} stderr={r.stderr!r}"
    assert "mismatch" not in r.stderr.lower(), f"unexpected mismatch complaint: {r.stderr!r}"


@pytest.mark.parametrize("case,expect_msg", [
    ("fewer_seqs", "reference sequences"),
    ("renamed",    "header mismatch at reference index"),
    ("relength",   "header mismatch for reference"),
])
def test_sq_mismatch_fails_loud(sq_bams, xtea2_scan_bin, case, expect_msg):
    """Each way @SQ can differ must exit non-zero, explain itself, and name the bad file.

    The filename assertion is the point: "header mismatch" alone, over a pooled input of five
    CRAMs from two centers, does not tell the user which file to fix.
    """
    r = _run_scan(xtea2_scan_bin, [sq_bams["base"], sq_bams[case]])
    assert r.returncode != 0, (
        f"@SQ mismatch ({case}) was ACCEPTED -- (tid,pos) keys would be silently wrong. "
        f"stdout={r.stdout[:200]!r}"
    )
    assert expect_msg in r.stderr, f"expected {expect_msg!r} in stderr, got {r.stderr!r}"
    assert os.path.basename(sq_bams[case]) in r.stderr, (
        f"error names no offending file, so it is not actionable: {r.stderr!r}"
    )


def test_sq_mismatch_detected_before_any_output(sq_bams, xtea2_scan_bin):
    """The guard must fire before scanning, not partway through emitting records.

    Validation is documented as running before any tile is built; if records had already been
    written the caller could consume a truncated stream as if it were complete.
    """
    r = _run_scan(xtea2_scan_bin, [sq_bams["base"], sq_bams["renamed"]])
    assert r.returncode != 0
    assert r.stdout.strip() == "", f"emitted records despite a bad header: {r.stdout[:200]!r}"
