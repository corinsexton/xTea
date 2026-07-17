"""Read-level disc-step tests: "the same discordant reads are found".

Mirrors test_clip_reads.py for the discordant step, using the `detect_disc_reads` oracle
(spec §2.2). Discordance is classified from each read's own record, so it is tile-safe and
directly comparable between BAM and CRAM.
"""
import pytest

from tests.utils.clipdisc_io import (
    detect_disc_reads, disc_read_keys, compare_read_sets, parse_disc_sites, run_xtea_scan,
)


def test_oracle_finds_disc_reads_in_mei_rich_region(test_bam, mei_rich_region):
    chrom, start, end = mei_rich_region
    reads = detect_disc_reads(test_bam, chrom, start, end)
    assert len(reads) > 50, "expected discordant reads in a 7 Mb MEI-rich window"


def test_disc_reads_bam_cram_equivalence(test_bam, test_cram, reference, mei_rich_region):
    """Same discordant reads whether the source is BAM or CRAM (dual-input requirement)."""
    chrom, start, end = mei_rich_region
    from_bam = disc_read_keys(detect_disc_reads(test_bam, chrom, start, end))
    from_cram = disc_read_keys(
        detect_disc_reads(test_cram, chrom, start, end, reference=reference)
    )
    diff = compare_read_sets(from_bam, from_cram)
    assert diff.equal, (
        f"BAM vs CRAM disc read-sets differ: {diff.summary()}\n"
        f"missing(sample)={list(diff.only_expected)[:5]}\n"
        f"spurious(sample)={list(diff.only_actual)[:5]}"
    )


def test_disc_mapq_gate_is_enforced(test_bam, mei_rich_region):
    """Raising the MAPQ floor can only shrink the discordant set (monotonic gate sanity)."""
    chrom, start, end = mei_rich_region
    lenient = disc_read_keys(detect_disc_reads(test_bam, chrom, start, end, min_mapq=20))
    strict = disc_read_keys(detect_disc_reads(test_bam, chrom, start, end, min_mapq=40))
    assert strict <= lenient
    assert len(strict) <= len(lenient)


def test_golden_disc_sites_parse(golden_disc_tables):
    if "ALU" not in golden_disc_tables:
        pytest.skip("no golden Alu disc table")
    sites = parse_disc_sites(golden_disc_tables["ALU"])
    assert len(sites) > 0
    # counts are non-negative ints
    for s in list(sites.values())[:50]:
        assert s.n_left >= 0 and s.n_right >= 0


def test_xtea2_scan_matches_oracle(test_bam, xtea2_scan_bin, mei_rich_region):
    """Read-for-read: the C++ xtea-scan disc set must equal the oracle set."""
    chrom, start, end = mei_rich_region
    _, scan_disc = run_xtea_scan(xtea2_scan_bin, test_bam, chrom, start, end)
    oracle = disc_read_keys(detect_disc_reads(test_bam, chrom, start, end))
    diff = compare_read_sets(oracle, scan_disc)
    assert diff.equal, (
        f"xtea-scan disc set != oracle: {diff.summary()}\n"
        f"missing(sample)={list(diff.only_expected)[:5]}\n"
        f"spurious(sample)={list(diff.only_actual)[:5]}"
    )


def test_xtea2_scan_matches_oracle_cram(test_cram, reference, xtea2_scan_bin, mei_rich_region):
    """Read-for-read on CRAM input: xtea-scan disc set == oracle set."""
    chrom, start, end = mei_rich_region
    _, scan_disc = run_xtea_scan(xtea2_scan_bin, test_cram, chrom, start, end, reference=reference)
    oracle = disc_read_keys(detect_disc_reads(test_cram, chrom, start, end, reference=reference))
    diff = compare_read_sets(oracle, scan_disc)
    assert diff.equal, f"xtea-scan (CRAM) disc set != oracle: {diff.summary()}"
