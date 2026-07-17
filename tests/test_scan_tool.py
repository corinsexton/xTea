"""Tool-level tests for the C++ xtea-scan Stage-1 collector.

Covers the capabilities beyond raw read detection:
  * deterministic output regardless of --threads (bp-tile parallelism);
  * multithreaded scan still matches the oracle read-for-read;
  * the aggregated --clip-sites table equals the oracle's per-position L/R aggregation.

All skip cleanly when XTEA2_SCAN_BIN is unset.
"""
from tests.utils.clipdisc_io import (
    detect_clip_reads, clip_read_keys, compare_read_sets,
    run_xtea_scan, run_xtea_scan_raw, run_xtea_scan_clip_sites,
    aggregate_clip_sites_from_oracle,
)


def test_scan_deterministic_across_threads(test_bam, xtea2_scan_bin, mei_rich_region):
    """threads=1 and threads=4 (small tiles → many tiles) must produce byte-identical output."""
    chrom, start, end = mei_rich_region
    out1 = run_xtea_scan_raw(xtea2_scan_bin, test_bam, chrom, start, end, threads=1, tile_size=1_000_000)
    out4 = run_xtea_scan_raw(xtea2_scan_bin, test_bam, chrom, start, end, threads=4, tile_size=1_000_000)
    assert out1 == out4, "xtea-scan output differs between 1 and 4 threads (non-deterministic tiling)"
    assert out1.strip(), "expected non-empty output"


def test_scan_multithread_matches_oracle(test_bam, xtea2_scan_bin, mei_rich_region):
    """Parallel scan (4 threads) must still match the oracle clip set read-for-read."""
    chrom, start, end = mei_rich_region
    scan_clip, _ = run_xtea_scan(xtea2_scan_bin, test_bam, chrom, start, end,
                                 threads=4, tile_size=1_000_000)
    oracle = clip_read_keys(detect_clip_reads(test_bam, chrom, start, end))
    diff = compare_read_sets(oracle, scan_clip)
    assert diff.equal, f"multithreaded scan != oracle: {diff.summary()}"


def test_clip_sites_aggregation_matches_oracle(test_bam, xtea2_scan_bin, mei_rich_region):
    """The aggregated --clip-sites table equals the oracle's per-position L/R counts, exactly."""
    chrom, start, end = mei_rich_region
    tool_sites = run_xtea_scan_clip_sites(xtea2_scan_bin, test_bam, chrom, start, end)
    oracle_sites = aggregate_clip_sites_from_oracle(detect_clip_reads(test_bam, chrom, start, end))
    assert tool_sites == oracle_sites, (
        f"clip-site aggregation mismatch: tool={len(tool_sites)} oracle={len(oracle_sites)}; "
        f"only_tool={list(set(tool_sites)-set(oracle_sites))[:5]} "
        f"only_oracle={list(set(oracle_sites)-set(tool_sites))[:5]}"
    )


def test_clip_sites_boundary_tiles_agree(test_bam, xtea2_scan_bin, mei_rich_region):
    """Aggregation is tile-size invariant: tiny tiles (many boundaries) == one big tile."""
    chrom, start, end = mei_rich_region
    small = run_xtea_scan_clip_sites(xtea2_scan_bin, test_bam, chrom, start, end,
                                     threads=4, tile_size=250_000)
    big = run_xtea_scan_clip_sites(xtea2_scan_bin, test_bam, chrom, start, end,
                                   threads=1, tile_size=1_000_000_000)
    assert small == big, "clip-site aggregation depends on tile size (boundary ownership bug?)"
