"""Read-level clip-step tests: "the same clipped reads are found".

The invariant (per user): as the clip step is refactored (Python → C xtea-scan), the *set of
clipped reads found* must not change. We express this against an independent pysam oracle
(`detect_clip_reads`, which encodes the spec) rather than master's opaque intermediate files.

Active now:
  * the oracle finds clipped reads in a MEI-rich region;
  * BAM and CRAM inputs yield the *identical* clip read-set (the dual-input requirement);
  * oracle-detected clip positions cover the golden Alu candidate sites.

Activates when built:
  * `test_xtea2_scan_matches_oracle` — the C tool's clip read-set == oracle read-set.
"""
import pytest

from tests.utils.clipdisc_io import (
    detect_clip_reads, clip_read_keys, compare_read_sets, parse_clip_sites, run_xtea_scan,
)


def test_oracle_finds_clips_in_mei_rich_region(test_bam, mei_rich_region):
    chrom, start, end = mei_rich_region
    reads = detect_clip_reads(test_bam, chrom, start, end)
    assert len(reads) > 100, "expected many clipped reads in a 7 Mb MEI-rich window"
    # both sides represented
    sides = {r.side for r in reads}
    assert sides == {"L", "R"}


def test_clip_reads_bam_cram_equivalence(test_bam, test_cram, reference, mei_rich_region):
    """Same clipped reads whether the source is BAM or CRAM (dual-input requirement)."""
    chrom, start, end = mei_rich_region
    from_bam = clip_read_keys(detect_clip_reads(test_bam, chrom, start, end))
    from_cram = clip_read_keys(
        detect_clip_reads(test_cram, chrom, start, end, reference=reference)
    )
    diff = compare_read_sets(from_bam, from_cram)
    assert diff.equal, (
        f"BAM vs CRAM clip read-sets differ: {diff.summary()}\n"
        f"missing(sample)={list(diff.only_expected)[:5]}\n"
        f"spurious(sample)={list(diff.only_actual)[:5]}"
    )


def test_oracle_clip_positions_cover_golden_alu_sites(test_bam, golden_clip_tables, mei_rich_region):
    """Oracle raw clip positions should cover master's golden Alu clip sites.

    This is a loose sanity check that the oracle finds real signal — NOT a strict invariant.
    Oracle positions are *raw* breakpoints (no cutoff/merge), while master's candidate positions
    are cutoff-filtered and slightly repositioned (nearby-region clustering + the L/R TSD offset),
    so exact coordinates differ by up to a few tens of bp for low-support sites. We therefore use
    our concordance tolerance (±25 bp): empirically ~94% covered at ±25 (97% at ±50). The strict,
    read-for-read invariant lives in `test_xtea2_scan_matches_oracle`.
    """
    import bisect
    if "ALU" not in golden_clip_tables:
        pytest.skip("no golden Alu clip table")
    chrom, start, end = mei_rich_region
    golden = parse_clip_sites(golden_clip_tables["ALU"])
    golden_in_region = sorted(pos for (c, pos) in golden if c == chrom and start <= pos < end)
    if not golden_in_region:
        pytest.skip("no golden Alu clip sites in the MEI-rich region")
    oracle_pos = sorted({r.clip_pos for r in detect_clip_reads(test_bam, chrom, start, end)})
    # Oracle is a raw superset: it must find far more candidate positions than the filtered golden.
    assert len(oracle_pos) > len(golden_in_region)
    window = 25
    covered = 0
    for g in golden_in_region:
        i = bisect.bisect_left(oracle_pos, g)
        nearest = min((abs(oracle_pos[j] - g) for j in (i - 1, i) if 0 <= j < len(oracle_pos)),
                      default=10**9)
        if nearest <= window:
            covered += 1
    frac = covered / len(golden_in_region)
    assert frac >= 0.90, f"only {frac:.2f} of golden Alu clip sites covered by oracle within ±{window}bp"


def test_xtea2_scan_matches_oracle(test_bam, xtea2_scan_bin, mei_rich_region):
    """Read-for-read: the C++ xtea-scan clip set must equal the oracle set."""
    chrom, start, end = mei_rich_region
    scan_clip, _ = run_xtea_scan(xtea2_scan_bin, test_bam, chrom, start, end)
    oracle = clip_read_keys(detect_clip_reads(test_bam, chrom, start, end))
    diff = compare_read_sets(oracle, scan_clip)
    assert diff.equal, (
        f"xtea-scan clip set != oracle: {diff.summary()}\n"
        f"missing(sample)={list(diff.only_expected)[:5]}\n"
        f"spurious(sample)={list(diff.only_actual)[:5]}"
    )


def test_xtea2_scan_matches_oracle_cram(test_cram, reference, xtea2_scan_bin, mei_rich_region):
    """Read-for-read on CRAM input: xtea-scan clip set == oracle set (dual-input via the C++ tool)."""
    chrom, start, end = mei_rich_region
    scan_clip, _ = run_xtea_scan(xtea2_scan_bin, test_cram, chrom, start, end, reference=reference)
    oracle = clip_read_keys(detect_clip_reads(test_cram, chrom, start, end, reference=reference))
    diff = compare_read_sets(oracle, scan_clip)
    assert diff.equal, f"xtea-scan (CRAM) clip set != oracle: {diff.summary()}"
