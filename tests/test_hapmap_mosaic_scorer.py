"""The HapMap mosaic scorer: the VAF derivation, and what does NOT count as a false positive.

Two things here would silently produce a wrong-but-plausible benchmark:

  * The VAF derivation. Only L1 ships a VAF file; Alu and SVA VAFs are derived from the donor
    columns times the mixture fractions. If that derivation is wrong the curve is wrong for
    two of three families and nothing complains -- so it is checked against all 422 shipped
    L1 values, which is the whole licence for using it elsewhere.
  * The false-positive definition. HG005 is 83.5% of the mixture, so its insertions are real
    and GERMLINE at ~42% VAF. Counting them as FPs would turn precision into a measure of how
    much germline signal the run saw, and a caller would look worse the better it worked.
"""
import os
import sys

import pytest

REPO_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(REPO_ROOT, "benchmark"))

import benchmark_hapmap_mosaic as bm  # noqa: E402

TRUTH_DIR = bm.DEFAULT_TRUTH_DIR
pytestmark = pytest.mark.skipif(not os.path.isdir(TRUTH_DIR),
                                reason=f"HapMap truth set not present at {TRUTH_DIR}")


class _Call:
    __slots__ = ("chrom", "pos", "family")

    def __init__(self, chrom, pos, family="L1"):
        self.chrom, self.pos, self.family = chrom, pos, family


# --------------------------------------------------------------------------------------
# the derivation
# --------------------------------------------------------------------------------------

def test_vaf_derivation_reproduces_every_shipped_l1_value():
    """422/422 exact. This is what licenses deriving Alu and SVA VAFs, which ship none."""
    assert bm.self_check(TRUTH_DIR) == 0


def test_mixture_is_parsed_from_the_file_not_the_hardcoded_fallback():
    parsed = bm.load_mixture(TRUTH_DIR)
    assert parsed["HG005"] == pytest.approx(0.835)
    assert parsed["HG00438"] == pytest.approx(0.005)
    assert sum(parsed.values()) == pytest.approx(1.0, abs=1e-9), \
        f"mixture fractions must sum to 1: {parsed}"


def test_vaf_is_half_the_carried_fraction():
    """Heterozygous carriers: a donor at 10% contributes 5% VAF, not 10%."""
    mix = {"A": 0.10, "B": 0.02}
    assert bm.derive_vaf([".", "x"], ["A", "B"], mix) == pytest.approx(0.01)
    assert bm.derive_vaf(["x", "."], ["A", "B"], mix) == pytest.approx(0.05)
    assert bm.derive_vaf(["x", "x"], ["A", "B"], mix) == pytest.approx(0.06)
    assert bm.derive_vaf([".", "."], ["A", "B"], mix) == 0.0


@pytest.mark.parametrize("fam,n", [("L1", 32), ("Alu", 196), ("SVA", 12)])
def test_chr1_truth_counts(fam, n):
    """Pinned so a truth-set swap that changes the denominator is visible, not silent."""
    assert len(bm.load_truth(TRUTH_DIR, fam, region="chr1")) == n


def test_region_filter_actually_restricts():
    """A chr1-only run scored against genome-wide truth reports a meaningless sensitivity."""
    all_l1 = bm.load_truth(TRUTH_DIR, "L1")
    chr1 = bm.load_truth(TRUTH_DIR, "L1", region="chr1")
    assert len(chr1) < len(all_l1)
    assert {t.chrom for t in chr1} == {"chr1"}


def test_every_truth_locus_lands_in_a_bin():
    """A locus whose VAF falls outside every bin edge would vanish from the curve."""
    for fam in ("L1", "Alu", "SVA"):
        for t in bm.load_truth(TRUTH_DIR, fam):
            assert any(lo <= t.vaf < hi for lo, hi, _ in bm.VAF_BINS), \
                f"{fam} {t.chrom}:{t.pos} vaf={t.vaf} falls outside VAF_BINS"


def test_bin_edges_do_not_split_a_mixture_rung():
    """The ladder's rungs are discrete. An edge landing ON one splits a single biological
    stratum across two rows and makes both look half-sized."""
    edges = {lo for lo, _, _ in bm.VAF_BINS} | {hi for _, hi, _ in bm.VAF_BINS}
    rungs = {round(t.vaf, 6) for fam in ("L1", "Alu", "SVA")
             for t in bm.load_truth(TRUTH_DIR, fam)}
    assert not (edges & rungs), f"bin edge coincides with a mixture rung: {edges & rungs}"


# --------------------------------------------------------------------------------------
# what counts as a false positive
# --------------------------------------------------------------------------------------

def test_hg005_calls_are_not_false_positives():
    """HG005 is 83.5% of the mix: its insertions are real germline at ~42% VAF. A mosaic
    caller finding them is working, and must not be penalised for it."""
    hg005, masked = bm.load_exclusions(TRUTH_DIR)
    chrom = next(c for c in hg005 if hg005[c])
    pos = hg005[chrom][0]
    _, overall, unm = bm.score([_Call(chrom, pos)], [], hg005, masked, 200)
    assert overall["fp"] == 0
    assert unm["germline"] == 1


def test_masked_calls_are_neither_credited_nor_penalised():
    """masked_regions.txt is tier3-not-tier1: lower-confidence truth the set does not claim."""
    hg005, masked = bm.load_exclusions(TRUTH_DIR)
    chrom = next(c for c in masked if masked[c])
    pos = masked[chrom][0]
    _, overall, unm = bm.score([_Call(chrom, pos)], [], hg005, masked, 200)
    assert overall["tp"] == 0 and overall["fp"] == 0
    assert unm["ambiguous"] == 1


def test_a_call_in_open_sequence_is_a_false_positive():
    """The exclusions must not swallow everything -- otherwise precision is always 100%."""
    hg005, masked = bm.load_exclusions(TRUTH_DIR)
    _, overall, unm = bm.score([_Call("chr21", 9_123_456)], [], hg005, masked, 200)
    assert overall["fp"] == 1 and unm["fp"] == 1


def test_a_matched_truth_locus_is_a_true_positive_and_bins_by_its_vaf():
    hg005, masked = bm.load_exclusions(TRUTH_DIR)
    truth = bm.load_truth(TRUTH_DIR, "L1", region="chr1")
    t = truth[0]
    bins, overall, _ = bm.score([_Call(t.chrom, t.pos + 50)], truth, hg005, masked, 200)
    assert overall["tp"] == 1 and overall["fn"] == len(truth) - 1
    hit = [b for b in bins if b["tp"] == 1]
    assert len(hit) == 1, "the TP landed in zero or several VAF bins"
    lo, hi, _lab = next((lo, hi, lab) for lo, hi, lab in bm.VAF_BINS if lab == hit[0]["label"])
    assert lo <= t.vaf < hi


def test_a_call_outside_the_window_does_not_match():
    hg005, masked = bm.load_exclusions(TRUTH_DIR)
    truth = bm.load_truth(TRUTH_DIR, "L1", region="chr1")
    t = truth[0]
    _, overall, _ = bm.score([_Call(t.chrom, t.pos + 5000)], truth, hg005, masked, 200)
    assert overall["tp"] == 0
