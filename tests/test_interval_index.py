"""Unit tests for the static interval index that replaced the vendored AVL interval tree.

**Why these tests do not diff against the old tree.** The one-off differential that validated the
replacement *did* run against the vendored `intervaltree` — over all 2,188,624 interval starts in
the real hg38 Alu and L1 annotations, at the real postF parameters — and it came back with zero
boolean differences and zero hit-set differences. But the whole point of step 2 is to delete that
tree, so a permanent test depending on it would block its own premise. The differential's results
are recorded in `docs/clipdisc_contract_spec.md` §4 P1-EVALUATED; what lives here is the
implementation-independent oracle: **brute force**.

The properties pinned below are exactly the ones the replacement had to preserve, and the ones a
future edit is most likely to break:

1. **Half-open semantics**, `begin <= p < end`. This was documented at the old `interval.py:46,48`
   and deliberately preserved through step 1's comment stripping so this module could be written
   against it. Off-by-one here silently moves calls.
2. **Every overlapping interval is found.** This is what the AVL tree was originally adopted to
   fix — the hand-rolled binary search it replaced (`XAnnotation.is_within_repeat_region`)
   returns at the first midpoint hit and can walk past a genuine hit when copies overlap. The
   sorted-array scan must not reintroduce that bug, which is why the randomized test below uses
   deliberately overlap-heavy inputs.
3. **Ascending hit order.** B8/B8b came from handing an unordered `set` to code that then picks or
   joins one element. Order is now part of the contract, not an accident.
4. **The boolean agrees with the hit list.** The hot path uses `contains()` and the cold path uses
   `hits()`; if they can disagree, B8's class of bug comes back through the side door.
"""

import random

import pytest

from xtea.interval_index import ChromIntervalIndex, GenomeIntervalIndex


def brute(ivs, p):
    """The oracle: every non-null interval containing p, ascending."""
    return sorted((s, e) for s, e in ivs if e > s and s <= p < e)


class TestHalfOpenSemantics:
    """begin <= p < end — the boundary convention the old tree used."""

    @pytest.fixture
    def idx(self):
        return ChromIntervalIndex([100], [200])

    @pytest.mark.parametrize("p,expected", [
        (99, False),    # just before
        (100, True),    # start is INCLUSIVE
        (150, True),
        (199, True),    # last contained position
        (200, False),   # end is EXCLUSIVE
        (201, False),
    ])
    def test_boundaries(self, idx, p, expected):
        assert idx.contains(p) is expected
        assert bool(idx.hits(p)) is expected

    def test_null_intervals_dropped(self):
        """The old tree refused begin >= end outright; so must we."""
        idx = ChromIntervalIndex([10, 20, 30], [10, 15, 40])  # first null, second inverted
        assert len(idx) == 1
        assert idx.contains(10) is False
        assert idx.contains(20) is False
        assert idx.contains(35) is True


class TestOverlap:
    """The property the AVL tree was adopted to guarantee."""

    def test_finds_every_overlapping_interval(self):
        # Nested, staggered and identical spans all overlapping p=150.
        ivs = [(100, 400), (120, 160), (140, 200), (150, 151), (300, 500)]
        idx = ChromIntervalIndex([s for s, _ in ivs], [e for _, e in ivs])
        assert idx.hits(150) == [(100, 400), (120, 160), (140, 200), (150, 151)]

    def test_hits_are_ascending_by_start(self):
        """Order is contract, not accident — B8/B8b came from unordered results."""
        ivs = [(300, 400), (100, 500), (200, 450)]
        idx = ChromIntervalIndex([s for s, _ in ivs], [e for _, e in ivs])
        assert idx.hits(350) == [(100, 500), (200, 450), (300, 400)]
        assert idx.hit_starts(350) == [100, 200, 300]

    def test_long_interval_does_not_hide_later_short_ones(self):
        """A long early interval keeps prefix-max-end high; the scan must still terminate
        correctly and still report the short intervals near the query."""
        ivs = [(0, 1_000_000), (999_000, 999_010), (999_500, 999_600)]
        idx = ChromIntervalIndex([s for s, _ in ivs], [e for _, e in ivs])
        assert idx.hits(999_550) == [(0, 1_000_000), (999_500, 999_600)]
        assert idx.hits(999_005) == [(0, 1_000_000), (999_000, 999_010)]


class TestAgainstBruteForce:
    """Randomized equivalence — the implementation-independent oracle."""

    def test_randomized(self):
        random.seed(11)
        for _ in range(300):
            ivs = []
            for _ in range(random.randint(0, 40)):
                s = random.randint(0, 200)
                # length 0 exercises null-dropping; 150 makes overlap the common case
                ivs.append((s, s + random.choice([0, 1, 2, 5, 50, 150])))
            idx = ChromIntervalIndex([s for s, _ in ivs], [e for _, e in ivs])
            for p in range(-5, 210):
                want = brute(ivs, p)
                assert idx.hits(p) == want
                assert idx.contains(p) is bool(want)

    def test_empty_index(self):
        idx = ChromIntervalIndex([], [])
        assert len(idx) == 0
        assert idx.contains(0) is False
        assert idx.hits(0) == []


class TestGenomeIndex:
    def test_per_chromosome_isolation(self):
        g = GenomeIntervalIndex()
        g.set_chrom("chr1", [100], [200])
        g.set_chrom("chr2", [500], [600])
        assert g.contains("chr1", 150) is True
        assert g.contains("chr2", 150) is False   # chr1's interval must not leak
        assert g.contains("chr1", 550) is False
        assert g.contains("chr2", 550) is True

    def test_unknown_chromosome_is_a_miss_not_an_error(self):
        """Callers pass whatever chrom name the BAM had; a miss must not raise."""
        g = GenomeIntervalIndex()
        g.set_chrom("chr1", [100], [200])
        assert g.contains("chrUnplaced", 150) is False
        assert g.hit_starts("chrUnplaced", 150) == []
