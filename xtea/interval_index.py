"""Static interval index — the replacement for the vendored AVL interval tree.

WHY THIS EXISTS (measured; full evidence in docs/clipdisc_contract_spec.md §4 P1-EVALUATED)
-------------------------------------------------------------------------------------------
`intervaltree.py` was a vendored copy of the PyPI `intervaltree` package: a general AVL
structure supporting dynamic insert/delete, range overlap, merging, chopping and slicing.
xTea used exactly ONE of those operations -- `tree[point]` -- after a build-once phase.
Nothing was ever inserted or deleted after construction, and no range query was ever issued.
A tree buys mutation that nothing uses, and charges 23.6 s (Alu) / 18.0 s (L1) to build.

The tree had itself replaced a sorted-array binary search (`XAnnotation.is_within_repeat_region`,
still present until step 4 deletes it) which is *incorrect under overlap*: it binary-searches on
start positions and returns at the first midpoint hit, so overlapping copies can hide a genuine
hit behind a miss. So the tree was adopted as a correctness fix. The correct-and-fast answer is
neither: sorted starts plus a prefix-max-end scan, which is exact under overlap and builds in
0.22 s (105x faster) while querying 12x faster.

WHAT MAKES THE SCAN CHEAP -- and the assumption to re-check if annotations change
--------------------------------------------------------------------------------
The backward scan below walks from `bisect_right(starts, p)` towards 0 and stops as soon as
`prefix_max_end[i] <= p`. That terminates within a handful of entries only because repeat
copies are short: Alu median 498 / max 1,013 bp, L1 median 1,181 / max 9,304 bp (extension
included). The `_max_len` bound gives a second, independent stopping rule so a single
pathologically long interval cannot force a long scan. If an annotation with genuinely long
intervals is ever introduced, re-measure -- correctness is unaffected, cost is not.

SEMANTICS -- must match the tree exactly
----------------------------------------
Half-open, `begin <= p < end`. That is what `Interval.contains_point` did (the convention was
documented at the old `interval.py:46,48` and deliberately preserved through step 1's
comment stripping precisely so this module could be written against it).

TWO QUERY MODES, AND WHY THEY ARE SEPARATE
------------------------------------------
`contains()` answers "does anything overlap" -- order-independent, so it is immune to the B8
family of bugs by construction. `hits()` returns every overlapping entry in ascending start
order, for the few callers that must then choose one. Handing an *unordered* container to code
that picks one element is the root cause of B8/B8b/B8c; this module never does that.
"""

from bisect import bisect_right


class ChromIntervalIndex:
    """Static index over one chromosome's intervals. Built once, queried read-only."""

    __slots__ = ("starts", "ends", "_prefix_max_end", "_max_len")

    def __init__(self, starts, ends):
        """`starts`/`ends` are parallel sequences; they need not be sorted.

        Intervals are half-open [start, end). Entries with end <= start are dropped, matching
        the old tree, which refused null intervals outright.
        """
        pairs = sorted(
            (s, e) for s, e in zip(starts, ends) if e > s
        )
        self.starts = [s for s, _ in pairs]
        self.ends = [e for _, e in pairs]

        # prefix_max_end[i] = max(ends[0..i]) -- lets the backward scan stop early.
        pmax = []
        running = None
        for e in self.ends:
            if running is None or e > running:
                running = e
            pmax.append(running)
        self._prefix_max_end = pmax
        self._max_len = max((e - s for s, e in pairs), default=0)

    def __len__(self):
        return len(self.starts)

    def _hit_indices(self, p):
        """Indices of every interval containing p, ascending by start."""
        starts = self.starts
        hi = bisect_right(starts, p)  # every i < hi has starts[i] <= p
        if hi == 0:
            return []
        ends = self.ends
        pmax = self._prefix_max_end
        # Nothing starting at or before `floor` can still be open at p.
        floor = bisect_right(starts, p - self._max_len)
        out = []
        i = hi - 1
        while i >= floor and pmax[i] > p:
            if ends[i] > p:
                out.append(i)
            i -= 1
        out.reverse()
        return out

    def contains(self, p):
        """True iff any interval contains p. Order-independent by construction."""
        starts = self.starts
        hi = bisect_right(starts, p)
        if hi == 0:
            return False
        ends = self.ends
        pmax = self._prefix_max_end
        floor = bisect_right(starts, p - self._max_len)
        i = hi - 1
        while i >= floor and pmax[i] > p:
            if ends[i] > p:
                return True
            i -= 1
        return False

    def hits(self, p):
        """Every (start, end) containing p, ascending by start. May be empty."""
        return [(self.starts[i], self.ends[i]) for i in self._hit_indices(p)]

    def hit_starts(self, p):
        """Every start position whose interval contains p, ascending."""
        return [self.starts[i] for i in self._hit_indices(p)]


class GenomeIntervalIndex:
    """A `ChromIntervalIndex` per chromosome, plus the whole-genome query surface."""

    __slots__ = ("_by_chrom",)

    def __init__(self, by_chrom=None):
        self._by_chrom = dict(by_chrom or {})

    def set_chrom(self, chrm, starts, ends):
        self._by_chrom[chrm] = ChromIntervalIndex(starts, ends)

    def __contains__(self, chrm):
        return chrm in self._by_chrom

    def __len__(self):
        return len(self._by_chrom)

    def chroms(self):
        return self._by_chrom.keys()

    def get(self, chrm):
        return self._by_chrom.get(chrm)

    def contains(self, chrm, p):
        idx = self._by_chrom.get(chrm)
        return False if idx is None else idx.contains(p)

    def hit_starts(self, chrm, p):
        idx = self._by_chrom.get(chrm)
        return [] if idx is None else idx.hit_starts(p)
