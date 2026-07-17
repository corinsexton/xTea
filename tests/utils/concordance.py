"""Call-set concordance: the primary reproducibility gate for xTea2.

Two MEI calls are "the same site" when they are on the same chromosome, of the same MEI family,
and within a position window (default ±25 bp). Genotype/quality are NOT required to match — the
gate is on *what is called where*, per the project decision.

`compare_callsets(truth, query)` returns sensitivity (fraction of truth recovered) and PPV
(fraction of query that is real), plus the matched/missed/extra breakdowns for debugging.

Matching is one-to-one and greedy-by-nearest, so a single query call cannot satisfy two truth
calls (and vice versa). Results are deterministic regardless of input order.
"""
from __future__ import annotations

from dataclasses import dataclass, field
from typing import List, Tuple

from .vcf_io import MEICall

DEFAULT_WINDOW = 25


@dataclass
class ConcordanceResult:
    n_truth: int
    n_query: int
    matched: List[Tuple[MEICall, MEICall]] = field(default_factory=list)
    missed: List[MEICall] = field(default_factory=list)   # truth calls with no match = false negatives
    extra: List[MEICall] = field(default_factory=list)     # query calls with no match = false positives
    window: int = DEFAULT_WINDOW

    @property
    def n_matched(self) -> int:
        return len(self.matched)

    @property
    def sensitivity(self) -> float:
        """Recall of the truth set. Undefined-as-1.0 when there is no truth to recover."""
        return self.n_matched / self.n_truth if self.n_truth else 1.0

    @property
    def ppv(self) -> float:
        """Precision of the query set. Undefined-as-1.0 when the query is empty."""
        return self.n_matched / self.n_query if self.n_query else 1.0

    def summary(self) -> str:
        return (
            f"truth={self.n_truth} query={self.n_query} matched={self.n_matched} "
            f"missed={len(self.missed)} extra={len(self.extra)} "
            f"sensitivity={self.sensitivity:.4f} ppv={self.ppv:.4f} (window=±{self.window})"
        )


def compare_callsets(
    truth: List[MEICall],
    query: List[MEICall],
    window: int = DEFAULT_WINDOW,
    match_type: bool = True,
) -> ConcordanceResult:
    """Compare a query call-set against a truth call-set.

    Args:
        truth: the golden/master calls we must recover.
        query: the calls produced by the pipeline under test.
        window: max |POS difference| (bp) for two calls to be the same site.
        match_type: require the MEI family to match (True for the real gate).
    """
    # Bucket query calls by (chrom, type) for efficient nearest-search, keeping original order.
    from collections import defaultdict

    buckets = defaultdict(list)
    for q in query:
        bkey = (q.chrom, q.mei_type if match_type else None)
        buckets[bkey].append(q)

    used = set()  # ids of query calls already matched
    result = ConcordanceResult(n_truth=len(truth), n_query=len(query), window=window)

    # Deterministic order: process truth sorted by (chrom, type, pos).
    for t in sorted(truth, key=lambda c: (c.chrom, c.mei_type, c.pos)):
        bkey = (t.chrom, t.mei_type if match_type else None)
        best = None
        best_d = window + 1
        for q in buckets.get(bkey, ()):
            if id(q) in used:
                continue
            d = abs(q.pos - t.pos)
            if d <= window and d < best_d:
                best, best_d = q, d
        if best is not None:
            used.add(id(best))
            result.matched.append((t, best))
        else:
            result.missed.append(t)

    for q in query:
        if id(q) not in used:
            result.extra.append(q)

    return result
