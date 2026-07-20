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

import os
import sys
from dataclasses import dataclass, field
from typing import List, Tuple

from .vcf_io import MEICall

DEFAULT_WINDOW = 25

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "..", "xtea"))
from locus_table import LocusTable, TYPE_PRIORITY, MERGE_WINDOW  # noqa: E402

# MEICall.mei_type (ALU/LINE1/SVA canonical tokens, tests/utils/vcf_io.py) <-> locus_table's
# L1/SVA/Alu tokens (xtea/locus_table.py TYPE_PRIORITY / TYPE_TO_RTYPE).
_CANON_TO_LOCUS_TYPE = {"LINE1": "L1", "SVA": "SVA", "ALU": "Alu"}


def dedupe_conflict_loci(calls: List[MEICall], window: int = MERGE_WINDOW) -> List[MEICall]:
    """Collapse calls that are the SAME physical locus reported under more than one MEI type.

    Master's per-type VCFs are generated independently and never re-merged, so a locus near a
    type boundary (e.g. an SVA insertion that also clips an Alu-family match) can legitimately
    appear once per type it was independently detected under. xTea2's move-A locus table
    (xtea/locus_table.py) resolves each such conflict to exactly ONE winning type before ever
    writing a VCF, by design (see docs/handoff/MOVE_A_PLAN.md) — so comparing xtea2's (correctly
    deduplicated) output against a golden built by naively unioning master's per-type VCFs makes
    every one of these ~18-20 known conflict loci genome-wide read as a false negative, even
    though xtea2 called the locus, just under the type the shared rule says should win.

    Applies the EXACT SAME rule xtea2 itself uses (reuses `LocusTable`, not a re-derivation of it):
    types resolve in priority order L1 -> SVA -> Alu, a call within +-window bp of an
    already-placed call is the same locus, and `two_side_tprt_both` (SUBTYPE=two_side_tprt_both
    in the VCF INFO field) overrides load-order priority. Only meaningful for a TRUTH set built
    by unioning multiple independent per-type VCFs (e.g. the golden) -- do not apply to xtea2's own
    query set, which is already deduplicated upstream of any VCF ever being written.
    """
    table = LocusTable(window=window)
    by_id = {}
    orig_index = {id(call): i for i, call in enumerate(calls)}
    for type_ in TYPE_PRIORITY:
        canon = {v: k for k, v in _CANON_TO_LOCUS_TYPE.items()}[type_]
        for call in calls:
            if call.mei_type != canon:
                continue
            tprt = call.info.get("SUBTYPE") == "two_side_tprt_both"
            loc = _IdLocus(call.chrom, call.pos, type_, tprt, str(id(call)))
            by_id[loc.line] = call
            table._place(loc)  # noqa: SLF001 -- reusing the production placement rule verbatim
    survivors = [loc for bucket in table._by_chrm.values() for loc in bucket]  # noqa: SLF001
    winners = [by_id[loc.line] for loc in survivors]
    # LocusTable groups by type (load-order) then chromosome; restore the caller's original
    # relative order so a no-conflict input is returned unchanged (order carries no semantic
    # meaning downstream, but a stable, predictable result is easier to reason about/test).
    winners.sort(key=lambda call: orig_index[id(call)])
    return winners


class _IdLocus:
    """A locus_table.Locus stand-in whose `line` slot carries an id-lookup key instead of a raw
    hc_pf text line — dedupe_conflict_loci works on parsed MEICall objects, not hc_pf files, but
    LocusTable._place only reads chrm/pos/type/tprt and copies `line` through unexamined."""
    __slots__ = ("chrm", "pos", "type", "tprt", "line")

    def __init__(self, chrm, pos, type_, tprt, line):
        self.chrm = chrm
        self.pos = pos
        self.type = type_
        self.tprt = tprt
        self.line = line


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
