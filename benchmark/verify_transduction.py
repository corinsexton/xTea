#!/usr/bin/env python3
"""Decide from the BAM ALONE whether a called transduction is real. No caller in the loop.

    python benchmark/verify_transduction.py \
        --pool-tsv <td_pool.tsv> --bam <sample.bam> --ref <hg38_no_alt.fa> \
        --rep-lib-dir /home/cos689/xTea/rep_lib_annotation --out <dir> \
        [--baseline <frozen verdicts.tsv>]     # sec 10.4's no-regression gate

WHAT THIS IS. `HANDOFF_2026-08-04_stage3_td_gate_spec.md` sec 3 calls this "the actual
deliverable": for each pooled `(insertion locus X, source Y)` pair, ask whether the insertion at X
carries sequence from downstream of Y, using the reads and not the caller's own bookkeeping. That
is what makes the answer truth rather than agreement -- sec 6.4's locked constraint is that
"agreement with master is explicitly not a requirement", so a gate that scored against master's
calls would contradict the plan it serves.

The pool comes from `docs/handoff/runs/probe_td_pool.py`; build it first. This program never reads
a VCF.

THE FIVE DECISIONS, and each one is either the codebase's own or reported as a sweep.

1. **The source's transduced window is `[Y_end, Y_end + W]`, `W = flklen = 3000`** -- the step's own
   flank length (`run_xtea`'s golden defaults `--nfclip 4 --nfdisc 5 --flklen 3000`), so the gate
   and the caller mean the same thing by "downstream". EVERY source class has such a window,
   including polymorphic: what a polymorphic source lacks in the reference is the element BODY, and
   the transduced sequence comes from downstream of its 3' end, which is ordinary genomic sequence.
   `[ins_pos, ins_pos + 3000]` is exactly the window `gnrt_flank_regions_of_polymerphic_insertions`
   cuts its flank from, which is why it writes `start == end`. A sibling's point already carries
   xTea's own `+-100` from the pool builder.

   **BOTH SIDES ARE SEARCHED BY DEFAULT, and this is spec sec 8's open question, not a shortcut.**
   `TD_SRC` does not say which strand the source element sits on, so "downstream" is only oriented
   once you know that. Searching `[end, end+W]` and `[start-W, start]` and REPORTING WHICH SIDE
   supported is strictly more informative than guessing, and it makes the orientation question
   answerable from output instead of from argument -- if support lands overwhelmingly on one side,
   sec 8 is answered. `--one-sided` restricts to downstream for anyone who wants the stricter claim.

2. **The reads at X** (sec 3.2): soft-clipped reads within the scorer's `+-200` bp whose clipped
   part is the candidate sequence, plus discordant reads whose MATE lands in the source window.
   The two evidence types are **counted and reported separately** because they are not the same
   claim and their (a)-side tests are not equally strong -- see 3.

3. **Two alignments per candidate, and for disc reads one of them is already done.** sec 3.3 asks
   for alignment to (a) the source window and (b) the family consensus, with support requiring (a)
   and NOT (b). For a clipped part both are run here, with **bwa -- the same aligner the step uses**,
   so "aligns to the source window" means what the caller means by it. For a discordant read, (a) is
   the BAM's OWN mate placement, which is a *stronger* test than this one: the mate was placed
   against the whole genome by the aligner, not forced into a 3 kb window. Only (b) is run.

   Reported separately rather than summed into one number precisely because of that asymmetry.

4. **A source window that is itself a same-family repeat is UNDECIDABLE, not REFUTED** (sec 3.4).
   `annotation1` is the index the transduction step is actually handed -- `hg38_L1.fa.out` for L1,
   `hg38_SVA.out` for SVA -- and a window overlapping a same-family copy below
   `TD_REP_DIVERGENT_CUTOFF = 5`% divergence is dropped, because there (a) and (b) are the same
   sequence and no read can distinguish them. **THIS IS THE EXACT TRAP sec 7.1 HIT**: 11 of 20 L1
   residuals >= 100 bp mapped INSIDE a reference L1 and every one was an older subfamily, so
   "residual" was divergent body. A verifier without this test confirms element body as transduction.
   The 5% is the step's own constant, and it is swept.

5. **Four states, and none of them is folded into another.** `VERIFIED` (support >= N), `REFUTED`
   (support == 0), `SUB_THRESHOLD` (0 < support < N), `UNDECIDABLE` (repeat window, or local
   coverage below the floor). Precision is `VERIFIED / (VERIFIED + REFUTED)` with the other two
   reported beside it and never in the denominator -- the same rule the decision ledger enforces
   between `no_ledger_row` and "unattributed". Getting that wrong is how a set difference produced
   184 "threshold rejections" where there were 38.

   **`SUB_THRESHOLD` IS A GAP IN THE SPEC, MADE VISIBLE RATHER THAN CLOSED.** sec 3.5 names
   `VERIFIED (>=N)` and `REFUTED (0)` and leaves `0 < support < N` unnamed. Folding it into REFUTED
   would let the choice of N manufacture false positives; folding it into VERIFIED would let N do
   nothing. It is its own state, and it is also the single clearest reason sec 3.1 is right that a
   point estimate relocates the problem: every SUB_THRESHOLD entry becomes VERIFIED at a lower N.

WHAT IT DOES NOT CLAIM (sec 6, and this belongs next to every number it prints):
  * **Not sensitivity.** The pool is bounded by what some caller found. It cannot say what all the
    arms missed and therefore cannot state how many transductions the sample has.
  * **Not master agreement.** Master enumerates; the reads decide.
  * **Not comparable to the paper** -- hs38DH vs no_alt (plan sec 7.1).
  * **Not a verdict on Alu**, which has no transduction flank, so the step never runs. Biology.

THE PRIMARY STRATUM IS `family x call class x shape x sidedness`, NOT the source class
(`HANDOFF_2026-08-04_td_gate_built.md` sec 10.2, Corinne 2026-08-05). The source-class axis is
defined for **7 of 115** pool entries, because `TD_SRC`'s BARE shape carries two discordant-pair
midpoints rather than a source element's interval (`benchmark_truth_set.td_src_shape`) -- so
`start < end` there says nothing about copy status and the "reference-copy source" column was
reading disc spans. Shape x sidedness is fully populated AND discriminating -- measured on the
3-arm HG002 pool, L1 regular bare **two_sided 85.2%** against **one_sided 56.2%**, siblings
50-56%. (Master's NA12878 pool gives 83.9% / 55.0%: two samples, same direction, same size.)
The source-class table is KEPT BENEATH IT with explicit
`EMPTY DENOMINATOR` markers rather than deleted, because it is the axis the SVA decoy story is
told in and a future arm could populate it -- deleting it loses the axis, silencing it lets a
reader trust an empty cell.

Both facts come from the pool row and neither is a new measurement: the shape from
`benchmark_truth_set.td_src_shape(td_src_raw)` and the sidedness from `src_start < 0`, the same
test `PoolRow.windows` already makes. NOTHING here reparses `TD_SRC` with its own code -- sec 1's
rule, because a probe with its own copy of a rule can agree with a scorer that disagrees with it.

`--baseline <verdicts.tsv>` IS THE NO-REGRESSION GATE (sec 10.4). It compares **absolute VERIFIED
counts per stratum** and fails if any stratum loses entries -- not a percentage tolerance, because
the strata that matter most hold 1-4 entries, where a percentage is noise and one reclassification
swings it by 100 points. Its margin prints `PROVISIONAL` and stays provisional until an arm exists
that is neither master nor xTea2: three-arm agreement (82 of 115) is a floor on independence and
not a measure of it, because two of the three arms differ by ONE FLAG.
"""
import argparse
import itertools
import os
import re
import subprocess
import sys
from collections import Counter, defaultdict, OrderedDict

HERE = os.path.dirname(os.path.abspath(__file__))
REPO_ROOT = os.path.dirname(HERE)
sys.path.insert(0, os.path.join(REPO_ROOT, "xtea"))
sys.path.insert(0, HERE)

import pysam  # noqa: E402

# The shape rule lives in ONE place. Importing it is sec 1's discipline, not convenience: this
# program already refuses to parse `TD_SRC` itself, and `td_src_shape` is the function that knows
# a bare interval is two disc-pair midpoints rather than a source element.
import benchmark_truth_set as bts  # noqa: E402

# ---------------------------------------------------------------------------------------
# constants: every one of these is the codebase's own value, cited, and swept where sec 3.1
# says to sweep it
# ---------------------------------------------------------------------------------------
FLKLEN = 3000                 # run_xtea's golden default --flklen; the step's "downstream"
POS_WINDOW = 200              # benchmark_truth_set.DEFAULT_POS_WINDOW: the scorer's locus window
MIN_CLIP_LEN = 20             # a clipped part shorter than this cannot be placed anywhere
DEFAULT_CLIP_MAPQ = 12        # global_values.MINIMUM_CLIP_MAPQ
DEFAULT_DISC_MAPQ = 20        # global_values.MINIMUM_DISC_MAPQ
DEFAULT_CNS_FRACTION = 0.80   # sec 3.3's "over >= 80%": explained-by-consensus threshold
DEFAULT_SRC_FRACTION = 0.65   # global_values.F_MIN_TRSDCT_DISC_MAP_RATION, the step's own
DEFAULT_REP_DIV = 5.0         # global_values.TD_REP_DIVERGENT_CUTOFF
DEFAULT_COV_FLOOR = 5         # below this many reads at X nothing can be concluded
DEFAULT_N = 2                 # support required for VERIFIED; swept

# annotation1, per family: the repeat index the transduction step is handed (plan sec 10).
ANNOTATION1_RELPATH = {
    "L1": "LINE/{g}/{g}_L1.fa.out",
    "SVA": "SVA/{g}/{g}_SVA.out",
}
# The family consensus, already bwa-indexed in rep_lib -- do NOT rebuild these indexes.
CONSENSUS_RELPATH = {"L1": "consensus/LINE1.fa", "SVA": "consensus/SVA.fa",
                     "Alu": "consensus/ALU.fa"}
# Same-family repeat classes in a RepeatMasker .out, matched on the class/family column.
FAMILY_REPEAT_CLASS = {"L1": "LINE/L1", "SVA": "Retroposon/SVA"}

VERIFIED, REFUTED, SUB_THRESHOLD, UNDECIDABLE = (
    "VERIFIED", "REFUTED", "SUB_THRESHOLD", "UNDECIDABLE")
STATES = (VERIFIED, REFUTED, SUB_THRESHOLD, UNDECIDABLE)


def log(msg):
    sys.stderr.write("[verify-td] {0}\n".format(msg))
    sys.stderr.flush()


# ---------------------------------------------------------------------------------------
# the pool
# ---------------------------------------------------------------------------------------

class PoolRow(object):
    __slots__ = ("idx", "family", "chrom", "pos", "call_class", "source_class",
                 "src_chrom", "src_start", "src_end", "arms", "raw")

    def __init__(self, idx, d):
        self.idx = idx
        self.family = d["family"]
        self.chrom = d["ins_chrom"]
        self.pos = int(d["ins_pos"])
        self.call_class = d["call_class"]
        self.source_class = None if d["source_class"] == "-" else d["source_class"]
        self.src_chrom = d["src_chrom"]
        self.src_start = int(d["src_start"])
        self.src_end = int(d["src_end"])
        self.arms = d["arms"].split(",")
        self.raw = d.get("td_src_raw", "")

    @property
    def shape(self):
        """Which WRITER produced this `TD_SRC`, which decides what its coordinates mean.

        Imported, never re-derived. `bare` is `x_orphan_transduction.py:371`'s two discordant-pair
        midpoints; `tailed` is `x_transduction.py:1252`'s flank-derived source element;
        `sibling_point` has no interval at all.
        """
        return bts.td_src_shape(self.raw)

    @property
    def sidedness(self):
        """`two_sided` | `one_sided` | `-`, and this is the axis that DISCRIMINATES.

        `src_start < 0` means the second discordant-pair midpoint never clustered **and its
        absence was EXCUSED because that side's disc mates fell predominantly inside repeat
        regions** -- `x_orphan_transduction.py:1120` is the only route by which a `-1` reaches
        surviving output (`_check_form_repeat_or_same_cluster` returning `True, -1`; the
        alternative at `:1126` returns a real position). Pinned to that line, not inferred.

        **So this axis is not just "less evidence", and that is why it belongs in the primary
        table**: a one-sided entry is one whose second side is REPEAT-DRIVEN, which is exactly
        what decision 4 below exists to distrust. It verifies ~30 points worse -- 56.2% against
        85.2% here, 55.0% against 83.9% on master's NA12878 -- and now with a mechanism rather
        than a correlation.

        `src_start < 0` is the same test `windows()` makes, so the report and the window cannot
        disagree about which coordinate is usable.

        A `sibling_point` gets `-`: the pool builder gave it xTea's own [POS-100, POS+100], so its
        `src_start` is an invented bound and reading sidedness off it would be reading noise.
        """
        if self.shape == bts.SHAPE_SIBLING_PT:
            return "-"
        return "one_sided" if self.src_start < 0 else "two_sided"

    @property
    def shape_stratum(self):
        """sec 10.2's PRIMARY stratum: fully populated, and it separates the verdicts."""
        return (self.family, self.call_class, self.shape, self.sidedness)

    @property
    def source_stratum(self):
        """The SECONDARY stratum, kept with empty cells rather than deleted (sec 10.2)."""
        return (self.family, self.call_class, self.source_class or "-")

    def windows(self, w, one_sided):
        """The reference intervals the transduced sequence could have come from.

        Downstream of the element's 3' end is `[end, end + w]`. Which genomic direction that is
        depends on the source element's strand, which `TD_SRC` does not carry -- so unless
        `--one-sided`, the upstream mirror `[start - w, start]` is searched too and the report
        says which side supported. See the module docstring, decision 1.
        """
        s = self.src_end if self.src_start < 0 else self.src_start
        out = [("D", self.src_end, self.src_end + w)]
        if not one_sided:
            out.append(("U", max(0, s - w), s))
        return out


def read_pool(path):
    rows = []
    with open(path) as fin:
        header = None
        for line in fin:
            if line.startswith("#") or not line.strip():
                continue
            f = line.rstrip("\n").split("\t")
            if header is None:
                header = f
                continue
            rows.append(PoolRow(len(rows), dict(zip(header, f))))
    return rows


# ---------------------------------------------------------------------------------------
# decision 4: is the source window itself a same-family repeat
# ---------------------------------------------------------------------------------------

def load_repeat_index(path, class_name):
    """chrom -> sorted [(start, end, divergence, name)] for one repeat class.

    A RepeatMasker `.out` row is whitespace-delimited: col 2 is % divergence, cols 5-7 are
    chrom/start/end, col 10 the repeat name, col 11 the class/family. The two header lines have
    no numeric divergence, so they fall out of the float() guard rather than needing a line count.
    """
    by_chrom = defaultdict(list)
    n = 0
    with open(path) as fin:
        for line in fin:
            f = line.split()
            if len(f) < 11:
                continue
            try:
                div = float(f[1])
                start, end = int(f[5]), int(f[6])
            except ValueError:
                continue
            if f[10] != class_name:
                continue
            by_chrom[f[4]].append((start, end, div, f[9]))
            n += 1
    for c in by_chrom:
        by_chrom[c].sort()
    return by_chrom, n


def window_min_divergence(rep_index, chrom, start, end):
    """(min divergence, name) over same-family copies overlapping the window; (None, None) if none.

    The MINIMUM, and the cutoff is applied later, because the cutoff is a swept parameter and not
    a property of the window. A low-divergence copy is the dangerous one: there the source
    window's sequence and the family consensus ARE the same sequence, so no read can tell "came
    from downstream of Y" from "is element body". `TD_REP_DIVERGENT_CUTOFF = 5` is the step's own
    tolerance for exactly this reason -- its comment says "as this is transduction region,
    tolerate the wrong discordant reads from repeats".

    Storing the number rather than a boolean is what lets the verdict be a function of the cutoff,
    which for SVA turns out to be the parameter that decides the answer: SVA is a young family, so
    its reference copies sit at low divergence and this test disqualifies most of them.
    """
    best_div, best_name = None, None
    for (s, e, div, name) in rep_index.get(chrom, ()):
        if s > end:
            break
        if e < start:
            continue
        if best_div is None or div < best_div:
            best_div, best_name = div, name
    return best_div, best_name


# ---------------------------------------------------------------------------------------
# decision 2: the reads at X
# ---------------------------------------------------------------------------------------

class Candidate(object):
    __slots__ = ("row_idx", "kind", "name", "seq", "side", "mate_chrom", "mate_pos", "mapq")

    def __init__(self, row_idx, kind, name, seq, side=None,
                 mate_chrom=None, mate_pos=None, mapq=0):
        self.row_idx, self.kind, self.name, self.seq = row_idx, kind, name, seq
        self.side, self.mate_chrom, self.mate_pos, self.mapq = side, mate_chrom, mate_pos, mapq

    @property
    def qname(self):
        """The FASTA id, which has to round-trip through bwa: entry, kind, side, read."""
        return "{0}|{1}|{2}|{3}".format(self.row_idx, self.kind, self.side or "-", self.name)


def _usable(read):
    return not (read.is_unmapped or read.is_duplicate or read.is_secondary
                or read.is_supplementary)


def gather_candidates(bam, row, wins, clip_mapq, disc_mapq, max_disc):
    """(candidates, local_depth, n_disc_capped) for one pool row.

    Clipped parts are taken unconditionally within the locus window -- they are cheap and their
    (a) test is a real alignment. Discordant reads are taken ONLY when their mate already lands
    in a source window, which is both the (a) test (decision 3) and the thing that keeps the
    mate-fetch count bounded: without it every disc read at every locus would need a random seek.
    """
    lo, hi = max(0, row.pos - POS_WINDOW), row.pos + POS_WINDOW
    cands, depth, n_capped = [], 0, 0
    n_disc = 0
    for read in bam.fetch(row.chrom, lo, hi):
        if not _usable(read):
            continue
        depth += 1
        cig = read.cigartuples
        if cig and read.mapping_quality >= clip_mapq and read.query_sequence:
            # 4 == BAM_CSOFT_CLIP. A clipped part IS by construction the part that did not map
            # at X, so sec 3.2's "whose clipped part does not map at X" needs no extra test.
            if cig[0][0] == 4 and cig[0][1] >= MIN_CLIP_LEN:
                cands.append(Candidate(row.idx, "clip", read.query_name + "/L",
                                       read.query_sequence[:cig[0][1]],
                                       mapq=read.mapping_quality))
            if cig[-1][0] == 4 and cig[-1][1] >= MIN_CLIP_LEN:
                cands.append(Candidate(row.idx, "clip", read.query_name + "/R",
                                       read.query_sequence[-cig[-1][1]:],
                                       mapq=read.mapping_quality))
        if (read.mapping_quality >= disc_mapq and read.is_paired
                and not read.mate_is_unmapped and read.next_reference_name):
            for (side, ws, we) in wins:
                if (read.next_reference_name == row.src_chrom
                        and ws <= read.next_reference_start <= we):
                    if n_disc >= max_disc:
                        n_capped += 1
                        break
                    n_disc += 1
                    cands.append(Candidate(row.idx, "disc", read.query_name, None, side=side,
                                           mate_chrom=read.next_reference_name,
                                           mate_pos=read.next_reference_start,
                                           mapq=read.mapping_quality))
                    break
    return cands, depth, n_capped


def fetch_mate_sequences(bam, cands):
    """Fill in `seq` for disc candidates. One random seek each, which is why decision 2 only
    admits disc reads whose mate ALREADY landed in a source window."""
    n_ok = n_fail = 0
    by_name = defaultdict(list)
    for c in cands:
        if c.kind == "disc":
            by_name[(c.mate_chrom, c.mate_pos, c.name)].append(c)
    for (chrom, pos, name), group in by_name.items():
        seq = None
        try:
            for read in bam.fetch(chrom, max(0, pos - 1), pos + 2):
                if read.query_name == name and _usable(read) and read.query_sequence:
                    seq = read.query_sequence
                    break
        except ValueError:
            seq = None
        for c in group:
            c.seq = seq
        if seq:
            n_ok += 1
        else:
            n_fail += 1
    return n_ok, n_fail


# ---------------------------------------------------------------------------------------
# decision 3: the two alignments, batched through bwa
# ---------------------------------------------------------------------------------------

def aligned_fraction(cigar, qlen):
    """Query bases consumed by M/=/X/I over the read's full length.

    `qlen` is the CANDIDATE's own length, taken from the FASTA we wrote and never from the SAM
    record. That is not fussiness: `bwa mem -a` emits secondary records with `SEQ = *`, so a
    candidate whose best hit is some other entry's window would have length 0 -- and its
    own-entry hit, the only one that counts, would score 0.0 and read as a refutation.
    """
    if not cigar or cigar == "*" or not qlen:
        return 0.0
    total = 0
    num = ""
    for ch in cigar:
        if ch.isdigit():
            num += ch
            continue
        if ch in "M=XI":
            total += int(num or 0)
        num = ""
    return float(total) / qlen


def write_fasta(path, records):
    with open(path, "w") as fout:
        for name, seq in records:
            fout.write(">{0}\n{1}\n".format(name, seq))


def bwa_align(fa_query, fa_ref, threads, out_sam, index=False):
    """bwa mem -a, so a query's alignment to ITS OWN entry is visible even when its best hit is
    somewhere else in the batched reference. Returns the SAM path."""
    if index:
        subprocess.check_call(["bwa", "index", fa_ref],
                              stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    with open(out_sam, "w") as fout:
        subprocess.check_call(["bwa", "mem", "-a", "-t", str(threads), fa_ref, fa_query],
                              stdout=fout, stderr=subprocess.DEVNULL)
    return out_sam


def best_fractions(sam_path, qlens, want=None):
    """qname -> (best aligned fraction, the reference it hit), over alignments that COUNT.

    `want` maps qname -> the set of reference names that count; a hit to any other reference is
    ignored, which is how one batched source-window reference serves every pool entry without
    cross-crediting a read to a neighbour's source. `qlens` comes from the candidates, for the
    reason `aligned_fraction`'s docstring gives -- and it must contain ONLY the queries this SAM
    was asked about, because every key in it is seeded to 0.0 and a caller that passes a superset
    silently zeroes results established elsewhere. (It did: disc candidates' (a) verdict, which is
    set by construction to 1.0, came back 0 and every disc read read as unsupporting.)

    The reference name comes back too, because for a clipped part it carries WHICH SIDE of the
    source window the support came from, which is spec sec 8's open question.
    """
    best = dict((q, (0.0, None)) for q in qlens)
    with open(sam_path) as fin:
        for line in fin:
            if line.startswith("@"):
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 6:
                continue
            qname, rname, cigar = f[0], f[2], f[5]
            if rname == "*":
                continue
            if want is not None and rname not in want.get(qname, ()):
                continue
            frac = aligned_fraction(cigar, qlens.get(qname, 0))
            if frac > best.get(qname, (0.0, None))[0]:
                best[qname] = (frac, rname)
    return best


# ---------------------------------------------------------------------------------------
# decision 5: the verdict, as a function of the thresholds
# ---------------------------------------------------------------------------------------

class RowResult(object):
    __slots__ = ("row", "rep_div", "rep_name", "rep_untested", "depth", "support")

    def __init__(self, row):
        self.row = row
        self.rep_div = None        # min divergence of an overlapping same-family copy
        self.rep_name = None
        self.rep_untested = False  # the annotation was missing: a GAP, not a pass
        self.depth = 0
        # per candidate: (kind, side, consensus_fraction, source_fraction)
        self.support = []

    def n_support(self, cns_frac, src_frac, kinds=("clip", "disc")):
        """Candidates that align to the source window and are NOT explained by the consensus."""
        return sum(1 for (kind, _side, f_cns, f_src) in self.support
                   if kind in kinds and f_src >= src_frac and f_cns < cns_frac)

    def sides(self, cns_frac, src_frac):
        return Counter(side for (_k, side, f_cns, f_src) in self.support
                       if f_src >= src_frac and f_cns < cns_frac)

    def verdict(self, n_req, cns_frac, src_frac, cov_floor, rep_div_cutoff):
        """The state, as a function of all five thresholds. Order matters: a window that cannot
        be interrogated is UNDECIDABLE regardless of how many reads support it, because those
        reads cannot distinguish transduction from element body."""
        if self.rep_div is not None and self.rep_div < rep_div_cutoff:
            return UNDECIDABLE, "repeat_window({0}@{1:.1f}%)".format(self.rep_name, self.rep_div)
        if self.depth < cov_floor:
            return UNDECIDABLE, "low_coverage(depth={0})".format(self.depth)
        n = self.n_support(cns_frac, src_frac)
        if n >= n_req:
            return VERIFIED, "support={0}".format(n)
        if n == 0:
            return REFUTED, "support=0"
        return SUB_THRESHOLD, "support={0}<{1}".format(n, n_req)


def _rate(num, den):
    return float("nan") if den == 0 else 100.0 * num / den


def tally_strata(results, point, which):
    """stratum tuple -> Counter of states, at one operating point. `which` names the property."""
    out = defaultdict(Counter)
    for res in results.values():
        out[getattr(res.row, which)][res.verdict(*point)[0]] += 1
    return out


def print_strata(strata, label_width):
    """One stratified table, with the empty-denominator marker under the rows that need it."""
    print()
    print("  {0:<{w}} {1:>8} {2:>8} {3:>8} {4:>8}  {5:>10}".format(
        "stratum", VERIFIED, REFUTED, "SUB_THR", "UNDEC", "precision", w=label_width))
    for key in sorted(strata, key=lambda t: tuple(str(x) for x in t)):
        t = strata[key]
        den = t[VERIFIED] + t[REFUTED]
        label = " ".join(str(x) for x in key)
        print("  {0:<{w}} {1:>8} {2:>8} {3:>8} {4:>8}  {5:>9.1f}%".format(
            label, t[VERIFIED], t[REFUTED], t[SUB_THRESHOLD], t[UNDECIDABLE],
            _rate(t[VERIFIED], den), w=label_width))
        if den == 0:
            print("  {0:<{w}} ^ EMPTY DENOMINATOR -- this cell says nothing. Do not read the "
                  "0.0% as a rate.".format("", w=label_width))


# ---------------------------------------------------------------------------------------
# sec 10.4: the no-regression gate
# ---------------------------------------------------------------------------------------

# `# per-entry verdicts at N=3 f_cns=0.8 f_src=0.65 cov>=5 rep_div<5.0.` -- the header this
# program writes. Parsed rather than trusted, so a baseline measured at another operating point
# is REFUSED instead of silently compared. See `read_baseline_verdicts`.
_OP_POINT_RE = re.compile(
    r"per-entry verdicts at N=(\S+) f_cns=(\S+) f_src=(\S+) cov>=(\S+) rep_div<([^\s.]+(?:\.\d+)?)")


def format_op_point(point):
    return "N={0} f_cns={1} f_src={2} cov>={3} rep_div<{4}".format(*point)


def read_baseline_verdicts(path):
    """(op_point_string | None, {shape_stratum: n_verified}, {source_stratum: n_verified}, n_rows).

    **THE STRATA ARE RE-DERIVED HERE, NOT READ FROM COLUMNS**, from `td_src_raw` and `src_start`.
    That is what lets the FROZEN baseline -- written before the shape axis existed, so it has no
    `shape` or `sidedness` column -- be compared against a run that has them, with exactly one
    implementation of the rule in the process. A file that carries the columns is not trusted over
    the derivation; if the two ever disagreed it would mean `td_src_shape` changed under a frozen
    file, which is a thing to be told about loudly rather than to paper over.
    """
    op_point = None
    shape_v, source_v = Counter(), Counter()
    n_rows = 0
    header = None
    with open(path) as fin:
        for line in fin:
            if line.startswith("#"):
                m = _OP_POINT_RE.search(line)
                if m:
                    op_point = format_op_point(m.groups())
                continue
            if not line.strip():
                continue
            f = line.rstrip("\n").split("\t")
            if header is None:
                header = f
                continue
            d = dict(zip(header, f))
            n_rows += 1
            raw = d.get("td_src_raw", "")
            shape = bts.td_src_shape(raw)
            if shape == bts.SHAPE_SIBLING_PT:
                sided = "-"
            else:
                sided = "one_sided" if int(d["src_start"]) < 0 else "two_sided"
            if d["verdict"] != VERIFIED:
                continue
            shape_v[(d["family"], d["call_class"], shape, sided)] += 1
            source_v[(d["family"], d["call_class"], d["source_class"])] += 1
    return op_point, shape_v, source_v, n_rows


def compare_no_regression(name, base_v, cur_strata):
    """PASS/FAIL on ABSOLUTE VERIFIED counts per stratum (sec 10.4). Returns True on PASS.

    Not a percentage tolerance, and the reason is the data: `SVA regular reference_copy` has ONE
    entry, where a percentage is noise and a single reclassification swings it by 100 points.
    Losses fail; gains are reported as growth and never offset a loss elsewhere -- a stratum that
    gained cannot pay for a stratum that lost, because they are different hypotheses about
    different biology.
    """
    cur_v = dict((k, t[VERIFIED]) for k, t in cur_strata.items())
    keys = sorted(set(base_v) | set(cur_v), key=lambda t: tuple(str(x) for x in t))
    lost, gained = [], []
    print()
    print("  --- {0} ---".format(name))
    print("  {0:<40} {1:>9} {2:>9} {3:>8}".format("stratum", "baseline", "current", "delta"))
    for k in keys:
        b, c = base_v.get(k, 0), cur_v.get(k, 0)
        if b == 0 and c == 0:
            continue
        mark = ""
        if c < b:
            lost.append((k, b, c))
            mark = "  <-- REGRESSION"
        elif c > b:
            gained.append((k, b, c))
            mark = "  growth"
        print("  {0:<40} {1:>9} {2:>9} {3:>+8}{4}".format(
            " ".join(str(x) for x in k), b, c, c - b, mark))
    print("  totals: baseline {0} VERIFIED, current {1} VERIFIED".format(
        sum(base_v.values()), sum(cur_v.values())))
    if lost:
        print("  {0}: FAIL -- {1} stratum/strata LOST verified entries".format(name, len(lost)))
        for k, b, c in lost:
            print("      {0}: {1} -> {2}".format(" ".join(str(x) for x in k), b, c))
        return False
    print("  {0}: PASS -- no stratum lost a verified entry ({1} gained)".format(name, len(gained)))
    return True


# ---------------------------------------------------------------------------------------
# main
# ---------------------------------------------------------------------------------------

def main():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--pool-tsv", required=True, help="from docs/handoff/runs/probe_td_pool.py")
    p.add_argument("--bam", required=True)
    p.add_argument("--ref", required=True, help="the reference the BAM's contigs live in")
    p.add_argument("--rep-lib-dir", required=True)
    p.add_argument("--genome", default="hg38")
    p.add_argument("--out", required=True, help="directory for the per-entry TSV and the SAMs")
    p.add_argument("--families", default="L1,SVA",
                   help="Alu has no transduction flank, so it is not in the default")
    p.add_argument("--flklen", type=int, default=FLKLEN)
    p.add_argument("--one-sided", action="store_true",
                   help="search only [end, end+W]; see decision 1 on orientation")
    p.add_argument("--threads", type=int, default=int(os.environ.get("SLURM_CPUS_PER_TASK", 4)))
    p.add_argument("--max-disc", type=int, default=200,
                   help="per-entry cap on discordant candidates; over-cap counts are REPORTED, "
                        "never silently dropped (default %(default)s)")
    # the sweep. Single values give a point estimate, which sec 3.1 says is not enough.
    p.add_argument("--sweep-n", default="1,2,3,4")
    p.add_argument("--sweep-cns-frac", default="0.70,0.80,0.90")
    p.add_argument("--sweep-src-frac", default="0.50,0.65,0.80")
    p.add_argument("--sweep-cov-floor", default="0,5,10")
    p.add_argument("--sweep-rep-div", default="2,5,10",
                   help="repeat-divergence cutoff; a window inside a same-family copy\n"
                        "BELOW this is UNDECIDABLE. The step's own value is 5 "
                        "(TD_REP_DIVERGENT_CUTOFF)")
    p.add_argument("--clip-mapq", type=int, default=DEFAULT_CLIP_MAPQ)
    p.add_argument("--disc-mapq", type=int, default=DEFAULT_DISC_MAPQ)
    p.add_argument("--baseline", metavar="VERDICTS_TSV",
                   help="a frozen verdicts.tsv to check for regression (sec 10.4). Fails if any "
                        "stratum LOSES verified entries. Exits 1 on regression, 3 if the "
                        "baseline's operating point does not match this run's midpoint")
    a = p.parse_args()

    families = [f.strip() for f in a.families.split(",") if f.strip()]
    os.makedirs(a.out, exist_ok=True)
    rows = [r for r in read_pool(a.pool_tsv) if r.family in families]

    print("=" * 78)
    print("VERIFY TRANSDUCTION -- from the BAM, no caller in the loop")
    print("  pool     : {0} ({1} entr(ies) in {2})".format(
        a.pool_tsv, len(rows), ",".join(families)))
    print("  bam      : {0}".format(a.bam))
    print("  window   : [Y_end, Y_end+{0}]{1}".format(
        a.flklen, "" if a.one_sided else " and its upstream mirror (orientation unknown)"))
    print("  NOT sensitivity, NOT master agreement, NOT comparable to the paper. See the")
    print("  module docstring's last block before quoting anything below.")
    print("=" * 78)

    # --- decision 4: the repeat context, as a NUMBER per entry, not a verdict -----------------
    # The cutoff is swept, so the divergence is measured for every entry and compared later. An
    # earlier version applied the cutoff here and skipped the read work for whatever it excluded,
    # which quietly made the sweep a no-op: an entry the 5% cutoff rejects has to be decidable at
    # 2% or the sweep is reporting the same answer four times.
    results = OrderedDict((r.idx, RowResult(r)) for r in rows)
    for fam in families:
        rel = ANNOTATION1_RELPATH.get(fam)
        cls = FAMILY_REPEAT_CLASS.get(fam)
        fam_rows = [r for r in rows if r.family == fam]
        if not rel or not cls or not os.path.isfile(
                os.path.join(a.rep_lib_dir, rel.format(g=a.genome))):
            log("{0}: no usable annotation1 -- the same-family-repeat test CANNOT RUN for its "
                "{1} entr(ies). They stay decidable, which is a GAP and not a pass."
                .format(fam, len(fam_rows)))
            for r in fam_rows:
                results[r.idx].rep_untested = True
            continue
        path = os.path.join(a.rep_lib_dir, rel.format(g=a.genome))
        idx, n = load_repeat_index(path, cls)
        log("{0}: {1} {2} copies from {3}".format(fam, n, cls, os.path.basename(path)))
        for r in fam_rows:
            for (_side, ws, we) in r.windows(a.flklen, a.one_sided):
                div, name = window_min_divergence(idx, r.src_chrom, ws, we)
                if div is not None and (results[r.idx].rep_div is None
                                        or div < results[r.idx].rep_div):
                    results[r.idx].rep_div, results[r.idx].rep_name = div, name

    n_untested = sum(1 for v in results.values() if v.rep_untested)
    print()
    print("  step 4  source-window repeat context measured for {0} of {1} entr(ies); the cutoff "
          "is applied in the sweep".format(len(rows) - n_untested, len(rows)))
    print("          A window inside a low-divergence same-family copy is UNDECIDABLE, never")
    print("          REFUTED: where (a) and (b) are the same sequence no read separates")
    print("          transduction from element body -- the exact trap sec 7.1 hit.")
    if n_untested:
        print("          *** {0} entr(ies) could not be tested at all (missing annotation1). They"
              " are NOT known-clean.".format(n_untested))

    # --- decision 2: gather the reads, for EVERY entry ---------------------------------------
    bam = pysam.AlignmentFile(a.bam, "rb")
    cands, n_capped_total = [], 0
    for i, r in enumerate(rows):
        if i and i % 50 == 0:
            log("gathered {0}/{1}".format(i, len(rows)))
        wins = r.windows(a.flklen, a.one_sided)
        c, depth, n_capped = gather_candidates(bam, r, wins, a.clip_mapq, a.disc_mapq, a.max_disc)
        results[r.idx].depth = depth
        cands.extend(c)
        n_capped_total += n_capped
    n_ok, n_fail = fetch_mate_sequences(bam, cands)
    bam.close()
    cands = [c for c in cands if c.seq]
    print()
    print("  step 2  {0} candidate read part(s) over {1} entr(ies): {2} clip, {3} disc"
          .format(len(cands), len(rows),
                  sum(1 for c in cands if c.kind == "clip"),
                  sum(1 for c in cands if c.kind == "disc")))
    print("          mate sequences resolved {0}, unresolved {1}".format(n_ok, n_fail))
    if n_capped_total:
        print("          *** {0} discordant candidate(s) OVER the per-entry cap of {1} were not "
              "gathered. Raise --max-disc to remove the cap from the answer."
              .format(n_capped_total, a.max_disc))

    if not cands:
        print("\n  no candidates: nothing to align. Every decidable entry will be REFUTED at any "
              "N, which is a statement about the reads and not about the verifier.")

    # --- decision 3: the two alignments ------------------------------------------------------
    ref = pysam.FastaFile(a.ref)
    win_records, want = [], defaultdict(set)
    for r in rows:
        for (side, ws, we) in r.windows(a.flklen, a.one_sided):
            name = "e{0}_{1}".format(r.idx, side)
            try:
                seq = ref.fetch(r.src_chrom, max(0, ws), we)
            except (KeyError, ValueError):
                seq = ""
            if len(seq) < MIN_CLIP_LEN:
                continue
            win_records.append((name, seq))
    ref.close()
    for c in cands:
        # a candidate may be credited by either side's window; the report says which won
        want[c.qname] = set("e{0}_{1}".format(c.row_idx, s)
                            for s in (("D",) if a.one_sided else ("D", "U")))

    # A DISCORDANT candidate's (a) test is already passed and is not re-run here. `gather_
    # candidates` only admitted it because the whole-genome aligner had ALREADY placed its mate
    # inside the source window -- a stronger statement than a forced local alignment to a 3 kb
    # reference, which can also fail spuriously when the mate straddles the window's edge.
    # Re-aligning it would only ever lose real support. Its (b) veto still runs, below, and its
    # side came from the window its mate landed in.
    frac_src = dict((c.qname, 1.0) for c in cands if c.kind == "disc")
    frac_cns = {}
    clip_side = {}
    clips = [c for c in cands if c.kind == "clip"]
    if clips and win_records:
        q_fa = os.path.join(a.out, "candidates.fa")
        w_fa = os.path.join(a.out, "source_windows.fa")
        write_fasta(q_fa, [(c.qname, c.seq) for c in clips])
        write_fasta(w_fa, win_records)
        log("bwa: {0} clipped part(s) vs {1} source window(s)".format(len(clips),
                                                                     len(win_records)))
        got = best_fractions(
            bwa_align(q_fa, w_fa, a.threads, os.path.join(a.out, "vs_source.sam"), index=True),
            dict((c.qname, len(c.seq)) for c in clips), want=want)
        for q, (frac, rname) in got.items():
            frac_src[q] = frac
            # `e{idx}_{side}` -- the window it actually aligned to, which is sec 8's answer for
            # the clipped evidence. Recovered from the reference name rather than assumed.
            clip_side[q] = rname.rsplit("_", 1)[1] if rname else "?"
    if cands:
        # (b) per family, against the pre-indexed consensus. Do NOT re-index rep_lib.
        by_fam = defaultdict(list)
        for c in cands:
            by_fam[results[c.row_idx].row.family].append(c)
        for fam, group in by_fam.items():
            cns = os.path.join(a.rep_lib_dir, CONSENSUS_RELPATH[fam])
            if not os.path.isfile(cns + ".bwt"):
                log("{0}: consensus {1} is not bwa-indexed; the (b) veto cannot run and every "
                    "candidate will count as unexplained. That INFLATES support.".format(fam, cns))
                continue
            fa = os.path.join(a.out, "candidates_{0}.fa".format(fam))
            write_fasta(fa, [(c.qname, c.seq) for c in group])
            for q, (frac, _r) in best_fractions(
                    bwa_align(fa, cns, a.threads,
                              os.path.join(a.out, "vs_cns_{0}.sam".format(fam))),
                    dict((c.qname, len(c.seq)) for c in group)).items():
                frac_cns[q] = frac

    for c in cands:
        side = c.side or clip_side.get(c.qname, "?")
        results[c.row_idx].support.append(
            (c.kind, side, frac_cns.get(c.qname, 0.0), frac_src.get(c.qname, 0.0)))

    # --- decision 5: the verdict, swept ------------------------------------------------------
    def parse_sweep(s, cast):
        return [cast(x) for x in s.split(",") if x.strip()]

    axes = OrderedDict((
        ("N", parse_sweep(a.sweep_n, int)),
        ("f_cns", parse_sweep(a.sweep_cns_frac, float)),
        ("f_src", parse_sweep(a.sweep_src_frac, float)),
        ("cov>=", parse_sweep(a.sweep_cov_floor, int)),
        ("rep_div<", parse_sweep(a.sweep_rep_div, float)),
    ))
    mid = tuple(v[len(v) // 2] for v in axes.values())

    def tally_at(point):
        t = Counter()
        for res in results.values():
            t[res.verdict(*point)[0]] += 1
        return t

    print()
    print("=" * 78)
    print("THE SWEEP -- sec 3.1: N, the consensus fraction, the source fraction, the coverage")
    print("floor and the repeat-divergence cutoff are as arbitrary as the constants this gate")
    print("exists to interrogate, so the verdict is reported AS A FUNCTION of them. A single")
    print("operating point relocates the problem rather than solving it.")
    print("=" * 78)
    n_points = 1
    for v in axes.values():
        n_points *= len(v)
    print("  ONE AXIS AT A TIME around the midpoint: the full cross product is {0} points and\n"
          "  printing it would bury the shape. Its precision RANGE is quoted below, so the honest\n"
          "  bounds still come from all {0} while this table is {1} rows.".format(
              n_points, sum(len(v) for v in axes.values())))
    print()
    print("  {0:<10} {1:>8}  {2:>8} {3:>8} {4:>8} {5:>8}  {6:>10}".format(
        "axis", "value", VERIFIED, REFUTED, "SUB_THR", "UNDEC", "precision"))
    for ai, (name, values) in enumerate(axes.items()):
        for v in values:
            point = list(mid)
            point[ai] = v
            t = tally_at(tuple(point))
            den = t[VERIFIED] + t[REFUTED]
            mark = " *" if v == mid[ai] else ""
            print("  {0:<10} {1:>8}  {2:>8} {3:>8} {4:>8} {5:>8}  {6:>9.1f}%{7}".format(
                name, v, t[VERIFIED], t[REFUTED], t[SUB_THRESHOLD], t[UNDECIDABLE],
                _rate(t[VERIFIED], den), mark))
    # the full cross product, computed for the range but not printed
    grid = []
    for point in itertools.product(*axes.values()):
        t = tally_at(point)
        den = t[VERIFIED] + t[REFUTED]
        grid.append((point, t, _rate(t[VERIFIED], den)))
    prec = [g[2] for g in grid if g[2] == g[2]]
    print()
    print("  * = the midpoint, which is what the stratified table below uses.")
    if prec:
        print("  precision over the FULL cross product ({0} points): {1:.1f}% - {2:.1f}%. QUOTE "
              "THE RANGE, not a point.".format(n_points, min(prec), max(prec)))
        lo = min(grid, key=lambda g: g[2] if g[2] == g[2] else 1e9)
        hi = max(grid, key=lambda g: g[2] if g[2] == g[2] else -1e9)
        print("    worst  {0} -> {1:.1f}%".format(
            dict(zip(axes, lo[0])), lo[2]))
        print("    best   {0} -> {1:.1f}%".format(
            dict(zip(axes, hi[0])), hi[2]))
    else:
        print("  every point has an EMPTY denominator: no entry is both decidable and read-backed."
              " There is no precision to quote, and that is the result.")
    print("  denominator is VERIFIED + REFUTED. SUB_THRESHOLD and UNDECIDABLE are reported")
    print("  beside it and are NEVER in it (decision 5).")

    # --- sec 4 + sec 10.2: stratified, at the sweep's midpoint --------------------------------
    shape_strata = tally_strata(results, mid, "shape_stratum")
    source_strata = tally_strata(results, mid, "source_stratum")

    print()
    print("=" * 78)
    print("PRIMARY STRATUM at N={0} f_cns={1} f_src={2} cov>={3} rep_div<{4}".format(*mid))
    print("  family x call class x SHAPE x SIDEDNESS -- sec 10.2, Corinne 2026-08-05")
    print("=" * 78)
    print("  Strata are not cosmetic: an `~orphan`-TAILED source has anchors bracketing the")
    print("  insertion point, so there is no element body at the site, and a sibling's source")
    print("  element is absent from the reference. sec 4.")
    print("  NOT the 0/37: that number is `SUBTYPE == orphan_or_sibling_transduction`, a")
    print("  population DISJOINT from this one by construction (the created-site writer emits a")
    print("  bare source with no tail, so all 37 are `regular` here). See benchmark_truth_set's")
    print("  CALL_CLASS_* block. Two fields, two words, and they agree on no call.")
    print("  SHAPE is which writer produced TD_SRC: `bare` is two discordant-pair midpoints")
    print("  (x_orphan_transduction.py:371), `tailed` a flank-derived source element")
    print("  (x_transduction.py:1252), `sibling_point` no interval at all. SIDEDNESS is whether")
    print("  both midpoints were established -- src_start = -1 says one was not.")
    print_strata(shape_strata, 40)

    print()
    print("=" * 78)
    print("SECONDARY STRATUM -- family x call class x source class")
    print("=" * 78)
    print("  KEPT, NOT PROMOTED, and kept for a stated reason (sec 10.2). The source class is")
    print("  only defined where the interval comes from a flank record, i.e. the `tailed` shape:")
    print("  7 of 115 pool entries. Every other row's `disc_span` means THE AXIS DOES NOT APPLY,")
    print("  not `unknown` -- so an empty cell here is the axis being absent, never a rate of 0.")
    print("  It stays because it is the axis the SVA decoy fix is measured in and a future arm")
    print("  could populate it.")
    print_strata(source_strata, 40)

    # which side the support came from: sec 8's orientation question, as data
    if not a.one_sided:
        sides = Counter()
        for res in results.values():
            for s, k in res.sides(mid[1], mid[2]).items():
                sides[s] += k
        print()
        print("  support by source-window side: {0}".format(dict(sides)))
        print("  D = [Y_end, Y_end+W], U = its upstream mirror. sec 8 asks whether a source's")
        print("  POS is a 3' end or an insertion point; a lopsided split here answers it from")
        print("  data. A near-even split means the window is catching something orientation-")
        print("  independent and the claim should be loosened, not the parameter tuned.")

    # --- per-entry output --------------------------------------------------------------------
    out_tsv = os.path.join(a.out, "verdicts.tsv")
    with open(out_tsv, "w") as fout:
        fout.write("# per-entry verdicts at N={0} f_cns={1} f_src={2} cov>={3} rep_div<{4}.\n"
                   "# THE SWEEP is the result; this file is one operating point's detail, and\n"
                   "# `rep_min_div` is here so a reader can see which rows the last column of\n"
                   "# the sweep moves.\n".format(*mid))
        # `shape` and `sidedness` are written for a HUMAN reader. `--baseline` re-derives both
        # from `td_src_raw` and `src_start` rather than reading these columns, which is what lets
        # a baseline frozen before this axis existed be compared without a migration.
        fout.write("\t".join(("family", "ins_chrom", "ins_pos", "call_class", "source_class",
                              "shape", "sidedness",
                              "src_chrom", "src_start", "src_end", "arms", "depth",
                              "rep_min_div", "n_clip_support", "n_disc_support", "verdict", "why",
                              "td_src_raw")) + "\n")
        for res in results.values():
            r = res.row
            v, why = res.verdict(*mid)
            fout.write("\t".join(str(x) for x in (
                r.family, r.chrom, r.pos, r.call_class, r.source_class or "-",
                r.shape, r.sidedness,
                r.src_chrom, r.src_start, r.src_end, ",".join(r.arms), res.depth,
                "untested" if res.rep_untested else
                ("none" if res.rep_div is None else "{0:.1f}".format(res.rep_div)),
                res.n_support(mid[1], mid[2], kinds=("clip",)),
                res.n_support(mid[1], mid[2], kinds=("disc",)),
                v, why, r.raw)) + "\n")
    print()
    print("  wrote {0}".format(out_tsv))

    # --- sec 10.4: the no-regression gate ----------------------------------------------------
    if not a.baseline:
        return 0

    print()
    print("=" * 78)
    print("NO-REGRESSION GATE -- absolute VERIFIED counts per stratum (sec 10.4)")
    print("=" * 78)
    print("  baseline : {0}".format(a.baseline))
    base_point, base_shape, base_source, n_base = read_baseline_verdicts(a.baseline)
    cur_point = format_op_point(tuple(str(x) for x in mid))
    print("  {0} baseline entr(ies) at {1}".format(n_base, base_point or "AN UNSTATED point"))
    print("  {0} current  entr(ies) at {1}".format(len(rows), cur_point))

    # A comparison across operating points is not a comparison. Every VERIFIED count is a
    # function of all five thresholds -- that is decision 5 and the whole reason the sweep
    # exists -- so a baseline measured at N=3 against a run at N=2 would report regressions and
    # growth that are entirely artifacts of N. REFUSED, not warned about.
    if base_point is None:
        print()
        print("  REFUSED: the baseline file states no operating point, so there is nothing to")
        print("  check this run's midpoint against. A VERIFIED count is meaningless without one.")
        return 3
    if base_point != cur_point:
        print()
        print("  REFUSED: the baseline was measured at a DIFFERENT operating point.")
        print("    baseline {0}".format(base_point))
        print("    current  {0}".format(cur_point))
        print("  Every VERIFIED count is a function of all five thresholds (decision 5), so the")
        print("  deltas would be artifacts of the thresholds and not of the caller. Re-run with")
        print("  sweeps whose midpoints match the baseline's, or re-freeze the baseline.")
        return 3

    ok = compare_no_regression("PRIMARY  family x call class x shape x sidedness",
                               base_shape, shape_strata)
    ok = compare_no_regression("SECONDARY  family x call class x source class",
                               base_source, source_strata) and ok

    print()
    print("  " + "-" * 74)
    print("  MARGIN: PROVISIONAL (sec 10.4). This gate is not validated until an arm exists that")
    print("  is neither master nor xTea2. Three-arm agreement -- 82 of 115 entries -- is a FLOOR")
    print("  on independence and not a measure of it, because two of the three arms differ by ONE")
    print("  FLAG. Do not quote this margin as validated.")
    print("  It is also not a precision claim: it says no stratum lost verified entries, which is")
    print("  a weaker and different statement than 'precision did not fall'.")
    print("  " + "-" * 74)
    print()
    print("  NO-REGRESSION: {0}".format("PASS" if ok else "FAIL"))
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
