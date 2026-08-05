#!/usr/bin/env python3
"""Score xTea calls against the HG002 truth sets, with transductions scored separately.

DOWNSTREAM_OVERHAUL_PLAN.md §7.1 ("Still to implement (Stage 3 prep)"). Stage 3 replaces
the six-phase veto cascade with an explicit evidence model, and its gate is per-family
sensitivity/PPV against these truth sets rather than agreement with master -- so this
harness has to exist and be trusted BEFORE Stage 3 starts.

Two subcommands, because they have different inputs and very different costs:

    derive-td   truth VCFs (+ bwa, consensus, reference, rmsk) -> truth_transductions.tsv
    score       calls VCF + truth VCFs (+ that tsv)            -> the report

`derive-td` depends only on the truth set and the reference, never on the calls, so its
output is a cacheable artifact you regenerate when the reference changes and not once per
scored run. `score` runs in seconds and needs no external tool.

WHAT MAKES THIS HARDER THAN "INTERSECT TWO VCFs"
------------------------------------------------
**The truth sets carry no transduction annotation at all.** Verified, not assumed: no
`transduct*`, `TD=`, `SOURCE=`, `ORPHAN` or `SIBLING` token appears in any of the five
VCFs, and `_tprt` carries no field `_all` lacks. Transductions are the point of the
rewrite, so the truth side has to be *derived* -- from the ALT sequence, which is the full
assembled insertion. `derive-td` does that: align the ALT to the family consensus, take the
3' soft-clip residual, and require the conjunction the plan specifies.

**`_tprt` is NOT "transduction-positive".** It is a canonical-TPRT-mechanism filter and a
strict subset of `_all` (0 `_tprt` positions absent from `_all`, in either family). The
sensitivity denominator is `_all` -- 197 L1 / 1355 Alu / 90 SVA. `_tprt` is reported here
as a secondary higher-confidence stratum, which separates "we miss hard calls" from "we
miss well-formed calls", and never as the headline. Do not re-litigate this from the
filename.

**Two denominators, never conflated.** Position sensitivity is over all truth records;
transduction sensitivity is over only the derived transduction-positive subset (~51-60 for
L1), with `TRANSDUCTION_INDETERMINATE` records excluded from that denominator rather than
silently counted as negatives. Quoting one number against the wrong denominator is the
specific failure the error taxonomy below exists to prevent.

**Alu has no transductions** -- non-autonomous SINE, no internal Pol II promoter, so it
cannot 3'-read-through (§6.2). Alu is scored on position only, and `derive-td` refuses to
run on it rather than producing an all-negative table that looks like a measurement.

**These numbers are not comparable to the paper's.** Our HG002 BAM is hs38DH-aligned, not
`no_alt`; xTea filters on MAPQ and alt-aware alignment lowers MAPQ in exactly the
repeat-homologous regions that matter. The report says so on every run, because the
warning is worthless if it lives only in a plan document.
"""
import argparse
import os
import re
import subprocess
import sys
import tempfile
from collections import Counter, defaultdict

# The paper's x_cmp.py defaults to --extnd 200 and matches over `range(-EXTND, EXTND)`,
# i.e. [-200, +199] -- asymmetric by one. This harness is symmetric; at 200 bp the
# difference cannot matter, but it is why a count here may differ from x_cmp's by one call.
DEFAULT_POS_WINDOW = 200
# How close a called source locus must be to the derived one to count as the same source.
# A transduced segment is typically 200-4600 bp and its mapped start is the segment's
# start, while xTea reports the source *region* it assembled -- so this is deliberately
# looser than the position window.
DEFAULT_SRC_WINDOW = 1000
# Below this, a 3' residual is polyA plus target-site duplication, not a transduction:
# the median residual over all 197 L1 truth records is 18 bp.
DEFAULT_MIN_RESIDUAL = 100
# A residual must map back to the reference unambiguously to name a source.
DEFAULT_MIN_MAPQ = 30
DEFAULT_MIN_ALIGNED_FRAC = 0.8
# A 3' transduction's source is immediately downstream of the source element, so the
# residual must land within this distance of a full-length reference copy's 3' end.
DEFAULT_ADJACENCY = 10000
# "Full length" per family, for the adjacency test.
FULL_LENGTH_MIN = {"L1": 5000, "SVA": 1000}

TRUTH_FILES = {
    "L1": {"all": "L1_all.vcf", "tprt": "L1_tprt.vcf"},
    "Alu": {"all": "Alu_all.vcf", "tprt": "Alu_tprt.vcf"},
    # SVA ships one combined file; it is the `_all` denominator AND the tprt stratum.
    "SVA": {"all": "SVA_tprt_all.vcf", "tprt": "SVA_tprt_all.vcf"},
}
# ALT tag in an xTea VCF -> family key.
VCF_ALT_TO_FAMILY = {"INS:ME:LINE1": "L1", "INS:ME:ALU": "Alu", "INS:ME:SVA": "SVA"}
CONSENSUS = {"L1": "consensus/LINE1.fa", "Alu": "consensus/ALU.fa", "SVA": "consensus/SVA.fa"}
# RepeatMasker .out per family, for the full-length-copy adjacency test.
RMSK = {"L1": "LINE/{g}/{g}_L1.fa.out", "SVA": "SVA/{g}/{g}_SVA.out"}

TD_POSITIVE = "TRANSDUCTION_POSITIVE"
TD_NEGATIVE = "TRANSDUCTION_NEGATIVE"
TD_INDETERMINATE = "TRANSDUCTION_INDETERMINATE"


# ---------------------------------------------------------------------------------------
# parsing
# ---------------------------------------------------------------------------------------

class TruthRecord:
    __slots__ = ("chrom", "pos", "ref", "alt", "info", "tprt")

    def __init__(self, chrom, pos, ref, alt, info):
        self.chrom = chrom
        self.pos = pos
        self.ref = ref
        self.alt = alt
        self.info = info
        self.tprt = False       # set by the loader when the position is in the _tprt file

    @property
    def key(self):
        return (self.chrom, self.pos)

    def inserted_sequence(self):
        """The insertion itself: ALT minus the anchoring reference base.

        VCF represents an insertion as REF=<one base> ALT=<that base + the insert>. Strip
        it only when it is actually there, so a truth file that ever stops using the
        convention loses one base rather than silently losing its first real base.
        """
        if len(self.ref) == 1 and self.alt[:1].upper() == self.ref[:1].upper():
            return self.alt[1:]
        return self.alt


class CallRecord:
    __slots__ = ("chrom", "pos", "family", "td_src", "subtype", "info")

    def __init__(self, chrom, pos, family, td_src, subtype, info):
        self.chrom = chrom
        self.pos = pos
        self.family = family
        self.td_src = td_src      # raw TD_SRC string, or None
        self.subtype = subtype
        self.info = info


def read_truth_vcf(path):
    """The truth VCFs are HEADERLESS -- no `##`, no `#CHROM`. Tolerate a header anyway,
    since nothing guarantees the next release keeps that quirk."""
    out = []
    with open(path) as fin:
        for line in fin:
            if not line.strip() or line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 8:
                raise ValueError("{0}: expected >=8 columns, got {1}: {2!r}".format(
                    path, len(f), line[:80]))
            out.append(TruthRecord(f[0], int(f[1]), f[3], f[4], f[7]))
    return out


def load_truth_family(truth_dir, family):
    """`_all` records, each flagged with whether its position is also in `_tprt`.

    `_tprt` is a strict subset of `_all` (measured: 0 positions absent), so it is applied
    as a flag rather than loaded as a separate set. If that ever stops holding, this raises
    instead of quietly dropping the extra records.
    """
    names = TRUTH_FILES[family]
    recs = read_truth_vcf(os.path.join(truth_dir, names["all"]))
    if names["tprt"] == names["all"]:
        for r in recs:
            r.tprt = True
        return recs
    tprt_keys = {(r.chrom, r.pos) for r in read_truth_vcf(os.path.join(truth_dir, names["tprt"]))}
    all_keys = {(r.chrom, r.pos) for r in recs}
    orphans = tprt_keys - all_keys
    if orphans:
        raise ValueError(
            "{0}: {1} _tprt position(s) absent from _all (e.g. {2}). The plan's measured "
            "subset property no longer holds; the denominator logic here assumes it."
            .format(family, len(orphans), sorted(orphans)[:3]))
    for r in recs:
        r.tprt = (r.chrom, r.pos) in tprt_keys
    return recs


_INFO_RE = re.compile(r"(?:^|;)([A-Za-z_0-9]+)=([^;]*)")


def parse_info(s):
    return dict(_INFO_RE.findall(s))


def read_calls_vcf(path, families=None):
    """xTea's own VCF. Families come from the ALT tag, so this reads a merged VCF and a
    per-family one identically."""
    out = []
    with open(path) as fin:
        for line in fin:
            if line.startswith("#") or not line.strip():
                continue
            f = line.rstrip("\n").split("\t")
            alt = f[4].strip("<>")
            fam = VCF_ALT_TO_FAMILY.get(alt)
            if fam is None or (families and fam not in families):
                continue
            info = parse_info(f[7])
            out.append(CallRecord(f[0], int(f[1]), fam,
                                  info.get("TD_SRC"), info.get("SUBTYPE"), info))
    return out


def parse_td_src(s):
    """xTea's TD_SRC field -> ("region", chrom, start, end) | ("sibling", chrom, pos) | None.

    FIVE shapes, all measured on real output rather than read off the writer. The counts are
    per WGS VCF, master's NA12878 serial run (job 48825977) and xTea2's BASELINE-0 on HG002:

        not_transduction                                  master 1423 / BASELINE-0 1403
        chr8:51817625-51818695            a source region        33 /           23
        chr10:-1-108181069                ... START UNKNOWN      23 /           14
        chrY:9591830~sibling              a sibling point        52 /           50
        chr6:32546845-32546845~0~orphan   an ORPHAN               7 /            3 *
        chr3:89460826-89466856~1~89460674-89460690              0 /            3 *

    *THE LAST TWO SHAPES USED TO RETURN None HERE, WHICH IS INDISTINGUISHABLE FROM
    `not_transduction`.* That silently dropped 7 of master's 115 transduction assertions and 3
    of BASELINE-0's, and 5 more in the `--td-decoy-structural` arm. It reached the scorer as
    `TP_POS_SOURCE_MISSED` for a call that named its source perfectly well. Fixed by stripping
    the tail (below) rather than by adding two more special cases.

    THE GRAMMAR, from the two writers. `transduction_evidence.call_candidates` writes
    `"{trsdct_source}~{trsdct_l}-{trsdct_r}"` for a regular call and `"{trsdct_source}~orphan"`
    for an orphan, and `trsdct_source` is itself `chrom:start-end~rc` -- so the source region is
    the FIRST `~`-delimited field after the colon and everything after it is provenance:
    the supporting disc cluster's bounds, or the literal `orphan`. Master's
    `x_orphan_transduction.py:1203` consumes the same grammar (`chr16:35608476-35614501~1~orphan`
    in its own comment). Master's regular path omits the tail, which is why both a bare region
    and a tailed one appear in one file.

    So: split off the tail, parse the head. The **call class** (regular/sibling/orphan) lives in
    the tail and is `td_src_call_class`'s job, not this function's -- this one answers "where did
    the sequence come from", which is the same question for all three classes.

    THE PARSING HAZARD is the third line: the region delimiter is `-` and an unknown start
    is `-1`, so splitting naively on `-` yields ("", "1", "108181069"). Splitting on the
    LAST `-` is correct for both region shapes because the end coordinate is never
    negative -- which is asserted here rather than assumed, so a future `-1` end fails
    loudly instead of parsing as an absurd region.
    """
    if not s or s == "not_transduction":
        return None
    if ":" not in s:
        return None
    chrom, rest = s.split(":", 1)
    # The head is the source locus; the tail is provenance (`sibling`, `orphan`, `rc`, or the
    # supporting cluster's bounds). Everything this function returns comes from the head.
    head, _sep, tail = rest.partition("~")
    if SIBLING_TAIL in tail.split("~"):
        try:
            return ("sibling", chrom, int(head))
        except ValueError:
            return None
    if "-" not in head:
        return None
    start_s, end_s = head.rsplit("-", 1)
    if start_s.endswith("-"):
        # `100--1`: the last `-` was the end's sign, not the delimiter. Splitting on the
        # last `-` is only correct while the END is non-negative, which is what was
        # measured. Refuse rather than return a nonsense region -- a silently wrong source
        # locus shows up downstream as a plausible TP_POS_WRONG_SOURCE rate, not as an error.
        #
        # THIS IS NOW SOURCED FROM THE WRITER AND NOT ONLY FROM THE MEASUREMENT. The bare
        # shape's writer takes `i_end = max(src_pos_start, src_pos_end)`
        # (`x_orphan_transduction.py:359-363`), and a record only reaches it if `b_consistent`
        # (`:596-619`), which requires at least one side's cluster to have FORMED --
        # `_position_of_same_chrm_form_cluster` (`:1086-1096`) returns the median of a non-empty
        # list, so a formed side is never -1. Hence `max(...) >= 0` always, and this branch is
        # unreachable for that writer rather than merely unobserved. Kept anyway: it is the
        # tripwire for a THIRD writer, and its message says what to re-measure.
        raise ValueError(
            "TD_SRC {0!r} has a negative END coordinate. This parser splits on the last "
            "'-' precisely because only the START was ever -1; that assumption is now "
            "broken and the region shapes need re-measuring.".format(s))
    try:
        start, end = int(start_s), int(end_s)
    except ValueError:
        return None
    return ("region", chrom, start, end)


# The two provenance words that name a call class. Both are literals in xTea's own writers
# (`transduction_evidence.call_candidates`) and consumers (`x_orphan_transduction.py:1180`,
# `:1203`); `global_values.SIBLING_LABEL` is the same string on the caller side.
SIBLING_TAIL = "sibling"
ORPHAN_TAIL = "orphan"

CALL_CLASS_REGULAR = "regular"
CALL_CLASS_SIBLING = "sibling"
CALL_CLASS_ORPHAN = "orphan"

SRC_CLASS_REFERENCE = "reference_copy"
SRC_CLASS_POLYMORPHIC = "polymorphic"
SRC_CLASS_UNKNOWN_START = "unknown_start"
# The fourth value, and it is not a subtype of the other three -- see `td_src_source_class`.
SRC_CLASS_DISC_SPAN = "disc_span"

SHAPE_NONE = "none"          # not_transduction, or unparseable
SHAPE_BARE = "bare"          # chrom:START-END, no tail -- the SIBLING module's writer
SHAPE_TAILED = "tailed"      # chrom:START-END~rc~tail -- the TRANSDUCTION step's writer
SHAPE_SIBLING_PT = "sibling_point"   # chrom:POS~sibling


def td_src_shape(s):
    """Which WRITER produced this TD_SRC, which decides what its coordinates mean.

    THIS DISTINCTION IS LOAD-BEARING AND IT IS NOT COSMETIC. There are two writers and they put
    different things in the same field:

      * `x_transduction.py:1252` -- the transduction step -- writes
        `s_source = "{chrm_fl_L1}:{source_start}-{source_end}~{sourc_rc}"` from a FLANK RECORD, so
        its interval is the source ELEMENT's. `call_candidates` then appends `~{l}-{r}` or
        `~orphan`, which is why this shape always has a tail. Byte-identical in master.
      * `x_orphan_transduction.py:371` -- the sibling / orphan-transduction module -- writes
        `s_src = "{src_chrm}:{i_start}-{i_end}"` with NO tail, and its own comments at `:352-353`
        say what those coordinates are: **"mid-position of left discordant pair"** and
        **"mid-position of right discordant pair"**. That is the span of the SUPPORTING EVIDENCE,
        not an element. Byte-identical in master at the same line numbers.

    So a bare shape's `start`/`end` are disc-pair midpoints. Measured consequences that made this
    visible before the code was read: bare spans are 1-2814 bp (median 28 master, 52 BASELINE-0)
    where a full-length L1 is ~6 kb; and **none of master's 33 bare intervals lie inside, or even
    within 5 kb of, any source region in `unique_trsdct_disc_only`**.

    WHAT `start == -1` MEANS, PINNED TO A LINE, and it is narrower than "a midpoint is missing".
    `-1` is `mid_pos`'s initialiser in `_form_unique_cluster`
    (`x_orphan_transduction.py:1064`), returned unchanged when that side has too few discordant
    pairs (`:1072`) or no dominant source chromosome (`:1082`). But a record with a failed side
    does not survive on that alone: `b_consistent` (`:596-619`) requires at least one side to have
    formed, and for the OTHER side it calls `_check_form_repeat_or_same_cluster`, which returns

        `True, -1`   (`:1120`)  the missing side's disc mates are PREDOMINANTLY INSIDE REPEATS
                                (`n_in_rep > n_half_disc` and the ratio clears `f_rep_ratio`)
        `True, <pos>` (`:1126`) they instead cluster with the other side -- a REAL position

    So in surviving output `start == -1` means **the second side never clustered and its absence
    was EXCUSED because its discordant mates fell inside repeat regions.** That is a stronger
    statement than "one-sided", and it predicts the verification deficit rather than merely
    accompanying it: those entries are exactly the ones whose second side is repeat-driven, and
    repeat-driven discordant support is what `verify_transduction`'s decision 4 exists to
    distrust. Measured: one-sided verifies at 56.2% against two-sided 85.2% on HG002's 3-arm pool
    (55.0% / 83.9% on master's NA12878) -- a ~30-point gap with a mechanism behind it.

    Byte-identical in master (`mid_pos = -1` at its `:1093`, the same two `-1` returns at
    `:1148`/`:1155`), so this is inherited, not an xTea2 behaviour.
    """
    if not s or s == "not_transduction" or ":" not in s:
        return SHAPE_NONE
    tail = s.split(":", 1)[1].split("~")[1:]
    if SIBLING_TAIL in tail:
        return SHAPE_SIBLING_PT
    return SHAPE_TAILED if tail else SHAPE_BARE


def td_src_call_class(s):
    """TD_SRC -> "regular" | "sibling" | "orphan" | None, the first of the gate's two strata.

    NOT cosmetic and not derivable from the source locus. An orphan transduction site is where
    transduced sequence landed WITHOUT the element body, so it is not an MEI insertion at all --
    BASELINE-0's L1 orphan calls are 0/37 on-truth against GIAB (plan sec 1). A sibling's source
    element is by definition absent from the reference. Scoring the three together hides a
    category error, which is why the gate spec strata by this.
    """
    if not s or s == "not_transduction" or ":" not in s:
        return None
    tail = s.split(":", 1)[1].split("~")[1:]
    if SIBLING_TAIL in tail:
        return CALL_CLASS_SIBLING
    if ORPHAN_TAIL in tail:
        return CALL_CLASS_ORPHAN
    return CALL_CLASS_REGULAR


def td_src_source_class(s):
    """TD_SRC -> "reference_copy" | "polymorphic" | "unknown_start" | "disc_span" | None.

    The gate's second stratum, and the one the SVA decoy defect is told in. Takes the RAW field,
    not `parse_td_src`'s result, because the classification depends on which writer produced it
    (`td_src_shape`) and the parsed tuple has thrown that away.

    THE RULE, and the scope it is valid in. `gnrt_flank_regions_of_polymerphic_insertions` writes a
    polymorphic source's flank with BOTH coordinates equal to `ins_pos`, so `start == end` means
    polymorphic and `start < end` means a reference copy. **That is a statement about FLANK RECORDS,
    so it is only valid for the shape whose interval comes from one** -- `SHAPE_TAILED`, written by
    the transduction step at `x_transduction.py:1252`.

    **A BARE shape gets `disc_span`, NOT `reference_copy`.** Its coordinates are the mid-positions
    of the left and right discordant pairs (`x_orphan_transduction.py:352-353`, its own comments),
    so `end - start` is how far apart two disc-pair midpoints are and `start < end` says nothing
    about whether the source is a reference copy. Reading it as one is a category error that
    silently relabels the majority of non-sibling entries: 33 of master's 36 NA12878 "L1 regular
    reference_copy" and 23 of BASELINE-0's 24 are this shape. `disc_span` is **not a subtype** of
    the other three and must not be folded into any of them or into None.

    `start == -1` is a source whose start the caller could not name -- kept separate because the
    verifier needs `end` and only `end`. A sibling POINT has no interval at all, so it has no
    source class: this returns None for it, and None is not "unknown_start".

    THIS AXIS IS INDEPENDENT OF THE CALL CLASS, measured: master's 4 SVA orphans are all
    polymorphic (`chr6:32546845-32546845~0~orphan`) while its 3 L1 orphans are reference copies
    spanning ~6 kb (`chr4:79966908-79972933~0~orphan`). An earlier reading of the gate spec had
    "orphan" implying polymorphic; it does not.
    """
    parsed = parse_td_src(s)
    if parsed is None or parsed[0] != "region":
        return None
    if td_src_shape(s) == SHAPE_BARE:
        return SRC_CLASS_DISC_SPAN
    _, _chrom, start, end = parsed
    if start < 0:
        return SRC_CLASS_UNKNOWN_START
    if start == end:
        return SRC_CLASS_POLYMORPHIC
    return SRC_CLASS_REFERENCE


def td_src_point(parsed):
    """One comparable coordinate for a called source. For a region with an unknown start
    (-1) that is the end, which is the only coordinate actually known."""
    if parsed is None:
        return None
    if parsed[0] == "sibling":
        return (parsed[1], parsed[2])
    _, chrom, start, end = parsed
    return (chrom, end if start < 0 else start)


# ---------------------------------------------------------------------------------------
# derive-td: the truth side's transduction annotation, built from the ALT sequence
# ---------------------------------------------------------------------------------------

def revcomp(s):
    return s.translate(str.maketrans("ACGTNacgtn", "TGCANtgcan"))[::-1]


_CIGAR_RE = re.compile(r"(\d+)([MIDNSHP=X])")


def parse_cigar(c):
    return [(int(n), op) for n, op in _CIGAR_RE.findall(c)]


def trim_polya(seq):
    """Drop a trailing polyA/polyT run. A 3' transduction's residual is
    [transduced genomic sequence][polyA], and the polyA is not part of the source."""
    i = len(seq)
    while i > 0 and seq[i - 1] in "Aa":
        i -= 1
    if i < len(seq):
        return seq[:i]
    i = len(seq)
    while i > 0 and seq[i - 1] in "Tt":
        i -= 1
    return seq[:i]


def run_bwa(bwa, ref, query_fa, threads, out_sam):
    with open(out_sam, "w") as fout:
        subprocess.run([bwa, "mem", "-t", str(threads), ref, query_fa],
                       stdout=fout, stderr=subprocess.DEVNULL, check=True)
    return out_sam


def primary_alignments(sam_path):
    """qname -> (flag, rname, pos, mapq, cigar, seq) for primary records only."""
    best = {}
    with open(sam_path) as fin:
        for line in fin:
            if line.startswith("@"):
                continue
            f = line.rstrip("\n").split("\t")
            flag = int(f[1])
            if flag & 0x100 or flag & 0x800:      # secondary / supplementary
                continue
            best[f[0]] = (flag, f[2], int(f[3]), int(f[4]), f[5], f[9])
    return best


def all_alignments(sam_path):
    """qname -> [(flag, rname, pos, mapq, cigar)], primary AND supplementary.

    Supplementary records matter here for one specific reason: bwa splits a long ALT
    rather than soft-clipping it whenever the far end also hits the consensus -- a 5'
    inversion, the commonest non-linear L1 structure, does exactly that. Reading only the
    primary's clip counts such a second consensus hit as "non-consensus residual", which
    is precisely what it is not. 38 of 189 aligned L1 truth ALTs carry one.
    """
    out = defaultdict(list)
    with open(sam_path) as fin:
        for line in fin:
            if line.startswith("@"):
                continue
            f = line.rstrip("\n").split("\t")
            flag = int(f[1])
            if flag & 0x100:                       # secondary: a duplicate placement
                continue
            if flag & 0x4:
                continue
            out[f[0]].append((flag, f[2], int(f[3]), int(f[4]), f[5]))
    return out


def query_span_in_alt(flag, cigar, alt_len):
    """The [start, end) of the ALT that this record aligns, in ALT coordinates.

    Clips at the ends give the offset -- counting H as well as S, because bwa hard-clips
    supplementary records. A reverse record's SEQ is revcomp(ALT), so its query
    coordinates mirror: i -> alt_len-1-i.
    """
    ops = parse_cigar(cigar)
    if not ops:
        return None
    lead = ops[0][0] if ops[0][1] in "SH" else 0
    trail = ops[-1][0] if len(ops) > 1 and ops[-1][1] in "SH" else 0
    aligned = sum(n for n, op in ops if op in "MI=X")
    if lead + aligned + trail != alt_len:
        # Only possible if the record belongs to a different sequence than we were told.
        return None
    if flag & 0x10:
        return (alt_len - (lead + aligned), alt_len - lead)
    return (lead, lead + aligned)


def three_prime_residual(records, alt_seq):
    """The part of the insertion at its 3' end that NO consensus alignment explains.

    Two things this gets right that "take the primary's trailing soft clip" does not:

      * **Supplementary hits count as explained.** The residual is the ALT minus the union
        of every alignment's query span, not minus the primary's alone.
      * **3' means the ELEMENT's 3', not the ALT's right end.** The ALT is written in
        reference-forward orientation, so a minus-strand insertion has its 3' end on the
        LEFT. The primary's strand decides which end, and the returned sequence is put
        into the element's own 5'->3' orientation so a downstream polyA trim works either
        way.

    Returns "" when the record is unmapped or the 3' end is fully explained.
    """
    if not records:
        return ""
    alt_len = len(alt_seq)
    spans = [s for s in (query_span_in_alt(f, c, alt_len) for f, _r, _p, _q, c in records)
             if s is not None]
    if not spans:
        return ""
    primary = records[0]
    for rec in records:
        if not rec[0] & 0x800:
            primary = rec
            break
    lo = min(s[0] for s in spans)
    hi = max(s[1] for s in spans)
    if primary[0] & 0x10:
        # element 3' is the ALT's left end; return it 5'->3' along the element
        return revcomp(alt_seq[:lo])
    return alt_seq[hi:]


def load_full_length_copies(rmsk_path, min_len):
    """chrom -> sorted list of (start, end) for reference copies at least min_len long.

    RepeatMasker .out layout, matching x_annotation's parser: fields[4]=chrom,
    fields[5]=start, fields[6]=end, with two header lines and a blank third.
    """
    by_chrom = defaultdict(list)
    with open(rmsk_path) as fin:
        for line in fin:
            f = line.split()
            if len(f) < 7:
                continue
            try:
                start, end = int(f[5]), int(f[6])
            except ValueError:
                continue                      # the two header lines
            if end - start >= min_len:
                by_chrom[f[4]].append((start, end))
    for c in by_chrom:
        by_chrom[c].sort()
    return by_chrom


def contains(copies, chrom, pos):
    """Does `pos` fall INSIDE any interval? Same linear scan as the adjacency test;
    the intervals are sorted so the walk can stop once one starts past pos."""
    for (start, end) in copies.get(chrom, ()):
        if start <= pos <= end:
            return True
        if start > pos:
            break
    return False


def adjacent_full_length_copy(copies, chrom, pos, adjacency):
    """Is `pos` within `adjacency` bp downstream of a full-length copy's 3' end?

    Checked on both sides: the reference copy may be on either strand, so "downstream of
    the source element" is `pos` shortly after the copy's end OR shortly before its start.
    Linear scan over one chromosome's copies -- a few tens of thousands of intervals and a
    few hundred queries, so an index would be premature.
    """
    for (start, end) in copies.get(chrom, ()):
        if start - adjacency <= pos <= end + adjacency:
            return (start, end)
        if start > pos + adjacency:
            break
    return None


def cmd_derive_td(a):
    if a.family == "Alu":
        sys.exit("Alu is non-autonomous and cannot 3'-read-through: there are no Alu "
                 "transductions to derive. Score Alu on position only.")
    truth = load_truth_family(a.truth_dir, a.family)
    cns = os.path.join(a.rep_lib_dir, CONSENSUS[a.family])
    rmsk = os.path.join(a.rep_lib_dir, RMSK[a.family].format(g=a.genome))
    for p in (cns, rmsk, a.reference):
        if not os.path.isfile(p):
            sys.exit("missing input: " + p)
    min_full = FULL_LENGTH_MIN[a.family]

    workdir = a.workdir or tempfile.mkdtemp(prefix="derive_td_")
    os.makedirs(workdir, exist_ok=True)
    log = lambda m: print("[derive-td] " + m, flush=True)
    log("{0}: {1} truth records; consensus {2}; rmsk {3}".format(
        a.family, len(truth), os.path.basename(cns), os.path.basename(rmsk)))

    # 1. every truth ALT against the family consensus
    alt_fa = os.path.join(workdir, "alt.fa")
    with open(alt_fa, "w") as fout:
        for r in truth:
            fout.write(">{0}:{1}\n{2}\n".format(r.chrom, r.pos, r.inserted_sequence()))
    sam1 = run_bwa(a.bwa, cns, alt_fa, a.threads, os.path.join(workdir, "alt_vs_cns.sam"))
    aln = all_alignments(sam1)
    n_supp = sum(1 for v in aln.values() if len(v) > 1)
    log("aligned {0}/{1} ALTs to the consensus ({2} with supplementary hits, whose spans "
        "count as EXPLAINED)".format(len(aln), len(truth), n_supp))

    # 2. the 3' residual, polyA trimmed
    residuals = {}
    unmapped = 0
    for r in truth:
        key = "{0}:{1}".format(r.chrom, r.pos)
        recs = aln.get(key)
        if not recs:
            unmapped += 1
            residuals[key] = None
            continue
        residuals[key] = trim_polya(three_prime_residual(recs, r.inserted_sequence()))
    lens = sorted(len(v) for v in residuals.values() if v is not None)
    if lens:
        log("3' residual length: median {0} bp, >={1} bp: {2}, >=200 bp: {3}, max {4}".format(
            lens[len(lens) // 2], a.min_residual,
            sum(1 for x in lens if x >= a.min_residual),
            sum(1 for x in lens if x >= 200), lens[-1]))
        log("  (this is the 3' residual the plan's criterion asks for. The plan's own "
            "'>=100 bp: 60, >=200 bp: 51' figures are TOTAL soft clip, both ends, and are "
            "not this number -- see the note in derive_truth_td.sbatch.)")

    # 3. the long residuals back to the reference genome
    res_fa = os.path.join(workdir, "residual.fa")
    n_cand = 0
    with open(res_fa, "w") as fout:
        for key, res in residuals.items():
            if res and len(res) >= a.min_residual:
                fout.write(">{0}\n{1}\n".format(key, res))
                n_cand += 1
    log("{0} residual(s) at or above {1} bp -> reference".format(n_cand, a.min_residual))
    res_aln = {}
    if n_cand:
        sam2 = run_bwa(a.bwa, a.reference, res_fa, a.threads,
                       os.path.join(workdir, "residual_vs_ref.sam"))
        res_aln = primary_alignments(sam2)

    # 4. the conjunction
    copies = load_full_length_copies(rmsk, min_full)
    log("{0} reference {1} copies >= {2} bp over {3} contigs".format(
        sum(len(v) for v in copies.values()), a.family, min_full, len(copies)))
    # ...and EVERY copy, at any length, for the divergent-body test below.
    all_copies = load_full_length_copies(rmsk, 0)

    rows, tally = [], Counter()
    for r in truth:
        key = "{0}:{1}".format(r.chrom, r.pos)
        res = residuals[key]
        reason, status, src_chrom, src_pos, src_copy = "", TD_NEGATIVE, "", -1, ""
        if res is None:
            status, reason = TD_INDETERMINATE, "alt_unmapped_to_consensus"
        elif len(res) < a.min_residual:
            status, reason = TD_NEGATIVE, "residual_{0}bp_below_{1}".format(
                len(res), a.min_residual)
        else:
            rec = res_aln.get(key)
            if rec is None or rec[0] & 0x4:
                status, reason = TD_INDETERMINATE, "residual_unmapped_to_reference"
            elif rec[3] < a.min_mapq:
                status, reason = TD_INDETERMINATE, "residual_mapq_{0}".format(rec[3])
            else:
                flag, rname, ppos, _mapq, cigar, _seq = rec
                matched = sum(n for n, op in parse_cigar(cigar) if op in "M=X")
                frac = matched / float(len(res)) if res else 0.0
                if frac < a.min_aligned_frac:
                    status = TD_INDETERMINATE
                    reason = "residual_aligned_frac_{0:.2f}".format(frac)
                elif contains(all_copies, rname, ppos):
                    # A residual that lands INSIDE a reference copy of its own family is
                    # far more likely to be divergent element body than transduced
                    # sequence: consensus/LINE1.fa has 4 L1HS-centric records, so an older
                    # subfamily's body partly fails to align and surfaces as "3' residual".
                    # Measured, not supposed: 11 of the 20 L1 residuals >=100 bp land in a
                    # reference L1, and EVERY one is an older subfamily -- L1PA6/7/8 at
                    # 4-7% divergence, L1MB/MD/ME at 16-25%. INDETERMINATE, not negative:
                    # the ALT genuinely carries non-consensus sequence, we just cannot tell
                    # what it is.
                    status = TD_INDETERMINATE
                    reason = "residual_inside_reference_{0}".format(a.family)
                    src_chrom, src_pos = rname, ppos
                else:
                    hit = adjacent_full_length_copy(copies, rname, ppos, a.adjacency)
                    if hit is None:
                        status = TD_INDETERMINATE
                        reason = "no_full_length_{0}_within_{1}bp".format(a.family, a.adjacency)
                        src_chrom, src_pos = rname, ppos
                    else:
                        status, reason = TD_POSITIVE, "residual_{0}bp".format(len(res))
                        src_chrom, src_pos = rname, ppos
                        src_copy = "{0}:{1}-{2}".format(rname, hit[0], hit[1])
        tally[status] += 1
        rows.append((r.chrom, r.pos, a.family, "1" if r.tprt else "0", status,
                     0 if res is None else len(res), src_chrom, src_pos, src_copy, reason))

    with open(a.out, "w") as fout:
        fout.write("#chrom\tpos\tfamily\tin_tprt\ttd_status\tresidual_bp\t"
                   "src_chrom\tsrc_pos\tsrc_copy\treason\n")
        for row in rows:
            fout.write("\t".join(str(x) for x in row) + "\n")
    log("wrote {0}".format(a.out))
    for k in (TD_POSITIVE, TD_NEGATIVE, TD_INDETERMINATE):
        log("  {0:28s} {1}".format(k, tally[k]))
    if unmapped:
        log("  ({0} ALT(s) did not align to the consensus at all)".format(unmapped))
    log("Reminder: INDETERMINATE records leave the transduction denominator. They are NOT "
        "negatives, and a run whose INDETERMINATE count rivals its POSITIVE count is not "
        "measuring transduction sensitivity -- it is measuring this derivation's limits.")


def load_td_table(path):
    """(chrom, pos) -> (status, src_chrom, src_pos, in_tprt)."""
    out = {}
    with open(path) as fin:
        for line in fin:
            if line.startswith("#") or not line.strip():
                continue
            f = line.rstrip("\n").split("\t")
            out[(f[0], int(f[1]))] = (f[4], f[6], int(f[7]), f[3] == "1")
    return out


# ---------------------------------------------------------------------------------------
# score
# ---------------------------------------------------------------------------------------

def match_one_to_one(truth, calls, window):
    """Greedy nearest-position matching, ONE call per truth record and vice versa.

    x_cmp.py instead asks, per truth record, "is any call within the window", which lets a
    single call satisfy several nearby truth records and counts them all as hits. At a
    200 bp window over clustered MEI loci that is not hypothetical. Matching one-to-one by
    increasing distance costs nothing here and makes TP + FN = |truth| and TP + FP =
    |calls| hold exactly, which is what makes the two rates comparable.
    """
    by_chrom = defaultdict(list)
    for i, c in enumerate(calls):
        by_chrom[c.chrom].append((c.pos, i))
    for c in by_chrom:
        by_chrom[c].sort()

    pairs = []
    for ti, t in enumerate(truth):
        for pos, ci in by_chrom.get(t.chrom, ()):
            d = abs(pos - t.pos)
            if d <= window:
                pairs.append((d, ti, ci))
    pairs.sort()
    t_taken, c_taken, matched = set(), set(), {}
    for d, ti, ci in pairs:
        if ti in t_taken or ci in c_taken:
            continue
        t_taken.add(ti)
        c_taken.add(ci)
        matched[ti] = (ci, d)
    return matched


def _rate(num, den):
    return float("nan") if den == 0 else 100.0 * num / den


def cmd_score(a):
    families = [f.strip() for f in a.families.split(",") if f.strip()]
    calls_all = read_calls_vcf(a.calls, families)
    td_table = load_td_table(a.td_table) if a.td_table else {}

    print("=" * 78)
    print("xTea vs HG002 truth set")
    print("  calls      : {0} ({1} record(s) in {2})".format(
        a.calls, len(calls_all), ",".join(families)))
    print("  truth      : {0}".format(a.truth_dir))
    print("  windows    : position +-{0} bp, source +-{1} bp".format(a.window, a.src_window))
    print("  truth TD   : {0}".format(a.td_table or "NOT SUPPLIED -- position metrics only"))
    print("  NOT COMPARABLE TO THE PAPER'S NUMBERS: our HG002 BAM is hs38DH-aligned, not")
    print("  no_alt, and xTea filters on MAPQ. Never quote these side by side without")
    print("  saying so (plan section 7.1).")
    print("=" * 78)

    grand = Counter()
    for fam in families:
        truth = load_truth_family(a.truth_dir, fam)
        calls = [c for c in calls_all if c.family == fam]
        matched = match_one_to_one(truth, calls, a.window)
        tp = len(matched)
        fn = len(truth) - tp
        fp = len(calls) - tp
        grand["tp"] += tp
        grand["fn"] += fn
        grand["fp"] += fp

        print()
        print("--- {0} ---".format(fam))
        print("  truth {0}   calls {1}".format(len(truth), len(calls)))
        print("  POSITION   TP {0}  FN {1}  FP {2}   sensitivity {3:.1f}%  precision {4:.1f}%"
              .format(tp, fn, fp, _rate(tp, len(truth)), _rate(tp, len(calls))))

        # An ORPHAN transduction site is where transduced sequence landed WITHOUT the
        # element body. It is not an MEI insertion, and the HG002 truth set enumerates MEI
        # insertions -- so every orphan call is an off-truth "false positive" BY
        # CONSTRUCTION, whatever the caller did. On BASELINE-0 this is not a small effect:
        # 37 of L1's 48 position-FPs are orphan calls, all 37 off-truth, and L1's precision
        # reads 74.5% with them and 92.7% without. Reported as its own stratum rather than
        # silently folded in either direction -- excluding them without saying so would
        # flatter the caller just as much as including them maligns it.
        orphan_idx = [i for i, c in enumerate(calls)
                      if c.subtype == "orphan_or_sibling_transduction"]
        if orphan_idx:
            matched_calls = {v[0] for v in matched.values()}
            orphan_on = sum(1 for i in orphan_idx if i in matched_calls)
            rest = len(calls) - len(orphan_idx)
            print("    of which {0} are orphan transduction SITES ({1} on-truth): these are "
                  "not MEI insertions".format(len(orphan_idx), orphan_on))
            print("    and the truth set does not enumerate them, so they cannot match. "
                  "Excluding them: precision {0:.1f}% over {1} call(s)."
                  .format(_rate(tp - orphan_on, rest), rest))
        # The _tprt stratum: same calls, denominator restricted to the well-formed subset.
        tprt_idx = [i for i, t in enumerate(truth) if t.tprt]
        if tprt_idx and len(tprt_idx) != len(truth):
            tprt_tp = sum(1 for i in tprt_idx if i in matched)
            print("  ...of which the _tprt stratum ({0} records): sensitivity {1:.1f}%  "
                  "(a canonical-TPRT-mechanism subset, NOT transduction-positive)"
                  .format(len(tprt_idx), _rate(tprt_tp, len(tprt_idx))))

        if fam == "Alu":
            print("  TRANSDUCTION  n/a -- Alu is non-autonomous and cannot 3'-read-through.")
            n_asserted = sum(1 for c in calls if parse_td_src(c.td_src))
            if n_asserted:
                print("    NOTE: {0} Alu call(s) nevertheless assert a TD_SRC. That is a "
                      "finding about the caller, not about the truth set.".format(n_asserted))
            continue
        if not td_table:
            print("  TRANSDUCTION  skipped (no --td-table; run `derive-td` first)")
            continue

        tax = Counter()
        details = []
        n_pos_truth = n_indet = 0
        for ti, t in enumerate(truth):
            entry = td_table.get((t.chrom, t.pos))
            if entry is None:
                continue
            status, src_chrom, src_pos, _ = entry
            if status == TD_INDETERMINATE:
                n_indet += 1
                continue        # excluded from the denominator, never counted as negative
            if status == TD_POSITIVE:
                n_pos_truth += 1
            hit = matched.get(ti)
            if hit is None:
                if status == TD_POSITIVE:
                    tax["FN"] += 1
                continue
            call = calls[hit[0]]
            parsed = parse_td_src(call.td_src)
            called_pt = td_src_point(parsed)
            if status == TD_POSITIVE:
                if called_pt is None:
                    tax["TP_POS_SOURCE_MISSED"] += 1
                    details.append(("TP_POS_SOURCE_MISSED", t.chrom, t.pos,
                                    "{0}:{1}".format(src_chrom, src_pos), "not_transduction"))
                elif (called_pt[0] == src_chrom
                      and abs(called_pt[1] - src_pos) <= a.src_window):
                    tax["TP_POS_AND_SOURCE"] += 1
                else:
                    tax["TP_POS_WRONG_SOURCE"] += 1
                    details.append(("TP_POS_WRONG_SOURCE", t.chrom, t.pos,
                                    "{0}:{1}".format(src_chrom, src_pos), call.td_src))
            else:   # TD_NEGATIVE
                if parsed is None:
                    pass
                elif parsed[0] == "sibling":
                    # NOT counted as a false positive. A sibling transduction's source
                    # element is by definition ABSENT from the reference, so the adjacency
                    # test at the heart of this derivation -- "the residual sits next to a
                    # full-length reference copy" -- cannot be satisfied even when the call
                    # is right. Folding these into FP_TRANSDUCTION would charge the caller
                    # for the derivation's blind spot. Reported, not scored.
                    tax["UNADJUDICATED_SIBLING"] += 1
                else:
                    tax["FP_TRANSDUCTION"] += 1
                    details.append(("FP_TRANSDUCTION", t.chrom, t.pos, "none", call.td_src))

        td_tp = tax["TP_POS_AND_SOURCE"]
        print("  TRANSDUCTION  denominator {0} derived-positive truth record(s); "
              "{1} INDETERMINATE excluded".format(n_pos_truth, n_indet))
        for k in ("TP_POS_AND_SOURCE", "TP_POS_WRONG_SOURCE", "TP_POS_SOURCE_MISSED",
                  "FP_TRANSDUCTION", "FN"):
            print("    {0:24s} {1}".format(k, tax[k]))
        if tax["UNADJUDICATED_SIBLING"]:
            print("    {0:24s} {1}   (sibling sources are non-reference by definition, so "
                  "this derivation can neither confirm nor refute them -- NOT scored)"
                  .format("UNADJUDICATED_SIBLING", tax["UNADJUDICATED_SIBLING"]))
        print("    source-resolved sensitivity {0:.1f}%   position-resolved {1:.1f}%".format(
            _rate(td_tp, n_pos_truth),
            _rate(td_tp + tax["TP_POS_WRONG_SOURCE"] + tax["TP_POS_SOURCE_MISSED"],
                  n_pos_truth)))
        if n_indet >= max(1, n_pos_truth):
            print("    WARNING: INDETERMINATE ({0}) is not smaller than POSITIVE ({1}). "
                  "These rates describe the derivation's reach as much as the caller's."
                  .format(n_indet, n_pos_truth))
        if a.details and details:
            print("    -- non-clean outcomes --")
            for row in details[:a.details]:
                print("      {0:24s} {1}:{2}  truth_src={3}  called={4}".format(*row))
            if len(details) > a.details:
                print("      ... {0} more".format(len(details) - a.details))

    print()
    print("--- all families ---")
    print("  POSITION   TP {0}  FN {1}  FP {2}   sensitivity {3:.1f}%  precision {4:.1f}%"
          .format(grand["tp"], grand["fn"], grand["fp"],
                  _rate(grand["tp"], grand["tp"] + grand["fn"]),
                  _rate(grand["tp"], grand["tp"] + grand["fp"])))


# ---------------------------------------------------------------------------------------

def main():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = p.add_subparsers(dest="cmd", required=True)

    d = sub.add_parser("derive-td", help="truth ALT sequences -> transduction status table")
    d.add_argument("--truth-dir", required=True)
    d.add_argument("--family", required=True, choices=["L1", "SVA", "Alu"])
    d.add_argument("--reference", required=True, help="bwa-indexed genome FASTA")
    d.add_argument("--rep-lib-dir", required=True)
    d.add_argument("--genome", default="hg38")
    d.add_argument("--bwa", default="bwa")
    d.add_argument("--threads", type=int, default=4)
    d.add_argument("--workdir", default="")
    d.add_argument("--out", required=True)
    d.add_argument("--min-residual", type=int, default=DEFAULT_MIN_RESIDUAL)
    d.add_argument("--min-mapq", type=int, default=DEFAULT_MIN_MAPQ)
    d.add_argument("--min-aligned-frac", type=float, default=DEFAULT_MIN_ALIGNED_FRAC)
    d.add_argument("--adjacency", type=int, default=DEFAULT_ADJACENCY)
    d.set_defaults(func=cmd_derive_td)

    s = sub.add_parser("score", help="score a calls VCF against the truth sets")
    s.add_argument("--calls", required=True, help="xTea VCF (merged or per-family)")
    s.add_argument("--truth-dir", required=True)
    s.add_argument("--families", default="L1,Alu,SVA")
    s.add_argument("--td-table", default="", help="output of `derive-td` (may be repeated "
                                                  "families concatenated)")
    s.add_argument("--window", type=int, default=DEFAULT_POS_WINDOW)
    s.add_argument("--src-window", type=int, default=DEFAULT_SRC_WINDOW)
    s.add_argument("--details", type=int, default=20,
                   help="print up to N non-clean transduction outcomes (0 to suppress)")
    s.set_defaults(func=cmd_score)

    a = p.parse_args()
    a.func(a)


if __name__ == "__main__":
    main()
