"""Phase C move A.1/A.2: explicit locus table for the per-type -> single-VCF merge.

Types load in order L1 -> SVA -> Alu (HERV unused -- xTea2 runs L1/Alu/SVA
only); a call within +-MERGE_WINDOW bp of an already-placed call is the same
locus, and the two calls are resolved by EVIDENCE (see evidence_score and
LocusTable._beats):

  1. higher evidence score wins: clusters on the family consensus (0-4)
     + 1 if poly-A reads + 1 if a TSD;
  2. tie -> a call that is NOT orphan_or_sibling_transduction wins;
  3. tie -> family priority L1 > SVA > Alu (the old load-order rule).

This replaced the original rule (x_vcf_merger's intended one: keep the
higher-priority family unless the incoming call carries two_side_tprt_both),
which in HG002 dropped true SVA insertions in favour of transduction-only L1
calls at the same spot (chr1:88507142, chr14:20707961). Losing calls are no
longer silently discarded: every winner keeps the calls it beat (Locus.beaten),
and assign_and_write can record all of them (merged/all_calls.tsv/.bed, and a
4th type-map column that becomes the VCF's OTHER_FAMILIES tag).

This intentionally does NOT reproduce a bug in master's merger (documented as
B5 in docs/clipdisc_contract_spec.md sec 4): master keys every locus by its own
exact VCF position in one dict, so a losing (non-TPRT-override) conflict whose
position differs from the winner's is written anyway (a missing `continue` on
the "keep existing" branch) -- a spurious duplicate call. When the two
positions are exactly equal, that same write instead silently clobbers the
dict key, so whichever type loaded LAST (always Alu) wins regardless of
priority or TPRT flags. This module keeps each chromosome's placed loci in a
list (not a position-keyed dict), so a resolved conflict always
replaces-or-drops the loser in place: no collision, no duplicate, independent
of whether the two positions happen to be equal.

Move A.2 moved the resolution point earlier: TWO_SIDE_TPRT_BOTH is known right
after post-filtering (x_post_filter.py), well before gene annotation,
genotyping, or VCF writing -- so this module now resolves loci from the three
types' raw `postF_hc` candidate-list files (the native xTea txt format, NOT
VCF) instead of their final VCFs. That lets gene/genotype/gVCF run ONCE over
the unified, already-deduplicated set instead of 3x + a post-hoc VCF merge.
The winning raw lines are written back out unchanged (same column format), so
gene annotation and genotyping need no changes at all; only the final gVCF
step needs each record's resolved type, via the `sf_type_map` sidecar this
module also writes.
"""
import os
import re
import subprocess

import global_values

MERGE_WINDOW = 50
TPRT = global_values.TWO_SIDE_TPRT_BOTH
TYPE_PRIORITY = ("L1", "SVA", "Alu")
# rtype bitmask codes x_TEA_main/x_rep_type use (mirrors run_xtea.py's REP_TYPES).
TYPE_TO_RTYPE = {"L1": 1, "Alu": 2, "SVA": 4}
RTYPE_TO_TYPE = {v: k for k, v in TYPE_TO_RTYPE.items()}
_PRIORITY_RANK = {t: i for i, t in enumerate(TYPE_PRIORITY)}  # lower = higher priority

# Raw-result columns the evidence score reads (0-based; x_clip_disc_filter's
# output format, unchanged through post-filtering).
COL_TSD = 4                       # TSD length, <= 0 when none
COLS_POLYA = (9, 10)              # left / right poly-A clip reads
COLS_CNS_CLUSTERS = (19, 20, 21, 22)  # clip/disc cluster ranges on the consensus, "-1:-1" = none
COL_SUBTYPE = 32                  # e.g. two_side_tprt_both, orphan_or_sibling_transduction
ORPHAN = "orphan_or_sibling_transduction"


def _int(fields, i):
    try:
        return int(fields[i])
    except (IndexError, ValueError):
        return 0


def evidence_components(fields):
    """(n_cns_clusters 0-4, has_polya 0/1, has_tsd 0/1) for one raw-result row.

    Missing columns count as no evidence, so a short row (a test fixture, or a
    file from before a column existed) scores 0 rather than raising."""
    n_cls = sum(1 for i in COLS_CNS_CLUSTERS if i < len(fields) and fields[i] != "-1:-1")
    polya = int(sum(_int(fields, i) for i in COLS_POLYA) > 0)
    tsd = int(_int(fields, COL_TSD) > 0)
    return n_cls, polya, tsd


def evidence_score(fields):
    """Clusters on consensus (0-4) + poly-A present + TSD present: 0..6."""
    return sum(evidence_components(fields))


def is_orphan(fields):
    """The call is a transduction-only (orphan/sibling) site. This label is
    written by the transduction step, so it never appears in cns_filter output
    -- the tie-break simply never fires at the early barrier."""
    return len(fields) > COL_SUBTYPE and fields[COL_SUBTYPE] == ORPHAN


class Locus:
    """One family's claim on one locus: the row, plus the fields the
    resolution rule reads (position, evidence score, orphan flag). `tprt`
    (two_side_tprt_both) no longer decides anything; it is kept for the
    conflict log.

    `idx` is the row's 0-based line number within its own family's file. It
    exists so a caller can map a resolution decision back to the exact input
    row -- `line` cannot serve as that key, since nothing forbids a file from
    carrying two identical rows. filter_rows_to_assigned() relies on it.

    `beaten` lists the calls this one won against (carried over to whichever
    call later beats this one), so nothing that lost is lost track of.
    """
    __slots__ = ("chrm", "pos", "type", "tprt", "line", "idx",
                 "evidence", "score", "orphan", "beaten")

    def __init__(self, chrm, pos, type_, tprt, line, idx=-1):
        self.chrm = chrm
        self.pos = pos
        self.type = type_
        self.tprt = tprt
        self.line = line
        self.idx = idx
        fields = line.split()
        self.evidence = evidence_components(fields)
        self.score = sum(self.evidence)
        self.orphan = is_orphan(fields)
        self.beaten = []

    def rank_key(self):
        """Higher key wins: score, then non-orphan, then family priority."""
        return (self.score, not self.orphan, -_PRIORITY_RANK.get(self.type, len(TYPE_PRIORITY)))


def call_position(fields):
    """The exact position key x_gvcf.cvt_raw_rslt_to_gvcf uses for a raw-result
    record: min(lclip_pos, rclip_pos), falling back to the other end if that's
    <= 0. MUST stay identical to x_gvcf.py's own derivation (x_gvcf.py:204-213)
    -- this is the key that ties a locus-table resolution to the record a later
    gVCF pass looks its type up for."""
    lclip_pos = int(fields[2])
    rclip_pos = int(fields[3])
    i_begin = lclip_pos
    i_end = rclip_pos
    if lclip_pos > rclip_pos:
        i_begin = rclip_pos
        i_end = lclip_pos
    if i_begin <= 0:
        i_begin = i_end
    return i_begin


def load_rows(sf_raw, type_):
    """Parse a native raw-result file into this family's Locus claims.

    The format is whitespace-separated, no header, one candidate per line, and
    two_side_tprt_both is a literal substring of the line when present.

    STAGE-AGNOSTIC ON PURPOSE, and that is the whole of Stage 2c.2's "lift the
    rule out of _place": every raw-result file the per-family chain produces
    shares this format, so the same loader and the same resolution rule run at
    whichever point in the chain the caller chooses. The two that matter:

      candidate_disc_filtered_cns.txt   -- cns_filter's output, the EARLIEST
          point the rule can run at all, because two_side_tprt_both is written
          by that step (x_clip_disc_filter.py:1845) and does not exist before
          it. Assigning any earlier means assigning without the override that
          decides 13 of BASELINE-0's 22 conflicts.
      ...high_confident.post_filtering.txt -- postF_hc, where the merge has
          always resolved family (run_shared_tail).

    Nothing here is specific to either; see assign_rows().
    """
    rows = []
    if not os.path.isfile(sf_raw):
        return rows
    with open(sf_raw) as fin:
        for idx, line in enumerate(fin):
            fields = line.split()
            if len(fields) < 4:
                continue
            line = line.rstrip("\n")
            rows.append(Locus(fields[0], call_position(fields), type_, TPRT in line, line, idx))
    return rows


# Pre-2c name, kept so nothing that imported it privately breaks.
_load_hc_pf_rows = load_rows


class Conflict:
    """One locus two families both called, and how it was resolved.

    Recorded rather than discarded because this population IS Stage 2's gate.
    Stage 2c moves family assignment from here -- after three independent
    per-family chains -- to a single up-front decision, and its gate is "calls may
    move, but only as family-assignment flips". You cannot check that against a
    population nobody counted: without this list, a Stage 2c run that loses a call
    for an unrelated reason is indistinguishable from one that reassigned it.

    `rule` is what decided it: `score` (higher evidence score), `tie_non_orphan`
    (equal scores; the other call was orphan_or_sibling_transduction), or
    `tie_priority` (equal on both; TYPE_PRIORITY, or the already-placed call for
    two calls of the same family).
    """

    __slots__ = ("chrm", "winner_type", "winner_pos", "winner_tprt",
                 "loser_type", "loser_pos", "loser_tprt", "rule",
                 "winner_score", "loser_score")

    def __init__(self, chrm, winner, loser, rule):
        self.chrm = chrm
        self.winner_type = winner.type
        self.winner_pos = winner.pos
        self.winner_tprt = winner.tprt
        self.loser_type = loser.type
        self.loser_pos = loser.pos
        self.loser_tprt = loser.tprt
        self.rule = rule
        self.winner_score = winner.score
        self.loser_score = loser.score

    def row(self):
        return "\t".join(str(x) for x in (
            self.chrm, self.winner_type, self.winner_pos, int(self.winner_tprt),
            self.loser_type, self.loser_pos, int(self.loser_tprt), self.rule,
            self.winner_score, self.loser_score))


# New columns go at the END so scripts that read the first 8 by index keep working.
CONFLICT_HEADER = ("#chrm\twinner_type\twinner_pos\twinner_tprt\t"
                   "loser_type\tloser_pos\tloser_tprt\trule\twinner_score\tloser_score")


class LocusTable:
    """Collapses per-type raw postF_hc calls to one call per locus.

    Types load in TYPE_PRIORITY order. A call within +-`window` bp of an
    already-placed call is the same locus; the one with the higher
    Locus.rank_key() keeps it (see the module docstring), and on a full tie the
    already-placed call stays.

    Every such collision is also appended to `self.conflicts` -- see Conflict --
    and every call ever placed to `self.all_loci`, in load order.
    """

    def __init__(self, window=MERGE_WINDOW):
        self.window = window
        self._by_chrm = {}  # chrm -> [Locus, ...], insertion order
        self.conflicts = []
        self.all_loci = []

    def assign_type(self, hc_pf_paths):
        """hc_pf_paths: {type_name: raw-result path}. Returns the surviving Locus records."""
        return self.assign_rows({t: load_rows(p, t) for t, p in hc_pf_paths.items() if p})

    def assign_rows(self, rows_by_type):
        """The rule, over already-parsed claims. rows_by_type: {type_name: [Locus]}.

        Split out from assign_type so the identical rule can be applied at a
        different point in the chain without a second copy of it existing --
        Stage 2c moves the decision from postF_hc up to cns_filter, and the
        thing that must NOT differ between the two is this."""
        for type_ in TYPE_PRIORITY:
            for locus in rows_by_type.get(type_) or ():
                self._place(locus)
        return [locus for bucket in self._by_chrm.values() for locus in bucket]

    def _place(self, locus):
        self.all_loci.append(locus)
        bucket = self._by_chrm.setdefault(locus.chrm, [])
        for i, placed in enumerate(bucket):
            if abs(placed.pos - locus.pos) <= self.window:
                if locus.rank_key() > placed.rank_key():
                    winner, loser = locus, placed
                    bucket[i] = locus
                else:
                    winner, loser = placed, locus
                self.conflicts.append(Conflict(locus.chrm, winner, loser, _deciding_rule(winner, loser)))
                winner.beaten.append(loser)
                winner.beaten.extend(loser.beaten)
                loser.beaten = []
                return  # conflict resolved either way -- never duplicate
        bucket.append(locus)

    def conflict_summary(self):
        """{(winner_type, loser_type, rule): count} -- the shape Stage 2c predicts against."""
        counts = {}
        for c in self.conflicts:
            key = (c.winner_type, c.loser_type, c.rule)
            counts[key] = counts.get(key, 0) + 1
        return counts

    def winner_of(self):
        """{id(loser): winning Locus} for every call that lost."""
        return {id(b): w for bucket in self._by_chrm.values() for w in bucket for b in w.beaten}


def _deciding_rule(winner, loser):
    if winner.score != loser.score:
        return "score"
    if winner.orphan != loser.orphan:
        return "tie_non_orphan"
    return "tie_priority"


def _natural_key(chrm):
    return [int(t) if t.isdigit() else t for t in re.split(r"(\d+)", chrm)]


def other_families_tag(locus):
    """`L1:88507142,...` for the calls `locus` beat, or `.` when it beat none."""
    if not locus.beaten:
        return "."
    return ",".join(f"{b.type}:{b.pos}" for b in sorted(locus.beaten, key=lambda b: (b.pos, b.type)))


ALL_CALLS_HEADER = ("#chrm\tpos\tfamily\tscore\tcns_clusters\tpolyA\tTSD\torphan\t"
                    "status\twinner_family\twinner_pos\tevidence_line")


def write_all_calls(table, sf_tsv=None, sf_bed=None):
    """Every call that reached the merge, won or lost, with its score.

    The TSV ends with the call's full raw-result line, so a lost call can be
    inspected (or rescued) without going back to the per-family directories.
    The BED is an IGV track: one 1-bp feature per call, named e.g.
    `SVA score=5 WON`.
    """
    winner_of = table.winner_of()
    rows = sorted(table.all_loci, key=lambda l: (_natural_key(l.chrm), l.pos, l.type))
    if sf_tsv:
        with open(sf_tsv, "w") as fout:
            fout.write(ALL_CALLS_HEADER + "\n")
            for l in rows:
                w = winner_of.get(id(l))
                status, wt, wp = ("LOST", w.type, w.pos) if w else ("WON", l.type, l.pos)
                n_cls, polya, tsd = l.evidence
                fout.write("\t".join(str(x) for x in (
                    l.chrm, l.pos, l.type, l.score, n_cls, polya, tsd, int(l.orphan),
                    status, wt, wp, l.line)) + "\n")
    if sf_bed:
        with open(sf_bed, "w") as fout:
            fout.write('track name="xTea2 all calls" description="every family call at merge; '
                       'WON = kept, LOST = beaten by another call"\n')
            for l in rows:
                status = "LOST" if id(l) in winner_of else "WON"
                fout.write(f"{l.chrm}\t{max(l.pos - 1, 0)}\t{l.pos}\t{l.type} score={l.score} {status}\n")


def assign_and_write(hc_pf_paths, sf_merged_hcpf, sf_type_map, window=MERGE_WINDOW,
                     sf_conflicts=None, sf_all_calls=None, sf_all_calls_bed=None):
    """Resolve hc_pf_paths ({type: postF_hc path}) into one deduplicated raw-result
    file at sf_merged_hcpf (sorted, native hc_pf column format -- unchanged, so
    the shared gene/genotype tail needs no changes) plus a sidecar type map at
    sf_type_map (`chrm\\tpos\\trtype_int\\tother_families` per surviving locus,
    pos = call_position) for the later per-record gVCF fan-out. other_families
    is `<family>:<pos>[,...]` for the calls that locus beat, or `.`.

    sf_all_calls / sf_all_calls_bed, when given, record every call that reached
    the merge, won or lost (see write_all_calls).

    sf_conflicts, when given, additionally writes every resolved family collision
    (see Conflict). Purely diagnostic -- it does not touch the two files above --
    and it is what Stage 2c's "only family-assignment flips" gate is checked
    against.

    Returns (loci, table) so the caller can report the conflict population.
    """
    table = LocusTable(window=window)
    loci = table.assign_type(hc_pf_paths)

    sf_unsorted = sf_merged_hcpf + ".unsorted"
    with open(sf_unsorted, "w") as fout:
        for locus in loci:
            fout.write(locus.line + "\n")
    subprocess.run(["sort", "-k1,1V", "-k2,2n", "-o", sf_merged_hcpf, sf_unsorted], check=True)
    os.remove(sf_unsorted)

    with open(sf_type_map, "w") as fout:
        for locus in loci:
            fout.write(f"{locus.chrm}\t{locus.pos}\t{TYPE_TO_RTYPE[locus.type]}\t"
                       f"{other_families_tag(locus)}\n")

    if sf_conflicts:
        with open(sf_conflicts, "w") as fout:
            fout.write(CONFLICT_HEADER + "\n")
            for c in table.conflicts:
                fout.write(c.row() + "\n")
    write_all_calls(table, sf_all_calls, sf_all_calls_bed)
    return loci, table


# Companion files cns_filter writes NEXT TO its candidate list, which a later step
# opens by DERIVING the path from its own `-i` rather than being handed it.
#
# `.high_confident` is the one that exists: x_TEA_main.py:551 builds
# `sf_candidate_list + HIGH_CONFIDENT_SUFFIX` and update_high_confident_callset
# (x_transduction.py:632) then looks every one of its rows up in the transduction
# output, keyed on (fields[0], fields[1]), with no membership test. So filtering the
# candidate list alone leaves a row in the companion that no longer exists downstream,
# and the step dies with a bare `KeyError: <pos>` (job 49498881, SVA:transduction).
#
# The other two companions -- `.before_calling_transduction` and `.before_filtering`
# -- are written AND consumed inside the cns_filter step (genotype features); no later
# step reads them, so they are deliberately not touched.
CNS_COMPANION_SUFFIXES = (".high_confident",)


def _companion_key(fields):
    """The key a companion file is joined on downstream: (chrm, fields[1]).

    NOT call_position(). update_high_confident_callset reads `int(fields[1])`
    directly, so that -- not the lclip/rclip derivation used to cluster loci -- is
    the identity under which the two files have to agree."""
    return fields[0], fields[1]


def assign_and_filter(raw_paths, window=MERGE_WINDOW, sf_conflicts=None,
                      backup_suffix=".pre_family_assign",
                      companion_suffixes=CNS_COMPANION_SUFFIXES):
    """Stage 2c.3: resolve family ONCE, mid-chain, and hand each family only the
    candidates it won.

    raw_paths is {type_name: raw-result path} -- in practice each family's
    `candidate_disc_filtered_cns.txt`, the earliest file the rule can read (see
    load_rows). Each file is moved aside to `<path><backup_suffix>` and rewritten
    with the losing families' rows removed, so every step after this one sees a
    partitioned candidate set instead of three overlapping ones.

    Two things this deliberately does NOT do:

    - It does not invent a rule. It is the same LocusTable resolution the merge
      has always applied, just applied sooner, so a difference between the two
      is a difference in the *population* being resolved and nothing else.
    - It does not drop a line the parser did not recognise as a claim (fewer
      than 4 fields). Those pass through untouched; only rows that lost a
      resolved conflict are removed.

    companion_suffixes are filtered alongside each candidate list, by key rather
    than by row index -- see CNS_COMPANION_SUFFIXES for why leaving them alone is
    not an option.

    Returns (stats, table): stats is {type_name: (n_rows_in, n_rows_kept)}.
    """
    rows_by_type = {}
    for type_, path in raw_paths.items():
        if path:
            rows_by_type[type_] = load_rows(path, type_)

    table = LocusTable(window=window)
    loci = table.assign_rows(rows_by_type)

    won = {}
    for locus in loci:
        won.setdefault(locus.type, set()).add(locus.idx)

    stats = {}
    for type_, path in raw_paths.items():
        if not path or not os.path.isfile(path):
            continue
        claimed = {r.idx for r in rows_by_type.get(type_, ())}
        keep = won.get(type_, set())
        backup = path + backup_suffix
        os.replace(path, backup)
        n_in = n_kept = 0
        kept_keys = set()
        with open(backup) as fin, open(path, "w") as fout:
            for idx, line in enumerate(fin):
                if idx in claimed:
                    n_in += 1
                    if idx not in keep:
                        continue
                    n_kept += 1
                    kept_keys.add(_companion_key(line.split()))
                fout.write(line)
        stats[type_] = (n_in, n_kept)

        # Companions carry a subset of the same rows and are joined by key, so they
        # are filtered by key -- a row index means nothing across the two files.
        for suffix in companion_suffixes:
            cpath = path + suffix
            if not os.path.isfile(cpath):
                continue
            cbackup = cpath + backup_suffix
            os.replace(cpath, cbackup)
            with open(cbackup) as fin, open(cpath, "w") as fout:
                for line in fin:
                    fields = line.split()
                    if len(fields) >= 4 and _companion_key(fields) not in kept_keys:
                        continue
                    fout.write(line)

    if sf_conflicts:
        with open(sf_conflicts, "w") as fout:
            fout.write(CONFLICT_HEADER + "\n")
            for c in table.conflicts:
                fout.write(c.row() + "\n")
    return stats, table


def load_other_families(sf_type_map):
    """{(chrm, pos): "L1:88507142,..."} for the loci in an assign_and_write
    sidecar that beat another call. A 3-column (pre-rule-B) map gives {}."""
    m = {}
    with open(sf_type_map) as fin:
        for line in fin:
            fields = line.split()
            if len(fields) >= 4 and fields[3] != ".":
                m[(fields[0], int(fields[1]))] = fields[3]
    return m


def load_type_map(sf_type_map):
    """Loads an assign_and_write sidecar into {(chrm, pos): rtype_int}.
    Reads the first 3 columns only, so old and new maps both load."""
    m = {}
    with open(sf_type_map) as fin:
        for line in fin:
            fields = line.split()
            if len(fields) < 3:
                continue
            m[(fields[0], int(fields[1]))] = int(fields[2])
    return m
