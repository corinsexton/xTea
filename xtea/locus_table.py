"""Phase C move A.1/A.2: explicit locus table for the per-type -> single-VCF merge.

Reimplements x_vcf_merger.SampleVCFMerger.merge_single_sample_ins's INTENDED
type-resolution rule (pinned in docs/handoff/MOVE_A_PLAN.md): types load in
priority order L1 -> SVA -> Alu (HERV unused -- xTea2 runs L1/Alu/SVA only); a
call within +-MERGE_WINDOW bp of an already-placed call is the same locus; the
placed call is kept UNLESS the incoming call carries two_side_tprt_both and the
placed one does not (no numeric QUAL comparison -- x_vcf_merger's own "check
quality first" comment is dead code).

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
import subprocess

import global_values

MERGE_WINDOW = 50
TPRT = global_values.TWO_SIDE_TPRT_BOTH
TYPE_PRIORITY = ("L1", "SVA", "Alu")
# rtype bitmask codes x_TEA_main/x_rep_type use (mirrors run_xtea.py's REP_TYPES).
TYPE_TO_RTYPE = {"L1": 1, "Alu": 2, "SVA": 4}


class Locus:
    """One family's claim on one locus: the row, plus the two fields the
    resolution rule actually reads (position and two_side_tprt_both).

    `idx` is the row's 0-based line number within its own family's file. It
    exists so a caller can map a resolution decision back to the exact input
    row -- `line` cannot serve as that key, since nothing forbids a file from
    carrying two identical rows. filter_rows_to_assigned() relies on it.
    """
    __slots__ = ("chrm", "pos", "type", "tprt", "line", "idx")

    def __init__(self, chrm, pos, type_, tprt, line, idx=-1):
        self.chrm = chrm
        self.pos = pos
        self.type = type_
        self.tprt = tprt
        self.line = line
        self.idx = idx


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
    """One locus two families both called, and how the tie was broken.

    Recorded rather than discarded because this population IS Stage 2's gate.
    Stage 2c moves family assignment from here -- after three independent
    per-family chains -- to a single up-front decision, and its gate is "calls may
    move, but only as family-assignment flips". You cannot check that against a
    population nobody counted: without this list, a Stage 2c run that loses a call
    for an unrelated reason is indistinguishable from one that reassigned it.

    `rule` is `priority` (TYPE_PRIORITY load order decided it) or `tprt_override`
    (the incoming call carried two_side_tprt_both and the placed one did not, so it
    displaced a higher-priority family).
    """

    __slots__ = ("chrm", "winner_type", "winner_pos", "winner_tprt",
                 "loser_type", "loser_pos", "loser_tprt", "rule")

    def __init__(self, chrm, winner, loser, rule):
        self.chrm = chrm
        self.winner_type = winner.type
        self.winner_pos = winner.pos
        self.winner_tprt = winner.tprt
        self.loser_type = loser.type
        self.loser_pos = loser.pos
        self.loser_tprt = loser.tprt
        self.rule = rule

    def row(self):
        return "{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\t{7}".format(
            self.chrm, self.winner_type, self.winner_pos, int(self.winner_tprt),
            self.loser_type, self.loser_pos, int(self.loser_tprt), self.rule)


CONFLICT_HEADER = ("#chrm\twinner_type\twinner_pos\twinner_tprt\t"
                   "loser_type\tloser_pos\tloser_tprt\trule")


class LocusTable:
    """Collapses per-type raw postF_hc calls to one call per locus.

    Types load in TYPE_PRIORITY order (priority = load order). A call within
    +-`window` bp of an already-placed call is the same locus; the placed call
    is kept UNLESS the incoming call carries two_side_tprt_both and the placed
    one does not.

    Every such collision is also appended to `self.conflicts` -- see Conflict.
    """

    def __init__(self, window=MERGE_WINDOW):
        self.window = window
        self._by_chrm = {}  # chrm -> [Locus, ...], insertion order
        self.conflicts = []

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
        bucket = self._by_chrm.setdefault(locus.chrm, [])
        for i, placed in enumerate(bucket):
            if abs(placed.pos - locus.pos) <= self.window:
                if locus.tprt and not placed.tprt:
                    # TPRT-both overrides load-order priority
                    self.conflicts.append(Conflict(locus.chrm, locus, placed, "tprt_override"))
                    bucket[i] = locus
                else:
                    self.conflicts.append(Conflict(locus.chrm, placed, locus, "priority"))
                return  # conflict resolved either way -- never duplicate
        bucket.append(locus)

    def conflict_summary(self):
        """{(winner_type, loser_type, rule): count} -- the shape Stage 2c predicts against."""
        counts = {}
        for c in self.conflicts:
            key = (c.winner_type, c.loser_type, c.rule)
            counts[key] = counts.get(key, 0) + 1
        return counts


def assign_and_write(hc_pf_paths, sf_merged_hcpf, sf_type_map, window=MERGE_WINDOW,
                     sf_conflicts=None):
    """Resolve hc_pf_paths ({type: postF_hc path}) into one deduplicated raw-result
    file at sf_merged_hcpf (sorted, native hc_pf column format -- unchanged, so
    the shared gene/genotype tail needs no changes) plus a sidecar type map at
    sf_type_map (`chrm\\tpos\\trtype_int` per surviving locus, pos = call_position)
    for the later per-record gVCF fan-out.

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
            fout.write(f"{locus.chrm}\t{locus.pos}\t{TYPE_TO_RTYPE[locus.type]}\n")

    if sf_conflicts:
        with open(sf_conflicts, "w") as fout:
            fout.write(CONFLICT_HEADER + "\n")
            for c in table.conflicts:
                fout.write(c.row() + "\n")
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


def load_type_map(sf_type_map):
    """Loads an assign_and_write sidecar into {(chrm, pos): rtype_int}."""
    m = {}
    with open(sf_type_map) as fin:
        for line in fin:
            fields = line.split()
            if len(fields) < 3:
                continue
            m[(fields[0], int(fields[1]))] = int(fields[2])
    return m
