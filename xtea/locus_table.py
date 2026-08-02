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
    """One surviving call at a resolved locus."""
    __slots__ = ("chrm", "pos", "type", "tprt", "line")

    def __init__(self, chrm, pos, type_, tprt, line):
        self.chrm = chrm
        self.pos = pos
        self.type = type_
        self.tprt = tprt
        self.line = line


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


def _load_hc_pf_rows(sf_hc_pf, type_):
    """Parse a postF_hc raw-result file (candidate_disc_filtered_cns.txt.high_
    confident.post_filtering.txt): whitespace-separated, no header, one call per
    line. two_side_tprt_both is a literal substring of the line when present."""
    rows = []
    if not os.path.isfile(sf_hc_pf):
        return rows
    with open(sf_hc_pf) as fin:
        for line in fin:
            fields = line.split()
            if len(fields) < 4:
                continue
            line = line.rstrip("\n")
            rows.append(Locus(fields[0], call_position(fields), type_, TPRT in line, line))
    return rows


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
        """hc_pf_paths: {type_name: postF_hc path}. Returns the surviving Locus records."""
        for type_ in TYPE_PRIORITY:
            sf_hc_pf = hc_pf_paths.get(type_)
            if not sf_hc_pf:
                continue
            for locus in _load_hc_pf_rows(sf_hc_pf, type_):
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
