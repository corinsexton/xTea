"""Phase C move A.1: explicit locus table for the per-type -> single-VCF merge.

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
priority or TPRT flags. Verified against the real merger (both by direct
invocation and by inspecting the WGS seeded merged golden, where all 5
non-TPRT-override conflicts -- all exact-position ties -- resolve to Alu
instead of the intended L1/SVA priority winner). This module keeps each
chromosome's placed loci in a list (not a position-keyed dict), so a resolved
conflict always replaces-or-drops the loser in place: no collision, no
duplicate, independent of whether the two positions happen to be equal.
"""
import os
import subprocess

import global_values

MERGE_WINDOW = 50
TPRT = global_values.TWO_SIDE_TPRT_BOTH
TYPE_PRIORITY = ("L1", "SVA", "Alu")


class Locus:
    """One surviving call at a resolved locus."""
    __slots__ = ("chrm", "pos", "type", "tprt", "line")

    def __init__(self, chrm, pos, type_, tprt, line):
        self.chrm = chrm
        self.pos = pos
        self.type = type_
        self.tprt = tprt
        self.line = line


def _load_vcf_rows(sf_vcf, type_):
    rows = []
    if not os.path.isfile(sf_vcf):
        return rows
    with open(sf_vcf) as fin:
        for line in fin:
            if not line or line[0] == "#":
                continue
            fields = line.split("\t")
            if len(fields) < 3:
                continue
            line = line.rstrip("\n")
            rows.append(Locus(fields[0], int(fields[1]), type_, TPRT in line, line))
    return rows


def _load_vcf_head(sf_vcf):
    """Mirrors x_vcf_merger.VCFMerger.load_vcf_head (contig-line filtering) --
    generic VCF-header plumbing, not part of the type-resolution rule."""
    head = []
    if not os.path.isfile(sf_vcf):
        return head
    with open(sf_vcf) as fin:
        for line in fin:
            if len(line) > 0 and line[0] == "#":
                if "contig" in line:
                    if "HLA" in line or "Un" in line or "_" in line:
                        continue
                head.append(line.rstrip())
    return head


class LocusTable:
    """Collapses per-type VCF calls to one call per locus.

    Types load in TYPE_PRIORITY order (priority = load order). A call within
    +-`window` bp of an already-placed call is the same locus; the placed call
    is kept UNLESS the incoming call carries two_side_tprt_both and the placed
    one does not.
    """

    def __init__(self, window=MERGE_WINDOW):
        self.window = window
        self._by_chrm = {}  # chrm -> [Locus, ...], insertion order

    def assign_type(self, type_vcfs):
        """type_vcfs: {type_name: vcf_path}. Returns the surviving Locus records."""
        for type_ in TYPE_PRIORITY:
            sf_vcf = type_vcfs.get(type_)
            if not sf_vcf:
                continue
            for locus in _load_vcf_rows(sf_vcf, type_):
                self._place(locus)
        return [locus for bucket in self._by_chrm.values() for locus in bucket]

    def _place(self, locus):
        bucket = self._by_chrm.setdefault(locus.chrm, [])
        for i, placed in enumerate(bucket):
            if abs(placed.pos - locus.pos) <= self.window:
                if locus.tprt and not placed.tprt:
                    bucket[i] = locus  # TPRT-both overrides load-order priority
                return  # conflict resolved either way -- never duplicate
        bucket.append(locus)


def write_merged_vcf(type_vcfs, sf_merged, window=MERGE_WINDOW):
    """Resolve type_vcfs ({type: path}) into one sorted VCF at sf_merged.
    Returns the list of surviving Locus records."""
    head = []
    for type_ in TYPE_PRIORITY:
        sf_vcf = type_vcfs.get(type_)
        if sf_vcf and os.path.isfile(sf_vcf):
            head = _load_vcf_head(sf_vcf)
            break

    loci = LocusTable(window=window).assign_type(type_vcfs)

    sf_unsorted = sf_merged + ".unsorted"
    with open(sf_unsorted, "w") as fout:
        for locus in loci:
            fout.write(locus.line + "\n")
    sf_sorted = sf_merged + ".sorted"
    subprocess.run(["sort", "-k1,1V", "-k2,2n", "-o", sf_sorted, sf_unsorted], check=True)

    with open(sf_merged, "w") as fout:
        for line in head:
            fout.write(line + "\n")
        with open(sf_sorted) as fin:
            for line in fin:
                fout.write(line.rstrip("\n") + "\n")
    os.remove(sf_unsorted)
    os.remove(sf_sorted)
    return loci
