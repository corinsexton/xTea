"""The transduction step's decoy consensus block: what it is, and how to rebuild it.

WHY THIS MODULE EXISTS -- a comment that is not behaviour, and a family asymmetry that
follows from it.

`x_transduction.construct_novel_flanks` builds the bwa reference that BOTH the clip and the
disc realignment pass use, and both of its call sites are preceded by this comment:

    ####Note, we add L1, Alu and SVA consensus to the flank regions to filter out those
    ####mapped to consensus

**It does not.** The function concatenates the family's packaged flank library with the
polymorphic flanks it just built, and nothing else. Whether any decoy consensus is in that
reference is therefore purely a property of how the build's library was PACKAGED -- and the
packaging is not uniform (measured 2026-08-05, `docs/handoff/runs/probe_decoy_library.py`):

    library                            records   decoy records
    LINE/hg38  hg38_FL_L1_flanks_3k        707      45
    LINE/hg19  hg19_FL_L1_flanks_3k        705      45
    LINE/chm13 chm13_FL_L1_flanks_3k       818       0
    SVA/hg38   hg38_FL_SVA_flanks_3k      5454       0
    SVA/hg19   hg19_FL_SVA_flanks_3k      5158      44   (all but HERV-K)
    SVA/chm13  chm13_FL_SVA_flanks_3k     6072       0

So on the path this project measures -- SVA on hg38 -- the decoy records the comment
promises are simply absent, and that absence is visible from two independent directions in
the step's own numbers:

  * the DISC side: `gate_decoy_cross_family`'s SVA arm has nothing to count, while L1's
    fires 144 times at WGS (plan sec 10);
  * the CLIP side: L1 counts **3,992 of 99,245** qualified clip records as support from a
    `>LINE1`/`>SVA`/`>ALU` consensus record and SVA counts **0 of 128,262**
    (`HANDOFF_2026-08-05_gate_a_measured.md` sec 2.4).

Both are the same packaging hole reached by different routes, and it is what leaves plan
sec 5's SVA/Alu mapping hazard -- SVA is SINE-VNTR-*Alu*, so SVA and Alu sequence really do
compete -- **untested** rather than measured-as-absent.

WHAT MAKES THE FIX CHEAP: the decoy block needs no new data and no library repackaging.

The 45 records are all reconstructible from `<rep_lib>/consensus/`, which the transduction
step already has in hand -- it is the directory of the `-r` consensus FASTA that x_TEA_main
is passed. And "reconstructible" is meant literally: `decoy_block_text` below emits the
byte-identical tail of both packaged L1 libraries, line wrapping and header comment fields
included (`tests/test_td_decoy_library.py`, and the same check on hg19's SVA library as an
unordered set, since that one packages the 44 in a different order and without HERV-K).

WHAT THIS MODULE DOES NOT DECIDE. Concatenating is CALL-MOVING on any library that lacks
the block, so the caller gates it behind `--td-decoy-concat` -- **ON by default since
2026-08-07**, once its two-scale arm was run and passed; `--no-td-decoy-concat` reproduces a
pre-flip baseline. This module only knows what the block IS; `construct_novel_flanks`
decides whether to append it, and appends only what the library does not already carry --
which is what makes the change a structural no-op on both L1 libraries rather than a source
of duplicate bwa reference names.
"""
import os

# The consensus files the block is assembled from, in the order the packaged L1 libraries
# carry it: the three transduction decoys, HERV-K, the six SVA subfamilies, then the 35 Alu
# subfamilies under an `ALU_` prefix.
#
# Records are taken BY NAME rather than by position. A positional cut would be shorter --
# the Alu 35 are exactly `ALU.fa`'s first 35 and the SVA 6 its first 6 -- but it would also
# be silently wrong if a consensus file were ever re-ordered or extended, and this block
# ends up in a bwa reference whose contents nothing downstream re-checks. Both files DO
# continue past the cut (`ALU.fa` has 47 records, `SVA.fa` 88), so the cut is real and
# naming it is what pins it.
#
# `SVA` and `ALU` come from `transduction_decoy.fa`, not from `SVA.fa`/`ALU.fa`, both of
# which also carry a record of that name. That is the packaged libraries' choice, and the
# byte-identity test is what proves it is theirs and not a guess.
DECOY_SOURCES = (
    ("transduction_decoy.fa", ("LINE1", "SVA", "ALU"), ""),
    ("HERV.fa", ("HERV-K",), ""),
    ("SVA.fa", ("SVA_A", "SVA_B", "SVA_C", "SVA_D", "SVA_E", "SVA_F"), ""),
    ("ALU.fa",
     ("AluJo", "AluJb", "AluSc", "AluSg", "AluSp", "AluSq", "AluSx", "AluSz", "AluY",
      "AluYa5", "AluYa8", "AluYb8", "AluYb9", "AluYc1", "AluYc2", "AluYd2", "AluYd3",
      "AluYd3a1", "AluYd8", "AluYe5", "AluYa1", "AluYa4", "AluYb3a1", "AluYb3a2", "AluYf1",
      "AluYg6", "AluYh9", "AluYi6", "AluYbc3a", "AluYe2", "AluYf2", "AluYk13", "AluYb11",
      "AluYb10", "AluYb8a1"),
     "ALU_"),
)

# 45. Named so a caller can assert the size without rebuilding the block, and so the number
# quoted in four handoffs has exactly one definition in the code.
N_DECOY_RECORDS = sum(len(names) for _, names, _ in DECOY_SOURCES)


def record_name(s_header):
    """The name bwa reports as `reference_name`: the header's first whitespace token.

    The distinction matters here because the packaged block keeps the source files' header
    COMMENT fields verbatim -- `>SVA_A  SINE    Homo sapiens`, `>AluYk13    SINE1/7SL...` --
    so the header and the name are not the same string, and the byte-identity this module
    claims is about the header while every name comparison is about the token.
    """
    return s_header.split()[0] if s_header.split() else ""


def read_fasta_raw(sf_fasta):
    """[(header without '>', [sequence lines])], preserving the file's own line wrapping.

    Deliberately not `pysam.FastaFile`: that needs a `.fai` beside a file in a rep_lib
    directory that is routinely read-only, and it returns one unwrapped string, which throws
    away exactly the property this module is checked on.
    """
    l_records = []
    s_header = None
    l_lines = []
    with open(sf_fasta) as fin:
        for line in fin:
            line = line.rstrip("\n")
            if line.startswith(">"):
                if s_header is not None:
                    l_records.append((s_header, l_lines))
                s_header = line[1:]
                l_lines = []
            else:
                l_lines.append(line)
    if s_header is not None:
        l_records.append((s_header, l_lines))
    return l_records


def build_decoy_block(s_consensus_dir):
    """The 45 decoy records as [(header, [sequence lines])], in the packaged order.

    Raises `ValueError` naming the file and the record if anything in `DECOY_SOURCES` is
    missing. A partial block is worse than no block: it would index cleanly, change some
    calls, and leave nothing in the run that says which records were actually there.
    """
    l_block = []
    for s_basename, l_names, s_prefix in DECOY_SOURCES:
        sf_src = os.path.join(s_consensus_dir, s_basename)
        if not os.path.isfile(sf_src):
            raise ValueError(
                "transduction decoy block: {0} is not in the consensus directory {1}".format(
                    s_basename, s_consensus_dir))
        m_by_name = {}
        for s_header, l_lines in read_fasta_raw(sf_src):
            m_by_name.setdefault(record_name(s_header), (s_header, l_lines))
        for s_name in l_names:
            if s_name not in m_by_name:
                raise ValueError(
                    "transduction decoy block: record {0} is not in {1}".format(
                        s_name, sf_src))
            s_header, l_lines = m_by_name[s_name]
            l_block.append((s_prefix + s_header, l_lines))
    return l_block


def decoy_block_text(l_block):
    """The block as FASTA text. Byte-identical to both packaged L1 libraries' tails."""
    l_out = []
    for s_header, l_lines in l_block:
        l_out.append(">" + s_header + "\n")
        for s_line in l_lines:
            l_out.append(s_line + "\n")
    return "".join(l_out)


def fasta_record_names(sf_fasta):
    """The set of bwa reference names in a FASTA, read from headers only."""
    s_names = set()
    with open(sf_fasta) as fin:
        for line in fin:
            if line.startswith(">"):
                s_names.add(record_name(line[1:].rstrip("\n")))
    return s_names


def write_missing_decoys(sf_flank, s_consensus_dir, sf_out):
    """Write to `sf_out` the decoy records `sf_flank` does not already carry, by NAME.

    Returns the list of names written -- empty on both packaged L1 libraries, which is the
    point: turning the concatenation on must be a structural no-op wherever the library was
    packaged with the block, and "no records appended" is a stronger statement of that than
    "the output happened to compare equal".

    By name and not by sequence, because a duplicate name in a bwa reference is the failure
    this is avoiding: bwa indexes both copies, `reference_name` becomes ambiguous for a
    reader, and the step's whole decoy classification is a test on that string.
    """
    s_have = fasta_record_names(sf_flank)
    l_block = build_decoy_block(s_consensus_dir)
    l_missing = [(s_header, l_lines) for s_header, l_lines in l_block
                 if record_name(s_header) not in s_have]
    with open(sf_out, "w") as fout:
        fout.write(decoy_block_text(l_missing))
    return [record_name(s_header) for s_header, _ in l_missing]
