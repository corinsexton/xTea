#!/usr/bin/env python3
"""Which flank records can be a transduction source, and which are silently classified as decoys.

Offline, seconds, no BAM. Two halves that answer one question --

    python docs/handoff/runs/probe_td_source_classes.py \
        --rep-lib-dir /home/cos689/xTea/rep_lib_annotation [--run-dir <out>/<sample>]

WHAT THIS IS FOR. `parse_transduction_disc_algnmt` (`x_transduction.py:1160`) decides whether a
disc read's flank hit is a real source or a decoy consensus by a **substring test on the record
name**, with `TD_DECOY_SVA = "SVA"`, `TD_DECOY_LINE = "LINE1"`, `TD_DECOY_ALU = "ALU"`,
`TD_DECOY_HERV = "HERV"`. When the branch fires it counts the hit into the cross-family tally and
`continue`s -- the record can never become a transduction source.

That test is exactly right for the L1 library and structurally wrong for the SVA one, and both
halves of this probe measure it rather than asserting it:

  * **the library half** reads each family's flank FASTA and asks, of every record that the
    substring test would call a decoy, whether it *looks like* a flank record -- a flank id is
    `chrom~start~end~subfamily~rc[LR]`, so it has at least three `~`-delimited fields and integer
    coordinates in fields 2 and 3, while a decoy consensus record (`>LINE1`, `>ALU_AluYb8`,
    `>SVA_D  SINE  Homo sapiens`) has none. A flank record the test calls a decoy is a COLLISION.
  * **the run half** reads a finished run's own `…unique_trsdct_disc_only` and classifies each
    picked source as a reference copy or a polymorphic one. `gnrt_flank_regions_of_polymerphic_
    insertions` writes polymorphic flanks with `start == end`, so the source id
    `chrom:start-end~rc` says which it was with no ambiguity.

MEASURED (job `49634290`, hg38 rep_lib): the L1 library's 45 decoy-substring matches are all
genuine consensus records and NONE of its 662 flank records collide, while **every one of the SVA
library's 5,454 flank records collides** -- they are named `chr13~51077120~51078923~SVA_D~0L`. And
the run agrees: of L1's 159 picked sources 114 are reference copies, of SVA's 105 **zero** are.

THIS PROBE DOES NOT FIX ANYTHING and the fix is not a tidy-up: making the test structural moves
calls and BASELINE-0 with them. What the probe is for is being run on both arms of that change.
The behaviour arrived from master byte-identical (`46520ab`), so it is not an xTea2 regression.
"""
import argparse
import os
import sys

REPO_ROOT = os.path.dirname(os.path.dirname(os.path.dirname(
    os.path.dirname(os.path.abspath(__file__)))))
sys.path.insert(0, os.path.join(REPO_ROOT, "xtea"))
import global_values  # noqa: E402

# The flank libraries, per family, exactly as run_xtea's path table names them.
FLANK_RELPATH = {
    "L1": "LINE/{g}/{g}_FL_L1_flanks_3k.fa",
    "SVA": "SVA/{g}/{g}_FL_SVA_flanks_3k.fa",
}
# The transduction step's per-family source file.
SOURCES = ("candidate_disc_filtered_cns2.txt.all_non_sibling_td.txt.unique_trsdct_disc_only")

DECOY_TAGS = (global_values.TD_DECOY_LINE, global_values.TD_DECOY_SVA,
              global_values.TD_DECOY_ALU, global_values.TD_DECOY_HERV)


def is_decoy_by_substring(name):
    """`parse_transduction_disc_algnmt`'s own test, verbatim in effect."""
    return [t for t in DECOY_TAGS if t in name]


def looks_like_a_flank_record(name):
    """True if the id has a flank's STRUCTURE: chrom~start~end~subfamily~rc, coords integral.

    This is what the substring test should have been, and it is deliberately written here rather
    than imported -- nothing in the shipping code does this yet, and a probe that imported a
    function that does not exist would be a plan, not a measurement.
    """
    f = name.split(global_values.S_DELIM)
    if len(f) < 3:
        return False
    try:
        int(f[1])
        int(f[2])
    except ValueError:
        return False
    return True


def scan_library(path):
    """(n_records, n_substring_decoy, [(name, tags)] for the collisions)."""
    n = 0
    n_decoy = 0
    collisions = []
    with open(path) as fin:
        for line in fin:
            if not line.startswith(">"):
                continue
            n += 1
            # The step sees pysam's reference_name, which is the id up to the first whitespace.
            name = line[1:].split()[0] if line[1:].split() else ""
            tags = is_decoy_by_substring(name)
            if not tags:
                continue
            n_decoy += 1
            if looks_like_a_flank_record(name):
                collisions.append((name, tags))
    return n, n_decoy, collisions


def classify_sources(path):
    """(n_reference_copy, n_polymorphic, n_unparsed) from a `unique_trsdct_disc_only` file.

    A source id is `chrom:start-end~rc`. Polymorphic flanks are written with start == end
    (`gnrt_flank_regions_of_polymerphic_insertions` uses `ins_pos` for both), reference copies
    carry the copy's real interval.
    """
    n_ref = n_poly = n_bad = 0
    with open(path) as fin:
        for line in fin:
            f = line.split()
            if len(f) < 3:
                continue
            sid = f[2].split(global_values.S_DELIM)[0]
            try:
                _chrm, span = sid.rsplit(":", 1)
                start, end = span.split("-")
                if int(start) == int(end):
                    n_poly += 1
                else:
                    n_ref += 1
            except ValueError:
                n_bad += 1
    return n_ref, n_poly, n_bad


def main():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--rep-lib-dir", required=True)
    p.add_argument("--genome", default="hg38")
    p.add_argument("--run-dir", help="a finished run's SAMPLE directory (optional second half)")
    a = p.parse_args()

    print("=" * 78)
    print("TRANSDUCTION SOURCES vs THE DECOY SUBSTRING TEST")
    print("  decoy tags: {0}".format(", ".join(repr(t) for t in DECOY_TAGS)))
    print("=" * 78)

    problems = []
    print("\n--- the flank libraries ---")
    for fam, rel in sorted(FLANK_RELPATH.items()):
        path = os.path.join(a.rep_lib_dir, rel.format(g=a.genome))
        if not os.path.isfile(path):
            print("  {0}: MISSING {1}".format(fam, path))
            problems.append("{0}'s flank library is missing".format(fam))
            continue
        n, n_decoy, collisions = scan_library(path)
        print("  {0:<4} {1}".format(fam, os.path.basename(path)))
        print("       records {0:>6}   called decoy by the substring test {1:>6}".format(n, n_decoy))
        print("       of those, records with a FLANK's structure (collisions): {0:>6}".format(
            len(collisions)))
        if collisions:
            frac = 100.0 * len(collisions) / n if n else 0.0
            print("       *** {0:.1f}% of this library can never be a transduction source ***".format(
                frac))
            for (name, tags) in collisions[:3]:
                print("         e.g. {0}   matched {1}".format(name, ",".join(tags)))
        # A library with NO decoy consensus records cannot support the cross-family tally at all.
        if n_decoy == len(collisions):
            print("       and it has NO genuine decoy records, so the cross-family tally for")
            print("       this family has nothing to count -- its arm of gate_decoy_cross_family")
            print("       is dead code, not a measurement of zero competition.")

    if not a.run_dir:
        print("\n(no --run-dir given; skipping the run half)")
        return

    print("\n--- a finished run's picked sources ---")
    for fam in sorted(FLANK_RELPATH):
        path = os.path.join(a.run_dir, fam, SOURCES)
        if not os.path.isfile(path):
            print("  {0}: no source file at {1}".format(fam, path))
            continue
        n_ref, n_poly, n_bad = classify_sources(path)
        tot = n_ref + n_poly + n_bad
        print("  {0:<4} picked sources {1:>5}   reference copies {2:>5}   polymorphic {3:>5}"
              "   unparsed {4}".format(fam, tot, n_ref, n_poly, n_bad))
        if tot and n_ref == 0:
            print("       *** not one reference copy: this family's transductions come only from")
            print("       polymorphic flanks, which is the collision above, end to end ***")

    print("\n  domain assumptions: {0}".format(
        "OK" if not problems else "FAILED -- " + "; ".join(problems)))


if __name__ == "__main__":
    main()
