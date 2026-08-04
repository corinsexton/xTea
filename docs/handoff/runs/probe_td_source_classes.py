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

MEASURED, and the per-build sweep separated TWO defects that look like one:

  * **the collision is on every build** -- SVA flank records collide at 99.1% (hg19, 5,114 of
    5,158), 100% (hg38, 5,454) and 100% (chm13, 6,072), while L1's collide nowhere. That is the
    code defect.
  * **whether a library has any genuine decoy records at all varies by build** -- hg19's SVA
    library carries all 44, hg38's and chm13's SVA carry none, and **chm13's L1 carries none
    either**, so on chm13 even L1's cross-family tally is dead. That is library packaging, not
    code, and it is why hg38's `decoy_cross_family` zero for SVA is a dead arm rather than a
    measurement.

And the run agrees with the collision end to end: of L1's 159 picked sources 114 are reference
copies, of SVA's 105 **zero** are (job `49634290`, hg38).

THE FIX NOW EXISTS AND THIS PROBE TESTS IT RATHER THAN DESCRIBING IT. `transduction_evidence.
is_flank_record` / `decoy_tags(name, b_structural)` are the shipping classification, behind
`run_xtea --td-decoy-structural` (`XTEA_TD_DECOY_STRUCTURAL`), off by default because turning it
on moves SVA calls and BASELINE-0 with them. Both halves below now IMPORT those functions
instead of writing their own copy, which is the whole point: what the probe measures and what
the step does are the same code, so a collision count here is a statement about the shipping
step. The behaviour it replaces arrived from master byte-identical (`46520ab`), so it was never
an xTea2 regression. Run this on both arms of the flip.
"""
import argparse
import os
import sys

REPO_ROOT = os.path.dirname(os.path.dirname(os.path.dirname(
    os.path.dirname(os.path.abspath(__file__)))))
sys.path.insert(0, os.path.join(REPO_ROOT, "xtea"))
import global_values  # noqa: E402
import transduction_evidence  # noqa: E402

# The flank libraries, per family, exactly as run_xtea's path table names them.
FLANK_RELPATH = {
    "L1": "LINE/{g}/{g}_FL_L1_flanks_3k.fa",
    "SVA": "SVA/{g}/{g}_FL_SVA_flanks_3k.fa",
}
# The transduction step's per-family source file.
SOURCES = ("candidate_disc_filtered_cns2.txt.all_non_sibling_td.txt.unique_trsdct_disc_only")

DECOY_TAGS = transduction_evidence.DECOY_TAGS


def is_decoy_by_substring(name):
    """The step's test with `--td-decoy-structural` OFF, which is the shipping default."""
    return transduction_evidence.decoy_tags(name, b_structural=False)


def looks_like_a_flank_record(name):
    """The step's test with `--td-decoy-structural` ON: does the id have a flank's STRUCTURE.

    Imported, not reimplemented. An earlier version of this probe wrote the rule itself because
    the shipping code had no such function; now that it does, a probe with its own copy could
    agree with a step that disagrees with it.
    """
    return transduction_evidence.is_flank_record(name)


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


def report_library(genome, fam, path):
    """One library's two independent verdicts, which the per-build sweep showed are separable.

    A library can be broken in either of two ways and they are NOT the same defect:

      * **collisions** -- flank records the substring test calls decoys, so they can never be a
        transduction source. This is a CODE defect and it is present on every build.
      * **no genuine decoy records** -- the consensus was never concatenated into this library, so
        the cross-family tally has nothing to count and that family's arm of
        `gate_decoy_cross_family` is dead. This is a LIBRARY PACKAGING difference and it varies by
        build: hg19's SVA library carries all 44 (`>ALU`, 35 `>ALU_Alu*`, `>LINE1`, `>SVA`, 6
        `>SVA_[A-F]`), hg38's and chm13's SVA carry none, and chm13's L1 carries none either.
    """
    n, n_decoy, collisions = scan_library(path)
    n_genuine = n_decoy - len(collisions)
    print("  {0:<6} {1:<4} {2}".format(genome, fam, os.path.basename(path)))
    print("         records {0:>6}   matched the substring test {1:>6}   of those, genuine decoy "
          "records {2:>4}".format(n, n_decoy, n_genuine))
    if collisions:
        frac = 100.0 * len(collisions) / n if n else 0.0
        print("         COLLISIONS {0:>6} ({1:.1f}% of the library) -- these can never be a "
              "transduction source".format(len(collisions), frac))
        for (name, tags) in collisions[:2]:
            print("           e.g. {0}   matched {1}".format(name, ",".join(tags)))
    if n_genuine == 0:
        print("         NO genuine decoy records, so the cross-family tally has nothing to count")
        print("         -- this family's arm of gate_decoy_cross_family is dead on this build,")
        print("         which is a library-packaging fact and not the collision above.")


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
    p.add_argument("--genome", default="hg38",
                   help="comma-separated; the collision is genome-independent, so sweeping every "
                        "build the rep_lib has is the honest check (hg19,hg38,chm13)")
    p.add_argument("--run-dir", help="a finished run's SAMPLE directory (optional second half)")
    a = p.parse_args()

    print("=" * 78)
    print("TRANSDUCTION SOURCES vs THE DECOY SUBSTRING TEST")
    print("  decoy tags: {0}".format(", ".join(repr(t) for t in DECOY_TAGS)))
    print("=" * 78)

    problems = []
    print("\n--- the flank libraries ---")
    genomes = [g for g in a.genome.split(",") if g]
    for genome in genomes:
        for fam, rel in sorted(FLANK_RELPATH.items()):
            path = os.path.join(a.rep_lib_dir, rel.format(g=genome))
            if not os.path.isfile(path):
                print("  {0}/{1}: MISSING {2}".format(genome, fam, path))
                if len(genomes) == 1:
                    problems.append("{0}'s flank library is missing".format(fam))
                continue
            report_library(genome, fam, path)

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
