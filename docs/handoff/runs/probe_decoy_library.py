#!/usr/bin/env python3
"""Which flank libraries carry the decoy consensus block, and whether it can be rebuilt.

    python docs/handoff/runs/probe_decoy_library.py [--rep-lib-dir ~/xTea/rep_lib_annotation]

OFFLINE, SECONDS, NO BAM, NO RUN, NO SLURM JOB. It reads six packaged FASTAs and four consensus
FASTAs and nothing else, which is why the packaging half of the SVA decoy defect never needed the
"call-moving, needs its own arm" treatment to be *measured* -- only to be turned on.

WHAT IT ESTABLISHES, and both halves are load-bearing for plan sec 5's packaging item.

1. **THE CENSUS.** `construct_novel_flanks` concatenates the family's packaged flank library with
   the polymorphic flanks it just built and nothing else, while the comment above both of its call
   sites claims the L1/Alu/SVA consensus is added. So whether the transduction step's bwa reference
   has a decoy channel at all is a property of the PACKAGING, and the packaging is not uniform:
   both L1 libraries ship 45 decoy records, hg19's SVA ships 44 (no HERV-K), and **hg38's SVA,
   chm13's SVA and chm13's L1 ship none**. hg38 SVA is the path every measurement in this project
   is taken on.

2. **THE RECIPE, AS BYTES.** The 45 records are rebuilt from `<rep_lib>/consensus/` --
   `transduction_decoy.fa`'s 3, `HERV.fa`'s 1, `SVA.fa`'s 6, `ALU.fa`'s 35 under an `ALU_` prefix
   -- and the rebuilt text is **byte-identical to the tail of both packaged L1 libraries**,
   including each source file's own line wrapping (60 / 120 / 70 / 71 columns) and the header
   comment fields. So making the SVA path symmetric needs no new data, no library repackaging and
   no build-specific file.

The recipe is **imported from `xtea/td_decoy_library.py`, not reimplemented here** -- plan sec 10's
`probe_leftmost_cns` entry is the reason: a probe carrying its own copy of a shipping rule read
2.6x wrong and was on its way into a handoff. What this file adds is the census over the six
libraries, which the module has no reason to know.

WHAT IT CANNOT SAY. Nothing here is evidence that concatenating the block IMPROVES anything. It
says the channel is absent on one family and reconstructible; whether restoring it moves SVA calls
in the right direction is what the two-scale `--td-decoy-concat` arm is for. Note also the one
asymmetry that survives the fix: the clip path does not skip decoy hits at all (the skip is
commented out at `x_transduction.py:872-875`), so on the clip side the block adds a source of
DEGENERATE-source support rather than a filter -- 3,992 of L1's 99,245 qualified clip records today
(HANDOFF_2026-08-05_gate_a_measured.md sec 2.4), and no L1 survivor rests on it alone.
"""
import argparse
import os
import sys

REPO_ROOT = os.path.dirname(os.path.dirname(os.path.dirname(
    os.path.dirname(os.path.abspath(__file__)))))
sys.path.insert(0, os.path.join(REPO_ROOT, "xtea"))
import td_decoy_library as tdl      # noqa: E402
import transduction_evidence as te  # noqa: E402

LIBRARIES = (
    ("LINE",  "hg38",  "LINE/hg38/hg38_FL_L1_flanks_3k.fa"),
    ("LINE",  "hg19",  "LINE/hg19/hg19_FL_L1_flanks_3k.fa"),
    ("LINE",  "chm13", "LINE/chm13/chm13_FL_L1_flanks_3k.fa"),
    ("SVA",   "hg38",  "SVA/hg38/hg38_FL_SVA_flanks_3k.fa"),
    ("SVA",   "hg19",  "SVA/hg19/hg19_FL_SVA_flanks_3k.fa"),
    ("SVA",   "chm13", "SVA/chm13/chm13_FL_SVA_flanks_3k.fa"),
)


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--rep-lib-dir", default=os.path.expanduser("~/xTea/rep_lib_annotation"))
    a = ap.parse_args()

    s_cns_dir = os.path.join(a.rep_lib_dir, "consensus")
    block = tdl.build_decoy_block(s_cns_dir)
    s_block = tdl.decoy_block_text(block)
    m_recipe = {tdl.record_name(h): "".join(l) for h, l in block}
    print("rep_lib: {0}".format(a.rep_lib_dir))
    print("recipe : {0} records, {1} bytes, from {2}".format(
        len(block), len(s_block), ", ".join(f for f, _, _ in tdl.DECOY_SOURCES)))
    print("")

    print("{0:<6} {1:<6} {2:>8} {3:>7}  {4}".format(
        "family", "build", "records", "decoy", "rebuild check"))
    n_fail = 0
    for s_fam, s_build, s_rel in LIBRARIES:
        sf_lib = os.path.join(a.rep_lib_dir, s_rel)
        if not os.path.isfile(sf_lib):
            print("{0:<6} {1:<6} {2:>8} {3:>7}  {4}".format(s_fam, s_build, "-", "-", "ABSENT"))
            continue
        l_recs = tdl.read_fasta_raw(sf_lib)
        l_decoy = [(h, l) for h, l in l_recs if not te.is_flank_record(tdl.record_name(h))]
        # The strongest available check, where it applies: the packaged block as raw TEXT.
        # Where the build packages a reordered subset (hg19 SVA), fall back to an unordered
        # (name, sequence) comparison and SAY SO, rather than quietly weakening the claim.
        s_verdict = "no decoy records packaged"
        if l_decoy:
            s_lib = open(sf_lib).read()
            i = s_lib.find(">" + l_decoy[0][0] + "\n")
            if len(l_decoy) == len(block) and i >= 0 and s_lib[i:] == s_block:
                s_verdict = "BYTE-IDENTICAL to the recipe (tail of file)"
            else:
                m_have = {tdl.record_name(h): "".join(l) for h, l in l_decoy}
                l_missing = sorted(set(m_recipe) - set(m_have))
                b_same = all(m_have[k] == m_recipe[k] for k in m_have if k in m_recipe)
                b_extra = sorted(set(m_have) - set(m_recipe))
                if b_same and not b_extra:
                    s_verdict = "same sequences, reordered; missing {0}".format(
                        ",".join(l_missing) or "nothing")
                else:
                    s_verdict = "DISAGREES with the recipe (extra {0})".format(b_extra)
                    n_fail += 1
        print("{0:<6} {1:<6} {2:>8} {3:>7}  {4}".format(
            s_fam, s_build, len(l_recs), len(l_decoy), s_verdict))

    print("")
    print("what --td-decoy-concat would APPEND, per library (by name, so a packaged block is a "
          "no-op):")
    for s_fam, s_build, s_rel in LIBRARIES:
        sf_lib = os.path.join(a.rep_lib_dir, s_rel)
        if not os.path.isfile(sf_lib):
            continue
        s_have = tdl.fasta_record_names(sf_lib)
        l_add = [tdl.record_name(h) for h, _ in block if tdl.record_name(h) not in s_have]
        print("  {0:<6} {1:<6} {2:>3} records{3}".format(
            s_fam, s_build, len(l_add),
            "" if len(l_add) > 4 else ("" if not l_add else ": " + ",".join(l_add))))
    return 1 if n_fail else 0


if __name__ == "__main__":
    sys.exit(main())
