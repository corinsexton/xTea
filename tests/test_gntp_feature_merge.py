"""Multi-BAM merge in `XGenotyper.call_genotype`, which produces the AF denominators.

WHY THIS GATE EXISTS
--------------------
The merge loop summed each site's features across BAMs, but only field 0 (`n_af_clip`)
was indexed by the loop variables `(ins_chrm, ins_pos)`. Fields 1-9 read
`m_info[chrm][pos]` -- `chrm`/`pos` are leftovers from the candidate-list loop ~40 lines
above, so they held the LAST line of `sf_before_filtering` and never varied inside the
merge. Every site seen in 2+ BAMs therefore got

    field_k = this_site_from_bam0 + ONE_FIXED_UNRELATED_SITE_from_bam1      (k = 1..9)

`n_full_map`, `n_disc_pairs` and `n_concd_pairs` are the denominators of the allele
fraction, so the AF that every mosaic gate and the VCF's own AF field read was wrong by
whatever that one site happened to carry.

WHAT IT COST, measured rather than reasoned about. On broad (2 BAMs, HapMap mixture) the
stale site for SVA was chrY:11531991, whose BAM-1 row is n_full_map=1102, n_disc=3471,
n_concod=3894. SVA's truth-matching records came out at medians 1147/3484/4070 -- about
45/13/176 of real signal sitting under that constant. Subtracting each family's own
stale constant collapses all three families onto the same real numbers (L1 ~37/12/162,
Alu ~39/10/170, SVA ~45/13/176), which is what one sample at one coverage should look
like. The visible symptom was SVA scoring 0/179 against the truth set: its clip AF read
~0.005 instead of ~0.14, below the mosaic filter's 0.02 floor, so every true positive
was rejected on an arithmetic artefact.

The bug is invisible in a single-BAM run (the merge branch is never entered), which is
why it survived: the chr1 harness and several goldens are one BAM. So the assertions
below deliberately exercise BOTH arities, and the multi-BAM case uses DISTINCT values
per site so a stale index cannot coincidentally produce the right answer.
"""
import os
import sys

import pytest

sys.path.insert(0, os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "xtea"))

from x_genotype_feature import XGenotyper  # noqa: E402


# field order written by collect_features_one_site / read by load_in_features_from_file:
# chrm pos n_af_clip n_full_map n_l_raw_clip n_r_raw_clip n_disc_pairs n_concd_pairs
#      n_disc_large_indel s_clip_lens n_polyA n_disc_chrms
def _feature_row(chrm, pos, base):
    """A features line whose 9 numeric fields are `base`+0 .. `base`+8, so a wrong index
    is arithmetically detectable rather than accidentally equal."""
    nums = [base + i for i in range(7)]          # n_af_clip .. n_disc_large_indel
    return "{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\t{7}\t{8}\t{9}\t{10}\t{11}\n".format(
        chrm, pos, nums[0], nums[1], nums[2], nums[3], nums[4], nums[5], nums[6],
        "20:30", base + 7, base + 8)


class _SerialPool:
    """Stand-in for multiprocessing.Pool so the stubbed per-site collector runs in THIS
    process (a forked worker would not see the monkeypatch on some platforms, and the
    real one needs BAMs we deliberately do not have)."""

    def __init__(self, n):
        pass

    def map(self, fn, iterable, chunksize=None):
        return [fn(x) for x in iterable]

    def close(self):
        pass

    def join(self):
        pass


def _run_merge(tmp_path, per_bam_rows, monkeypatch, clip_lens=None):
    """Drive the real call_genotype, stubbing ONLY the per-site BAM reader.

    per_bam_rows[i] is [(chrm, pos, base), ...] for BAM i: the features that BAM's pass
    should produce. Everything after that -- the per-site file naming, the per-BAM
    concatenation, and the merge under test -- is the shipped code path.
    """
    import x_genotype_feature as xgf

    sf_out = str(tmp_path / "feat")
    wf = str(tmp_path) + "/"

    sites = sorted({(c, p) for rows in per_bam_rows for c, p, _ in rows})
    cand = str(tmp_path / "cands.txt")
    with open(cand, "w") as fout:
        for c, p in sites:
            fout.write("{0}\t{1}\n".format(c, p))

    bam_list = str(tmp_path / "bams.txt")
    with open(bam_list, "w") as fout:
        for i in range(len(per_bam_rows)):
            fout.write("/nonexistent/bam{0}.bam\n".format(i))

    # (bam path, chrm, pos) -> base value for that BAM's pass over that site
    lookup = {}
    for i, rows in enumerate(per_bam_rows):
        for c, p, base in rows:
            lookup[("/nonexistent/bam{0}.bam".format(i), c, p)] = base

    def _fake_collect(self, record):
        (chrm, ins_pos, _extnd), sf_bam, working_folder = record
        base = lookup.get((sf_bam, chrm, ins_pos))
        if base is None:                       # this BAM saw nothing at this site
            return
        path = working_folder + "{0}_{1}".format(chrm, ins_pos) + xgf.GNTP_FEATURE_SUFFIX
        with open(path, "w") as fout:
            if clip_lens is None:
                fout.write(_feature_row(chrm, ins_pos, base))
            else:
                fout.write(clip_lens(chrm, ins_pos, base))

    monkeypatch.setattr(xgf, "Pool", _SerialPool)
    monkeypatch.setattr(XGenotyper, "collect_features_one_site", _fake_collect, raising=True)

    XGenotyper("ref.fa", wf, 1).call_genotype(bam_list, cand, 450, sf_out)

    merged = {}
    with open(sf_out) as fin:
        for line in fin:
            f = line.rstrip("\n").split("\t")
            merged[(f[0], int(f[1]))] = f[2:]
    return merged


def test_single_bam_features_pass_through_unchanged(tmp_path, monkeypatch):
    """One BAM never enters the merge branch -- the arity where the bug was invisible."""
    merged = _run_merge(tmp_path, [[("chr1", 1000, 10), ("chr1", 2000, 100)]], monkeypatch)
    assert merged[("chr1", 1000)][:7] == [str(10 + i) for i in range(7)]
    assert merged[("chr1", 2000)][:7] == [str(100 + i) for i in range(7)]


def test_two_bams_sum_each_site_with_its_own_values(tmp_path, monkeypatch):
    """THE REGRESSION. Two sites per BAM with very different magnitudes; the second site
    is the 'last line' whose values used to leak into every other site's fields 1-9."""
    merged = _run_merge(tmp_path, [
        [("chr1", 1000, 10), ("chrY", 9999, 1000)],
        [("chr1", 1000, 20), ("chrY", 9999, 2000)],
    ], monkeypatch)
    # chr1:1000 must be its own values summed: (10+i) + (20+i)
    got = [int(v) for v in merged[("chr1", 1000)][:7]]
    assert got == [(10 + i) + (20 + i) for i in range(7)], (
        "fields 1-9 of chr1:1000 must sum chr1:1000 across BAMs. Got {0}; with the stale "
        "index they were chr1:1000(bam0) + chrY:9999(bam1) = {1}".format(
            got, [(10 + 0) + (20 + 0)] + [(10 + i) + (2000 + i) for i in range(1, 7)]))
    # and the high-magnitude site is unaffected by the low-magnitude one
    assert [int(v) for v in merged[("chrY", 9999)][:7]] == \
           [(1000 + i) + (2000 + i) for i in range(7)]


def test_three_bams_accumulate_per_site(tmp_path, monkeypatch):
    """The merge folds left across BAMs; with >2 the stale index compounds."""
    merged = _run_merge(tmp_path, [
        [("chr1", 1000, 1), ("chrY", 9999, 500)],
        [("chr1", 1000, 2), ("chrY", 9999, 600)],
        [("chr1", 1000, 4), ("chrY", 9999, 700)],
    ], monkeypatch)
    assert [int(v) for v in merged[("chr1", 1000)][:7]] == \
           [(1 + i) + (2 + i) + (4 + i) for i in range(7)]


def test_clip_len_string_concatenates_the_right_site(tmp_path, monkeypatch):
    """s_clip_lens (field 7) is a string join, so a stale index corrupts it too --
    silently, since any concatenation still 'looks like' a clip-length vector."""
    tag = {10: "AAA", 20: "BBB", 1000: "ZZZ", 2000: "YYY"}

    def _row_with_tag(chrm, pos, base):
        return "{0}\t{1}\t1\t2\t3\t4\t5\t6\t7\t{2}\t8\t9\n".format(chrm, pos, tag[base])

    merged = _run_merge(tmp_path, [
        [("chr1", 1000, 10), ("chrY", 9999, 1000)],
        [("chr1", 1000, 20), ("chrY", 9999, 2000)],
    ], monkeypatch, clip_lens=_row_with_tag)

    lens = merged[("chr1", 1000)][7]
    assert "ZZZ" not in lens and "YYY" not in lens, (
        "chr1:1000 joined chrY:9999's clip lengths: {0!r}".format(lens))
    assert set(lens.split(":")) == {"AAA", "BBB"}


if __name__ == "__main__":
    sys.exit(pytest.main([__file__, "-v"]))
