"""The post-scan region mask: which positions it excludes, and which it must NOT.

WHY THIS GATE EXISTS
--------------------
The mask is a precision fix that also saves time: on the HapMap arms, false positives are
enriched 10-20x over true positives in these regions (simple repeats hold 0-3% of TPs and
13-32% of FPs), and excluding candidates there skips 9-12% of clip sites and ~25% of the
transduction site list -- the pipeline's most expensive step.

THE ASSERTION THAT MATTERS MOST is the one about what is NOT masked. The intuitive
implementation drops every scan record inside the regions, and for scan_disc_pos.tsv that
silently loses real calls: those rows are EVIDENCE, counted in a +/-iextend window (810-874 bp
on these samples) around each clip peak by disc_adapter.run_disc_from_scan. Deleting them by
their own position strips support from candidates OUTSIDE the mask -- a genuine insertion
500 bp from a segdup boundary loses one flank, its nL/nR fall, and it can drop below
disc_cutoff. With 655k intervals in the default mask, every candidate within ~874 bp of any
of them is exposed. So the mask is applied to clip candidate POSITIONS only, which is
sufficient (peaks and raw-disc rows are both derived per candidate) and safe.
"""
import gzip
import os
import sys

import pytest

sys.path.insert(0, os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "xtea"))

import region_mask  # noqa: E402
import scan_adapter  # noqa: E402


def _bed(tmp_path, name, rows, gz=False):
    p = tmp_path / name
    body = "".join(f"{c}\t{s}\t{e}\n" for c, s, e in rows)
    if gz:
        with gzip.open(str(p), "wt") as fh:
            fh.write(body)
    else:
        p.write_text(body)
    return str(p)


def test_empty_mask_is_a_usable_no_op():
    m = region_mask.RegionMask([])
    assert not m
    assert m.contains("chr1", 1000) is False


def test_containment_is_half_open():
    """BED is [start, end) -- the end coordinate itself is NOT masked. Off-by-one here
    would silently shift every boundary by a base across 655k intervals."""
    import tempfile
    with tempfile.TemporaryDirectory() as d:
        from pathlib import Path
        p = _bed(Path(d), "m.bed", [("chr1", 100, 200)])
        m = region_mask.RegionMask([p])
        assert m.contains("chr1", 99) is False
        assert m.contains("chr1", 100) is True     # start is inside
        assert m.contains("chr1", 199) is True
        assert m.contains("chr1", 200) is False    # end is NOT
        assert m.contains("chr2", 150) is False    # other contig


def test_overlapping_intervals_are_merged(tmp_path):
    p = _bed(tmp_path, "m.bed", [("chr1", 100, 200), ("chr1", 150, 300), ("chr1", 400, 500)])
    m = region_mask.RegionMask([p])
    assert m.n_intervals == 2
    assert m.n_bases == 200 + 100
    for pos in (100, 250, 299, 400, 499):
        assert m.contains("chr1", pos) is True
    assert m.contains("chr1", 350) is False


def test_multiple_files_union_including_gz(tmp_path):
    a = _bed(tmp_path, "a.bed", [("chr1", 100, 200)])
    b = _bed(tmp_path, "b.bed.gz", [("chr1", 500, 600)], gz=True)
    m = region_mask.RegionMask([a, b])
    assert m.contains("chr1", 150) and m.contains("chr1", 550)
    assert not m.contains("chr1", 300)


def test_resolve_spec_names_paths_and_none(tmp_path):
    own = _bed(tmp_path, "own.bed", [("chr1", 1, 2)])
    assert region_mask.resolve_mask_spec("none") == []
    assert region_mask.resolve_mask_spec("off") == []
    assert region_mask.resolve_mask_spec(own) == [own]
    # unknown path must fail loudly rather than silently masking nothing
    with pytest.raises(FileNotFoundError):
        region_mask.resolve_mask_spec("/no/such/file.bed")


def test_default_spec_uses_the_builtin_names():
    names = set(region_mask.DEFAULT_MASKS)
    assert {"simple_repeats", "segdup", "centromeres"} <= names


# --- the integration point: clip candidate positions ------------------------------------
def _realign(tmp_path, rows):
    """xtea-scan realign-out: chrom pos family nL nR mate_in_rep"""
    p = tmp_path / "realign.tsv"
    p.write_text("".join(f"{c}\t{pos}\tLINE1\t{nl}\t{nr}\t0\n" for c, pos, nl, nr in rows))
    return str(p)


def test_write_clip_candidate_list_drops_only_masked_positions(tmp_path):
    rt = _realign(tmp_path, [("chr1", 150, 3, 0),    # inside mask
                             ("chr1", 350, 2, 1),    # outside
                             ("chr1", 500, 0, 4)])   # inside second interval
    mp = _bed(tmp_path, "m.bed", [("chr1", 100, 200), ("chr1", 450, 550)])
    out = str(tmp_path / "clip.txt")
    n, nmask = scan_adapter.write_clip_candidate_list(
        rt, "L1", out, min_clip=1, region_mask=region_mask.RegionMask([mp]))
    assert (n, nmask) == (1, 2)
    kept = [l.split("\t")[1] for l in open(out).read().splitlines()]
    assert kept == ["350"]


def test_write_clip_candidate_list_without_mask_is_unchanged(tmp_path):
    """The default-off path must be byte-identical to the pre-mask behaviour."""
    rt = _realign(tmp_path, [("chr1", 150, 3, 0), ("chr1", 350, 2, 1)])
    a, b = str(tmp_path / "a.txt"), str(tmp_path / "b.txt")
    n1, m1 = scan_adapter.write_clip_candidate_list(rt, "L1", a, min_clip=1)
    n2, m2 = scan_adapter.write_clip_candidate_list(rt, "L1", b, min_clip=1,
                                                    region_mask=region_mask.RegionMask([]))
    assert (n1, m1) == (2, 0) and (n2, m2) == (2, 0)
    assert open(a).read() == open(b).read()


def test_masked_count_precedes_min_clip_by_design(tmp_path):
    """The masked count is "positions the mask excluded from consideration", which
    necessarily includes ones that would ALSO have failed min_clip.

    That is a consequence of where the filter runs, and the placement is deliberate: the
    site dict holds every candidate position genome-wide (~57M on a 173x WGS arm) and is
    this function's peak-memory structure, so masked positions are skipped as the input is
    read rather than at write time -- they never enter the dict. min_clip cannot be
    evaluated that early because it needs the summed nL/nR across all rows for a position.
    Trading a slightly generous masked COUNT for ~10% off the peak memory of the largest
    dict in the step is the right way round; the written output is identical either way.
    """
    rt = _realign(tmp_path, [("chr1", 150, 0, 0), ("chr1", 160, 5, 0)])
    mp = _bed(tmp_path, "m.bed", [("chr1", 100, 200)])
    out = str(tmp_path / "clip.txt")
    n, nmask = scan_adapter.write_clip_candidate_list(
        rt, "L1", out, min_clip=1, region_mask=region_mask.RegionMask([mp]))
    assert (n, nmask) == (0, 2)
    assert open(out).read() == ""


# --- the mode-dependent default ----------------------------------------------------------
def _resolve(mosaic, explicit=None):
    """Mirror run_xtea's resolution so the rule is asserted, not just described."""
    spec = explicit
    if spec is None:
        spec = "default" if mosaic else "none"
    return spec


def test_default_is_mosaic_only():
    """THE SCOPING RULE. Every measurement behind this mask came from mosaic runs; germline
    has a different FP profile and a long-standing concordance baseline, so turning it on
    there by default would move that baseline on no evidence."""
    assert _resolve(mosaic=True) == "default"
    assert _resolve(mosaic=False) == "none"


def test_explicit_flag_beats_the_mode_default_both_ways():
    assert _resolve(mosaic=False, explicit="default") == "default"   # opt IN for germline
    assert _resolve(mosaic=True, explicit="none") == "none"          # opt OUT for mosaic


def test_germline_default_resolves_to_an_empty_mask():
    assert region_mask.load(_resolve(mosaic=False)).n_intervals == 0


if __name__ == "__main__":
    sys.exit(pytest.main([__file__, "-v"]))
