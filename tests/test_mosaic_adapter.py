"""M.1: the mosaic candidate list must match master's contract, including its edge cases.

Three behaviours here are easy to get backwards, and each one silently produces a plausible
candidate list rather than an error:

  * the cutoff is AND, not OR. A mosaic insertion at 1-5% VAF often has clip support on ONE side
    only; an OR would drop exactly the sites this mode exists to find.
  * a site with raw clip support but no realigned clip is KEPT with zeros. master's own miss
    branch writes "0\t0\t0\t0" and the post-filter expects those rows to exist.
  * a 6-column (germline) realign-out is a hard error. Zero-filling it instead would give every
    site npL=npR=0, which reads as "no polyA anywhere" -- indistinguishable from a real result.
"""
import os
import sys

import pytest

REPO_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(REPO_ROOT, "xtea"))

import mosaic_adapter as ma  # noqa: E402


def _write(p, rows):
    p.write_text("".join("\t".join(str(x) for x in r) + "\n" for r in rows))
    return str(p)


@pytest.fixture
def inputs(tmp_path):
    clip = _write(tmp_path / "clip_sites.txt", [
        ("chr1", 100, 5, 4),     # both sides strong
        ("chr1", 200, 3, 0),     # LEFT only  -> must survive the AND rule
        ("chr1", 300, 0, 3),     # RIGHT only -> must survive the AND rule
        ("chr1", 400, 1, 1),     # both below cutoff -> dropped
        ("chr1", 500, 900, 900), # pileup artefact -> dropped by max_cov
        ("chr1", 600, 4, 4),     # no realigned row -> kept with zeros
    ])
    realign = _write(tmp_path / "realign.txt", [
        ("chr1", 100, "LINE1", 3, 2, 7, 2, 1),
        ("chr1", 200, "LINE1", 2, 0, 5, 2, 0),
        ("chr1", 300, "LINE1", 0, 2, 4, 0, 2),
        ("chr1", 100, "ALU",   9, 9, 9, 9, 9),   # different family -> must not leak into L1
    ])
    return clip, realign, tmp_path


def _rows(path):
    return [ln.split("\t") for ln in open(path).read().splitlines() if ln]


def test_one_sided_support_survives(inputs):
    """The cutoff is AND: one strong side is enough. An OR here loses the mosaic case."""
    clip, realign, tmp = inputs
    out = str(tmp / "cand.txt")
    ma.write_mosaic_candidate_list(clip, realign, "L1", out,
                                   cutoff_left_clip=2, cutoff_right_clip=2, max_cov_cutoff=1000)
    kept = {r[1] for r in _rows(out)}
    assert "200" in kept, "left-only site dropped -- cutoff applied as OR"
    assert "300" in kept, "right-only site dropped -- cutoff applied as OR"


def test_below_cutoff_and_above_coverage_are_dropped(inputs):
    clip, realign, tmp = inputs
    out = str(tmp / "cand.txt")
    n, low, high = ma.write_mosaic_candidate_list(
        clip, realign, "L1", out, cutoff_left_clip=2, cutoff_right_clip=2, max_cov_cutoff=1000)
    kept = {r[1] for r in _rows(out)}
    assert "400" not in kept and low == 1
    assert "500" not in kept and high == 1
    assert n == 4


def test_site_without_a_realigned_clip_is_kept_with_zeros(inputs):
    """master's miss branch writes 0\t0\t0\t0 -- the row must exist, not vanish."""
    clip, realign, tmp = inputs
    out = str(tmp / "cand.txt")
    ma.write_mosaic_candidate_list(clip, realign, "L1", out)
    row = [r for r in _rows(out) if r[1] == "600"]
    assert row, "site with raw support but no realignment was dropped"
    assert row[0][4:] == ["0", "0", "0", "0"]


def test_other_families_do_not_leak(inputs):
    """chr1:100 has an ALU row with 9s; the L1 list must show only the LINE1 numbers."""
    clip, realign, tmp = inputs
    out = str(tmp / "cand.txt")
    ma.write_mosaic_candidate_list(clip, realign, "L1", out)
    row = [r for r in _rows(out) if r[1] == "100"][0]
    assert row[4:] == ["3", "2", "2", "1"], f"family leak: {row}"


def test_polya_columns_are_carried_through(inputs):
    """npL/npR are the point of the whole mode; they must reach the candidate list intact."""
    clip, realign, tmp = inputs
    out = str(tmp / "cand.txt")
    ma.write_mosaic_candidate_list(clip, realign, "L1", out)
    by_pos = {r[1]: r for r in _rows(out)}
    assert by_pos["200"][6:] == ["2", "0"]
    assert by_pos["300"][6:] == ["0", "2"]


def test_six_column_realign_out_is_a_hard_error(tmp_path):
    """A germline realign-out must not be silently accepted as 'no polyA anywhere'."""
    clip = _write(tmp_path / "clip.txt", [("chr1", 100, 5, 4)])
    germline = _write(tmp_path / "realign6.txt", [("chr1", 100, "LINE1", 3, 2, 7)])
    with pytest.raises(ma.MosaicPolyAMissing) as e:
        ma.write_mosaic_candidate_list(clip, germline, "L1", str(tmp_path / "o.txt"))
    assert "--mosaic-polya" in str(e.value)


def test_output_is_position_sorted(inputs):
    clip, realign, tmp = inputs
    out = str(tmp / "cand.txt")
    ma.write_mosaic_candidate_list(clip, realign, "L1", out)
    pos = [(r[0], int(r[1])) for r in _rows(out)]
    assert pos == sorted(pos)


# ---------------------------------------------------------------------------------------------
# The mosaic AF-cutoff table must be TOTAL over insertion subtypes.
#
# Found by running it, not by reading it: job 50434768 reached L1:postF_mosaic after 2h45m and
# died with KeyError: 'orphan_or_sibling_transduction' (x_post_filter.py:1125). The table is
# built from AFConflictFilter.get_rep_type(), which enumerates 13 subtypes and omits that one --
# though global_values.ORPHAN_TRANSDUCTION defines it. It was 20 of 32 L1 records on the chr1
# HapMap run, i.e. the DOMINANT subtype, so mosaic crashed on essentially the first record.
# ---------------------------------------------------------------------------------------------
import importlib  # noqa: E402

_mos = importlib.import_module("x_mosaic_calling")


def test_unenumerated_subtype_does_not_raise():
    """The exact crash: a subtype get_rep_type() does not list must still resolve."""
    c = _mos.MosaicCaller("/tmp", 1)
    m = c.get_AF_cutoff(["two_side", "one_side"])
    assert m["orphan_or_sibling_transduction"] == _mos.MosaicCaller.MOSAIC_AF_BOUNDS


def test_unknown_subtype_gets_the_same_bounds_as_every_known_one():
    """Defaulting is only sound because every branch of master's chain assigns one pair.
    If that ever stops being true, this test fails and the default must be revisited."""
    c = _mos.MosaicCaller("/tmp", 1)
    known = ["two_side", "one_side", "one_half_side", "something_else"]
    m = c.get_AF_cutoff(known)
    assert len({m[k] for k in known}) == 1, "per-type AF bounds now differ; the default is unsafe"
    assert m["a_subtype_nobody_has_seen"] == m["two_side"]


def test_the_gap_is_reported_not_swallowed(capsys):
    """A silent default would hide a real enumeration gap; it must warn, once per subtype."""
    c = _mos.MosaicCaller("/tmp", 1)
    m = c.get_AF_cutoff(["two_side"])
    _ = m["orphan_or_sibling_transduction"]
    _ = m["orphan_or_sibling_transduction"]
    err = capsys.readouterr().err
    assert err.count("orphan_or_sibling_transduction") == 1, f"expected one warning, got: {err}"
    assert "not enumerated" in err


def test_orphan_transduction_really_is_absent_from_get_rep_type():
    """Pins the root cause. If it is ever added upstream this fails, and the default can go."""
    from x_post_filter import AFConflictFilter
    import global_values
    assert global_values.ORPHAN_TRANSDUCTION not in AFConflictFilter("/tmp", 1).get_rep_type()


# ---------------------------------------------------------------------------------------------
# The mosaic blacklist step must honour the "null" sentinel.
#
# "null" is the DEFAULT for both entry points (run_xtea.py:1282, x_TEA_main.py:211) and every
# other post-filter step tolerates it -- postF_hc ran with `--blacklist null` in the same chain
# and succeeded. MosaicCaller.filter_by_blacklist passed it straight to open():
#
#     FileNotFoundError: [Errno 2] No such file or directory: 'null'
#
# job 50435713 (broad WGS), 4h07m in, on the very LAST step. Any mosaic run that does not
# explicitly name a blacklist dies after all the expensive work is already done.
# ---------------------------------------------------------------------------------------------

@pytest.mark.parametrize("sentinel", ["null", "NULL", ""])
def test_blacklist_sentinel_passes_records_through(tmp_path, sentinel):
    c = _mos.MosaicCaller(str(tmp_path), 1)
    src = tmp_path / "in.txt"
    src.write_text("chr1\t1000\tx\ty\nchr2\t2000\ta\tb\n")
    out = tmp_path / "out.txt"
    c.filter_by_blacklist(str(src), sentinel, 100, str(out))
    assert out.exists(), f"sentinel {sentinel!r} produced no output"
    assert len(out.read_text().splitlines()) == 2, "records were dropped, not passed through"


def test_missing_blacklist_file_warns_but_does_not_crash(tmp_path, capsys):
    """A named-but-absent blacklist is a caller error; losing 4 hours of work to it is worse."""
    c = _mos.MosaicCaller(str(tmp_path), 1)
    src = tmp_path / "in.txt"
    src.write_text("chr1\t1000\tx\ty\n")
    out = tmp_path / "out.txt"
    c.filter_by_blacklist(str(src), str(tmp_path / "nope.bed"), 100, str(out))
    assert out.exists()
    assert "does not exist" in capsys.readouterr().err


def test_sentinel_does_not_substitute_a_default_blacklist(tmp_path):
    """Passing through is not the same as quietly filtering against some default file:
    that would change the call set unasked, and change what precision measures."""
    c = _mos.MosaicCaller(str(tmp_path), 1)
    src = tmp_path / "in.txt"
    src.write_text("chr1\t1000\tx\ty\nchr1\t1001\tx\ty\nchr2\t2000\ta\tb\n")
    out = tmp_path / "out.txt"
    c.filter_by_blacklist(str(src), "null", 100, str(out))
    assert len(out.read_text().splitlines()) == 3, "sentinel path filtered something"
