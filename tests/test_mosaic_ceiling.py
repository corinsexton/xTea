"""Per-sample mosaic AF_PE ceiling (x_mosaic_calling.per_sample_disc_upper, 2026-11-13).

disc_upper = MOSAIC_DISC_UPPER_K x the median AF_PE of the sample's germline het calls (the same
family's postF_hc list: raw SR and PE fractions >= 0.2, PE < 0.8), with the fixed 0.35 as the
fallback when there is no list or too few het calls. The step must log which one it used.
"""
import os
import sys

import pytest

REPO_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(REPO_ROOT, "xtea"))

import x_mosaic_calling as xm  # noqa: E402

MC = xm.MosaicCaller


def _row(clip, fmap, disc, conc, chrom="chr1", pos=1000):
    f = ["0"] * 49
    f[0], f[1] = chrom, str(pos)
    f[35], f[36], f[39], f[40] = str(clip), str(fmap), str(disc), str(conc)
    return "\t".join(f) + "\n"


def _hc(tmp_path, rows, name="hc.txt"):
    p = tmp_path / name
    p.write_text("".join(rows))
    return str(p)


def _hets(n, pe, sr=0.4):
    """n het-like rows with AF_PE pe (disc of 100 pairs) and AF_SR sr (of 100 reads)."""
    return [_row(int(sr * 100), 100 - int(sr * 100), int(pe * 100), 100 - int(pe * 100), pos=1000 + i)
            for i in range(n)]


def test_ceiling_is_k_times_the_median_het_afpe(tmp_path):
    rows = _hets(15, 0.40) + _hets(16, 0.50)          # median 0.50
    c, why = MC.per_sample_disc_upper(_hc(tmp_path, rows))
    assert c == pytest.approx(MC.MOSAIC_DISC_UPPER_K * 0.50)
    assert why.startswith("per-sample") and "31 germline het calls" in why


def test_hom_and_non_germline_calls_are_not_counted(tmp_path):
    """Hom calls (PE >= 0.8), low-SR calls (mosaic or background) and low-PE calls are left
    out, so they cannot pull the median."""
    rows = (_hets(25, 0.46)
            + _hets(30, 0.90)                         # hom
            + _hets(30, 0.45, sr=0.05)                # SR below 0.2
            + _hets(30, 0.10)                         # PE below 0.2
            + [_row(0, 0, 5, 5), _row(5, 5, 0, 0)])   # no SR or no PE pairs
    assert MC.germline_het_afpe(_hc(tmp_path, rows)) == pytest.approx([0.46] * 25)
    c, _ = MC.per_sample_disc_upper(_hc(tmp_path, rows, "b.txt"))
    assert c == pytest.approx(MC.MOSAIC_DISC_UPPER_K * 0.46)


@pytest.mark.parametrize("arg", [None, "null", ""])
def test_fallback_without_a_germline_list(arg):
    c, why = MC.per_sample_disc_upper(arg)
    assert c == MC.MOSAIC_AF_BOUNDS[2] == 0.35
    assert why.startswith("fixed fallback") and "no germline hc list" in why


def test_fallback_on_a_missing_file(tmp_path):
    c, why = MC.per_sample_disc_upper(str(tmp_path / "absent.txt"))
    assert c == 0.35
    assert "0 germline het calls" in why


def test_fallback_on_too_few_het_calls(tmp_path):
    n = MC.MOSAIC_MIN_GERMLINE_HET - 1
    c, why = MC.per_sample_disc_upper(_hc(tmp_path, _hets(n, 0.5) + _hets(50, 0.95)))
    assert c == 0.35
    assert "fixed fallback: %d germline het calls" % n in why
    c, _ = MC.per_sample_disc_upper(_hc(tmp_path, _hets(n + 1, 0.5), "c.txt"))
    assert c == pytest.approx(MC.MOSAIC_DISC_UPPER_K * 0.5)


def test_ceiling_reaches_the_cutoff_table():
    """The disc ceiling goes into every subtype's bounds, including unenumerated subtypes;
    clip bounds stay as they were."""
    c = MC("/tmp", 1)
    m = c.get_AF_cutoff(["two_side", "one_side"], 2, (0.35, 0.0, 0.29))
    assert m["two_side"] == (0.35, 0.0, 0.29)
    assert m["orphan_or_sibling_transduction"] == (0.35, 0.0, 0.29)
    af = xm.AFConflictFilter("/tmp", 1)
    assert af.is_ins_pass_mosaic_cutoff(m, "two_side", 0, 0, 0.05, 0.29)
    assert not af.is_ins_pass_mosaic_cutoff(m, "two_side", 0, 0, 0.05, 0.30)
    assert c.get_AF_cutoff(["two_side"])["two_side"] == MC.MOSAIC_AF_BOUNDS


def _run_stubbed(monkeypatch, sf_hc, tmp_path):
    """run_call_mosaic with the record loading and filtering stubbed out; returns the m_cutoff
    the family branch received."""
    seen = {}
    monkeypatch.setattr(xm.XTEARsltParser, "load_in_xTEA_rslt", lambda self, sf: [])
    monkeypatch.setattr(MC, "construct_interval_tree", lambda self, *a: None)
    monkeypatch.setattr(MC, "call_mosaic_Alu_from_bulk",
                        lambda self, l, p, x, af, ann, m, out: seen.setdefault("m", m))
    monkeypatch.setattr(MC, "filter_by_blacklist", lambda self, *a: None)
    MC(str(tmp_path), 1).run_call_mosaic("in.txt", "rmsk", 225, 2, "null", str(tmp_path / "o.txt"), sf_hc)
    return seen["m"]


def test_run_call_mosaic_uses_and_logs_the_per_sample_ceiling(monkeypatch, tmp_path, capsys):
    hc = _hc(tmp_path, _hets(40, 0.45))
    m = _run_stubbed(monkeypatch, hc, tmp_path)
    assert m["two_side"][2] == pytest.approx(MC.MOSAIC_DISC_UPPER_K * 0.45)
    out = capsys.readouterr().out
    assert "[mosaic] disc_upper (AF_PE ceiling) = 0.270 (per-sample" in out


def test_run_call_mosaic_logs_the_fallback(monkeypatch, tmp_path, capsys):
    m = _run_stubbed(monkeypatch, None, tmp_path)
    assert m["two_side"] == MC.MOSAIC_AF_BOUNDS
    assert "disc_upper (AF_PE ceiling) = 0.350 (fixed fallback" in capsys.readouterr().out
