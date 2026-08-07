"""Stage 3 check 2b, the claim handoff sec 5.2 asks its reader to CHECK rather than take.

The sequence sec 5.2 flags out loud is: the gate failed, the rule was amended, the gate passed --
and the amendment was written by the same session that ran the gate. What makes it a rule change
rather than tuning is that the hazard the old rule protected is preserved exactly: **a dropped
entry that the OFF arm's reads VERIFIED must still FAIL.** That is untestable while the logic
lives in a `python - <<EOF` heredoc inside an sbatch, which is why `benchmark/adjudicate_dropped.py`
exists. These are the tests that make the claim checkable.

Three branches, and all three matter:
  1. VERIFIED drop           -> FAIL. The preserved hazard.
  2. UNDECIDABLE/REFUTED/... -> PASS, with the count still printed. The amendment itself.
  3. not adjudicated at all  -> FAIL. It cannot be waived on evidence nobody has.
"""
import subprocess
import sys
from pathlib import Path

import pytest

REPO = Path(__file__).resolve().parents[1]
SCRIPT = REPO / "benchmark" / "adjudicate_dropped.py"

sys.path.insert(0, str(REPO / "benchmark"))
import adjudicate_dropped as ad  # noqa: E402


# The current writer's layout. Deliberately NOT minimal: `shape` and `sidedness` sit between
# `td_src_raw` and `verdict` exactly as the restratification left them, so a reader that took
# columns by index would land on the wrong field and these tests would catch it.
HEADER = ["family", "ins_chrom", "ins_pos", "td_src_raw", "n_clip_support", "n_disc_support",
          "shape", "sidedness", "verdict", "why"]


def _verdicts(tmp_path, rows, header=HEADER, name="verdicts.tsv"):
    p = tmp_path / name
    with open(p, "w") as fh:
        fh.write("# a comment line the reader must skip\n")
        fh.write("\t".join(header) + "\n")
        for r in rows:
            fh.write("\t".join(r) + "\n")
    return str(p)


def _row(family, chrom, pos, src, verdict, why="", n_clip="143", n_disc="9",
         shape="regular", sided="two_sided"):
    return [family, chrom, pos, src, n_clip, n_disc, shape, sided, verdict, why]


def _lost(tmp_path, keys, name="lost_keys.txt"):
    """check 2a writes `cut -f1-8,11` of the pool: 9 fields, td_src_raw last."""
    p = tmp_path / name
    with open(p, "w") as fh:
        for (family, chrom, pos, src) in keys:
            fh.write("\t".join([family, chrom, pos, "x", "x", "x", "x", "x", src]) + "\n")
    return str(p)


def _run(lost, verdicts):
    r = subprocess.run([sys.executable, str(SCRIPT), lost, verdicts],
                       capture_output=True, text=True)
    return r.returncode, r.stdout + r.stderr


# --------------------------------------------------------------------------------------
# 1. The preserved hazard: a dropped entry the OFF arm VERIFIED must FAIL.
# --------------------------------------------------------------------------------------

def test_a_dropped_VERIFIED_entry_still_FAILS_the_gate(tmp_path):
    """THE claim sec 5.2 asks the reader to check. If this ever passes, the 2026-08-07
    amendment stopped being a rule change and became a way to make the gate green."""
    v = _verdicts(tmp_path, [_row("SVA", "chr17", "45585847", "chr17:45583245~1~regular",
                                  "VERIFIED", "9 disc reads")])
    lost = _lost(tmp_path, [("SVA", "chr17", "45585847", "chr17:45583245~1~regular")])
    rc, out = _run(lost, v)
    assert rc == 1, "a dropped VERIFIED entry MUST fail the gate:\n" + out
    assert "FAILS: the OFF arm's reads VERIFIED this" in out
    assert "FAIL:" in out


def test_one_VERIFIED_drop_among_waivable_ones_still_FAILS(tmp_path):
    """The realistic shape of the hazard: the VERIFIED entry does not arrive alone, and a
    check that only looked at the first or the count would wave this through."""
    v = _verdicts(tmp_path, [
        _row("SVA", "chr17", "45585847", "srcA", "UNDECIDABLE", "self_window(D@+2602)"),
        _row("SVA", "chr6", "31329020", "srcB", "VERIFIED", "76 disc reads"),
        _row("L1", "chr1", "1000", "srcC", "REFUTED", "no support"),
    ])
    lost = _lost(tmp_path, [("SVA", "chr17", "45585847", "srcA"),
                            ("SVA", "chr6", "31329020", "srcB"),
                            ("L1", "chr1", "1000", "srcC")])
    rc, out = _run(lost, v)
    assert rc == 1, out
    assert "3 dropped hypothes(es)" in out, "every drop must still be listed:\n" + out
    assert out.count("<-- FAILS") == 1, "exactly the VERIFIED one fails:\n" + out


# --------------------------------------------------------------------------------------
# 2. The amendment: un-verified drops are reported and waived.
# --------------------------------------------------------------------------------------

@pytest.mark.parametrize("verdict", ["UNDECIDABLE", "REFUTED", "SUB_THRESHOLD"])
def test_an_unverified_drop_PASSES_and_is_still_counted(tmp_path, verdict):
    """The 2026-08-07 rule. `--td-decoy-concat`'s real case is the UNDECIDABLE one:
    chr17:45585847, which the truth set counts a false positive."""
    v = _verdicts(tmp_path, [_row("SVA", "chr17", "45585847", "srcA", verdict,
                                  "self_window(D@+2602)")])
    lost = _lost(tmp_path, [("SVA", "chr17", "45585847", "srcA")])
    rc, out = _run(lost, v)
    assert rc == 0, "an un-VERIFIED drop must be waived:\n" + out
    assert "1 dropped hypothes(es)" in out, (
        "THE COUNT PRINTS ON PASS TOO -- it is the only signal for a change that drops many "
        "un-verified but real hypotheses:\n" + out)
    assert verdict in out
    assert "PASS:" in out


def test_nothing_dropped_is_a_pass_with_no_adjudication(tmp_path):
    """check 2a's purely-additive case. The file exists and holds a blank line, because `comm`
    produced nothing and the caller `printf`ed it -- so "empty" is not an empty file."""
    v = _verdicts(tmp_path, [_row("SVA", "chr17", "45585847", "srcA", "VERIFIED")])
    p = tmp_path / "lost_keys.txt"
    p.write_text("\n")
    rc, out = _run(str(p), v)
    assert rc == 0, out
    assert "nothing was dropped" in out


# --------------------------------------------------------------------------------------
# 3. Silence is not a pass.
# --------------------------------------------------------------------------------------

def test_a_drop_absent_from_the_verdicts_FAILS(tmp_path):
    """It cannot be waived on evidence nobody has. This is the branch that separates "the reads
    say it is fine" from "we never asked the reads"."""
    v = _verdicts(tmp_path, [_row("L1", "chr1", "1000", "srcC", "REFUTED")])
    lost = _lost(tmp_path, [("SVA", "chr17", "45585847", "srcMISSING")])
    rc, out = _run(lost, v)
    assert rc == 1, "an unadjudicable drop MUST fail:\n" + out
    assert "cannot adjudicate, so it cannot be waived" in out


# --------------------------------------------------------------------------------------
# The sec 4.1 lesson: by NAME, never by index -- and a wrong file is not a gate failure.
# --------------------------------------------------------------------------------------

def test_the_frozen_16_column_layout_reads_identically(tmp_path):
    """`verdict` is field 13 in the frozen verdicts.tsv and 15 in the current writer's, because
    the restratification added `shape` and `sidedness` after the freeze. Same verdict, different
    position: an index-based reader compares n_clip_support against verdict and can report
    anything at all. Both layouts must give the same answer."""
    frozen_hdr = [c for c in HEADER if c not in ("shape", "sidedness")]
    frozen_row = [c for c, name in zip(_row("SVA", "chr17", "45585847", "srcA", "VERIFIED"),
                                       HEADER) if name not in ("shape", "sidedness")]
    v_frozen = _verdicts(tmp_path, [frozen_row], header=frozen_hdr, name="frozen.tsv")
    v_now = _verdicts(tmp_path, [_row("SVA", "chr17", "45585847", "srcA", "VERIFIED")])
    lost = _lost(tmp_path, [("SVA", "chr17", "45585847", "srcA")])
    assert _run(lost, v_frozen)[0] == _run(lost, v_now)[0] == 1


def test_a_verdicts_file_missing_a_required_column_is_LOUD_and_not_a_gate_failure(tmp_path):
    """Exit 2, not 1. "The harness is reading the wrong file" and "the gate failed" must not
    look alike -- that confusion is precisely what cost check 1 a spurious FAIL on 2026-08-06."""
    bad = [c for c in HEADER if c != "verdict"]
    v = _verdicts(tmp_path, [[c for c, n in zip(_row("SVA", "chr17", "1", "s", "VERIFIED"),
                                                HEADER) if n != "verdict"]], header=bad)
    lost = _lost(tmp_path, [("SVA", "chr17", "1", "s")])
    rc, out = _run(lost, v)
    assert rc == 2, out
    assert "NOT a gate failure" in out


def test_an_empty_verdicts_file_is_LOUD(tmp_path):
    """The OFF arm did not adjudicate. Every drop is then unwaivable, and reporting that as a
    plain gate FAIL would send the reader hunting for a regression that is not there."""
    p = tmp_path / "empty.tsv"
    p.write_text("")
    rc, out = _run(_lost(tmp_path, [("SVA", "chr17", "1", "s")]), str(p))
    assert rc == 2, out
    assert "did not adjudicate" in out


def test_the_fatal_verdict_string_has_exactly_one_definition():
    """The module's only public contract with verify_transduction.py. A second copy of the
    string is the private-copy-of-a-shipping-rule mistake plan sec 10 records three times."""
    assert ad.FATAL_VERDICT == "VERIFIED"
    src = (REPO / "benchmark" / "adjudicate_dropped.py").read_text()
    body = src.split('"""', 2)[2]  # skip the module docstring
    assert body.count('"VERIFIED"') == 1, "VERIFIED is hardcoded more than once outside the docs"
