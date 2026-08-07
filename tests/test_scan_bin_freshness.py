"""The xtea-scan binary must be FOUND when it is built, and REFUSED when it is stale.

WHY THIS FILE EXISTS. Two defects in how the suite reached the C++ tool, both of which made a
green run mean less than it looked:

  1. `xtea2_scan_bin` resolved only `$XTEA2_SCAN_BIN`, with no default. A binary built in-tree by
     the shipped Makefile was invisible, so a plain `pytest` skipped all 30 tests of the tool that
     owns `scan_pass` -- the largest step in a WGS run -- and said only "not built yet".
  2. Nothing compared the binary's mtime to its sources. The in-tree binary was found five days
     older than `xtea_scan.cpp`; several sbatch harnesses default `XTEA2_SCAN_BIN` to a *different
     commit's* worktree build. A stale binary is the one failure mode that lies: the suite passes
     while testing code that is not in the tree.

So the contract under test is three-valued, not two: present-and-fresh RUNS, absent SKIPS (not
built is legitimate, and the fast unit tests must still run), stale FAILS.

The staleness helper is checked against a synthetic `src/csrc/`, never the real one, so these
tests neither depend on nor perturb whatever is actually built in the working tree.
"""
import os

import pytest

from tests import conftest as ct


#: The Makefile's `SRC` + `HDR`, written out rather than read from `ct.SCAN_BIN_SOURCES` ON
#: PURPOSE. Deriving the expectation from the constant under test makes the test agree with any
#: mutation of it: shrinking SCAN_BIN_SOURCES to just the .cpp left a loop over one element that
#: passed happily. Measured -- that mutation survived until this list was hardcoded.
EXPECTED_SOURCES = ("xtea_scan.cpp", "realign.hpp", "cns_emit.hpp", "td_gather.hpp")


# --- helpers -------------------------------------------------------------------------------------

def _fake_csrc(tmp_path, bin_age=100.0, source_ages=None):
    """A synthetic repo root with src/csrc/. Ages are SECONDS BEFORE NOW, so larger = older.

    Returns (repo_root, binary_path). Sources default to older than the binary, i.e. fresh.
    """
    csrc = tmp_path / "src" / "csrc"
    csrc.mkdir(parents=True)
    now = 1_000_000_000.0  # fixed epoch: these tests compare mtimes, never wall-clock
    binary = csrc / "xtea_scan"
    binary.write_text("#!/bin/sh\n")
    os.utime(binary, (now - bin_age, now - bin_age))
    ages = source_ages if source_ages is not None else {n: 200.0 for n in EXPECTED_SOURCES}
    for name, age in ages.items():
        p = csrc / name
        p.write_text("// %s\n" % name)
        os.utime(p, (now - age, now - age))
    return str(tmp_path), str(binary)


def _call_fixture():
    """Invoke the session fixture's body directly, bypassing pytest's fixture machinery."""
    return ct.xtea2_scan_bin.__wrapped__()


# `pytest.fail` and `pytest.skip` raise OutcomeExceptions, which derive from BaseException and
# NOT from Exception -- so `pytest.raises(Exception)` would let them straight through and the
# assertion below it would never run. Name the two outcomes explicitly, which also makes
# "stale must FAIL, absent must SKIP" an assertion about the type rather than about a string.
SKIPPED = pytest.skip.Exception
FAILED = pytest.fail.Exception


def _expect_stale_failure():
    """Call the fixture, require a FAILED outcome, return its message.

    Catches BaseException rather than FAILED, and asserts the type afterwards, BECAUSE A LEAKED
    `Skipped` IS NOT A TEST FAILURE. With `pytest.raises(FAILED)` a fixture that skipped instead
    of failing would let `Skipped` propagate, pytest would report this test as *skipped*, and the
    run would still exit 0 -- so the mutation "stale skips instead of fails" survived. Measured,
    not hypothesised: it survived exactly that way before this helper existed.
    """
    with pytest.raises(BaseException) as e:      # noqa: PT011 -- the type IS the assertion below
        _call_fixture()
    assert not isinstance(e.value, SKIPPED), \
        "stale was reported as a SKIP; a skip is indistinguishable from `not built` and exits 0"
    assert isinstance(e.value, FAILED), \
        "expected a pytest FAILURE for a stale binary, got %s" % type(e.value).__name__
    return str(e.value)


# --- the staleness helper ------------------------------------------------------------------------

def test_fresh_binary_reports_no_stale_sources(tmp_path, monkeypatch):
    root, binary = _fake_csrc(tmp_path)
    monkeypatch.setattr(ct, "REPO_ROOT", root)
    assert ct._scan_bin_stale_against(binary) == []


def test_a_newer_header_alone_makes_the_binary_stale(tmp_path, monkeypatch):
    """The Makefile's own hazard: headers carry most of the code and are #included, so a
    header-only edit leaves `make` with nothing to do and the old binary in place."""
    ages = {"xtea_scan.cpp": 200.0, "realign.hpp": 10.0,
            "cns_emit.hpp": 200.0, "td_gather.hpp": 200.0}
    root, binary = _fake_csrc(tmp_path, bin_age=100.0, source_ages=ages)
    monkeypatch.setattr(ct, "REPO_ROOT", root)
    assert ct._scan_bin_stale_against(binary) == ["realign.hpp"]


def test_the_watched_source_list_is_the_makefiles():
    assert tuple(ct.SCAN_BIN_SOURCES) == EXPECTED_SOURCES


def test_every_source_is_checked_not_just_the_cpp(tmp_path, monkeypatch):
    """A check that watched only xtea_scan.cpp would pass three of these four edits."""
    for name in EXPECTED_SOURCES:
        sub = tmp_path / name.replace(".", "_")
        ages = {n: (10.0 if n == name else 200.0) for n in EXPECTED_SOURCES}
        root, binary = _fake_csrc(sub, bin_age=100.0, source_ages=ages)
        monkeypatch.setattr(ct, "REPO_ROOT", root)
        assert ct._scan_bin_stale_against(binary) == [name], name


def test_stale_sources_come_back_newest_first(tmp_path, monkeypatch):
    ages = {"xtea_scan.cpp": 30.0, "realign.hpp": 10.0,
            "cns_emit.hpp": 20.0, "td_gather.hpp": 200.0}
    root, binary = _fake_csrc(tmp_path, bin_age=100.0, source_ages=ages)
    monkeypatch.setattr(ct, "REPO_ROOT", root)
    assert ct._scan_bin_stale_against(binary) == ["realign.hpp", "cns_emit.hpp", "xtea_scan.cpp"]


def test_absent_source_tree_is_not_stale(tmp_path, monkeypatch):
    """An installed binary with no sources beside it cannot be checked; refusing to run it
    would be worse than running it."""
    binary = tmp_path / "xtea_scan"
    binary.write_text("#!/bin/sh\n")
    monkeypatch.setattr(ct, "REPO_ROOT", str(tmp_path / "nonexistent"))
    assert ct._scan_bin_stale_against(str(binary)) == []


def test_a_missing_individual_source_is_skipped_not_crashed(tmp_path, monkeypatch):
    ages = {"xtea_scan.cpp": 200.0, "realign.hpp": 10.0}  # cns_emit/td_gather absent
    root, binary = _fake_csrc(tmp_path, bin_age=100.0, source_ages=ages)
    monkeypatch.setattr(ct, "REPO_ROOT", root)
    assert ct._scan_bin_stale_against(binary) == ["realign.hpp"]


def test_equal_mtimes_are_fresh(tmp_path, monkeypatch):
    """Strictly-newer, not newer-or-equal: `make` can leave both at the same second, and
    failing on a tie would make a correct build unusable."""
    ages = {n: 100.0 for n in EXPECTED_SOURCES}
    root, binary = _fake_csrc(tmp_path, bin_age=100.0, source_ages=ages)
    monkeypatch.setattr(ct, "REPO_ROOT", root)
    assert ct._scan_bin_stale_against(binary) == []


# --- the fixture's three-valued contract ---------------------------------------------------------

def test_the_default_is_the_in_tree_build_path():
    """Defect 1: with no env var set, the module must still name the shipped build's output.
    This is what the old `os.environ.get("XTEA2_SCAN_BIN", "")` could never do."""
    assert os.environ.get("XTEA2_SCAN_BIN") or \
        ct.XTEA2_SCAN_BIN == os.path.join(ct.REPO_ROOT, "src", "csrc", "xtea_scan")


def test_a_built_and_fresh_binary_is_returned(tmp_path, monkeypatch):
    root, binary = _fake_csrc(tmp_path)
    monkeypatch.setattr(ct, "REPO_ROOT", root)
    monkeypatch.setattr(ct, "XTEA2_SCAN_BIN", binary)
    assert _call_fixture() == binary


def test_an_absent_binary_skips_and_names_the_build_command(tmp_path, monkeypatch):
    missing = str(tmp_path / "src" / "csrc" / "xtea_scan")
    monkeypatch.setattr(ct, "XTEA2_SCAN_BIN", missing)
    with pytest.raises(SKIPPED) as e:
        _call_fixture()
    msg = str(e.value)
    assert "make -C src/csrc" in msg and missing in msg


def test_a_stale_binary_FAILS_rather_than_skipping(tmp_path, monkeypatch):
    """The point of the whole file. A skip here would be indistinguishable from `not built`,
    and letting it run would be worse: a green suite against code not in the tree."""
    ages = {"xtea_scan.cpp": 10.0, "realign.hpp": 200.0,
            "cns_emit.hpp": 200.0, "td_gather.hpp": 200.0}
    root, binary = _fake_csrc(tmp_path, bin_age=100.0, source_ages=ages)
    monkeypatch.setattr(ct, "REPO_ROOT", root)
    monkeypatch.setattr(ct, "XTEA2_SCAN_BIN", binary)
    monkeypatch.delenv("XTEA2_ALLOW_STALE_SCAN_BIN", raising=False)
    msg = _expect_stale_failure()               # FAILED, not SKIPPED -- that is the whole claim
    assert "STALE" in msg and "xtea_scan.cpp" in msg


def test_the_stale_failure_names_the_rebuild(tmp_path, monkeypatch):
    ages = {"xtea_scan.cpp": 10.0, "realign.hpp": 200.0,
            "cns_emit.hpp": 200.0, "td_gather.hpp": 200.0}
    root, binary = _fake_csrc(tmp_path, bin_age=100.0, source_ages=ages)
    monkeypatch.setattr(ct, "REPO_ROOT", root)
    monkeypatch.setattr(ct, "XTEA2_SCAN_BIN", binary)
    monkeypatch.delenv("XTEA2_ALLOW_STALE_SCAN_BIN", raising=False)
    assert "make -C src/csrc" in _expect_stale_failure()


@pytest.mark.parametrize("value", ["1", "yes", "true"])
def test_the_escape_hatch_downgrades_stale_to_a_warning(tmp_path, monkeypatch, value):
    """For the deliberate case: a frozen binary against a moved tree, where the mismatch IS
    the intent. It must still say so out loud."""
    ages = {"xtea_scan.cpp": 10.0, "realign.hpp": 200.0,
            "cns_emit.hpp": 200.0, "td_gather.hpp": 200.0}
    root, binary = _fake_csrc(tmp_path, bin_age=100.0, source_ages=ages)
    monkeypatch.setattr(ct, "REPO_ROOT", root)
    monkeypatch.setattr(ct, "XTEA2_SCAN_BIN", binary)
    monkeypatch.setenv("XTEA2_ALLOW_STALE_SCAN_BIN", value)
    with pytest.warns(UserWarning, match="STALE"):
        assert _call_fixture() == binary


@pytest.mark.parametrize("value", ["", "0"])
def test_an_empty_or_zero_escape_hatch_still_fails(tmp_path, monkeypatch, value):
    """`XTEA2_ALLOW_STALE_SCAN_BIN=0` must not read as permission -- the same
    absent-or-"0"-means-off rule the driver's own env flags use."""
    ages = {"xtea_scan.cpp": 10.0, "realign.hpp": 200.0,
            "cns_emit.hpp": 200.0, "td_gather.hpp": 200.0}
    root, binary = _fake_csrc(tmp_path, bin_age=100.0, source_ages=ages)
    monkeypatch.setattr(ct, "REPO_ROOT", root)
    monkeypatch.setattr(ct, "XTEA2_SCAN_BIN", binary)
    monkeypatch.setenv("XTEA2_ALLOW_STALE_SCAN_BIN", value)
    assert "STALE" in _expect_stale_failure()


def test_an_explicit_env_var_still_wins_over_the_default(monkeypatch):
    """The sbatch harnesses depend on this: they point at other worktrees' builds."""
    import importlib
    monkeypatch.setenv("XTEA2_SCAN_BIN", "/somewhere/else/xtea_scan")
    reloaded = importlib.reload(ct)
    try:
        assert reloaded.XTEA2_SCAN_BIN == "/somewhere/else/xtea_scan"
    finally:
        monkeypatch.delenv("XTEA2_SCAN_BIN", raising=False)
        importlib.reload(ct)
