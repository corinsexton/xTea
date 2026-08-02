"""Stage 2c.3: the family-assignment barrier's effect on the step chain.

The barrier splits the per-family downstream into two waves -- cns_filter, then
transduction onward -- with one cross-family decision in between. What could go
wrong is not the decision (tests/test_locus_table.py pins that rule); it is the
splitting:

  1. **The off path must be untouched.** `--family-assign-early` is opt-in
     precisely because it is the first change in the overhaul expected to move
     calls. If turning it OFF produced anything but the pre-2c chain, every
     agreement gate in docs/handoff/runs/ would be measuring the wrong thing.
  2. **The waves must partition the chain exactly.** A step run twice would be
     run against its own output; a step run zero times would vanish silently,
     and postF's inputs exist either way, so nothing downstream would raise.
  3. **A single --types entry must not take the split at all** -- there is no
     other family to contend with, and `--types L1` is what the Stage 0 gate
     compares against to attribute the 12 Alu-over-L1 conflicts.
"""
import os
import sys

REPO_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(REPO_ROOT, "xtea"))
import run_xtea  # noqa: E402


class _Opt:
    cores = 4
    nfclip = 3
    nfdisc = 5
    flklen = 3000
    tumor = False
    purity = 0.45
    case_control = False
    sample = "s"
    family_assign_early = False


def _steps(case_control=False, **kw):
    """Step names for one family's chain, with build_type_steps' wave arguments."""
    opt = _Opt()
    opt.case_control = case_control
    cfg = {"reference": "/ref.fa", "blacklist": "/bl.bed", "gene_gff3": "/g.gff3",
           "model": "/m.onnx"}
    paths = {"rtype": 1, "flank": "/flank.fa", "cns": "/cns.fa",
             "annotation": "/a.out", "annotation1": "/a1.out"}
    dirs = {"type": "/out/L1", "bam_list": "/bl.txt", "tmp_cns": "/tmp_cns",
            "tmp_td": "/tmp_td", "control_bam_list": "/cbl.txt"}
    return [name for name, _ in run_xtea.build_type_steps("L1", cfg, paths, dirs, opt, **kw)]


def test_waves_partition_the_chain_exactly():
    """wave1 ++ wave2 == the single-wave chain, with no step repeated or dropped."""
    single = _steps(stop_after="postF_hc")
    wave1 = _steps(stop_after="cns_filter")
    wave2 = _steps(start_from="transduction", stop_after="postF_hc")
    assert wave1 + wave2 == single
    assert wave1 == ["cns_filter"]
    assert "cns_filter" not in wave2


def test_waves_partition_the_case_control_chain_too():
    """--case-control stops one step later; the split must follow it rather than
    silently truncating the two somatic steps off the end of wave 2."""
    single = _steps(case_control=True, stop_after="case_control_hc")
    wave1 = _steps(case_control=True, stop_after="cns_filter")
    wave2 = _steps(case_control=True, start_from="transduction",
                   stop_after="case_control_hc")
    assert wave1 + wave2 == single
    assert wave2[-1] == "case_control_hc"


def test_off_is_one_wave_and_is_the_pre_2c_chain():
    assert run_xtea.downstream_waves(False, "postF_hc", 3) == [(None, "postF_hc")]


def test_on_is_two_waves_split_at_cns_filter():
    assert run_xtea.downstream_waves(True, "postF_hc", 3) == [
        (None, "cns_filter"), ("transduction", "postF_hc")]


def test_single_type_never_takes_the_barrier():
    """--types L1 alone has no other family to contend with, and it is the arm the
    Stage 0 gate diffs to attribute the 12 Alu-over-L1 conflicts."""
    assert run_xtea.downstream_waves(True, "postF_hc", 1) == [(None, "postF_hc")]


def test_default_is_off():
    """Opt-in until a WGS run has scored it against the truth set -- this is the first
    change in the overhaul that is EXPECTED to move calls."""
    src = open(os.path.join(REPO_ROOT, "xtea", "run_xtea.py")).read()
    i = src.index('"--family-assign-early"')
    assert "default=False" in src[i:i + 300]


def test_intermediate_wave_returns_no_artifact():
    """run_type_downstream returns the path the merge consumes. A wave that stops at
    cns_filter has no such path, and must return None rather than warn about a VCF
    nobody asked it to produce."""
    import inspect
    src = inspect.getsource(run_xtea.run_type_downstream)
    assert 'stop_after not in ("case_control_hc", "postF_hc")' in src, (
        "the intermediate-wave guard is gone; a cns_filter-only wave will fall through "
        "to the VCF branch and log a spurious WARNING every run")
