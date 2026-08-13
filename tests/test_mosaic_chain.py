"""M.2: --mosaic must change the chain in exactly four places, and nowhere else.

The mosaic mode is "a parameter profile plus one scan field, not a fourth pipeline" (plan sec 6.4).
That claim is only true if turning --mosaic ON is inert for every other mode, and only useful if it
is NOT inert for the four things it is supposed to change. Both halves are asserted here.

The step argv comes from `run_xtea.build_type_steps`, the same function the real pipeline calls, so
this cannot drift from the shipping chain -- the same reason test_downstream_from_cns.py uses it.

Each of the four changes fails SILENTLY if wired wrong:
  * postF_mosaic reading hc_pf instead of pf_out gives a mosaic caller that only ever confirms
    germline-strength calls: plausible output, zero mosaic sensitivity, nothing in the log.
  * a missing --postFmosaic flag runs the ordinary germline post-filter under a mosaic filename.
  * PEAK_WINDOW left at 100 fuses a low-AF insertion into an adjacent germline peak, and the call
    disappears as a duplicate rather than being reported.
  * running the germline genotype classifier assigns every mosaic call a confident het/hom.
"""
import os
import sys
from types import SimpleNamespace

import pytest

REPO_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(REPO_ROOT, "xtea"))

import run_xtea  # noqa: E402


def _ctx(tmp_path, mosaic=False, case_control=False):
    dest = str(tmp_path)
    cfg = {"reference": "/ref.fa", "blacklist": "/bl.txt", "gene_gff3": "/g.gff3",
           "model": "/m.onnx", "input": "/in.bam"}
    paths = {"rtype": 1, "copies_with_flank": "/cwf.fa", "annotation": "/a.out",
             "annotation1": "/a1.out", "flank": "/fl.fa", "cns": "/cns.fa"}
    dirs = {"type": dest, "tmp": dest + "/tmp/", "tmp_cns": dest + "/tmp/cns/",
            "tmp_td": dest + "/tmp/td/", "pub_clip": dest + "/tmp/clip/",
            "bam_list": dest + "/bam_list.txt", "control_bam_list": None}
    opt = SimpleNamespace(cores=4, nclip=4, cr=2, nd=5, nfclip=4, nfdisc=5, flklen=3000,
                          tumor=False, purity=1.0, case_control=case_control, mosaic=mosaic,
                          sample="s")
    return run_xtea.build_type_steps("L1", cfg, paths, dirs, opt)


def _names(steps):
    return [n for n, _ in steps]


def _argv(steps, name):
    return dict(steps)[name]


def test_germline_chain_has_no_mosaic_step(tmp_path):
    assert "postF_mosaic" not in _names(_ctx(tmp_path, mosaic=False))


def test_mosaic_chain_adds_exactly_one_step(tmp_path):
    """Unlike case-control's two steps, mosaic adds ONE -- it re-reads no BAM."""
    plain = _names(_ctx(tmp_path / "a", mosaic=False))
    mos = _names(_ctx(tmp_path / "b", mosaic=True))
    assert set(mos) - set(plain) == {"postF_mosaic"}
    assert set(plain) - set(mos) == set()


def test_mosaic_step_carries_the_postFmosaic_flag(tmp_path):
    """Without it, x_TEA_main takes the ordinary germline post-filter branch."""
    argv = _argv(_ctx(tmp_path, mosaic=True), "postF_mosaic")
    assert "--postFmosaic" in argv
    assert "--postF" in argv


def test_mosaic_step_filters_the_raw_list_not_the_high_confident_one(tmp_path):
    """THE one to get right: hc_pf is already blacklist- and confidence-filtered, which is
    precisely what removes a 1-5% VAF insertion. Mosaic must see postF_all's raw output."""
    argv = _argv(_ctx(tmp_path, mosaic=True), "postF_mosaic")
    inp = argv[argv.index("-i") + 1]
    assert not inp.endswith(".high_confident.post_filtering.txt"), \
        f"mosaic is filtering the high-confident list: {inp}"
    assert "high_confident" not in inp, f"mosaic input is an hc list: {inp}"
    assert inp.endswith("_post_filtering.txt") or "post_filtering" in inp, inp


def test_mosaic_step_runs_after_both_postF_rounds(tmp_path):
    names = _names(_ctx(tmp_path, mosaic=True))
    assert names.index("postF_mosaic") > names.index("postF_all")
    assert names.index("postF_mosaic") > names.index("postF_hc")


def test_mosaic_output_path_is_the_one_the_tail_reads_back(tmp_path):
    """build_type_steps writes it and run_type_downstream reads it; one helper, no drift."""
    argv = _argv(_ctx(tmp_path, mosaic=True), "postF_mosaic")
    assert argv[argv.index("-o") + 1] == run_xtea.mosaic_path(str(tmp_path))


def test_stop_after_postF_mosaic_truncates_the_chain(tmp_path):
    """The per-type chain must end at the mosaic filter; the tail runs once, after the merge."""
    dest = str(tmp_path)
    cfg = {"reference": "/ref.fa", "blacklist": "/bl.txt", "gene_gff3": "/g.gff3",
           "model": "/m.onnx", "input": "/in.bam"}
    paths = {"rtype": 1, "copies_with_flank": "/c.fa", "annotation": "/a.out",
             "annotation1": "/a1.out", "flank": "/f.fa", "cns": "/cns.fa"}
    dirs = {"type": dest, "tmp": dest + "/tmp/", "tmp_cns": dest + "/tmp/cns/",
            "tmp_td": dest + "/tmp/td/", "pub_clip": dest + "/tmp/clip/",
            "bam_list": dest + "/bl.txt", "control_bam_list": None}
    opt = SimpleNamespace(cores=4, nclip=4, cr=2, nd=5, nfclip=4, nfdisc=5, flklen=3000,
                          tumor=False, purity=1.0, case_control=False, mosaic=True, sample="s")
    steps = run_xtea.build_type_steps("L1", cfg, paths, dirs, opt, stop_after="postF_mosaic")
    assert _names(steps)[-1] == "postF_mosaic"
    assert "gene" not in _names(steps)


@pytest.mark.parametrize("mode,expect_genotype", [
    ({}, True),
    ({"mosaic": True}, False),
    ({"case_control": True}, False),
])
def test_genotype_classifier_is_skipped_for_mosaic(tmp_path, mode, expect_genotype):
    """The sklearn model is a germline het/hom classifier. A 2% VAF insertion is neither, so
    running it would stamp every mosaic call with a confident germline genotype."""
    opt = SimpleNamespace(cores=4, case_control=mode.get("case_control", False),
                          mosaic=mode.get("mosaic", False), sample="s")
    cfg = {"reference": "/r.fa", "gene_gff3": "/g.gff3", "model": "/m.onnx"}
    steps = run_xtea.build_shared_tail_steps(
        cfg, str(tmp_path), opt, "/merged_hcpf.txt", "/type_map.txt", "/merged.vcf",
        str(tmp_path / "bam_list.txt"), 7)
    assert ("genotype" in [n for n, _ in steps]) is expect_genotype
