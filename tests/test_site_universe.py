"""The site-creation rule, lifted out of the transduction step -- pinned against the code
it replaced.

WHY THIS GATE EXISTS
--------------------
`XTransduction.re_slct_with_clip_raw_disc_sites` decides which raw-disc rows become NEW
candidate sites, and it is the only thing that creates the 52 of BASELINE-0's 1,515 final
per-family records that have no `candidate_disc_filtered_cns.txt` row (50 L1, 2 SVA -- see
`docs/handoff/runs/probe_site_creation.py`). It moved to `xtea/site_universe.py` because it
needs no BAM, no bwa and no cns evidence, so a driver can run it for every family before
the per-family chains start -- the site-enumeration half of "unified site creation".

The whole safety argument for that move is "moved, not rewritten", so this file's centre is
`test_library_reproduces_the_replaced_implementation`: `_legacy_re_slct` below is a VERBATIM
copy of the loop that was deleted from `x_transduction.py`, including its use of
`XIntemediateSites.are_sites_close`, and the two are compared line-for-line over randomized
inputs. That is what makes the module's claim about the proximity window (half-open,
asymmetric, because `are_sites_close` walks `range(pos-iwin, pos+iwin)`) a verified property
of the original rather than a reading of it -- the oracle inherits the window from the real
helper, so if the arithmetic replacement got it wrong the comparison fails.

The remaining tests pin each rejection separately, so a failure says which one broke.
"""
import os
import random
import sys

REPO_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(REPO_ROOT, "xtea"))
import global_values  # noqa: E402
import site_universe  # noqa: E402
from x_intermediate_sites import XIntemediateSites  # noqa: E402
from x_transduction import XTransduction  # noqa: E402

L1_RTYPE = 1
SVA_RTYPE = 4


class FakeAnnotation:
    """The two `XAnnotation` methods the rule calls, over a plain {chrom: {pos: rcd}} map.

    Return shapes match `x_annotation.py:656` and `:702` exactly, including the misses:
    `(False, -1)` and `(-1, "", None, None, None)`. Each entry is
    (start, end, div_rate, sub_family); a query position hits the entry whose [start, end]
    it lands in.
    """

    def __init__(self, entries):
        self.entries = entries  # {chrom: [(start, end, div, subfamily)]}

    def is_within_repeat_region_interval_tree(self, chrm, pos):
        for (start, end, _div, _sf) in self.entries.get(chrm, ()):
            if start <= pos <= end:
                return True, start
        return False, -1

    def get_div_subfamily_with_min_div(self, chrm, pos):
        for (start, end, div, sf) in self.entries.get(chrm, ()):
            if start == pos:
                return div, sf, None, start, end
        return -1, "", None, None, None


def _legacy_re_slct(sf_raw_disc, sf_candidates, n_disc_cutoff, xannotation,
                    i_rep_type, b_tumor, sf_out, b_with_candidate=True):
    """VERBATIM copy of the body deleted from x_transduction.py at the lift.

    Do not tidy this. Its value is being byte-for-byte the previous behaviour, including
    the O(iwin) `are_sites_close` loop and the `is_matched_rep_type` name test (called on a
    real XTransduction so the delegation is exercised too).
    """
    xintmdt = XIntemediateSites()
    xtd = XTransduction("/tmp", 1, None)
    m_te_candidates = xintmdt.load_in_candidate_list_str_version(sf_candidates)
    with open(sf_raw_disc) as fin_raw, open(sf_out, "w") as fout_slct:
        if b_with_candidate == True:
            for ins_chrm in m_te_candidates:
                for ins_pos in m_te_candidates[ins_chrm]:
                    fout_slct.write(ins_chrm + "\t" + str(ins_pos) + "\n")
        for line in fin_raw:
            fields = line.split()
            ins_chrm = fields[0]
            ins_pos = int(fields[1])
            if (ins_chrm in m_te_candidates):
                b_close = False
                for tmp_pos in m_te_candidates[ins_chrm]:
                    if xintmdt.are_sites_close(ins_pos, tmp_pos, global_values.TWO_SITES_CLOSE_DIST) == True:
                        b_close = True
                        break
                if b_close == True:
                    continue
            n_ldisc = int(fields[2])
            n_rdisc = int(fields[3])
            n_half_cutoff = int(n_disc_cutoff / 2)
            if b_tumor == False:
                n_half_cutoff = int(n_disc_cutoff * 3 / 4)

            if n_ldisc < n_half_cutoff and n_rdisc < n_half_cutoff:
                continue

            b_in_rep, i_pos = xannotation.is_within_repeat_region_interval_tree(ins_chrm, int(ins_pos))
            div_rate, sub_family, family, pos_start, pos_end = xannotation.get_div_subfamily_with_min_div(
                ins_chrm, i_pos)
            b_match_rep = xtd.is_matched_rep_type(i_rep_type, sub_family)
            if (b_in_rep is True) and (b_match_rep is True) and (div_rate < global_values.TD_REP_DIVERGENT_CUTOFF):
                continue
            fout_slct.write(line.rstrip() + "\n")


def _write(path, rows):
    with open(path, "w") as f:
        for r in rows:
            f.write("\t".join(str(x) for x in r) + "\n")
    return path


def _raw(chrm, pos, nl, nr):
    """A raw_disc row as disc_adapter writes it: chrm pos nrawL nrawR nL nR."""
    return (chrm, pos, nl, nr, nl, nr)


def _cns(chrm, pos):
    """A cns candidate row -- load_in_candidate_list_str_version needs >= 3 fields."""
    return (chrm, pos, pos, "-1", "-1")


def _run(tmp_path, raw_rows, cns_rows, annotation, i_rep_type=L1_RTYPE,
         n_disc_cutoff=4, b_tumor=False, b_with_candidate=True):
    """Run the library form and return the output file's lines."""
    raw = _write(str(tmp_path / "raw_disc.txt"), raw_rows)
    cns = _write(str(tmp_path / "cns.txt"), cns_rows)
    out = str(tmp_path / "out.slct")
    site_universe.write_slct_file(raw, cns, n_disc_cutoff, annotation, i_rep_type,
                                 b_tumor, out, b_with_candidate=b_with_candidate,
                                 log=lambda *a, **k: None)
    with open(out) as f:
        return [ln.rstrip("\n") for ln in f]


# --- the load-bearing one --------------------------------------------------------------------

def test_library_reproduces_the_replaced_implementation(tmp_path):
    """Randomized differential test: library output == the deleted loop's output.

    Positions are drawn densely around a small number of cns candidates so the proximity
    window's edges (including exactly +-TWO_SITES_CLOSE_DIST) get hit often, and disc counts
    straddle the germline and tumor half-cutoffs. Both rtypes and both b_tumor values run.
    """
    rnd = random.Random(20260803)
    iwin = global_values.TWO_SITES_CLOSE_DIST
    for trial in range(12):
        cns_rows = []
        anchors = {}
        for chrm in ("chr1", "chr2", "chrX"):
            anchors[chrm] = sorted(rnd.randrange(10_000, 200_000) for _ in range(4))
            for p in anchors[chrm]:
                cns_rows.append(_cns(chrm, p))
        raw_rows = []
        for chrm in ("chr1", "chr2", "chrX", "chr3"):  # chr3 has no cns row at all
            base = anchors.get(chrm) or [50_000]
            for a in base:
                for delta in (-iwin - 1, -iwin, -iwin + 1, -1, 0, 1,
                              iwin - 1, iwin, iwin + 1, 5_000):
                    raw_rows.append(_raw(chrm, a + delta, rnd.randrange(0, 6), rnd.randrange(0, 6)))
        rnd.shuffle(raw_rows)
        # repeat entries covering roughly half the raw positions, mixed family and divergence
        entries = {}
        for (chrm, pos, _nl, _nr, _a, _b) in raw_rows:
            if rnd.random() < 0.5:
                sf = rnd.choice(["L1HS", "L1PA7", "SVA_D", "AluY", "MIR3"])
                div = rnd.choice([0.0, 2.5, 4.9, 5.0, 12.0])
                entries.setdefault(chrm, []).append((pos, pos + 300, div, sf))
        annotation = FakeAnnotation(entries)

        i_rep_type = L1_RTYPE if trial % 2 == 0 else SVA_RTYPE
        b_tumor = (trial % 4) >= 2
        n_disc_cutoff = 4 + (trial % 3)

        raw = _write(str(tmp_path / f"raw_{trial}.txt"), raw_rows)
        cns = _write(str(tmp_path / f"cns_{trial}.txt"), cns_rows)
        got = str(tmp_path / f"got_{trial}.slct")
        want = str(tmp_path / f"want_{trial}.slct")
        site_universe.write_slct_file(raw, cns, n_disc_cutoff, annotation, i_rep_type,
                                     b_tumor, got, log=lambda *a, **k: None)
        _legacy_re_slct(raw, cns, n_disc_cutoff, annotation, i_rep_type, b_tumor, want)
        with open(got) as fg, open(want) as fw:
            assert fg.read() == fw.read(), f"trial {trial}: lifted rule diverged from the original"


def test_step_method_delegates_to_the_library(tmp_path):
    """XTransduction.re_slct_with_clip_raw_disc_sites and site_universe.write_slct_file
    must produce the identical file -- the method is a delegation and nothing else."""
    annotation = FakeAnnotation({"chr1": [(1_000, 1_300, 1.0, "L1HS")]})
    raw_rows = [_raw("chr1", 1_050, 4, 4), _raw("chr1", 9_000, 4, 0), _raw("chr2", 7_000, 0, 5)]
    cns_rows = [_cns("chr1", 500)]
    raw = _write(str(tmp_path / "raw.txt"), raw_rows)
    cns = _write(str(tmp_path / "cns.txt"), cns_rows)
    via_method = str(tmp_path / "method.slct")
    via_lib = str(tmp_path / "lib.slct")
    XTransduction("/tmp", 1, None).re_slct_with_clip_raw_disc_sites(
        raw, cns, 4, annotation, L1_RTYPE, False, via_method)
    site_universe.write_slct_file(raw, cns, 4, annotation, L1_RTYPE, False, via_lib,
                                 log=lambda *a, **k: None)
    with open(via_method) as fm, open(via_lib) as fl:
        assert fm.read() == fl.read()


# --- one test per rejection ------------------------------------------------------------------

def test_candidate_proximity_window_is_half_open(tmp_path):
    """`are_sites_close(raw_pos, cand_pos, iwin)` walks range(raw_pos-iwin, raw_pos+iwin),
    so the exclusion window on the raw-disc site's offset from the candidate is
    `(-iwin, +iwin]`: exactly TWO_SITES_CLOSE_DIST to the LEFT of a candidate is KEPT,
    exactly that far to the RIGHT is suppressed.

    This asymmetry is not a nicety: it decides whether a site becomes a transduction
    candidate, so "fixing" it to a symmetric window changes calls.
    """
    annotation = FakeAnnotation({})
    iwin = global_values.TWO_SITES_CLOSE_DIST
    cand = 100_000
    rows = [_raw("chr1", cand + d, 4, 4) for d in (-iwin - 1, -iwin, -1, 0, iwin - 1, iwin, iwin + 1)]
    lines = _run(tmp_path, rows, [_cns("chr1", cand)], annotation)
    kept = {int(ln.split()[1]) for ln in lines[1:]}   # lines[0] is the candidate row itself
    assert cand - iwin - 1 in kept
    assert cand - iwin in kept              # offset -iwin is outside the half-open window
    assert cand not in kept
    assert cand - iwin + 1 not in kept and cand + iwin - 1 not in kept
    assert cand + iwin not in kept          # offset +iwin is inside it
    assert cand + iwin + 1 in kept


def test_proximity_exclusion_is_per_chromosome(tmp_path):
    """A candidate at the same coordinate on another contig must not suppress anything."""
    annotation = FakeAnnotation({})
    lines = _run(tmp_path, [_raw("chr2", 100_000, 4, 4)], [_cns("chr1", 100_000)], annotation)
    assert any(ln.split()[0] == "chr2" for ln in lines[1:])


def test_half_disc_cutoff_is_higher_for_germline_than_tumor(tmp_path):
    """germline int(nd*3/4) vs tumor int(nd/2); a site needs ONE side at or above it."""
    annotation = FakeAnnotation({})
    rows = [_raw("chr1", 10_000, 2, 2), _raw("chr1", 50_000, 3, 0), _raw("chr1", 90_000, 0, 3)]
    germline = {int(ln.split()[1]) for ln in _run(tmp_path, rows, [], annotation,
                                                  n_disc_cutoff=4, b_tumor=False)}
    tumor = {int(ln.split()[1]) for ln in _run(tmp_path, rows, [], annotation,
                                               n_disc_cutoff=4, b_tumor=True)}
    assert germline == {50_000, 90_000}          # 3 >= int(4*3/4)=3, but 2 < 3
    assert tumor == {10_000, 50_000, 90_000}     # 2 >= int(4/2)=2


def test_own_family_low_divergence_site_is_dropped_others_kept(tmp_path):
    """Rejection 3: the insertion point sitting inside a fresh copy of the SAME family is
    indistinguishable from reads mis-mapped out of that copy. A divergent copy of the same
    family, or any copy of a different family, is kept."""
    annotation = FakeAnnotation({"chr1": [
        (10_000, 10_300, 1.0, "L1HS"),     # own family, fresh -> dropped for rtype L1
        (20_000, 20_300, 12.0, "L1PA7"),   # own family, divergent -> kept
        (30_000, 30_300, 1.0, "AluY"),     # other family -> kept
        (40_000, 40_300, 1.0, "SVA_D"),    # other family for L1, own for SVA
    ]})
    rows = [_raw("chr1", p, 4, 4) for p in (10_000, 20_000, 30_000, 40_000)]
    as_l1 = {int(ln.split()[1]) for ln in _run(tmp_path, rows, [], annotation, i_rep_type=L1_RTYPE)}
    as_sva = {int(ln.split()[1]) for ln in _run(tmp_path, rows, [], annotation, i_rep_type=SVA_RTYPE)}
    assert as_l1 == {20_000, 30_000, 40_000}
    assert as_sva == {10_000, 20_000, 30_000}


def test_divergence_cutoff_boundary(tmp_path):
    """`div_rate < TD_REP_DIVERGENT_CUTOFF` is strict, so exactly 5% is kept."""
    annotation = FakeAnnotation({"chr1": [(10_000, 10_300, 4.999, "L1HS"),
                                          (20_000, 20_300, 5.0, "L1HS")]})
    rows = [_raw("chr1", 10_000, 4, 4), _raw("chr1", 20_000, 4, 4)]
    kept = {int(ln.split()[1]) for ln in _run(tmp_path, rows, [], annotation)}
    assert kept == {20_000}


def test_layout_candidates_first_then_raw_disc_order(tmp_path):
    """Candidates come first as bare chrom/pos pairs; new sites follow in raw_disc file
    order with all six fields. The td gather emits reads in site-file order and bwa's
    tie-break rides on it, so this ordering is not cosmetic (src/csrc/td_gather.hpp)."""
    annotation = FakeAnnotation({})
    rows = [_raw("chr1", 90_000, 4, 4), _raw("chr1", 10_000, 4, 4), _raw("chr2", 50_000, 4, 4)]
    lines = _run(tmp_path, rows, [_cns("chr1", 500), _cns("chr2", 700)], annotation)
    assert lines[0].split() == ["chr1", "500"]
    assert lines[1].split() == ["chr2", "700"]
    assert [ln.split()[:2] for ln in lines[2:]] == [
        ["chr1", "90000"], ["chr1", "10000"], ["chr2", "50000"]]
    assert len(lines[2].split()) == 6


def test_b_with_candidate_false_omits_the_candidate_block(tmp_path):
    """The sibling step calls this with the previous step's output as candidates and still
    wants only the new sites (x_TEA_main.py:622 -> `.slct2`)."""
    annotation = FakeAnnotation({})
    lines = _run(tmp_path, [_raw("chr1", 90_000, 4, 4)], [_cns("chr1", 500)], annotation,
                 b_with_candidate=False)
    assert [ln.split()[:2] for ln in lines] == [["chr1", "90000"]]


def test_is_matched_rep_type_has_no_alu_arm():
    """Alu (bit 2) matches nothing, deliberately: it has no transduction flank, so the
    transduction step never runs for it."""
    assert site_universe.is_matched_rep_type(1, "L1HS")
    assert site_universe.is_matched_rep_type(1, "LINE/L1")
    assert site_universe.is_matched_rep_type(4, "SVA_D")
    assert not site_universe.is_matched_rep_type(2, "AluY")
    assert not site_universe.is_matched_rep_type(1, "AluY")
    assert not site_universe.is_matched_rep_type(4, "L1HS")
    # a combined mask still matches either arm
    assert site_universe.is_matched_rep_type(7, "L1HS")
    assert site_universe.is_matched_rep_type(7, "SVA_F")


def test_read_positions_skips_unparseable_rows(tmp_path):
    p = str(tmp_path / "sites.txt")
    with open(p, "w") as f:
        f.write("chr1\t300\nchr1\t100\tx\nchr1\n\nchr2\tnotanint\nchr2\t7\n")
    assert site_universe.read_positions(p) == {"chr1": [100, 300], "chr2": [7]}
