"""The driver-side locus universe -- the per-locus loop's iteration domain.

WHY THIS GATE EXISTS
--------------------
`xtea/downstream_unified.py` re-derives, in the driver, the site list each family's
transduction step derives for itself. That replica is only worth having if it is EXACTLY the
step's list, so the load-bearing test here is
`test_driver_enumeration_reproduces_the_step`: it runs `enumerate_family_sites` and
`XTransduction.re_slct_with_clip_raw_disc_sites` over the same inputs and compares the files
byte for byte. The production version of that same check is a `cmp` after a real run (see
`run_xtea.run_locus_universe`'s docstring); this is its unit form, and it is what would catch
the driver drifting from the step -- a wrong disc cutoff, a wrong annotation, a wrong
boundary extend.

The rest pin the parts that only exist driver-side: the union across families with its
per-family provenance (which family reached a locus by cns row and which by creation, the
distinction the 52-of-1,515 finding is about), the flankless family's cns-only contribution,
and the `iextnd` derivation the `--td-manifest` caller needs.

`FakeAnnotation` is imported from tests/test_site_universe.py rather than re-written, so the
two files cannot disagree about what the annotation returns on a miss.
"""
import os
import sys

REPO_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(REPO_ROOT, "xtea"))
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import downstream_unified  # noqa: E402
from x_transduction import XTransduction  # noqa: E402
from test_site_universe import FakeAnnotation, _raw, _write  # noqa: E402

L1_RTYPE = 1
SVA_RTYPE = 4


def _slct(path, cns=(), new=()):
    """A `.slct` file in the two layouts write_slct_file uses."""
    with open(path, "w") as f:
        for (chrm, pos) in cns:
            f.write("{0}\t{1}\n".format(chrm, pos))
        for row in new:
            f.write("\t".join(str(x) for x in row) + "\n")
    return path


# --------------------------------------------------------------------------- the replica


def test_driver_enumeration_reproduces_the_step(tmp_path):
    """The whole safety argument for enumerating driver-side, in one comparison.

    Same raw_disc, same cns list, same annotation, same cutoff -> the same bytes as the
    step's own `re_slct_with_clip_raw_disc_sites`. Rows are placed to exercise all three
    rejections and both sides of the half-open proximity window, so the comparison is not
    vacuous: a run where every row is trivially kept would pass with the cutoff wrong.
    """
    raw = _write(str(tmp_path / "raw_disc.txt"), [
        _raw("chr1", 1000, 9, 9),      # kept
        _raw("chr1", 5000, 9, 9),      # suppressed: a candidate sits at 5000
        _raw("chr1", 5100, 9, 9),      # exactly +100 from the candidate -> suppressed
        _raw("chr1", 4900, 9, 9),      # exactly -100 from the candidate -> kept
        _raw("chr1", 2000, 1, 1),      # rejected: neither side reaches the half cutoff
        _raw("chr2", 700, 9, 9),       # rejected: low-divergence L1 at that position
        _raw("chr3", 300, 0, 9),       # kept on the right side alone
    ])
    cns = _write(str(tmp_path / "cns.txt"), [("chr1", 5000, "x"), ("chr9", 42, "x")])
    ann = FakeAnnotation({"chr2": [(600, 800, 1.0, "L1HS")]})

    driver_out = str(tmp_path / "universe" / "L1.slct")
    step_out = str(tmp_path / "step.slct")
    downstream_unified.enumerate_family_sites(
        raw, cns, None, L1_RTYPE, 8, False, driver_out, xannotation=ann, log=lambda _m: None)
    XTransduction(str(tmp_path), 1, None).re_slct_with_clip_raw_disc_sites(
        raw, cns, 8, ann, L1_RTYPE, False, step_out)

    with open(driver_out) as a, open(step_out) as b:
        assert a.read() == b.read()
    # and it is not the trivial "everything kept" file: name the survivors, so a change
    # that relaxed a rejection could not pass by moving both sides together.
    assert [l.split()[:2] for l in open(driver_out)] == [
        ["chr1", "5000"], ["chr9", "42"],                       # the cns candidates
        ["chr1", "1000"], ["chr1", "4900"], ["chr3", "300"]]    # created, in raw_disc order


def test_flankless_family_contributes_its_cns_list_only(tmp_path):
    """Alu has no transduction flank, so its step creates nothing and its universe is cns.

    Measured from the other direction by probe_site_creation.py: all 1,230 of BASELINE-0's
    Alu records first exist at cns. The file is written in cns FILE order, and a raw_disc
    that would have contributed rows must not be read at all.
    """
    raw = _write(str(tmp_path / "raw_disc.txt"), [_raw("chr1", 77, 9, 9)])
    cns = _write(str(tmp_path / "cns.txt"),
                 [("chr1", 300, "x"), ("chr1", 100, "x"), ("chr2", 5, "x")])
    out = str(tmp_path / "universe" / "Alu.slct")
    n_cand, n_new = downstream_unified.enumerate_family_sites(
        raw, cns, None, 2, 8, False, out, b_flank=False, xannotation=None)
    assert (n_cand, n_new) == (3, 0)
    assert [l.split() for l in open(out)] == [["chr1", "300"], ["chr1", "100"], ["chr2", "5"]]


# --------------------------------------------------------------------------- provenance


def test_slct_layouts_are_the_provenance(tmp_path):
    """A 2-field line is a cns row the barrier can already see; a longer one is created."""
    p = _slct(str(tmp_path / "L1.slct"),
              cns=[("chr1", 100), ("chr2", 200)],
              new=[_raw("chr1", 900, 9, 9), _raw("chr3", 5, 9, 9)])
    claims = downstream_unified.read_slct_claims(p)
    assert claims == {
        "chr1": {100: "cns", 900: "new"},
        "chr2": {200: "cns"},
        "chr3": {5: "new"},
    }


def test_a_position_claimed_both_ways_records_cns(tmp_path):
    """Only reachable if the cns list itself repeats a position, but the tie must not be
    order-dependent: the cns row is the one a pre-transduction decision point can see."""
    p = _slct(str(tmp_path / "L1.slct"), cns=[("chr1", 100)], new=[_raw("chr1", 100, 9, 9)])
    assert downstream_unified.read_slct_claims(p) == {"chr1": {100: "cns"}}


# --------------------------------------------------------------------------- the union


def _universe(tmp_path, **per_family):
    paths = {}
    for tag, (cns, new) in per_family.items():
        paths[tag] = _slct(str(tmp_path / (tag + ".slct")), cns=cns, new=new)
    return downstream_unified.build_locus_universe(paths)


def test_union_is_exact_position_with_per_family_counts(tmp_path):
    u = _universe(
        tmp_path,
        L1=([("chr1", 100)], [_raw("chr1", 900, 9, 9), _raw("chr1", 950, 9, 9)]),
        SVA=([("chr1", 100)], [_raw("chr1", 900, 9, 9)]),
        Alu=([("chr1", 100), ("chr2", 7)], []),
    )
    assert u.n_loci == 4                      # chr1:100, chr1:900, chr1:950, chr2:7
    assert u.by_n_families() == {1: 2, 2: 1, 3: 1}
    assert u.per_family() == {"L1": 3, "SVA": 2, "Alu": 2}


def test_union_does_not_merge_neighbouring_positions(tmp_path):
    """Exact position, no window. The post-cns chain does not re-anchor coordinates, so two
    families naming one locus name the same integer; folding 1 bp apart together here would
    hide distinct loci inside the count. Cross-family reconciliation is locus_table's job,
    at its own window, after evidence exists."""
    u = _universe(tmp_path, L1=([("chr1", 100)], []), SVA=([("chr1", 101)], []))
    assert u.n_loci == 2
    assert u.by_n_families() == {1: 2}


def test_created_only_loci_are_the_pre_transduction_ceiling(tmp_path):
    """A locus no family has a cns row for cannot be resolved by any decision point placed
    before transduction -- the universe-side form of plan sec 4's 7-of-22."""
    u = _universe(
        tmp_path,
        L1=([("chr1", 100)], [_raw("chr1", 900, 9, 9)]),
        SVA=([], [_raw("chr1", 900, 9, 9), _raw("chr2", 5, 9, 9)]),
    )
    assert u.n_created_only() == 2             # chr1:900 (both created), chr2:5
    assert u.n_loci == 3


def test_missing_family_file_is_skipped_not_an_error(tmp_path):
    u = downstream_unified.build_locus_universe({
        "L1": _slct(str(tmp_path / "L1.slct"), cns=[("chr1", 1)]),
        "SVA": str(tmp_path / "absent.slct"),
    })
    assert u.n_loci == 1 and u.families == ["L1"]


def test_written_table_is_sorted_and_carries_both_columns(tmp_path):
    u = _universe(
        tmp_path,
        L1=([("chr2", 5)], [_raw("chr1", 900, 9, 9)]),
        SVA=([("chr1", 900)], []),
    )
    out = u.write(str(tmp_path / "out" / "locus_universe.tsv"))
    lines = [l.rstrip("\n") for l in open(out)]
    assert lines[0].startswith("#chrom\tpos\tn_families")
    assert lines[1:] == ["chr1\t900\t2\tL1,SVA\tnew,cns", "chr2\t5\t1\tL1\tcns"]


# --------------------------------------------------------------------------- iextnd


def test_transduction_extnd_is_the_step_formula():
    """x_TEA_main.py:473-499: iextnd = 400, raised to int(mean_is + 3*std) when larger.

    The driver's cns_extnd is that second term already, so the floor is the only thing it
    was missing -- and at WGS insert sizes the floor never binds (HG002: 954).
    """
    assert downstream_unified.transduction_extnd(954) == 954
    assert downstream_unified.transduction_extnd(120) == 400
    assert downstream_unified.transduction_extnd(400) == 400
    assert downstream_unified.transduction_extnd(401) == 401


# --------------------------------------------------------------------------- td manifest


def test_manifest_layout_is_what_the_cli_parses(tmp_path):
    """`tag sites extnd clip_fastq disc_fa [disc_pos]`, `-` for an unwanted output.

    The C++ side treats a short line as a HARD error rather than a skip (src/csrc/
    xtea_scan.cpp), because a dropped field means a family silently does not run -- which on
    a WGS gate reads as "that family found nothing". So the field count is the contract.
    """
    out = str(tmp_path / "td_manifest.txt")
    n = downstream_unified.write_td_manifest([
        ("L1", "/x/L1.slct", 954, "/x/L1/clip.fq", "/x/L1/disc.fa"),
        ("SVA", "/x/SVA.slct", 954.0, "/x/SVA/clip.fq", "/x/SVA/disc.fa"),
    ], out)
    assert n == 2
    rows = [l.rstrip("\n").split("\t") for l in open(out)]
    assert all(len(r) == 6 for r in rows)
    assert rows[0] == ["L1", "/x/L1.slct", "954", "/x/L1/clip.fq", "/x/L1/disc.fa", "-"]
    # a float extnd must not reach the CLI as "954.0"
    assert rows[1][2] == "954"


def test_gather_paths_are_the_names_the_step_writes():
    """The driver and the step agree on two file names across a PROCESS boundary.

    `call_candidate_transduction_v3` builds them as `self.working_folder + "<name>"` and the
    driver has to land on exactly those or the step will not find them and will silently
    gather again. There is no shared constant to lean on -- the step hard-codes them -- so
    this pins the literals in the step's own source. If a rename breaks this test, the fix is
    to change BOTH, not to loosen the test.
    """
    src = os.path.join(REPO_ROOT, "xtea", "x_transduction.py")
    with open(src) as fin:
        body = fin.read()
    for name in (downstream_unified.TD_CLIP_FQ_NAME, downstream_unified.TD_DISC_FA_NAME):
        assert 'self.working_folder + "{0}"'.format(name) in body, name
    clip_fq, disc_fa = downstream_unified.td_gather_paths("/w/tmp/transduction/")
    assert clip_fq == "/w/tmp/transduction/" + downstream_unified.TD_CLIP_FQ_NAME
    assert disc_fa == "/w/tmp/transduction/" + downstream_unified.TD_DISC_FA_NAME


def test_the_step_honours_the_precomputed_flag_by_name():
    """The bypass is one constant read on both sides, not two string literals."""
    src = os.path.join(REPO_ROOT, "xtea", "x_transduction.py")
    with open(src) as fin:
        body = fin.read()
    assert "downstream_unified.ENV_TD_PRECOMPUTED" in body
    assert downstream_unified.ENV_TD_PRECOMPUTED == "XTEA_TD_GATHER_PRECOMPUTED"


# --------------------------------------------------------------------------- the waves


def test_locus_universe_opens_the_barrier_split_on_its_own():
    """It needs the same cns_filter -> transduction point --family-assign-early uses, but
    must not imply the assignment, and must not inherit its single-family exemption."""
    import run_xtea
    assert run_xtea.downstream_waves(False, "postF_hc", 3) == [(None, "postF_hc")]
    assert run_xtea.downstream_waves(False, "postF_hc", 3, locus_universe=True) == [
        (None, "cns_filter"), ("transduction", "postF_hc")]
    # one family: no assignment barrier, but still a universe
    assert run_xtea.downstream_waves(True, "postF_hc", 1) == [(None, "postF_hc")]
    assert run_xtea.downstream_waves(False, "postF_hc", 1, locus_universe=True) == [
        (None, "cns_filter"), ("transduction", "postF_hc")]
