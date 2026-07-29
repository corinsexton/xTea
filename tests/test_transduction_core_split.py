"""The transduction step gets its own core split, because Alu has no transduction.

WHY THIS GATE EXISTS
--------------------
`--parallel-downstream` gives each MEI type an even `cores // ntypes` slice. That is right for
every step except transduction: `REP_TYPES["Alu"]["flank"] is None`, so Alu's `--transduction`
takes the empty branch in x_TEA_main and finishes in ~1.5s, while L1's and SVA's take ~25 and ~30
minutes. Measured on the NA12878 WGS cpp-cns run (49070030, step_timings.tsv): Alu transduction
1.5s vs L1 1481s vs SVA 1772s, and transduction is ~93% of the whole downstream phase. An even
three-way split therefore parks a third of the machine on a no-op for ~25 minutes.

The fix is static -- which types have a flank is known from REP_TYPES before anything runs -- so
what this gate pins is the *derivation*, not a timing. Two things must hold, and both are easy to
break with an innocent-looking edit to build_type_steps:
  1. the transduction step's `-n` uses the flanked-types-only divisor, and
  2. every OTHER step keeps the even all-types divisor (widening those too would oversubscribe).
"""
import os
import sys
from types import SimpleNamespace

REPO_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(REPO_ROOT, "xtea"))
import run_xtea  # noqa: E402  (needs the path insert above)


def _opt(cores=24):
    return SimpleNamespace(
        cores=cores, nfclip=4, nfdisc=5, flklen=3000, sample="s",
        purity=None, tumor=False, case_control=False, user_specified=None,
    )


def _paths(t):
    spec = run_xtea.REP_TYPES[t]
    return {"rtype": spec["rtype"], "flank": spec["flank"] or "null",
            "cns": "cns.fa", "annotation": "a.out", "annotation1": "a1.out",
            "copies_with_flank": "cwf.fa"}


def _dirs():
    return {"type": "/tmp/x", "bam_list": "/tmp/bams.txt",
            "tmp_cns": "/tmp/x/tmp/cns/", "tmp_td": "/tmp/x/tmp/td/",
            "tmp": "/tmp/x/tmp/"}


def _n_of(argv):
    """The value of this step's -n flag."""
    return argv[argv.index("-n") + 1]


def _steps(t, n_cores, n_cores_transduction):
    # build_type_steps assembles the whole chain and only then applies stop_after, so the
    # cfg keys the later steps read must be present even though they get truncated away.
    cfg = {"reference": "ref.fa", "blacklist": "bl.txt", "gene_gff3": "genes.gff3",
           "model": "gt.onnx", "input": "in.txt"}
    return dict(run_xtea.build_type_steps(
        t, cfg, _paths(t), _dirs(), _opt(),
        n_cores=n_cores, stop_after="postF_hc",
        n_cores_transduction=n_cores_transduction))


def test_alu_has_no_transduction_flank():
    """The premise the whole split rests on. If this ever changes, the split must be revisited."""
    assert run_xtea.REP_TYPES["Alu"]["flank"] is None
    assert run_xtea.REP_TYPES["L1"]["flank"]
    assert run_xtea.REP_TYPES["SVA"]["flank"]


def test_flanked_types_get_the_wider_transduction_split():
    """24 cores / 3 types = 8 per type, but transduction splits over the 2 flanked types = 12."""
    for t in ("L1", "SVA"):
        steps = _steps(t, n_cores=8, n_cores_transduction=12)
        assert _n_of(steps["transduction"]) == "12", t


def test_every_other_step_keeps_the_even_split():
    """Only transduction is widened -- widening the rest would oversubscribe the allocation."""
    for t in ("L1", "SVA"):
        steps = _steps(t, n_cores=8, n_cores_transduction=12)
        for name, argv in steps.items():
            if name == "transduction":
                continue
            assert _n_of(argv) == "8", f"{t}:{name} should keep the even split"


def test_flankless_type_gets_one_worker_for_its_noop():
    """Alu's transduction is a no-op; it does not need -- and must not hold -- a third of the box."""
    steps = _steps("Alu", n_cores=8, n_cores_transduction=1)
    assert _n_of(steps["transduction"]) == "1"
    assert _n_of(steps["cns_filter"]) == "8"   # its real work keeps the even slice


def test_omitting_the_override_is_unchanged_behaviour():
    """Back-compat: no n_cores_transduction means every step uses n_cores, as before."""
    steps = _steps("L1", n_cores=8, n_cores_transduction=None)
    assert all(_n_of(argv) == "8" for argv in steps.values())


def test_the_divisor_the_caller_derives():
    """Reproduce the caller's arithmetic: flanked types only, floor-divided, min 1."""
    type_ctx = [("L1",), ("Alu",), ("SVA",)]
    td_types = [c for c in type_ctx if run_xtea.REP_TYPES[c[0]].get("flank")]
    assert [c[0] for c in td_types] == ["L1", "SVA"]
    assert max(1, 24 // len(td_types)) == 12
    # An Alu-only run has no flanked type at all; the divisor must not blow up.
    alu_only = [c for c in [("Alu",)] if run_xtea.REP_TYPES[c[0]].get("flank")]
    assert alu_only == []
