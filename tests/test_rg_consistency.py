"""@RG SM-tag consistency guard for multi-file input.

The invariant: pooling several files under one `-i`/`--normal-bam` is the intended use of
multi-file input, but pooling files that turn out to be *different samples* by accident is a
real bug class — so differing SM tags must fail loud by default, and `--ignore-rg` is the
explicit opt-out (SMaHT-style data where centers legitimately use different SM tags for the
same donor/tissue).

Two layers, because they can fail independently:
  * `check_rg_consistency` branch matrix — outcome *and* message. A guard that exits without
    naming which file carried which tag has failed at its job, so the text is asserted too.
  * the real `run_xtea.py` CLI — that the flag is wired to *both* call sites (`-i` and
    `--normal-bam`), and that the guard fires before any pipeline step, so a mismatched pair
    costs nothing.

Everything here is hermetic: the fixtures are ~250-byte pysam-built BAMs differing only in
their @RG block, and the CLI runs under `--dry-run` with placeholder reference/rep-lib paths
(nothing on the dry-run path dereferences them). No golden data, no reference FASTA, no binary.
"""
import os
import subprocess
import sys

import pytest

pysam = pytest.importorskip("pysam")

# run_xtea.py uses flat imports (scan_adapter, disc_adapter, ...), so the package dir itself
# has to be importable — same reason run_step sets cwd to x_TEA_main.py's dir.
REPO_ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
sys.path.insert(0, os.path.join(REPO_ROOT, "xtea"))
try:
    # Importing run_xtea pulls the whole pipeline import graph (disc_adapter → x_annotation →
    # interval_index). Skip rather than error out during collection if that graph is
    # unsatisfiable — the harness's contract is that the fast units always run.
    from run_xtea import check_rg_consistency  # noqa: E402
except ImportError as e:  # pragma: no cover - environment-dependent
    pytest.skip(f"run_xtea import graph unsatisfiable ({e}); "
                f"install the pipeline deps to run the @RG tests",
                allow_module_level=True)

RUN_XTEA = os.path.join(REPO_ROOT, "xtea", "run_xtea.py")

# name -> @RG list (None = omit the @RG block entirely)
_RG_CASES = {
    "sampA_lib1": [{"ID": "rg1", "SM": "SAMP_A", "PL": "ILLUMINA", "LB": "lib1"}],
    "sampA_lib2": [{"ID": "rg2", "SM": "SAMP_A", "PL": "ILLUMINA", "LB": "lib2"}],
    "sampB_lib1": [{"ID": "rg3", "SM": "SAMP_B", "PL": "ILLUMINA", "LB": "lib1"}],
    "no_rg": None,
    "rg_no_sm": [{"ID": "rg4", "PL": "ILLUMINA", "LB": "lib1"}],
    "two_sm": [{"ID": "rg5", "SM": "SAMP_A", "PL": "ILLUMINA", "LB": "lib1"},
               {"ID": "rg6", "SM": "SAMP_B", "PL": "ILLUMINA", "LB": "lib2"}],
}


@pytest.fixture(scope="module")
def rg_bams(tmp_path_factory):
    """{name: path} for tiny indexed BAMs that differ only in their @RG block."""
    d = tmp_path_factory.mktemp("rg_bams")
    out = {}
    for name, rgs in _RG_CASES.items():
        path = str(d / f"{name}.bam")
        header = {"HD": {"VN": "1.6", "SO": "coordinate"},
                  "SQ": [{"SN": "chr1", "LN": 248956422}]}
        if rgs is not None:
            header["RG"] = rgs
        with pysam.AlignmentFile(path, "wb", header=header) as bf:
            for i in range(2):
                a = pysam.AlignedSegment(bf.header)
                a.query_name = f"read{i}"
                a.query_sequence = "A" * 50
                a.query_qualities = pysam.qualitystring_to_array("I" * 50)
                a.flag = 0
                a.reference_id = 0
                a.reference_start = 1_000_000 + i * 100
                a.mapping_quality = 60
                a.cigarstring = "50M"
                if rgs is not None:
                    a.set_tag("RG", rgs[0]["ID"])
                bf.write(a)
        pysam.index(path)
        out[name] = path
    return out


def _run_check(paths, ignore_rg, capsys):
    """Call check_rg_consistency, returning (exited, combined_output)."""
    exited = False
    try:
        check_rg_consistency(paths, None, ignore_rg)
    except SystemExit as e:
        exited = True
        print(e.code)
    return exited, capsys.readouterr().out


# --- branch matrix ---------------------------------------------------------------------------
# (label, fixture names, ignore_rg, expect_exit, substrings the output must contain)
MATRIX = [
    ("single file with an SM tag proceeds",
     ["sampA_lib1"], False, False, ["@RG check: single SM tag", "SM:SAMP_A"]),
    ("two files sharing one SM tag proceed (the legitimate pooling case)",
     ["sampA_lib1", "sampA_lib2"], False, False, ["single SM tag", "SM:SAMP_A"]),
    ("two files with different SM tags fail loud, naming both files and tags",
     ["sampA_lib1", "sampB_lib1"], False, True,
     ["2 distinct SM tags", "SM:SAMP_A", "SM:SAMP_B",
      "sampA_lib1.bam", "sampB_lib1.bam", "--ignore-rg"]),
    ("--ignore-rg proceeds past differing SM tags, still reporting what it saw",
     ["sampA_lib1", "sampB_lib1"], True, False,
     ["--ignore-rg: proceeding", "SAMP_A", "SAMP_B"]),
    ("a missing @RG block fails loud",
     ["no_rg"], False, True, ["no @RG header", "--ignore-rg"]),
    ("--ignore-rg downgrades a missing @RG block to a warning",
     ["no_rg"], True, False, ["WARNING (--ignore-rg): no @RG header", "proceeding"]),
    ("@RG present but no SM tag fails loud",
     ["rg_no_sm"], False, True, ["@RG present but no SM tag", "--ignore-rg"]),
    ("--ignore-rg downgrades a missing SM tag to a warning",
     ["rg_no_sm"], True, False,
     ["WARNING (--ignore-rg): @RG present but no SM tag", "proceeding"]),
    ("one file carrying two SM tags fails loud",
     ["two_sm"], False, True,
     ["itself has multiple SM tags", "SAMP_A", "SAMP_B", "--ignore-rg"]),
    ("--ignore-rg downgrades an internally-mixed file to a warning",
     ["two_sm"], True, False,
     ["WARNING (--ignore-rg): file", "multiple SM tags", "proceeding"]),
    ("a bad file alongside a good one still fails",
     ["no_rg", "sampA_lib1"], False, True, ["no @RG header"]),
]


@pytest.mark.parametrize("label,names,ignore_rg,expect_exit,needles", MATRIX,
                         ids=[c[0] for c in MATRIX])
def test_rg_check_branch(rg_bams, capsys, label, names, ignore_rg, expect_exit, needles):
    paths = [rg_bams[n] for n in names]
    exited, out = _run_check(paths, ignore_rg, capsys)
    assert exited == expect_exit, (
        f"{label}: expected {'exit' if expect_exit else 'proceed'}, "
        f"got {'exit' if exited else 'proceed'}\noutput:\n{out}"
    )
    missing = [n for n in needles if n not in out]
    assert not missing, f"{label}: message missing {missing}\noutput:\n{out}"


# --- the real CLI ----------------------------------------------------------------------------
def _run_cli(tmp_path, extra_args):
    """Run run_xtea.py --dry-run with placeholder inputs; returns CompletedProcess.

    --dry-run prints the per-step commands and exits without executing them, and nothing on
    that path opens the reference / rep-lib / gff3 — so placeholders keep this hermetic.
    """
    fake = tmp_path / "placeholder"
    (fake / "replib").mkdir(parents=True, exist_ok=True)
    ref, gff3 = fake / "ref.fa", fake / "genes.gff3"
    ref.touch()
    gff3.touch()
    argv = [sys.executable, RUN_XTEA,
            "-s", "rgtest", "-o", str(tmp_path / "out"),
            "--reference", str(ref), "--rep-lib-dir", str(fake / "replib"),
            "--gene-gff3", str(gff3), "--blacklist", "null",
            "--types", "Alu", "--dry-run"] + extra_args
    return subprocess.run(argv, capture_output=True, text=True)


def test_cli_input_mismatch_fails_before_any_step(rg_bams, tmp_path):
    """Default: two SM tags in -i is a hard error, raised before any pipeline step runs."""
    r = _run_cli(tmp_path, ["-i", f"{rg_bams['sampA_lib1']} {rg_bams['sampB_lib1']}"])
    both = r.stdout + r.stderr
    assert r.returncode != 0, f"expected nonzero exit\n{both}"
    assert "2 distinct SM tags across input" in both, both
    assert "sampA_lib1.bam" in both and "sampB_lib1.bam" in both, \
        f"error must say which file carried which tag\n{both}"
    assert "pass --ignore-rg" in both, f"error must name the opt-out\n{both}"
    # fail-fast is the point: a mismatched pair must not cost a single pipeline step
    assert "=== STEP" not in both, f"guard fired too late — a step already ran\n{both}"


def test_cli_ignore_rg_proceeds_on_input(rg_bams, tmp_path):
    """--ignore-rg lets the same input through, and does so audibly."""
    r = _run_cli(tmp_path, ["-i", f"{rg_bams['sampA_lib1']} {rg_bams['sampB_lib1']}",
                            "--ignore-rg"])
    both = r.stdout + r.stderr
    assert "--ignore-rg: proceeding" in both, both
    assert "SAMP_A" in both and "SAMP_B" in both, \
        f"opting out must not be silent about what was pooled\n{both}"
    assert "2 distinct SM tags" not in both, both
    assert r.returncode == 0, f"dry run should complete\n{both}"


def test_cli_control_bam_is_checked_too(rg_bams, tmp_path):
    """The --normal-bam call site is guarded as well, not just -i."""
    r = _run_cli(tmp_path, ["-i", rg_bams["sampA_lib1"], "--case-control",
                            "--normal-bam", f"{rg_bams['sampA_lib1']} {rg_bams['sampB_lib1']}",
                            "--tumor", "--purity", "0.86"])
    both = r.stdout + r.stderr
    assert r.returncode != 0, f"expected nonzero exit\n{both}"
    assert "2 distinct SM tags across input" in both, both
    # the -i check logging success first is what proves we reached the *second* call site
    assert "single SM tag across all input files" in both, \
        f"failure should come from the control check, not from -i\n{both}"
    assert "=== STEP" not in both, f"guard fired too late — a step already ran\n{both}"


def test_cli_ignore_rg_proceeds_on_control_bam(rg_bams, tmp_path):
    r = _run_cli(tmp_path, ["-i", rg_bams["sampA_lib1"], "--case-control",
                            "--normal-bam", f"{rg_bams['sampA_lib1']} {rg_bams['sampB_lib1']}",
                            "--tumor", "--purity", "0.86", "--ignore-rg"])
    both = r.stdout + r.stderr
    assert "--ignore-rg: proceeding" in both, both
    assert "2 distinct SM tags" not in both, both
    assert r.returncode == 0, f"dry run should complete\n{both}"
