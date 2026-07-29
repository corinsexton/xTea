"""Retention policy for completed xTea2 runs."""
import os
import sys

import pytest

REPO_ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
sys.path.insert(0, os.path.join(REPO_ROOT, "xtea"))
try:
    from run_xtea import cleanup_successful_run  # noqa: E402
except ImportError as e:  # pragma: no cover - environment-dependent import graph
    pytest.skip(f"run_xtea import graph unsatisfiable ({e})", allow_module_level=True)


def test_cleanup_retains_only_final_merged_vcf(tmp_path):
    root = tmp_path / "sample"
    nested = root / "L1" / "tmp"
    nested.mkdir(parents=True)
    final_vcf = root / "sample.vcf"
    final_vcf.write_text("##fileformat=VCFv4.2\n")
    (nested / "temporary.sam").write_text("intermediate\n")
    (root / "sample_xtea2.log").write_text("log\n")
    (root / "bam_list.txt").write_text("input.bam\tillumina\n")

    cleanup_successful_run(root, "sample")

    assert final_vcf.read_text() == "##fileformat=VCFv4.2\n"
    assert list(root.iterdir()) == [final_vcf]


def test_cleanup_refuses_when_final_vcf_is_missing(tmp_path):
    root = tmp_path / "sample"
    root.mkdir()
    intermediate = root / "candidate_list.txt"
    intermediate.write_text("diagnostics\n")

    with pytest.raises(RuntimeError, match="final merged VCF is missing"):
        cleanup_successful_run(root, "sample")

    assert intermediate.exists()
