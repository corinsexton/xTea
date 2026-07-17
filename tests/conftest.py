"""Shared fixtures + path resolution for the xTea2 test harness.

Locations are overridable via environment variables so the suite is not hard-wired to one
machine; they default to the frozen NA12878 chr1 golden on O2. Tests that need a missing input
are skipped (not failed) so the fast unit tests always run.
"""
import os
import pytest

# --- resolvable locations --------------------------------------------------------------------
TESTDATA_DIR = os.environ.get("XTEA2_TESTDATA", "/home/cos689/scratch/xtea2_testdata")
GOLDEN_DIR = os.environ.get("XTEA2_GOLDEN_DIR", os.path.join(TESTDATA_DIR, "master_golden", "NA12878"))
REFERENCE = os.environ.get("XTEA2_REFERENCE", "/home/cos689/hg38_no_alt.fa")
TEST_BAM = os.environ.get("XTEA2_TEST_BAM", os.path.join(TESTDATA_DIR, "NA12878.chr1.bam"))
# Full-genome CRAM (same reads as the chr1 BAM over chr1); used for BAM↔CRAM equivalence tests.
TEST_CRAM = os.environ.get("XTEA2_TEST_CRAM", os.path.join(TESTDATA_DIR, "NA12878.final.cram"))
# A MEI-rich chr1 window (0-based half-open) for fast read-level checks — spans several golden Alu.
MEI_RICH_REGION = ("chr1", 23_000_000, 30_000_000)
# The xtea-scan binary under test (set once it is built) — enables oracle-vs-scan tests.
XTEA2_SCAN_BIN = os.environ.get("XTEA2_SCAN_BIN", "")
# Output of the pipeline under test (set once xtea2 can run); enables the real e2e gate.
XTEA2_OUTPUT_DIR = os.environ.get("XTEA2_OUTPUT_DIR", "")

# Golden VCF per MEI type (relative to GOLDEN_DIR).
GOLDEN_VCF_RELPATHS = {
    "ALU": os.path.join("Alu", "NA12878.chr1_ALU.vcf"),
    "LINE1": os.path.join("L1", "NA12878.chr1_LINE1.vcf"),
    "SVA": os.path.join("SVA", "NA12878.chr1_SVA.vcf"),
}
# Golden clip/disc site tables per type.
GOLDEN_CLIP_RELPATHS = {
    "ALU": os.path.join("Alu", "candidate_list_from_clip.txt"),
    "LINE1": os.path.join("L1", "candidate_list_from_clip.txt"),
    "SVA": os.path.join("SVA", "candidate_list_from_clip.txt"),
}
GOLDEN_DISC_RELPATHS = {
    "ALU": os.path.join("Alu", "candidate_list_from_disc.txt"),
    "LINE1": os.path.join("L1", "candidate_list_from_disc.txt"),
    "SVA": os.path.join("SVA", "candidate_list_from_disc.txt"),
}


def _require(path, what):
    if not os.path.exists(path):
        pytest.skip(f"{what} not found: {path} (set the matching XTEA2_* env var)")
    return path


@pytest.fixture(scope="session")
def golden_dir():
    return _require(GOLDEN_DIR, "golden dir")


@pytest.fixture(scope="session")
def golden_vcfs(golden_dir):
    """{type: path} for the golden VCFs — glob per-type subdir so it works for the chr1 OR WGS golden.

    Point XTEA2_GOLDEN_DIR at master_golden_wgs/NA12878 to use the whole-genome golden (the primary
    baseline — chr1-only cripples L1; see docs/clipdisc_contract_spec.md V4).
    """
    import glob
    patterns = {"ALU": ("Alu", "*_ALU.vcf"), "LINE1": ("L1", "*_LINE1.vcf"), "SVA": ("SVA", "*_SVA.vcf")}
    out = {}
    for t, (sub, pat) in patterns.items():
        hits = sorted(glob.glob(os.path.join(golden_dir, sub, pat))) or \
               sorted(glob.glob(os.path.join(golden_dir, os.path.basename(GOLDEN_VCF_RELPATHS[t]))))
        if hits:
            out[t] = hits[0]
    if not out:
        pytest.skip(f"no golden VCFs under {golden_dir}")
    return out


@pytest.fixture(scope="session")
def golden_clip_tables(golden_dir):
    return {t: os.path.join(golden_dir, rel) for t, rel in GOLDEN_CLIP_RELPATHS.items()
            if os.path.exists(os.path.join(golden_dir, rel))}


@pytest.fixture(scope="session")
def golden_disc_tables(golden_dir):
    return {t: os.path.join(golden_dir, rel) for t, rel in GOLDEN_DISC_RELPATHS.items()
            if os.path.exists(os.path.join(golden_dir, rel))}


@pytest.fixture(scope="session")
def reference():
    return _require(REFERENCE, "reference FASTA")


@pytest.fixture(scope="session")
def test_bam():
    return _require(TEST_BAM, "test BAM")


@pytest.fixture(scope="session")
def test_cram():
    return _require(TEST_CRAM, "test CRAM")


@pytest.fixture(scope="session")
def mei_rich_region():
    return MEI_RICH_REGION


@pytest.fixture(scope="session")
def xtea2_scan_bin():
    if not XTEA2_SCAN_BIN or not os.path.exists(XTEA2_SCAN_BIN):
        pytest.skip("XTEA2_SCAN_BIN not set — the C xtea-scan tool is not built yet")
    return XTEA2_SCAN_BIN


@pytest.fixture(scope="session")
def xtea2_output_dir():
    if not XTEA2_OUTPUT_DIR or not os.path.isdir(XTEA2_OUTPUT_DIR):
        pytest.skip("XTEA2_OUTPUT_DIR not set — pipeline-under-test output unavailable "
                    "(expected until xtea2 can run end-to-end)")
    return XTEA2_OUTPUT_DIR
