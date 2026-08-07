"""Shared fixtures + path resolution for the xTea2 test harness.

Locations are overridable via environment variables so the suite is not hard-wired to one
machine; they default to the frozen NA12878 chr1 golden layout under $HOME. Tests that need a
missing input are skipped (not failed) so the fast unit tests always run.

Defaults are $HOME-relative on purpose — no absolute per-user path belongs in a committed
harness. Point the env vars at your own data (see tests/README.md) if your layout differs.
"""
import os
import warnings

import pytest

# --- resolvable locations --------------------------------------------------------------------
TESTDATA_DIR = os.environ.get("XTEA2_TESTDATA",
                              os.path.expanduser("~/scratch/xtea2_testdata"))
GOLDEN_DIR = os.environ.get("XTEA2_GOLDEN_DIR", os.path.join(TESTDATA_DIR, "master_golden", "NA12878"))
REFERENCE = os.environ.get("XTEA2_REFERENCE", os.path.expanduser("~/hg38_no_alt.fa"))
TEST_BAM = os.environ.get("XTEA2_TEST_BAM", os.path.join(TESTDATA_DIR, "NA12878.chr1.bam"))
# Full-genome CRAM (same reads as the chr1 BAM over chr1); used for BAM↔CRAM equivalence tests.
TEST_CRAM = os.environ.get("XTEA2_TEST_CRAM", os.path.join(TESTDATA_DIR, "NA12878.final.cram"))
# A MEI-rich chr1 window (0-based half-open) for fast read-level checks — spans several golden Alu.
MEI_RICH_REGION = ("chr1", 23_000_000, 30_000_000)
# The xtea-scan binary under test — enables the oracle-vs-scan tests.
#
# DEFAULTS TO THE IN-TREE BUILD, and that is the whole point: for a long time this was
# `os.environ.get("XTEA2_SCAN_BIN", "")`, so a `make -C src/csrc` sitting right there in the
# working tree was invisible to `pytest` and 30 tests — every test of the C++ tool that owns
# `scan_pass`, the largest step in a run — skipped in silence on the default path. Repo-relative,
# so this stays free of per-user absolute paths (see the module docstring).
REPO_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
XTEA2_SCAN_BIN = os.environ.get("XTEA2_SCAN_BIN",
                                os.path.join(REPO_ROOT, "src", "csrc", "xtea_scan"))
# The sources the binary must be newer than. This is the Makefile's `SRC` + `HDR` — and the
# Makefile lists the headers for exactly the reason the staleness check exists: most of the code
# lives in the .hpp files, they are #included rather than compiled separately, so a header-only
# edit leaves `make` seeing nothing to do and the old binary in place.
SCAN_BIN_SOURCES = ("xtea_scan.cpp", "realign.hpp", "cns_emit.hpp", "td_gather.hpp")
# Output of the pipeline under test (set once xtea2 can run); enables the real e2e gate.
XTEA2_OUTPUT_DIR = os.environ.get("XTEA2_OUTPUT_DIR", "")
# Annotation inputs the post-cns downstream steps need (test_downstream_from_cns.py).
REP_LIB_DIR = os.environ.get("XTEA2_REP_LIB_DIR", os.path.expanduser("~/xTea/rep_lib_annotation"))
GENE_GFF3 = os.environ.get("XTEA2_GENE_GFF3",
                           os.path.expanduser("~/data1/corinne/ref/gencode.v45.basic.annotation.gff3"))
BLACKLIST = os.environ.get("XTEA2_BLACKLIST",
                           os.path.join(REP_LIB_DIR, "blacklist", "hg38", "centromere.bed"))

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


def _require(path, what, env_var):
    """Skip (never fail) when an external input is absent, naming the var that would fix it."""
    if not os.path.exists(path):
        pytest.skip(f"{what} not found: {path} (set ${env_var} to your own copy)")
    return path


@pytest.fixture(scope="session")
def golden_dir():
    return _require(GOLDEN_DIR, "golden dir", "XTEA2_GOLDEN_DIR")


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
    return _require(REFERENCE, "reference FASTA", "XTEA2_REFERENCE")


@pytest.fixture(scope="session")
def test_bam():
    return _require(TEST_BAM, "test BAM", "XTEA2_TEST_BAM")


@pytest.fixture(scope="session")
def test_cram():
    return _require(TEST_CRAM, "test CRAM", "XTEA2_TEST_CRAM")


@pytest.fixture(scope="session")
def mei_rich_region():
    return MEI_RICH_REGION


def _scan_bin_stale_against(binary):
    """Sources under src/csrc/ that are NEWER than `binary`, in mtime order. Empty when fresh.

    Returns [] rather than skipping when the source tree is absent — an installed binary with no
    sources beside it cannot be checked, and refusing to run it would be worse than running it.
    """
    csrc = os.path.join(REPO_ROOT, "src", "csrc")
    if not os.path.isdir(csrc):
        return []
    bin_mtime = os.path.getmtime(binary)
    newer = [(os.path.getmtime(p), name)
             for name, p in ((n, os.path.join(csrc, n)) for n in SCAN_BIN_SOURCES)
             if os.path.exists(p) and os.path.getmtime(p) > bin_mtime]
    return [name for _, name in sorted(newer, reverse=True)]


@pytest.fixture(scope="session")
def xtea2_scan_bin():
    """The xtea-scan binary, and a hard FAIL if it is older than its own sources.

    Absent is a SKIP: not built is a legitimate state, and the fast unit tests must still run.
    Stale is a FAILURE, not a skip, because it is the one state that lies — the suite goes green
    while testing code that is not in the tree. This has already happened here: the in-tree binary
    was found five days older than `xtea_scan.cpp`, and several sbatch harnesses default
    XTEA2_SCAN_BIN to a *different commit's* worktree build, so pointing pytest at a binary that
    never matched the sources under test takes only inheriting someone else's default.

    Set XTEA2_ALLOW_STALE_SCAN_BIN=1 to downgrade the failure to a warning — for the deliberate
    case of testing a frozen binary against a moved tree, where the mismatch is the intent.
    """
    if not XTEA2_SCAN_BIN or not os.path.exists(XTEA2_SCAN_BIN):
        pytest.skip(
            "xtea-scan not built at {0} — run `PKG_CONFIG_PATH=<htslib>/lib/pkgconfig "
            "make -C src/csrc`, or set $XTEA2_SCAN_BIN".format(XTEA2_SCAN_BIN))
    stale = _scan_bin_stale_against(XTEA2_SCAN_BIN)
    if stale:
        msg = ("xtea-scan is STALE: {0} is older than {1}. A green run against a stale binary "
               "tests code that is not in this tree. Rebuild with `make -C src/csrc` (the "
               "Makefile's `clean all` if only headers changed), or set "
               "XTEA2_ALLOW_STALE_SCAN_BIN=1 if the mismatch is deliberate."
               .format(XTEA2_SCAN_BIN, ", ".join(stale)))
        if os.environ.get("XTEA2_ALLOW_STALE_SCAN_BIN", "") not in ("", "0"):
            warnings.warn(msg)
        else:
            pytest.fail(msg, pytrace=False)
    return XTEA2_SCAN_BIN


@pytest.fixture(scope="session")
def rep_lib_dir():
    return _require(REP_LIB_DIR, "rep_lib_annotation dir", "XTEA2_REP_LIB_DIR")


@pytest.fixture(scope="session")
def gene_gff3():
    return _require(GENE_GFF3, "gene GFF3", "XTEA2_GENE_GFF3")


@pytest.fixture(scope="session")
def blacklist():
    return _require(BLACKLIST, "blacklist BED", "XTEA2_BLACKLIST")


@pytest.fixture(scope="session")
def xtea2_output_dir():
    if not XTEA2_OUTPUT_DIR or not os.path.isdir(XTEA2_OUTPUT_DIR):
        pytest.skip("XTEA2_OUTPUT_DIR not set — pipeline-under-test output unavailable "
                    "(expected until xtea2 can run end-to-end)")
    return XTEA2_OUTPUT_DIR
