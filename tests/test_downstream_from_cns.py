"""Fast downstream regression gate — run the post-cns chain against pre-computed cns outputs.

**Why this exists.** The full chr1 concordance gate takes ~18 minutes because it re-runs the whole
funnel (scan -> clip -> disc -> cns) before it reaches anything downstream. But the downstream steps
(transduction -> sibling -> postF x2 -> gene -> genotype -> gVCF) only ever read the *files* the cns
step wrote. So we stage the golden's own cns outputs as the input and run just the tail, comparing
byte-for-byte against the golden's own downstream outputs.

**Measured cost: ~11 min for all three types** (2026-07-25, two runs: 11m48s and 11m02s). That is a modest saving on the 18-minute
full gate — the point is not raw speed but *isolation*: a downstream edit is gated without re-running
scan/cns, and a failure names the exact step and file that moved. Gate one type while iterating
(`pytest -k Alu tests/test_downstream_from_cns.py`) and the loop is a few minutes.

**The interval tree has since been fixed (2026-07-25, trim step 2) and this gate now runs in
8m25s**, down from ~11 min. That is less than the "~2-3 min" the plan predicted, and the reason is
worth keeping: the AVL tree was 96% of *postF*, not 96% of the whole downstream chain, so the
saving was always bounded by postF's share of it. Annotation construction is now ~100x cheaper
(23.6s -> 0.22s on Alu) and cached, so what remains is genuinely other work. Profile before
optimizing it further — that is how the interval tree was found in the first place.

**⚠ The goldens predate the B8 fix.** Step 2 deliberately changed which repeat copy the divergence
filter consults, always in the more-conservative direction. On the chr1 golden no final call moved
(this gate passes), but four L1 transduction intermediates did, and those diffs are expected rather
than regressions until the goldens are re-baselined. See docs/handoff/TRIM_PLAN.md
"Step 2's carry-forward".

**Contract: the gate is VCF columns 1-8** (user decision, 2026-07-25). The call set is what must
not move; intermediate files are allowed to churn during the trim without failing the build.

The byte-level comparison of every intermediate still runs and prints what moved — it is
diagnostics, not a verdict. That is deliberate: B7 (the SVA transduction divergence) was found only
because intermediates were compared, and a trim that silently changes an intermediate is exactly
what you want to see in the log before it becomes a call-set bug. Two classes are recognised and
labelled rather than reported as churn: the genotype column (B6, DF21 golden vs ONNX) and the
enumerated known divergences below.

**What this does NOT cover.** Nothing before the cns step. The scan/clip/disc/cns adapters are
gated by the full e2e concordance run (test_e2e_concordance.py) — this harness deliberately trusts
its staged inputs so it can be fast.

The step argv comes from `run_xtea.build_type_steps(..., start_from="transduction")` — the same
function the real pipeline uses — so this harness cannot drift from the real chain.
"""
import os
import shutil
import subprocess
import sys
from types import SimpleNamespace

import pytest

from tests.conftest import _require

REPO_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(REPO_ROOT, "xtea"))
import run_xtea  # noqa: E402  (needs the path insert above)

# The genotype column in the post-gene genotype table. master's golden used DF21; xTea2 uses ONNX.
GENOTYPE_COL = 50
# Files carrying that genotype value, and therefore allowed to differ in it and nothing else.
GENOTYPE_BEARING_SUFFIXES = ("post_filtering_with_gene_gntp.txt",
                             "post_filtering_with_gene_gntp.txt.sorted")

# Known intermediate divergences from the golden — explicitly enumerated so that everything NOT
# listed here stays a byte-exact gate. Each entry must say why it is tolerated; an entry whose
# downstream effect is not "none" does not belong here.
#
# SVA transduction: our x_TEA_main.py emits one extra site (chr1:182305180) in cns2.txt, and its
# high_confident sibling gains that site while losing chr1:50381610 and chr1:182510697. Verified
# deterministic (two back-to-back runs byte-identical to each other) and verified to WASH OUT
# completely: candidate_disc_filtered_cns_post_filtering.txt, the high_confident post_filtering
# table and the final SVA VCF are all byte-identical to the golden, 5 calls either way. Root cause
# not yet attributed — candidates are the xTea2 edits to x_clip_disc_filter.py and the py3.8->3.10
# interpreter change. Tracked as B7 in docs/clipdisc_contract_spec.md.
KNOWN_INTERMEDIATE_DIVERGENCES = {
    ("SVA", "candidate_disc_filtered_cns2.txt"): "B7 — washes out by postF; VCF identical",
    ("SVA", "candidate_disc_filtered_cns2.txt.high_confident"): "B7 — washes out by postF; VCF identical",
}

# Everything the cns step (`x_TEA_main.py -N`) writes, identified by mtime in the golden tree.
# These are the harness's *inputs*; staging exactly this set is what makes the tail runnable
# without re-running the funnel. A file absent for a given type is skipped (L1/SVA emit fewer).
CNS_PRODUCTS = [
    "candidate_disc_filtered_cns.txt",
    "candidate_disc_filtered_cns.txt.high_confident",
    "candidate_disc_filtered_cns.txt.before_filtering",
    "candidate_disc_filtered_cns.txt.before_calling_transduction",
    "candidate_disc_filtered_cns.txt.before_calling_transduction.sites_cov",
    "candidate_disc_filtered_cns.txt.gntp.features",
    "candidate_disc_filtered_cns.txt.gntp.features0.out",
]
# transduction/sibling read the raw disc stream the *disc* step wrote.
EXTRA_INPUTS = ["candidate_list_from_disc.txt.clip_sites_raw_disc.txt"]
# Cached coverage/insert-size estimate. Staging it pins the value instead of re-deriving it from
# the BAM per run — which is both slow and the suspected source of B3's nondeterminism.
COV_CACHE = "basic_cov_is_rlth_info.txt"

TYPES = ["L1", "Alu", "SVA"]


def _stage(golden_type_dir, dest, test_bam):
    """Copy the golden's post-cns state into `dest`, and point bam_list at the test BAM."""
    os.makedirs(os.path.join(dest, "tmp", "cns"), exist_ok=True)
    os.makedirs(os.path.join(dest, "tmp", "transduction"), exist_ok=True)
    staged = 0
    for name in CNS_PRODUCTS + EXTRA_INPUTS:
        src = os.path.join(golden_type_dir, name)
        if os.path.exists(src):
            shutil.copy2(src, os.path.join(dest, name))
            staged += 1
    if not staged:
        pytest.skip(f"no cns outputs staged from {golden_type_dir}")
    for sub in ("cns", "transduction"):
        src = os.path.join(golden_type_dir, "tmp", sub, COV_CACHE)
        if os.path.exists(src):
            shutil.copy2(src, os.path.join(dest, "tmp", sub, COV_CACHE))
    # Rewrite rather than copy: the golden's bam_list holds an absolute path from the machine that
    # produced it, and the whole point of the env-var indirection is not to depend on that.
    with open(os.path.join(dest, "bam_list.txt"), "w") as fh:
        fh.write(f"{test_bam}\tillumina\n")


def _build_steps(t, dest, rep_lib_dir, reference, gene_gff3, blacklist_path, cores):
    """Real argv, from the real builder — see module docstring."""
    spec = run_xtea.REP_TYPES[t]
    cfg = {
        "reference": reference,
        "blacklist": blacklist_path,
        "gene_gff3": gene_gff3,
        "model": os.path.join(REPO_ROOT, "xtea", "genotyping", "genotype_model_6_7_2024.onnx"),
    }
    paths = {
        "rtype": spec["rtype"],
        "copies_with_flank": run_xtea.resolve_annot(rep_lib_dir, "hg38", spec["copies_with_flank"]),
        "annotation": run_xtea.resolve_annot(rep_lib_dir, "hg38", spec["annotation"]),
        "annotation1": run_xtea.resolve_annot(rep_lib_dir, "hg38", spec["annotation1"]),
        "flank": run_xtea.resolve_annot(rep_lib_dir, "hg38", spec["flank"]),
        "cns": run_xtea.resolve_annot(rep_lib_dir, "hg38", spec["cns"]),
    }
    dirs = {
        "type": dest,
        "tmp": os.path.join(dest, "tmp") + "/",
        "tmp_cns": os.path.join(dest, "tmp", "cns") + "/",
        "tmp_td": os.path.join(dest, "tmp", "transduction") + "/",
        "pub_clip": os.path.join(dest, "tmp", "clip") + "/",
        "bam_list": os.path.join(dest, "bam_list.txt"),
        "control_bam_list": None,
    }
    # Defaults pinned to the golden run (run_xtea.py's own argparse defaults).
    opt = SimpleNamespace(cores=cores, nclip=4, cr=2, nd=5, nfclip=4, nfdisc=5, flklen=3000,
                          tumor=False, purity=1.0, case_control=False, sample="NA12878.chr1")
    return run_xtea.build_type_steps(t, cfg, paths, dirs, opt,
                                     start_from="transduction")


def _run_steps(steps, dest):
    xtea_main = os.path.join(REPO_ROOT, "xtea", "x_TEA_main.py")
    for name, argv in steps:
        proc = subprocess.run([sys.executable, xtea_main] + argv,
                              cwd=dest, capture_output=True, text=True)
        if proc.returncode != 0:
            pytest.fail(f"downstream step '{name}' failed (exit {proc.returncode})\n"
                        f"--- stderr ---\n{proc.stderr[-3000:]}")


def _compare(mei_type, golden_type_dir, dest):
    """Compare every file we produced that the golden also has.

    Self-maintaining on purpose: a new downstream output starts being checked the moment both
    sides have it, with no list here to keep in sync.

    Returns diagnostics, NOT verdicts — the gate is VCF cols 1-8 (see the test docstring). The
    fourth return value is "intermediates that moved", which is printed, not asserted on.
    """
    identical, genotype_diffs, tolerated, byte_diffs = [], {}, [], []
    for name in sorted(os.listdir(dest)):
        ours = os.path.join(dest, name)
        theirs = os.path.join(golden_type_dir, name)
        if not os.path.isfile(ours) or not os.path.isfile(theirs):
            continue
        if name in CNS_PRODUCTS or name in EXTRA_INPUTS or name == "bam_list.txt":
            continue  # staged inputs, not results
        with open(ours, "rb") as a, open(theirs, "rb") as b:
            if a.read() == b.read():
                identical.append(name)
                continue
        if name.endswith(GENOTYPE_BEARING_SUFFIXES):
            genotype_diffs[name] = _genotype_only_diff(theirs, ours, byte_diffs, name)
        elif name.endswith(".vcf"):
            _vcf_genotype_only_diff(theirs, ours, byte_diffs, name)
        elif (mei_type, name) in KNOWN_INTERMEDIATE_DIVERGENCES:
            tolerated.append(f"{name}: {KNOWN_INTERMEDIATE_DIVERGENCES[(mei_type, name)]}")
        else:
            byte_diffs.append(f"{name}: differs from golden (expected byte-identical)\n"
                            f"    {_first_diff(theirs, ours)}")
    return identical, genotype_diffs, tolerated, byte_diffs


def _vcf_genotype_only_diff(golden_path, ours_path, byte_diffs, name):
    """A VCF may differ only in the sample column's GT field (col 10) — cols 1-9 are the call set."""
    g = [l.rstrip("\n").split("\t") for l in open(golden_path) if not l.startswith("#")]
    w = [l.rstrip("\n").split("\t") for l in open(ours_path) if not l.startswith("#")]
    if len(g) != len(w):
        byte_diffs.append(f"{name}: call count changed ({len(g)} golden vs {len(w)} ours)")
        return
    for a, b in zip(g, w):
        if a[:9] != b[:9]:
            byte_diffs.append(f"{name}: {a[0]}:{a[1]} differs in VCF cols 1-9 (the call set)\n"
                            f"      golden: {a[:9]}\n      ours:   {b[:9]}")


def _genotype_only_diff(golden_path, ours_path, byte_diffs, name):
    """The ONE tolerated difference: same rows/positions, only the genotype column may move."""
    g = [l.rstrip("\n").split("\t") for l in open(golden_path)]
    w = [l.rstrip("\n").split("\t") for l in open(ours_path)]
    if len(g) != len(w):
        byte_diffs.append(f"{name}: row count changed ({len(g)} golden vs {len(w)} ours) — "
                        f"the genotype step must not add or drop calls")
        return []
    moved = []
    for a, b in zip(g, w):
        bad = [j for j, (x, y) in enumerate(zip(a, b)) if x != y and j != GENOTYPE_COL]
        if bad:
            byte_diffs.append(f"{name}: {a[0]}:{a[1]} differs outside the genotype column "
                            f"(cols {bad}) — only col{GENOTYPE_COL} may move")
        if len(a) > GENOTYPE_COL and a[GENOTYPE_COL] != b[GENOTYPE_COL]:
            moved.append((a[0], a[1], a[GENOTYPE_COL], b[GENOTYPE_COL]))
    return moved


def _first_diff(golden_path, ours_path):
    g = open(golden_path, errors="replace").read().splitlines()
    w = open(ours_path, errors="replace").read().splitlines()
    if len(g) != len(w):
        return f"line count {len(g)} golden vs {len(w)} ours"
    for i, (a, b) in enumerate(zip(g, w), 1):
        if a != b:
            return f"first diff at line {i}:\n      golden: {a[:200]}\n      ours:   {b[:200]}"
    return "content differs only in trailing bytes"


@pytest.fixture(scope="session")
def downstream_run(golden_dir, test_bam, reference, rep_lib_dir, gene_gff3, blacklist,
                   tmp_path_factory):
    """Stage + run each type's post-cns chain at most ONCE per session.

    Session-scoped and memoized because the chain costs ~40-100s per type; without this the two
    assertions below would each pay for their own run and double the gate's wall-clock.
    """
    cache = {}

    def run(mei_type):
        if mei_type not in cache:
            golden_type_dir = os.path.join(golden_dir, mei_type)
            if not os.path.isdir(golden_type_dir):
                pytest.skip(f"golden has no {mei_type} dir: {golden_type_dir}")
            _require(os.path.join(golden_type_dir, "candidate_disc_filtered_cns.txt"),
                     f"{mei_type} cns output", "XTEA2_GOLDEN_DIR")
            dest = str(tmp_path_factory.mktemp(f"ds_{mei_type}") / mei_type)
            os.makedirs(dest, exist_ok=True)
            _stage(golden_type_dir, dest, test_bam)
            steps = _build_steps(mei_type, dest, rep_lib_dir, reference, gene_gff3, blacklist,
                                 cores=8)
            _run_steps(steps, dest)
            cache[mei_type] = (golden_type_dir, dest)
        return cache[mei_type]

    return run


@pytest.mark.slow
@pytest.mark.parametrize("mei_type,vcf_alt", [("L1", "LINE1"), ("Alu", "ALU"), ("SVA", "SVA")])
def test_downstream_call_set_matches_golden(mei_type, vcf_alt, downstream_run):
    """THE GATE: VCF columns 1-8 must be identical to the golden.

    Scope is a deliberate choice (user, 2026-07-25): gate on the call set, not on every intermediate
    file, so that trimming work is not blocked by intermediate churn that never reaches the VCF.

    The byte-level comparison of every intermediate still RUNS, and prints what moved — it just
    does not fail the build. Keeping it is not ceremony: B7 (the SVA transduction divergence) was
    found only because intermediates were compared, and a trim that silently changes an
    intermediate is exactly the thing worth seeing in the log before it becomes a call-set bug.
    """
    golden_type_dir, dest = downstream_run(mei_type)
    golden_vcf = os.path.join(golden_type_dir, f"NA12878.chr1_{vcf_alt}.vcf")
    if not os.path.exists(golden_vcf):
        pytest.skip(f"no golden VCF for {mei_type}: {golden_vcf}")

    # --- diagnostics (never fail) -------------------------------------------------------------
    identical, genotype_diffs, tolerated, byte_diffs = _compare(mei_type, golden_type_dir, dest)
    print(f"\n[{mei_type}] {len(identical)} intermediate files byte-identical to golden")
    for t in tolerated:
        print(f"[{mei_type}] known divergence: {t}")
    for d in byte_diffs:
        print(f"[{mei_type}] NOTE — intermediate moved (not gated): {d}")
    for name, moved in genotype_diffs.items():
        print(f"[{mei_type}] {name}: {len(moved)} genotype flips (DF21 golden -> ONNX)")
        for chrom, pos, was, now in moved:
            print(f"    {chrom}:{pos}  {was} -> {now}")

    # --- the gate ------------------------------------------------------------------------------
    ours_vcf = os.path.join(dest, f"NA12878.chr1_{vcf_alt}.vcf")
    assert os.path.exists(ours_vcf), f"gVCF step produced no {ours_vcf}"
    cols = lambda p: [l.split("\t")[:8] for l in open(p) if not l.startswith("#")]
    g, w = cols(golden_vcf), cols(ours_vcf)
    assert len(g) == len(w), f"{mei_type}: call count {len(g)} golden vs {len(w)} ours"
    for a, b in zip(g, w):
        assert a == b, f"{mei_type}: VCF cols 1-8 differ at {a[0]}:{a[1]}\n  golden={a}\n  ours={b}"
