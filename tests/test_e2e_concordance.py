"""End-to-end concordance gate: pipeline-under-test vs the master golden.

Two layers:

* ``test_golden_self_concordance`` — parse the real golden VCFs and compare them to themselves.
  This proves the VCF parser + comparator work on real xTea output and that the golden is
  internally consistent. Runs now (golden exists).

* ``test_xtea2_vs_golden`` — the real gate: the xtea2 output VCFs must recover ~100% of the
  master high-confidence calls (±25 bp + type). Skipped until ``XTEA2_OUTPUT_DIR`` is set (i.e.
  once xtea2 can run end-to-end). Wired now so it activates automatically.
"""
import glob
import os

import pytest

from tests.utils.vcf_io import parse_vcf, parse_vcfs
from tests.utils.concordance import compare_callsets, dedupe_conflict_loci

# Fraction of master calls the new pipeline must recover. 1.0 = the agreed 100% HC gate.
CONCORDANCE_TARGET = 1.0


def test_golden_parses_nonempty(golden_vcfs):
    total = sum(len(parse_vcf(p)) for p in golden_vcfs.values())
    # Alu(27)+SVA(5)+L1(0) = 32 on the chr1 golden; assert the parser sees the expected order.
    assert total >= 1, "golden VCFs parsed to zero calls — parser or golden problem"


def test_golden_self_concordance(golden_vcfs):
    truth = parse_vcfs(golden_vcfs.values())
    r = compare_callsets(truth, list(truth))
    assert r.sensitivity == 1.0, r.summary()
    assert r.ppv == 1.0, r.summary()
    assert not r.missed and not r.extra


@pytest.mark.parametrize("mei_type", ["ALU", "SVA", "LINE1"])
def test_golden_per_type_self_concordance(golden_vcfs, mei_type):
    if mei_type not in golden_vcfs:
        pytest.skip(f"no golden VCF for {mei_type}")
    truth = parse_vcf(golden_vcfs[mei_type])
    r = compare_callsets(truth, list(truth))
    assert r.sensitivity == 1.0 and r.ppv == 1.0, r.summary()


def _find_xtea2_vcfs(out_dir):
    # Only the PER-TYPE VCFs (mirrors the golden_vcfs per-type globbing). Excludes the
    # C.1 merged VCF ({sample}.vcf at the tree root), which otherwise double-counts every
    # call (per-type + merged) and halves PPV. Match the per-type suffixes emitted by the
    # gVCF step: *_ALU.vcf / *_LINE1.vcf / *_SVA.vcf.
    hits = []
    for suf in ("*_ALU.vcf", "*_LINE1.vcf", "*_SVA.vcf"):
        hits += glob.glob(os.path.join(out_dir, "**", suf), recursive=True)
    return sorted(set(hits))


def test_xtea2_vs_golden(golden_vcfs, xtea2_output_dir):
    """The real reproducibility gate (active once xtea2 produces VCFs).

    `truth` is deduplicated with the same locus-table rule xtea2 itself applies
    (dedupe_conflict_loci) before comparing. Master's per-type VCFs are generated independently
    and never re-merged, so a small number of loci near a type boundary (~18-20 genome-wide at
    WGS scale; see docs/handoff/MOVE_A_PLAN.md) are listed under MORE THAN ONE type in the raw
    golden. xtea2 correctly resolves each such locus to exactly one winning type before ever
    writing a VCF (that's the whole point of the move-A locus table) -- so a query call under the
    resolved type can only ever satisfy ONE of the golden's duplicate entries; comparing against
    the undeduplicated union makes every remaining duplicate an unconditional false negative,
    regardless of how correct xtea2's output is. Deduping truth with the identical rule first
    makes this an apples-to-apples comparison.
    """
    truth = dedupe_conflict_loci(parse_vcfs(golden_vcfs.values()))
    query_paths = _find_xtea2_vcfs(xtea2_output_dir)
    if not query_paths:
        pytest.skip(f"no VCFs found under XTEA2_OUTPUT_DIR={xtea2_output_dir}")
    query = parse_vcfs(query_paths)
    r = compare_callsets(truth, query, window=25)
    assert r.sensitivity >= CONCORDANCE_TARGET, (
        f"xtea2 recovered only {r.sensitivity:.3f} of master calls. {r.summary()}\n"
        f"missed: {[(c.chrom, c.pos, c.mei_type) for c in r.missed][:20]}"
    )
