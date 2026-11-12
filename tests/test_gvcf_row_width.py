"""Row-width contract between the post-filter steps and x_gvcf._cvt_to_vcf_rcd.

Layout: fields 0-48 are the core (48 = postF_all's repeat-copy label), 49 = gene
annotation, 50.. = genotype FORMAT values. The mosaic callers used to append a second
repeat-copy label, which pushed the gene annotation to field 50; gVCF then wrote the
repeat label as GENE_INFO and the gene string as the sample's GT. These tests pin
(1) the mosaic writer passing rows through unchanged, and (2) gVCF reading each
legal width correctly and refusing anything else instead of mis-reading it.
"""
import os
import sys

import pytest

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "xtea"))

import x_genotype_lik  # noqa: E402
from x_gvcf import gVCF  # noqa: E402
from x_mosaic_calling import MosaicCaller  # noqa: E402
sys.path.insert(0, os.path.dirname(__file__))
from test_gvcf_multi import _row  # noqa: E402

N_FMT = len(x_genotype_lik.FORMAT_KEYS)
GENE = "intron:ENSG00000078369.19:GNB1"
FMT_VALS = ["0/1", "35", "0,35,400", "0.42", "0.38", "0.40", "0.30-0.50"]
assert len(FMT_VALS) == N_FMT


def _fields(extra=()):
    return _row("chr1", 1000, 1015).split("\t") + list(extra)


def _cvt(l_fields):
    s = gVCF()._cvt_to_vcf_rcd(l_fields, "LINE1", "AAAAAAAAAAAAAAA", ".", "A")
    cols = s.rstrip("\n").split("\t")
    info = dict(kv.split("=", 1) for kv in cols[7].split(";") if "=" in kv)
    return info, cols[8], cols[9]


# ---------------------------------------------------------------- mosaic writer
class _Parser:
    def get_ins_sub_type(self, rcd):
        return "two_side"


class _Tprt:
    def has_polyA_signal(self, rcd):
        return True

    def is_polyA_dominant_two_side(self, rcd):
        return False


class _AF:
    def is_qualified_mosaic_rcd(self, line, cutoff):
        return True


class _Annot:
    def is_within_repeat_region_bool(self, chrm, pos):
        return False


@pytest.mark.parametrize("method", ["call_mosaic_L1_from_bulk",
                                    "call_mosaic_Alu_from_bulk",
                                    "call_mosaic_SVA_from_bulk"])
def test_mosaic_writer_keeps_row_unchanged(tmp_path, method):
    line = _row("chr1", 1000, 1015) + "\n"
    rcd = line.split() + [line]  # load_in_xTEA_rslt's shape: fields + raw line
    sf_out = str(tmp_path / "out.txt")
    caller = MosaicCaller(str(tmp_path), 1)
    getattr(caller, method)([rcd], _Parser(), _Tprt(), _AF(), _Annot(), 0.0, sf_out)
    with open(sf_out) as fin:
        out = fin.read()
    assert out == line
    assert len(out.split()) == 49


# ---------------------------------------------------------------- gVCF reader
def test_core_only_row():
    info, fmt, gt = _cvt(_fields())
    assert info["GENE_INFO"] == "."
    assert (fmt, gt) == ("GT", "./.")


def test_gene_annotated_row():
    info, fmt, gt = _cvt(_fields([GENE]))
    assert info["GENE_INFO"] == GENE
    assert (fmt, gt) == ("GT", "./.")


def test_gt_only_row():
    info, fmt, gt = _cvt(_fields([GENE, "1/1"]))
    assert info["GENE_INFO"] == GENE
    assert (fmt, gt) == ("GT", "1/1")


def test_full_genotype_row():
    info, fmt, gt = _cvt(_fields([GENE] + FMT_VALS))
    assert info["GENE_INFO"] == GENE
    assert fmt == ":".join(x_genotype_lik.FORMAT_KEYS)
    assert gt == ":".join(FMT_VALS)


def test_old_mosaic_shape_is_rejected():
    # what the old mosaic writer produced: extra repeat label at 49, gene at 50
    with pytest.raises(ValueError, match="not a genotype"):
        _cvt(_fields(["Fall_in_Alu_copy", GENE]))


@pytest.mark.parametrize("n_extra", [2, N_FMT - 1, N_FMT + 1])
def test_unknown_width_is_rejected(n_extra):
    with pytest.raises(ValueError, match="columns"):
        _cvt(_fields([GENE] + ["0/1"] * n_extra))


def test_full_genotype_row_with_bad_gt_is_rejected():
    with pytest.raises(ValueError, match="not a genotype"):
        _cvt(_fields([GENE, "Fall_in_Alu_copy"] + FMT_VALS[1:]))
