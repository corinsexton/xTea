"""Shared loaders for the GQ-calibration analysis."""
import sys
import numpy as np
import pandas as pd
REV = "/n/scratch/users/c/cos689/xtea2_review/gq_calib"
sys.path.insert(0, "/n/scratch/users/c/cos689/xtea2_wt/gq_calib/xtea")
import x_genotype_lik as L  # noqa
CNT = ["AF_CLIP", "AF_FMAP", "AF_DISC", "AF_CONCORDNT"]
ORDER = ["HG002_NovaSeq30x", "HG002_HiSeq30x", "HG002_ds20x", "HG002_ds15x", "HG002_ds10x",
         "HG00438", "HG02257", "HG02622"]
BINS = [0, 10, 20, 30, 50, 60, 1000]


def load(rec=REV + "/records_b15.tsv.gz", q=REV + "/q100_hg002.tsv.gz"):
    """Records + labels. truth = score_final label (xTea-paper set for HG002, HPRC tier1 for 1KG).
    lab_q = Q100 label for HG002 autosomal sites inside the Q100 benchmark bed (het/hom only),
    the HPRC tier1 label for the 1KG samples."""
    d = pd.read_csv(rec, sep="\t", keep_default_na=False)
    qq = pd.read_csv(q, sep="\t", keep_default_na=False)
    m = d.merge(qq, on=["dataset", "idx"], how="left")
    hg = m.dataset.str.startswith("HG002")
    auto = ~m.chrom.isin(["chrX", "chrY"])
    ok = (m.inbed == 1) & m.q100.isin(["het", "hom"]) & auto
    m["lab_q"] = np.where(hg, np.where(ok, m.q100, "NA"), m.truth)
    return m


def counts(x):
    return x[CNT].astype(float).values
