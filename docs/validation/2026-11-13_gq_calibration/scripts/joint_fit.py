#!/usr/bin/env python
"""Labeled fit of (rho_sr, rho_pe, eps) on the training sets (HG002 NovaSeq 30x, HG00438):
labeled log-likelihood of the called-GT-correct outcome, P(correct) = 1 - perr', with
perr' = eps + (1-2eps) perr. Also reports held-out log-likelihood and GT errors.
usage: joint_fit.py LABEL OUT.tsv     (LABEL = lab_q or truth)
"""
import sys
import numpy as np
from common import load
from evaluate import run, TRAIN

lab, outp = sys.argv[1:3]
GR = [0.02, 0.03, 0.04, 0.05, 0.06, 0.075, 0.1, 0.125, 0.15]
EPS = np.concatenate([[0], np.logspace(-5, -1.5, 71)])
m = load()
fo = open(outp, "w")
fo.write("rho_sr\trho_pe\teps\ttrain_ll\theldout_ll\ttrain_gt_err\theldout_gt_err\n")
for rs in GR:
    for rp in GR:
        x = run(m, rs, rp)
        k = x[lab].isin(["het", "hom"])
        x = x[k]
        tr = x.dataset.isin(TRAIN).values
        y = (x[lab] != x.gt_new).values
        pe = x.perr.values

        def ll(e, s):
            p = e + (1 - 2 * e) * pe[s]
            return float(np.sum(np.where(y[s], np.log(np.maximum(p, 1e-300)), np.log1p(-p))))
        lls = [ll(e, tr) for e in EPS]
        i = int(np.argmax(lls))
        fo.write("%g\t%g\t%.6f\t%.3f\t%.3f\t%d\t%d\n" % (rs, rp, EPS[i], lls[i], ll(EPS[i], ~tr),
                                                        y[tr].sum(), y[~tr].sum()))
        fo.flush()
