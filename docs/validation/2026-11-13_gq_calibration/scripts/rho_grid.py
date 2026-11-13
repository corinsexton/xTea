#!/usr/bin/env python
"""Grid over (rho_sr, rho_pe): per-sample EM fit at each rho pair (rho passed through to the
beta-binomial), then marginal log-likelihood (unsupervised), labeled log loss, GT errors and
GQ calibration per dataset. Uses the module's own EM with a per-type rho (monkeypatched
site_loglik) so the grid does not depend on the code change.
usage: rho_grid.py OUT.tsv
"""
import sys, math, itertools
import numpy as np
from common import L, load, counts, ORDER

GRID = [0.0, 0.005, 0.01, 0.02, 0.03, 0.05, 0.075, 0.1, 0.15, 0.2, 0.3]
_orig = L.site_loglik


def make_ll(rsr, rpe):
    def site_loglik(c, params, rho=None):
        c = np.asarray(c, float).reshape(-1, 4)
        out = np.empty((len(c), 2))
        for j, g in enumerate(("het", "hom")):
            out[:, j] = (L.betabinom_logpmf(c[:, 0], c[:, 0] + c[:, 1], params["sr_" + g], rsr)
                         + L.betabinom_logpmf(c[:, 2], c[:, 2] + c[:, 3], params["pe_" + g], rpe))
        return out
    return site_loglik


m = load()
out = open(sys.argv[1], "w")
out.write("rho_sr\trho_pe\tdataset\tsource\tmarg_ll\tn_lab\tlogloss\tgt_err\tgt_err_q\t"
          "n30\terr30\tpred30\tn60\terr60\n")
for rsr, rpe in itertools.product(GRID, GRID):
    L.site_loglik = make_ll(rsr, rpe)
    for n in ORDER:
        x = m[m.dataset == n]
        c = counts(x)
        P = L.fit_params(c, rho=0.05)
        ll = L.site_loglik(c, P)
        lp = ll + np.log([1 - P["hom_prior"], P["hom_prior"]])
        marg = float(np.logaddexp(lp[:, 0], lp[:, 1]).sum())
        r = L.genotype(c, P)
        gt = np.where(r["gt"] == "1/1", "hom", "het")
        lab = x.lab_q.values
        k = np.isin(lab, ["het", "hom"])
        ph = np.clip(r["post_hom"][k], 1e-12, 1 - 1e-12)
        y = lab[k] == "hom"
        logloss = float(-np.mean(np.where(y, np.log(ph), np.log(1 - ph))))
        lt = x.truth.values
        kt = np.isin(lt, ["het", "hom"])
        gq = r["gq"][k]
        err = gt[k] != lab[k]
        s30 = gq >= 30
        out.write("%g\t%g\t%s\t%s\t%.2f\t%d\t%.5f\t%d\t%d\t%d\t%d\t%.5f\t%d\t%d\n" % (
            rsr, rpe, n, P["source"], marg, k.sum(), logloss, (gt[kt] != lt[kt]).sum(),
            err.sum(), s30.sum(), err[s30].sum(), np.mean(10 ** (-gq[s30] / 10)),
            (gq >= 60).sum(), err[gq >= 60].sum()))
    out.flush()
