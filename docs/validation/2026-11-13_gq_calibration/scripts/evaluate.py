#!/usr/bin/env python
"""Re-genotype every dataset offline under a configuration (rho_sr, rho_pe, eps) and write
per-record GQ/GT plus calibration tables.

perr' = eps + (1 - 2 eps) * perr     (error floor; argmax unchanged)
usage: evaluate.py  (writes into ../eval/)
"""
import os, sys, math
import numpy as np, pandas as pd
from common import L, load, counts, ORDER, BINS, REV

TRAIN = ["HG002_NovaSeq30x", "HG00438"]


def site_ll(c, P, rsr, rpe):
    out = np.empty((len(c), 2))
    for j, g in enumerate(("het", "hom")):
        out[:, j] = (L.betabinom_logpmf(c[:, 0], c[:, 0] + c[:, 1], P["sr_" + g], rsr)
                     + L.betabinom_logpmf(c[:, 2], c[:, 2] + c[:, 3], P["pe_" + g], rpe))
    return out


def run(m, rsr, rpe):
    """per-record called GT and raw perr (no floor), per-sample EM at (rsr, rpe)."""
    orig = L.site_loglik
    L.site_loglik = lambda c, P, rho=None: site_ll(np.asarray(c, float).reshape(-1, 4), P, rsr, rpe)
    try:
        res = []
        for n in ORDER:
            x = m[m.dataset == n]
            c = counts(x)
            P = L.fit_params(c)
            r = L.genotype(c, P)
            post = r["post_hom"]
            perr = np.where(post >= 0.5, 1 - post, post)
            res.append(pd.DataFrame({"dataset": n, "idx": x.idx.values,
                                     "gt_new": np.where(post >= 0.5, "hom", "het"),
                                     "perr": perr}))
    finally:
        L.site_loglik = orig
    return m.merge(pd.concat(res), on=["dataset", "idx"])


def gq_of(perr, eps, cap=60):
    p = eps + (1 - 2 * eps) * perr
    return np.floor(np.minimum(cap, -10 * np.log10(np.maximum(p, 1e-300))) + 1e-9)


def fit_eps(x, lab):
    k = x[lab].isin(["het", "hom"]) & x.dataset.isin(TRAIN)
    y = (x[lab] != x.gt_new)[k].values
    pe = x.perr[k].values
    grid = np.concatenate([[0], np.logspace(-5, -1.5, 141)])
    ll = [np.sum(np.where(y, np.log(np.maximum(e + (1 - 2 * e) * pe, 1e-300)),
                          np.log1p(-(e + (1 - 2 * e) * pe)))) for e in grid]
    i = int(np.argmax(ll))
    return grid[i], ll[i] - ll[0], grid, np.array(ll)


def calib(x, lab, gqcol, group):
    k = x[lab].isin(["het", "hom"])
    y = x[k].copy()
    y["err"] = (y[lab] != y.gt_new).astype(int)
    y["pred"] = 10 ** (-y[gqcol] / 10)
    y["bin"] = pd.cut(y[gqcol], BINS, right=False, labels=["0-10", "10-20", "20-30", "30-50",
                                                            "50-60", "60"])
    return y.groupby(group + ["bin"], observed=True).agg(n=("err", "size"), errors=("err", "sum"),
                                                          obs_err=("err", "mean"),
                                                          pred_err=("pred", "mean")).reset_index()


if __name__ == "__main__":
    out = os.path.join(REV, "eval"); os.makedirs(out, exist_ok=True)
    m = load()
    cfgs = [(0.05, 0.05), (0.02, 0.02), (0.02, 0.03)]
    for rsr, rpe in cfgs:
        x = run(m, rsr, rpe)
        for lab in ("truth", "lab_q"):
            e, dll, grid, ll = fit_eps(x, lab)
            print("rho %.3f/%.3f lab %s: eps_hat %.5f (GQ max %.1f) dLL %.2f; 95%% CI eps %s" % (
                rsr, rpe, lab, e, -10 * math.log10(e) if e else 60, dll,
                (grid[ll >= ll.max() - 1.92].min(), grid[ll >= ll.max() - 1.92].max())))
        x.to_csv(os.path.join(out, "perr_rho%g_%g.tsv.gz" % (rsr, rpe)), sep="\t", index=False)
