#!/usr/bin/env python
"""Before/after GQ calibration. Before = GQ in the runs_b15 VCFs; after = the gq-calib module
(x_genotype_lik.fit_params + genotype, i.e. the pipeline code path) rerun on the same INFO counts.
Two label sets: 'truth' (score_final.sh labels) and 'lab_q' (GIAB Q100 v5.0q for HG002, HPRC
tier1 for the 1KG samples). Writes calib_tables.tsv, gt_errors.tsv, records_after.tsv.gz and
reliability.png into OUT.
usage: calib_tables.py OUT_DIR
"""
import os, sys
import numpy as np, pandas as pd
from common import L, load, counts, ORDER

TRAIN = ["HG002_NovaSeq30x", "HG00438"]
BINS = {"before": ([0, 10, 20, 30, 50, 60, 1000], ["0-10", "10-20", "20-30", "30-50", "50-60", "60"]),
        "after": ([0, 10, 20, 25, 29, 1000], ["0-10", "10-20", "20-25", "25-29", "29-30"])}
out = sys.argv[1]
os.makedirs(out, exist_ok=True)
m = load()
parts = []
for n in ORDER:
    x = m[m.dataset == n].copy()
    c = counts(x)
    P = L.fit_params(c)
    r = L.genotype(c, P)
    x["gt_after"] = np.where(r["gt"] == "1/1", "hom", "het")
    x["gq_after"] = r["gq"]
    parts.append(x)
m = pd.concat(parts)
m["set"] = np.where(m.dataset.isin(TRAIN), "train", "heldout")
m[["dataset", "idx", "chrom", "pos", "family", "AF_CLIP", "AF_FMAP", "AF_DISC", "AF_CONCORDNT",
   "gt", "gq", "gt_after", "gq_after", "truth", "q100", "lab_q"]].to_csv(
    os.path.join(out, "records_after.tsv.gz"), sep="\t", index=False)


def table(lab, arm):
    gtc, gqc = ("gt", "gq") if arm == "before" else ("gt_after", "gq_after")
    y = m[m[lab].isin(["het", "hom"])].copy()
    y["err"] = (y[lab] != y[gtc]).astype(int)
    y["pred"] = 10 ** (-y[gqc] / 10)
    b, names = BINS[arm]
    y["gq_bin"] = pd.cut(y[gqc], b, right=False, labels=names)
    rows = []
    for grp, sub in [("pooled", y), ("heldout", y[y.set == "heldout"]), ("train", y[y.set == "train"])] + \
            [(n, y[y.dataset == n]) for n in ORDER]:
        t = sub.groupby("gq_bin", observed=True).agg(n=("err", "size"), errors=("err", "sum"),
                                                     obs_err=("err", "mean"), pred_err=("pred", "mean"))
        t = t.reset_index()
        t.insert(0, "group", grp); t.insert(0, "arm", arm); t.insert(0, "labels", lab)
        rows.append(t)
    return pd.concat(rows)


T = pd.concat([table(lab, arm) for lab in ("truth", "lab_q") for arm in ("before", "after")])
T.to_csv(os.path.join(out, "calib_tables.tsv"), sep="\t", index=False, float_format="%.5f")

g = []
for n in ORDER:
    x = m[m.dataset == n]
    row = {"dataset": n, "set": x.set.iloc[0], "gt_changed": int((x["gt"] != x.gt_after).sum())}
    for lab in ("truth", "lab_q"):
        k = x[lab].isin(["het", "hom"])
        row["n_" + lab] = int(k.sum())
        row["err_before_" + lab] = int((x[k][lab] != x[k]["gt"]).sum())
        row["err_after_" + lab] = int((x[k][lab] != x[k].gt_after).sum())
    g.append(row)
G = pd.DataFrame(g)
G.to_csv(os.path.join(out, "gt_errors.tsv"), sep="\t", index=False)
pd.set_option("display.width", 250)
print(G.to_string(index=False))
for lab in ("truth", "lab_q"):
    for grp in ("pooled", "heldout"):
        print("==", lab, grp)
        print(T[(T.labels == lab) & (T.group == grp)].drop(columns=["labels", "group"]).round(4).to_string(index=False))

# reliability plot
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib import font_manager
plt.rcParams["font.family"] = "Nimbus Sans"
from scipy.stats import beta as B


def ci(k, n):
    lo = B.ppf(0.025, k, n - k + 1) if k > 0 else 0.0
    hi = B.ppf(0.975, k + 1, n - k)
    return lo, hi


fig, axes = plt.subplots(1, 2, figsize=(10, 4.6), sharey=True)
col = {"before": "#9a9a9a", "after": "#2a6fbb"}
FLOOR = 1e-4
for ax, lab, title in zip(axes, ("lab_q", "truth"),
                          ("Labels: GIAB Q100 (HG002), HPRC tier1 (1KG)",
                           "Labels: score_final.sh (xTea-paper HG002 set)")):
    ax.plot([FLOOR, 1], [FLOOR, 1], color="#444444", lw=0.8, ls="--", zorder=0)
    for arm, mk in (("before", "o"), ("after", "s")):
        t = T[(T.labels == lab) & (T.group == "heldout") & (T.arm == arm)]
        j = 0
        for _, r in t.iterrows():
            lo, hi = ci(r.errors, r.n)
            xo, yo = max(r.pred_err, FLOOR), max(r.obs_err, FLOOR)
            if r.pred_err < FLOOR:          # spread the clamped bins so they stay readable
                xo = FLOOR / 2.2 ** j; j += 1
            ax.errorbar(xo, yo, yerr=[[yo - max(lo, FLOOR)], [max(hi - yo, 0)]], fmt=mk,
                        color=col[arm], ms=6, capsize=2, lw=1,
                        label=arm if _ == t.index[0] else None)
            ax.annotate("%s\n%d/%d" % (r.gq_bin, r.errors, r.n), (xo, yo), fontsize=7,
                        xytext=(5, -10 if arm == "before" else 4), textcoords="offset points",
                        color=col[arm])
    ax.set_xscale("log"); ax.set_yscale("log")
    ax.set_xlim(FLOOR / 4, 0.6); ax.set_ylim(FLOOR * 0.7, 0.6)
    ax.set_xlabel("Predicted GT error (mean 10^(-GQ/10))")
    ax.set_title(title, fontsize=10)
    ax.text(FLOOR / 3.5, FLOOR * 0.75, "predicted <1e-4: drawn from 1e-4, spread left; 0 errors drawn at 1e-4", fontsize=7,
            color="#666666")
axes[0].set_ylabel("Observed GT error (95% CI)")
handles = [plt.Line2D([], [], marker="o", color=col["before"], ls="", label="before (GQ cap 60)"),
           plt.Line2D([], [], marker="s", color=col["after"], ls="", label="after (error floor 0.001)")]
axes[0].legend(handles=handles, loc="upper left", frameon=False, fontsize=9)
fig.suptitle("GQ reliability, held-out datasets (HiSeq 30x, 20x/15x/10x, HG02257, HG02622)",
             fontsize=11)
fig.tight_layout()
fig.savefig(os.path.join(out, "reliability.png"), dpi=150)
