#!/usr/bin/env python
"""Median AF_SR = AF_CLIP/(AF_CLIP+AF_FMAP) of truth-matched final calls by family and truth GT,
control vs counting-v2 arm, for the 8 validation datasets; writes sr_af_medians.tsv and fig_sr_af.png.
usage: sr_af_summary.py OLD_RUNS NEW_RUNS OUTDIR"""
import os, sys
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from compare_arms import ORDER, G, H, S
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
plt.rcParams["font.family"] = "Nimbus Sans"
old, new, out = sys.argv[1:4]
def info(v):
    m = {}
    for l in open(v):
        if l.startswith("#"): continue
        x = l.split("\t"); d = dict(kv.split("=", 1) for kv in x[7].split(";") if "=" in kv)
        m[(x[0], int(x[1]))] = (int(d["AF_CLIP"]), int(d["AF_FMAP"]))
    return m
rows = []; vals = {}
for n in ORDER:
    truth = S.load_truth(H) if n.startswith("HG002") else S.load_truth(os.path.join(G, "truth", "truth_%s.tsv" % n), True, True)
    for arm, r in (("control", old), ("counting_v2", new)):
        v = os.path.join(r, n, n, n + ".vcf"); I = info(v)
        res = S.score(truth, S.load_calls([v]), 200)
        for fam in S.FAMILIES:
            for gt in ("het", "hom"):
                a = [c_ / (c_ + f_) for t, c in res[fam]["pairs"] if t.gt == gt
                     for c_, f_ in [I[(c.chrom, c.pos)]] if c_ + f_ > 0]
                vals[(n, arm, fam, gt)] = a
                rows.append((n, arm, fam, gt, len(a), np.median(a) if a else float("nan")))
with open(os.path.join(out, "sr_af_medians.tsv"), "w") as fo:
    fo.write("dataset\tarm\tfamily\ttruth_gt\tn\tmedian_af_sr\n")
    for r in rows: fo.write("%s\t%s\t%s\t%s\t%d\t%.3f\n" % r)
fig, axes = plt.subplots(2, 4, figsize=(13, 6.5), sharey=True)
for ax, n in zip(axes.flat, ORDER):
    pos = 0; ticks = []; labs = []
    for gt in ("het", "hom"):
        for arm, col in (("control", "#8c8c8c"), ("counting_v2", "#2a6fbb")):
            a = sum((vals[(n, arm, f, gt)] for f in S.FAMILIES), [])
            b = ax.boxplot([a], positions=[pos], widths=0.7, patch_artist=True, showfliers=False)
            for p in b["boxes"]: p.set_facecolor(col); p.set_alpha(0.75)
            for p in b["medians"]: p.set_color("black")
            pos += 1
        ticks.append(pos - 1.5); labs.append("truth %s" % gt); pos += 0.8
    ax.set_xticks(ticks); ax.set_xticklabels(labs); ax.set_title(n, fontsize=10)
    ax.axhline(0.5, ls=":", c="0.6", lw=0.8); ax.set_ylim(0, 1.02)
for ax in axes[:, 0]: ax.set_ylabel("AF_SR")
from matplotlib.patches import Patch
fig.legend(handles=[Patch(color="#8c8c8c", label="control (xTea2 0553edb)"),
                    Patch(color="#2a6fbb", label="counting-v2 (B15 + B16)")], loc="lower center", ncol=2, frameon=False)
fig.suptitle("Split-read AF (AF_CLIP / (AF_CLIP + AF_FMAP)) of truth-matched calls, all families pooled", fontsize=12)
fig.tight_layout(rect=(0, 0.05, 1, 0.96)); fig.savefig(os.path.join(out, "fig_sr_af.png"), dpi=150)
