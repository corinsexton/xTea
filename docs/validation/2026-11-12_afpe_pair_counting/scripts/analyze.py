#!/usr/bin/env python
"""Summarize the cutoff sweeps (sweeps/*.tsv) and draw the figures.

Germline (afpe, td): choose on HG002_NovaSeq30x + HG00438 pooled (TP and calls summed across the two
datasets and the three families), take the middle of the plateau that holds the maximum F1 (values
whose F1 is within PLATEAU of the max, contiguous with it), round to 2 decimals, then compare the
chosen value with the converted one on each held-out dataset.
Mosaic (disc): leave one run out across BCM_300X, UW_100X and WashU_100X (choose on the other two,
evaluate on the held-out one); sensitivity by truth VAF bin.
Writes summary_*.tsv and fig_*.png in this directory and prints a report.
"""
import csv, glob, os, sys
from collections import defaultdict
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

T = os.path.dirname(os.path.abspath(__file__))
FONT = "Nimbus Sans"
plt.rcParams.update({"font.family": FONT, "font.size": 9})
PLATEAU = 0.0005
TUNE_SETS = ["HG002_NovaSeq30x", "HG00438"]
HELD = ["HG002_HiSeq30x", "HG002_ds10x", "HG002_ds15x", "HG002_ds20x", "HG02257", "HG02622"]
CONV = {"afpe": 0.14, "td": 0.148, "disc": 0.52}
MASTER = {"afpe": 0.075, "td": 0.08, "disc": 0.35}
OFF = {"afpe": -1.0, "td": 0.0, "disc": 1.0}
REC = {"afpe": 0.14, "td": 0.15, "disc": 0.30}   # recommended (see summary.txt / loro_disc.txt)
MOSAIC = ["BCM_300X", "UW_100X", "WashU_100X"]
FAMS = ["Alu", "L1", "SVA", "pooled"]


def load(pattern):
    d = defaultdict(dict)            # d[(dataset, family)][value] = row
    for f in sorted(glob.glob(os.path.join(T, "sweeps", pattern))):
        for r in csv.DictReader(open(f), delimiter="\t"):
            for k in ("n_truth", "calls", "TP", "FP"):
                r[k] = int(r[k])
            d[(r["dataset"], r["family"])][float(r["value"])] = r
    return d


def f1(tp, calls, nt):
    p = tp / calls if calls else 0.0
    s = tp / nt if nt else 0.0
    return 2 * p * s / (p + s) if p + s else 0.0


def combo(d, sets, fam, v):
    rs = [d[(s, fam)][v] for s in sets]
    tp = sum(r["TP"] for r in rs); c = sum(r["calls"] for r in rs); nt = sum(r["n_truth"] for r in rs)
    return tp, c - tp, f1(tp, c, nt)


def choose(d, sets, fam, grid):
    sc = [(v, combo(d, sets, fam, v)[2]) for v in grid]
    best = max(s for _, s in sc)
    i = max(range(len(sc)), key=lambda k: sc[k][1])
    lo = hi = i
    while lo > 0 and sc[lo - 1][1] >= best - PLATEAU:
        lo -= 1
    while hi < len(sc) - 1 and sc[hi + 1][1] >= best - PLATEAU:
        hi += 1
    mid = round((sc[lo][0] + sc[hi][0]) / 2, 2)
    mid = min(grid, key=lambda g: abs(g - mid))
    return mid, best, sc[lo][0], sc[hi][0]


def fig(d, param, datasets, chosen, fname, title):
    fig_, ax = plt.subplots(1, 3, figsize=(11, 3.4))
    cmap = plt.get_cmap("tab10")
    for k, ds in enumerate(datasets):
        if (ds, "pooled") not in d:
            continue
        g = sorted(v for v in d[(ds, "pooled")] if v != OFF["afpe"] or param != "afpe")
        ref = d[(ds, "pooled")][min(g, key=lambda x: abs(x - CONV[param]))]
        rows = [d[(ds, "pooled")][v] for v in g]
        ls = "-" if ds in TUNE_SETS or param == "disc" else "--"
        ax[0].plot(g, [float(r["F1"]) for r in rows], ls, color=cmap(k), lw=1.2, label=ds)
        ax[1].plot(g, [r["TP"] - ref["TP"] for r in rows], ls, color=cmap(k), lw=1.2)
        ax[2].plot(g, [r["FP"] - ref["FP"] for r in rows], ls, color=cmap(k), lw=1.2)
    for a_, lab in zip(ax, ("F1 (pooled families)", "TP minus TP at converted value",
                            "FP minus FP at converted value")):
        a_.axvline(CONV[param], color="0.3", ls=":", lw=1)
        a_.axvline(MASTER[param], color="0.6", ls=":", lw=1)
        if chosen is not None:
            a_.axvline(chosen, color="C3", lw=1, ls="--")
        a_.axvline(REC[param], color="k", lw=1.2)
        a_.set_ylabel(lab); a_.set_xlabel(title)
        a_.spines[["top", "right"]].set_visible(False)
    fig_.legend(*ax[0].get_legend_handles_labels(), loc="lower center", ncol=8, fontsize=7,
                frameon=False)
    fig_.suptitle("%s: black = recommended (%s); red dashed = F1-plateau middle, %s (%s); "
                  "dark dotted = converted (%s); light dotted = master (%s).%s"
                  % (title, REC[param], "pooled over runs" if param == "disc" else "tuning sets", chosen,
                     CONV[param], MASTER[param],
                     "" if param == "disc" else " Solid line = tuning set, dashed = held out."), fontsize=8)
    fig_.tight_layout(rect=(0, 0.07, 1, 0.95))
    fig_.savefig(os.path.join(T, fname), dpi=200)
    plt.close(fig_)


out = open(os.path.join(T, "summary.txt"), "w")


def P(*x):
    print(*x); print(*x, file=out)


for param in ("afpe", "td"):
    d = load("*_%s.tsv" % param)
    if not all((s, "pooled") in d for s in TUNE_SETS):
        P("## %s: tuning sets incomplete, skipped" % param); continue
    grid = sorted(v for v in d[(TUNE_SETS[0], "pooled")] if v >= 0)
    P("\n## %s  (master %s, converted %s, off %s)" % (param, MASTER[param], CONV[param], OFF[param]))
    ch = {}
    for fam in FAMS:
        mid, best, lo, hi = choose(d, TUNE_SETS, fam, grid)
        ch[fam] = mid
        P("tune %-6s max F1 %.4f plateau [%.2f, %.2f] -> middle %.2f; F1 at converted %.4f, at master %.4f"
          % (fam, best, lo, hi, mid, combo(d, TUNE_SETS, fam, CONV[param])[2],
             combo(d, TUNE_SETS, fam, MASTER[param])[2]))
    with open(os.path.join(T, "summary_%s.tsv" % param), "w") as fo:
        fo.write("dataset\tset\tfamily\tvalue\tlabel\tcalls\tTP\tFP\tF1\n")
        for ds in TUNE_SETS + HELD:
            for fam in FAMS:
                if (ds, fam) not in d:
                    continue
                for lab, v in (("off", OFF[param]), ("master", MASTER[param]),
                               ("converted", CONV[param]), ("chosen", ch["pooled"])):
                    if v not in d[(ds, fam)]:
                        continue
                    r = d[(ds, fam)][v]
                    fo.write("\t".join(map(str, [ds, "tune" if ds in TUNE_SETS else "held_out", fam, v, lab,
                                                 r["calls"], r["TP"], r["FP"], r["F1"]])) + "\n")
    P("held-out, chosen %.2f vs converted %.2f (pooled families):" % (ch["pooled"], CONV[param]))
    for ds in TUNE_SETS + HELD:
        if (ds, "pooled") not in d:
            continue
        a_ = d[(ds, "pooled")][CONV[param]]; b = d[(ds, "pooled")][ch["pooled"]]
        o = d[(ds, "pooled")].get(OFF[param]); m = d[(ds, "pooled")].get(MASTER[param])
        P("  %-17s %-8s conv TP %4d FP %3d F1 %s | chosen dTP %+d dFP %+d (per 1000 calls %+.2f/%+.2f) F1 %s"
          " | master dTP %+d dFP %+d | off dTP %+d dFP %+d"
          % (ds, "tune" if ds in TUNE_SETS else "held", a_["TP"], a_["FP"], a_["F1"],
             b["TP"] - a_["TP"], b["FP"] - a_["FP"], 1000 * (b["TP"] - a_["TP"]) / a_["calls"],
             1000 * (b["FP"] - a_["FP"]) / a_["calls"], b["F1"], m["TP"] - a_["TP"], m["FP"] - a_["FP"],
             o["TP"] - a_["TP"], o["FP"] - a_["FP"]))
    fig(d, param, TUNE_SETS + HELD, ch["pooled"], "fig_%s.png" % param,
        {"afpe": "germline AF_PE floor", "td": "transduction AF_PE floor"}[param])

# mosaic
d = load("*_disc_new*.tsv")
dold = load("*_disc_old.tsv")
if all((s, "pooled") in d for s in MOSAIC):
    grid = sorted(d[(MOSAIC[0], "pooled")])
    P("\n## disc_upper (mosaic; master 0.35, converted 0.52, off 1.0)")
    chosen_all = choose(d, MOSAIC, "pooled", grid)
    P("all three runs: max F1 %.4f plateau [%.2f, %.2f] -> %.2f" % (chosen_all[1], chosen_all[2],
                                                                     chosen_all[3], chosen_all[0]))
    for fam in FAMS[:3]:
        P("  family %s: %s" % (fam, choose(d, MOSAIC, fam, grid)))
    for h in MOSAIC:
        rest = [m for m in MOSAIC if m != h]
        mid = choose(d, rest, "pooled", grid)
        a_ = d[(h, "pooled")][CONV["disc"]]; b = d[(h, "pooled")][mid[0]]
        P("LORO hold %-10s chosen on others %.2f (plateau %.2f-%.2f): held-out dTP %+d dFP %+d vs 0.52"
          % (h, mid[0], mid[2], mid[3], b["TP"] - a_["TP"], b["FP"] - a_["FP"]))
    bins = [k for k in d[(MOSAIC[0], "pooled")][grid[0]] if k.startswith("sens_")]
    with open(os.path.join(T, "summary_disc.tsv"), "w") as fo:
        fo.write("dataset\tinputs\tvalue\tcalls\tTP\tFP\thg005\tF1\t" + "\t".join(bins) + "\n")
        for h in MOSAIC:
            rows = [("old", 0.35, dold.get((h, "pooled"), {}).get(0.35)),
                    ("old", 1.0, dold.get((h, "pooled"), {}).get(1.0))]
            rows += [("new", v, d[(h, "pooled")].get(v)) for v in sorted({0.36, 0.44, 0.52, chosen_all[0], 0.7, 1.0})]
            for inp, v, r in rows:
                if r is None:
                    continue
                fo.write("\t".join(map(str, [h, inp, v, r["calls"], r["TP"], r["FP"], r["hg005"], r["F1"]]
                                       + [r[b] for b in bins])) + "\n")
                P("  %-10s %s %.2f calls %5d TP %4d FP %5d %s" % (h, inp, v, r["calls"], r["TP"], r["FP"],
                                                               " ".join(r[b] for b in bins)))
    fig(d, "disc", MOSAIC, chosen_all[0], "fig_disc.png", "mosaic disc_upper (new-scale inputs)")
else:
    P("## disc: mosaic sweeps incomplete")
