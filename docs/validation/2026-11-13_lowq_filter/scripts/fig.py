#!/usr/bin/env python
"""Figure: what the low-MAPQ clip filter removes, current (all alignments) vs primary-only.
(a) clip-check stage: share of candidates near a truth site vs others that the filter drops
    (sim_truth.txt, offline replay). (b) final calls: TP and FP the filter removes, i.e. present
    with the filter off and absent with it on (arms_cmp.metrics.tsv)."""
import re, sys
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib import font_manager as fm
FONT = "Nimbus Sans"
fm.findfont(FONT, fallback_to_default=False)
plt.rcParams["font.family"] = FONT
L = "/n/scratch/users/c/cos689/xtea2_review/lowq_filter"
ORDER = ["HG002_NovaSeq30x", "HG002_HiSeq30x", "HG002_ds20x", "HG002_ds15x", "HG002_ds10x",
         "HG00438", "HG02257", "HG02622"]
LAB = ["HG002\nNova30x", "HG002\nHiSeq30x", "HG002\n20x", "HG002\n15x", "HG002\n10x",
       "HG00438", "HG02257", "HG02622"]
st = {}
cur = None
for ln in open(f"{L}/sim_truth.txt"):
    m = re.match(r"(\S+): pass-off", ln)
    if m: cur = m.group(1); continue
    m = re.match(r"\s+(\S+)\s+drops truth\s+\d+ \(\s*([\d.]+)%\)\s+other\s+\d+ \(\s*([\d.]+)%\)", ln)
    if m: st[(cur, m.group(1))] = (float(m.group(2)), float(m.group(3)))
met = {}
for ln in list(open(f"{L}/arms_cmp.metrics.tsv"))[1:]:
    x = ln.split("\t"); met[(x[0], x[1])] = (int(x[5]), int(x[6]))
fig, ax = plt.subplots(1, 2, figsize=(11, 4.2))
xs = range(len(ORDER)); w = 0.2
C = {"cur": "#7f7f7f", "prim": "#2a6fdb"}
for j, (arm, col, nm) in enumerate((("cur", "scanner", "current"), ("prim", "est_prim", "primary only"))):
    t = [st[(n, col)][0] for n in ORDER]; o = [st[(n, col)][1] for n in ORDER]
    ax[0].bar([i + (2 * j - 1.5) * w for i in xs], t, w, color=C[arm], label=f"{nm}: near truth")
    ax[0].bar([i + (2 * j - 0.5) * w for i in xs], o, w, color=C[arm], alpha=0.45, label=f"{nm}: other")
ax[0].set_xticks(list(xs)); ax[0].set_xticklabels(LAB, fontsize=8)
ax[0].set_ylabel("Candidates dropped at clip check (%)")
ax[0].set_title("a  Clip-check stage (offline replay)", loc="left", fontsize=10)
ax[0].set_ylim(0, 15); ax[0].legend(fontsize=7, frameon=False, ncol=2, loc="upper center")
for j, (arm, nm) in enumerate((("cur", "current"), ("prim", "primary only"))):
    tp = [met[(n, "off")][0] - met[(n, arm)][0] for n in ORDER]
    fp = [met[(n, "off")][1] - met[(n, arm)][1] for n in ORDER]
    ax[1].bar([i + (2 * j - 1.5) * w for i in xs], tp, w, color=C[arm], label=f"{nm}: TP removed")
    ax[1].bar([i + (2 * j - 0.5) * w for i in xs], fp, w, color=C[arm], alpha=0.45, label=f"{nm}: FP removed")
ax[1].set_xticks(list(xs)); ax[1].set_xticklabels(LAB, fontsize=8)
ax[1].set_ylabel("Final calls removed by the filter (vs off)")
ax[1].set_title("b  Final calls", loc="left", fontsize=10)
ax[1].legend(fontsize=7, frameon=False)
for a in ax:
    a.spines[["top", "right"]].set_visible(False)
fig.tight_layout()
fig.savefig(f"{L}/fig/lowq_filter_arms.png", dpi=200)
