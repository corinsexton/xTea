"""Validation of the per-site VAF model on per-pair (B14/B15) counts, production x_vaf_site.

Mixtures: mosaic estimate (control background from random positions x CONTEXT_FACTOR; default
calibration and self-calibration on germline-like calls) at somatic truth sites, VAF/truth by bin
and family. Germline: 7 genomes, het/hom VAF and 95% CI coverage; het VAF by family and SVLEN.
Part C ablation: the spanning-pair term E[n_S] added to C_alt (l_PE refit on HG002 NovaSeq 30x;
production leaves it out), and an empirical insert-size histogram in place of the normal.
usage: vaf_eval.py OUTDIR"""
import os, sys, numpy as np, pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib import font_manager
for f in font_manager.findSystemFonts():
    if "NimbusSans-" in f and f.endswith(".otf"):
        font_manager.fontManager.addfont(f)
plt.rcParams.update({"font.family": "Nimbus Sans", "font.size": 9, "axes.spines.top": False,
                     "axes.spines.right": False})
sys.path.insert(0, "/n/scratch/users/c/cos689/xtea2_wt/vaf_site_v2/xtea")
import x_vaf_site as xv
C = "/n/scratch/users/c/cos689/xtea2_review/mosaic_v2/counts"
O = sys.argv[1]
os.makedirs(O, exist_ok=True)
BINS = [0, 0.004, 0.0075, 0.015, 0.035, 0.075, 0.2, 1]
LAB = ["0.25%", "0.5%", "1%", "2-3%", "5-6%", "7-20%", ">20%"]
FAM = {"Alu": "ALU", "ALU": "ALU", "L1": "LINE1", "LINE1": "LINE1", "SVA": "SVA"}
MIX = ["BCM_300X", "WashU_100X", "UW_100X", "NYGC_100X"]
GERM = ["HG002_NovaSeq30x", "HG002_HiSeq30x", "HG002_ds10x", "HG002_2x250", "HG00438", "HG02257", "HG02622"]
_yields = xv.yields


def isize_of(name):
    lines = open(f"{C}/{name}.tsv.isize").read().split("\n")
    h = dict(x.split("=") for x in lines[0][1:].split())
    hist = np.array([[float(v) for v in l.split()] for l in lines[1:] if l.strip()])
    return h, hist


def lib_of(name, emp=False):
    h, hist = isize_of(name)
    lib = xv.Library(float(h["rlth"]), float(h["mean_is"]), float(h["sd"]))
    if emp:
        F, n = hist[:, 0], hist[:, 1]
        keep = F <= np.percentile(np.repeat(F, n.astype(int)), 99.9)
        lib.F, lib.w = F[keep], n[keep] / n[keep].sum()
    return lib


def inputs(d):
    S = np.array([xv.parse_svlen(s, FAM.get(f, f)) for s, f in zip(d.svlen, d.family)])
    t = np.array([xv.parse_tsdlen(x) for x in d.tsd_est])
    return S, t, d[["clip", "fmap", "disc", "conc"]].values.astype(float)


def yields_span(lib, S, t, cal):
    """Production yields plus the spanning-pair term E[n_S], n_S(F) = max(0, F - 2R - S + 2c)."""
    y = _yields(lib, S, t, cal)
    S = np.atleast_1d(np.asarray(S, float))
    n_s = np.maximum(0.0, lib.F[None, :] - 2 * lib.R - S[:, None] + 2 * xv.CLIP_C)
    return dict(y, C_alt=y["C_alt"] + n_s @ lib.w)


VARIANTS = {"model": (_yields, False), "with_span": (yields_span, False), "emp_isize": (_yields, True)}


def set_variant(name):
    xv.yields = VARIANTS[name][0]
    d, lib0 = pd.read_csv(f"{C}/HG002_NovaSeq30x.tsv", sep="\t", na_values=["NA"]), None
    d = d[d["gt"].isin(["het", "hom"])]
    S, t, c = inputs(d)
    fit = xv.fit_calibration(c, lib_of("HG002_NovaSeq30x", VARIANTS[name][1]), S, t)
    return {k: fit[k] for k in xv.FIT_KEYS}, fit["w"]


def run_variant(var):
    cal0, w0 = set_variant(var)
    xv.DEFAULT_CAL.update(cal0)
    emp = VARIANTS[var][1]
    sites = []
    for n in MIX:
        d = pd.read_csv(f"{C}/{n}.tsv", sep="\t", na_values=["NA"])
        lib = lib_of(n, emp)
        r = pd.read_csv(f"{C}/rnd_{n}.tsv", sep="\t")
        bg = xv.control_background(r[["clip", "fmap", "disc", "conc"]].values, lib)
        S, t, c = inputs(d)
        none = (np.zeros((0, 4)), np.zeros(0), np.zeros(0))
        tags = (("control", bg, none), ("constant", xv.MOSAIC_BACKGROUND, none), ("control_selfcal", bg, None))
        if var != "model":
            tags = tags[:1]
        for tag, b, calib in tags:
            vaf, lo, hi, P = xv.estimate(c, lib, S, t, mosaic=True, background=dict(b), calib=calib)
            s = d.assign(vaf=vaf, lo=lo, hi=hi, bam=n, bg=tag, b_SR=b["b_SR"], b_PE=b["b_PE"],
                         mean_is=lib.mean_is, cal_src=P["source"], n_calib=P["n_calib"], e_SR=P["e_SR"],
                         e_PE=P["e_PE"], l_PE=P["l_PE"], S=S, variant=var)
            sites.append(s[s.kind.isin(["somatic", "hg005"])])
    ms = pd.concat(sites)
    ms["fam"] = ms.family.map(FAM)
    ms["bin"] = pd.cut(ms.expected_vaf, BINS, labels=LAB)
    ms["ratio"] = ms.vaf / ms.expected_vaf
    ms["cov"] = (ms.lo <= ms.expected_vaf) & (ms.expected_vaf <= ms.hi)
    gsites = []
    for n in GERM:
        d = pd.read_csv(f"{C}/{n}.tsv", sep="\t", na_values=["NA"])
        lib = lib_of(n, emp)
        alld = d.copy()
        d = d[d["gt"].isin(["het", "hom"])].reset_index(drop=True)
        S, t, c = inputs(d)
        Sa, ta, ca = inputs(alld)
        vaf, lo, hi, P = xv.estimate(c, lib, S, t, calib=(ca, Sa, ta))
        g = d.assign(sample=n, fam=d.family.map(FAM), S=S, vaf=vaf, lo=lo, hi=hi, mean_is=lib.mean_is,
                     R=lib.R, cal_src=P["source"], e_SR=P["e_SR"], e_PE=P["e_PE"], l_SR=P["l_SR"],
                     l_PE=P["l_PE"], hom_frac=P["hom_frac"], variant=var)
        g["truth"] = np.where(g["gt"] == "hom", 1.0, 0.5)
        g["cov"] = (g.lo <= g.truth) & (g.truth <= g.hi)
        gsites.append(g)
    xv.yields = _yields
    return cal0, w0, ms, pd.concat(gsites)


def summ(g):
    return pd.Series({"n": len(g), "med_ratio": g.ratio.median(), "q25": g.ratio.quantile(.25),
                      "q75": g.ratio.quantile(.75), "cov": g["cov"].mean(),
                      "frac0": (g.vaf < 1e-3).mean(), "ci_hi_med": g.hi.median()})


if __name__ == "__main__":
    pd.set_option("display.width", 250)
    import copy
    base_cal = copy.deepcopy(xv.DEFAULT_CAL)
    res = {}
    for var in VARIANTS:
        xv.DEFAULT_CAL.clear(); xv.DEFAULT_CAL.update(base_cal)
        res[var] = run_variant(var)
        print("== variant", var, "HG002 Nova fit", {k: round(v, 3) for k, v in res[var][0].items()}, "w", round(res[var][1], 3))
    xv.DEFAULT_CAL.clear(); xv.DEFAULT_CAL.update(base_cal)
    _, _, ms, gs = res["model"]
    som = ms[ms.kind == "somatic"]
    tab_bin = som.groupby(["bg", "bam", "bin"], observed=True).apply(summ).reset_index()
    tab_fam = som[som.bg == "control"].groupby(["bam", "fam", "bin"], observed=True).apply(summ).reset_index()
    tab_bin.to_csv(O + "/mosaic_calib_by_bin.tsv", sep="\t", index=False)
    tab_fam.to_csv(O + "/mosaic_calib_by_family.tsv", sep="\t", index=False)
    print("== mosaic by bin")
    print(tab_bin.pivot_table(index=["bam", "bin"], columns="bg", values=["med_ratio", "cov", "frac0"], observed=True).round(2).to_string())
    print("== mosaic by family (control), bins >= 1%")
    print(tab_fam[tab_fam.bin.isin(["1%", "2-3%", "5-6%", "7-20%"])].round(2).to_string(index=False))
    print("== backgrounds / calibration", ms.groupby(["bam", "bg"])[["b_SR", "b_PE", "mean_is", "cal_src", "n_calib", "e_SR", "e_PE", "l_PE"]].first().round(4).to_string())
    h5 = ms[(ms.kind == "hg005") & (ms.bg == "control_selfcal")]
    print("== HG005 germline sites in mixtures (selfcal), VAF median (truth het 0.4175 / hom 0.835)")
    print(h5.groupby("bam").vaf.describe().round(3).to_string())
    gt = gs.groupby(["sample", "gt"]).apply(lambda g: pd.Series({
        "n": len(g), "R": g.R.iloc[0], "mean_is": g.mean_is.iloc[0], "cal_src": g.cal_src.iloc[0],
        "e_SR": g.e_SR.iloc[0], "e_PE": g.e_PE.iloc[0], "l_PE": g.l_PE.iloc[0], "hom_frac": g.hom_frac.iloc[0],
        "med_vaf": g.vaf.median(), "cov": g["cov"].mean(), "ciw": (g.hi - g.lo).median()})).reset_index()
    print("== germline"); print(gt.round(3).to_string(index=False))
    gt.to_csv(O + "/germline_vaf.tsv", sep="\t", index=False)
    SB = [0, 200, 280, 330, 1000, 5500, 7000, 20001]
    SL = ["<200", "200-280", "280-330", "330-1k", "1k-5.5k", "5.5k-7k", ">7k"]
    # Part C ablation tables
    rows = []
    for var, (cal0, w0, m2, g2) in res.items():
        s2 = m2[(m2.kind == "somatic") & (m2.bg == "control")]
        for (bam, fam, b), g in s2.groupby(["bam", "fam", "bin"], observed=True):
            rows.append(dict(variant=var, set="mixture", sample=bam, fam=fam, bin=b, n=len(g),
                             med_ratio=g.ratio.median(), cov=g["cov"].mean()))
        g2 = g2.assign(Sbin=pd.cut(g2.S, SB, labels=SL))
        for (smp, fam), g in g2[g2["gt"] == "het"].groupby(["sample", "fam"]):
            rows.append(dict(variant=var, set="germline_het_by_sample", sample=smp, fam=fam, bin="het", n=len(g),
                             med_ratio=g.vaf.median() / 0.5, cov=g["cov"].mean()))
        for (fam, gtp), g in g2.groupby(["fam", "gt"]):
            rows.append(dict(variant=var, set="germline", sample="all7", fam=fam, bin=gtp, n=len(g),
                             med_ratio=g.vaf.median() / g.truth.iloc[0], cov=g["cov"].mean()))
        for (fam, sb), g in g2[g2["gt"] == "het"].groupby(["fam", "Sbin"], observed=True):
            rows.append(dict(variant=var, set="germline_het_svlen", sample="all7", fam=fam, bin=sb, n=len(g),
                             med_ratio=g.vaf.median() / 0.5, cov=g["cov"].mean()))
    ab = pd.DataFrame(rows)
    ab.to_csv(O + "/spanning_pair_ablation.tsv", sep="\t", index=False)
    print("== Part C ablation (median VAF / truth)")
    print(ab.pivot_table(index=["set", "sample", "fam", "bin"], columns="variant", values="med_ratio", observed=True).round(3).to_string())
    print(ab.pivot_table(index=["set", "sample", "fam", "bin"], columns="variant", values="cov", observed=True).round(3).to_string())
    gs["Sbin"] = pd.cut(gs.S, SB, labels=SL)
    fe = gs[gs["gt"] == "het"].groupby(["fam", "Sbin"], observed=True).agg(
        n=("vaf", "size"), med=("vaf", "median"), cov=("cov", "mean")).reset_index()
    fe.to_csv(O + "/germline_het_by_family_svlen.tsv", sep="\t", index=False)
    gs.to_csv(O + "/germline_persite.tsv.gz", sep="\t", index=False)
    ms.to_csv(O + "/mixture_persite.tsv.gz", sep="\t", index=False)

    # figure 1: observed vs expected
    s = som[som.bg == "control"]
    fams = ["ALU", "LINE1", "SVA"]
    fig, ax = plt.subplots(len(MIX), 3, figsize=(7.2, 2.2 * len(MIX)), sharex=True, sharey=True)
    rng = np.random.default_rng(0)
    FLOOR = 5e-4
    for i, n in enumerate(MIX):
        for j, f in enumerate(fams):
            a = ax[i, j]
            g = s[(s.bam == n) & (s.fam == f)]
            x = g.expected_vaf.values * np.exp(rng.normal(0, 0.08, len(g)))
            a.scatter(x, np.maximum(g.vaf.values, FLOOR), s=3, alpha=0.25, color="#3b6ea8", lw=0)
            gb = g.groupby("expected_vaf").vaf
            m, lo, hi = gb.median(), gb.quantile(0.25), gb.quantile(0.75)
            a.errorbar(m.index, np.maximum(m.values, FLOOR), yerr=[np.maximum(m - lo, 0), np.maximum(hi - m, 0)],
                       fmt="o", color="#c0392b", ms=3.5, lw=1, capsize=0)
            a.plot([1e-3, 1], [1e-3, 1], color="0.4", lw=0.8, ls="--")
            a.axhline(FLOOR, color="0.8", lw=0.6)
            a.set_xscale("log"); a.set_yscale("log"); a.set_xlim(1.2e-3, 0.5); a.set_ylim(3e-4, 1.2)
            if i == 0:
                a.set_title(f)
            if j == 0:
                a.set_ylabel(f"{n}\nestimated VAF")
            if i == len(MIX) - 1:
                a.set_xlabel("expected VAF")
            a.text(0.03, 0.97, "n=%d" % len(g), transform=a.transAxes, va="top", fontsize=7)
    fig.suptitle("Mosaic VAF at HapMap-mixture truth sites, per-pair counts (estimates of 0 drawn at 5e-4;\n"
                 "red: median and IQR; background measured at random positions)", fontsize=9)
    fig.tight_layout(); fig.savefig(O + "/mixture_obs_vs_exp.png", dpi=200)
    # figure 2: germline het/hom per sample
    fig, a = plt.subplots(1, 1, figsize=(5.2, 3.2))
    data, pos, col = [], [], []
    for i, n in enumerate(GERM):
        for k, (g, c) in enumerate((("het", "#3b6ea8"), ("hom", "#c0392b"))):
            data.append(gs[(gs["sample"] == n) & (gs["gt"] == g)].vaf.dropna().values); pos.append(i * 3 + k); col.append(c)
    bp = a.boxplot(data, positions=pos, widths=0.8, showfliers=False, patch_artist=True, medianprops=dict(color="k", lw=1))
    for b, c in zip(bp["boxes"], col):
        b.set_facecolor(c); b.set_alpha(0.6)
    a.axhline(0.5, color="#3b6ea8", ls="--", lw=0.8); a.axhline(1.0, color="#c0392b", ls="--", lw=0.8)
    a.set_xticks([i * 3 + 0.5 for i in range(len(GERM))])
    a.set_xticklabels([n.replace("HG002_", "HG002\n") for n in GERM], fontsize=7)
    a.set_ylim(0, 1.3); a.set_ylabel("VAF (box: IQR; blue het, red hom)")
    a.set_title("Germline VAF, per-pair counts")
    fig.tight_layout(); fig.savefig(O + "/germline_vaf.png", dpi=200)
    # figure 3: het VAF vs SVLEN, model vs no-span
    fig, ax = plt.subplots(1, 2, figsize=(7.2, 3.0), sharey=True)
    for a, var in zip(ax, ("with_span", "model")):
        h = res[var][3]; h = h[h["gt"] == "het"]
        for f, c in (("ALU", "#3b6ea8"), ("LINE1", "#c0392b"), ("SVA", "#2e8b57")):
            g = h[h.fam == f]
            e = np.unique(np.quantile(g.S, np.linspace(0, 1, 9)))
            b = pd.cut(g.S, e, include_lowest=True)
            a.scatter(g.S, g.vaf, s=2, alpha=0.08, color=c, lw=0)
            a.plot(g.groupby(b, observed=True).S.median(), g.groupby(b, observed=True).vaf.median(), "o-", color=c, ms=3, lw=1.2, label=f)
        a.axhline(0.5, color="0.4", ls="--", lw=0.8); a.set_xscale("log"); a.set_ylim(0, 1)
        a.set_xlabel("insertion length (SVLEN, bp)")
        a.set_title({"with_span": "with spanning-pair term (l_PE refit)", "model": "without (production)"}[var])
    ax[0].set_ylabel("het VAF (line: binned median)"); ax[1].legend(frameon=False, fontsize=8)
    fig.tight_layout(); fig.savefig(O + "/het_vaf_vs_svlen_span_ablation.png", dpi=200)
    print("outputs in", O)
