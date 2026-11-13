"""Refit the per-site VAF model constants on per-pair (B14/B15) counts.

DEFAULT_CAL / DEFAULT_HOM_FRAC: fit_calibration on HG002 NovaSeq 30x truth het+hom sites.
MOSAIC_BACKGROUND: ML b at the BCM_300X somatic truth sites <= 0.5% VAF (default calibration).
CONTEXT_FACTOR: that fit / ML b at 1000 random clip peaks of BCM_300X.
Also: per-genome calibration fits (range), and the raw SR/PE fractions at HG005 germline sites in
the mixtures (for the germline_like threshold).
usage: vaf_refit.py [XTEA_DIR]"""
import sys, numpy as np, pandas as pd
XT = sys.argv[1] if len(sys.argv) > 1 else "/n/scratch/users/c/cos689/xtea2_wt/vaf_site_v2/xtea"
sys.path.insert(0, XT)
import x_vaf_site as xv
from scipy.optimize import minimize_scalar
C = "/n/scratch/users/c/cos689/xtea2_review/mosaic_v2/counts"
FAM = {"Alu": "ALU", "ALU": "ALU", "L1": "LINE1", "LINE1": "LINE1", "SVA": "SVA"}


def lib_of(name):
    h = dict(x.split("=") for x in open(f"{C}/{name}.tsv.isize").readline()[1:].split())
    return xv.Library(float(h["rlth"]), float(h["mean_is"]), float(h["sd"]))


def load(name):
    d = pd.read_csv(f"{C}/{name}.tsv", sep="\t", na_values=["NA"])
    S = np.array([xv.parse_svlen(s, FAM.get(f, f)) for s, f in zip(d.svlen, d.family)]) if "svlen" in d else np.full(len(d), 300.)
    t = np.array([xv.parse_tsdlen(x) for x in d.tsd_est])
    return d, lib_of(name), S, t, d[["clip", "fmap", "disc", "conc"]].values.astype(float)


def fit_b(c, v, y, cal, rho):
    out = {}
    for j, key in enumerate(("b_SR", "b_PE")):
        k, n = c[:, 2 * j], c[:, 2 * j] + c[:, 2 * j + 1]

        def f(lb):
            p = xv.p_alt(v, y, dict(cal, **{key: np.exp(lb)}))[j]
            return -xv.betabinom_logpmf(k, n, p, rho).sum()
        out[key] = float(np.exp(minimize_scalar(f, bounds=(-12, 0), method="bounded").x))
    return out


if __name__ == "__main__":
    print("== calibration fits (germline truth het+hom sites)")
    for n in ("HG002_NovaSeq30x", "HG002_HiSeq30x", "HG002_ds10x", "HG002_2x250", "HG00438", "HG02257", "HG02622"):
        d, lib, S, t, c = load(n)
        m = d["gt"].isin(["het", "hom"]).values
        fit = xv.fit_calibration(c[m], lib, S[m], t[m])
        print(n, lib, {k: (round(v, 3) if isinstance(v, float) else v) for k, v in fit.items()})
    d, lib, S, t, c = load("BCM_300X")
    cal = dict(xv.DEFAULT_CAL)
    m = ((d.kind == "somatic") & (d.expected_vaf <= 0.005)).values
    bm = fit_b(c[m], d.expected_vaf.values[m], xv.yields(lib, S[m], t[m], cal), cal, xv.MOSAIC_RHO)
    r = pd.read_csv(f"{C}/rnd_BCM_300X.tsv", sep="\t")
    rc = r[["clip", "fmap", "disc", "conc"]].values.astype(float)
    br = xv.fit_background(rc, lib, cal=cal)
    print("== BCM_300X n_low=%d b_MEI %s  b_random %s  factor %s" % (
        m.sum(), {k: round(v, 4) for k, v in bm.items()}, {k: round(v, 4) for k, v in br.items()},
        {k: round(bm[k] / br[k], 3) for k in bm}))
    print("== raw fractions at mixture sites (germline_like threshold)")
    for n in ("BCM_300X", "WashU_100X", "UW_100X", "NYGC_100X"):
        d, lib, S, t, c = load(n)
        with np.errstate(divide="ignore", invalid="ignore"):
            sr = c[:, 0] / (c[:, 0] + c[:, 1]); pe = c[:, 2] / (c[:, 2] + c[:, 3])
        h = (d.kind == "hg005").values
        s = ((d.kind == "somatic") & (d.expected_vaf >= 0.05)).values
        for thr in (0.2, 0.25, 0.3, 0.333):
            print(n, "thr", thr, "HG005 pass SR&PE %.3f  PE-only %.3f | somatic>=5%% pass %.3f" % (
                np.mean((sr[h] >= 0.2) & (pe[h] >= thr)), np.mean(pe[h] >= thr),
                np.mean((sr[s] >= 0.2) & (pe[s] >= thr))))
        print(n, "HG005 PE q05/q25/q50 %s  SR q05/q50 %s" % (np.round(np.nanpercentile(pe[h], [5, 25, 50]), 3),
                                                          np.round(np.nanpercentile(sr[h], [5, 50]), 3)))
