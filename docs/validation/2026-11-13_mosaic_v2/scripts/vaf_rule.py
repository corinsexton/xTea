"""Offline comparison of a VAF-model germline rule with AF_PE ceilings, on postF_mosaic output
with disc_upper opened to 1.0 (pf_open, B15 counts). The VAF rule drops calls whose mosaic VAF
estimate (x_vaf_site.estimate, mosaic=True, calibrated on the run's germline-like postF_hc calls
of that family, background from 1000 random positions x CONTEXT_FACTOR) exceeds v_max. Scored
with offline_ceiling.score (no cross-family locus merge, so FP runs a little above the real
sweep). Fit v_max on BCM_300X + WashU_100X, test UW_100X, and leave one run out."""
import sys, numpy as np, pandas as pd
B = "/n/scratch/users/c/cos689/xtea2_review/mosaic_v2"
sys.path.insert(0, B + "/scripts")
sys.path.insert(0, "/n/scratch/users/c/cos689/xtea2_wt/vaf_site_v2/xtea")
import offline_ceiling as oc, x_vaf_site as xv
R = "/n/scratch/users/c/cos689/xtea2_testdata/vafmix/runs"
I = "/n/scratch/users/c/cos689/xtea2_review/counting_v2/mosaic/inputs"
RT = {"L1": 1, "Alu": 2, "SVA": 4}
RUNS = ["BCM_300X", "UW_100X", "WashU_100X"]
gs = pd.read_csv(f"{B}/ceiling/germ_stats.tsv", sep="\t").set_index(["run", "fam"])


def lib_of(n):
    h = dict(x.split("=") for x in open(f"{B}/counts/{n}.tsv.isize").readline()[1:].split())
    return xv.Library(float(h["rlth"]), float(h["mean_is"]), float(h["sd"]))


D = {}
for n in RUNS:
    lib = lib_of(n)
    r = pd.read_csv(f"{B}/counts/rnd_{n}.tsv", sep="\t")
    bg = xv.control_background(r[["clip", "fmap", "disc", "conc"]].values, lib)
    d = oc.records(n)
    d["vaf"] = np.nan
    for f in oc.FAMS:
        inp = {}
        for ln in open(f"{I}/{n}_{f}.input_cv2.txt"):
            x = ln.split(); inp[(x[0], x[1])] = x
        hc = [inp[(x[0], x[1])] for x in (l.split() for l in open(
            f"{R}/{n}/{n}/{f}/candidate_disc_filtered_cns.txt.high_confident.post_filtering.txt")) if (x[0], x[1]) in inp]
        calib = xv.raw_site_inputs(hc, RT[f])
        m = (d.fam == f).values
        c, S, t = xv.raw_site_inputs([x.split() for x in d.raw[m]], RT[f])
        v, lo, hi, P = xv.estimate(c, lib, S, t, mosaic=True, calib=calib, background=dict(bg))
        d.loc[m, "vaf"] = v
        print(n, f, "calib", P["source"], P["n_calib"], "e_PE %.2f l_PE %.1f" % (P["e_PE"], P["l_PE"]), file=sys.stderr)
    d["ceil_het"] = [0.60 * gs.loc[(n, f), "med_het"] for f in d.fam]
    D[n] = d


def f1(n, keep):
    return oc.score(D[n], keep)["pooled"]


VG = np.round(np.concatenate([np.arange(0.05, 0.30, 0.01), np.arange(0.30, 0.61, 0.02)]), 2)
CG = np.round(np.arange(0.20, 0.51, 0.01), 2)
RULES = {"vaf": (VG, lambda d, p: d.vaf <= p), "afpe_fixed": (CG, lambda d, p: d.af_pe <= p)}
rows = []
for name, (grid, kf) in RULES.items():
    for test in RUNS:
        train = [r for r in RUNS if r != test]
        p = max(grid, key=lambda p: np.mean([f1(r, kf(D[r], p))["F1"] for r in train]))
        s = f1(test, kf(D[test], p))
        rows.append(dict(rule=name, held_out=test, param=p, **{k: s[k] for k in ("TP", "FP", "hg005", "F1")}))
for test in RUNS:
    for lab, keep in (("afpe_0.35", D[test].af_pe <= 0.35), ("afpe_0.60xmed_het", D[test].af_pe <= D[test].ceil_het)):
        s = f1(test, keep)
        rows.append(dict(rule=lab, held_out=test, param=np.nan, **{k: s[k] for k in ("TP", "FP", "hg005", "F1")}))
t = pd.DataFrame(rows)
t.to_csv(f"{B}/ceiling/vaf_rule_offline.tsv", sep="\t", index=False, float_format="%.4f")
print(t.round(4).to_string(index=False))
# VAF of TP vs FP vs hg005 calls at 0.35
for n in RUNS:
    d = D[n]
    print(n, "median VAF of calls with AF_PE > 0.30:", round(d.vaf[d.af_pe > 0.30].median(), 3),
          "<= 0.30:", round(d.vaf[d.af_pe <= 0.30].median(), 3))
