import sys, numpy as np
sys.path.insert(0, "/n/scratch/users/c/cos689/xtea2_wt/vaf_site_v2/xtea")
import x_genotype_lik as gl
R = "/n/scratch/users/c/cos689/xtea2_testdata/vafmix/runs"
I = "/n/scratch/users/c/cos689/xtea2_review/counting_v2/mosaic/inputs"
for n in ("BCM_300X", "UW_100X", "WashU_100X"):
    for f in ("Alu", "L1", "SVA"):
        inp = {}
        for ln in open(f"{I}/{n}_{f}.input_cv2.txt"):
            x = ln.split(); inp[(x[0], x[1])] = [float(x[i]) for i in (35, 36, 39, 40)]
        hc = [ln.split() for ln in open(f"{R}/{n}/{n}/{f}/candidate_disc_filtered_cns.txt.high_confident.post_filtering.txt")]
        c = np.array([inp[(x[0], x[1])] for x in hc if (x[0], x[1]) in inp])
        sr = c[:, 0] / np.maximum(c[:, 0] + c[:, 1], 1); pe = c[:, 2] / np.maximum(c[:, 2] + c[:, 3], 1)
        out = [n, f, len(c)]
        for t in (0.0, 0.1, 0.2, 0.25, 0.3):
            m = (sr >= t) & (c[:, 2] + c[:, 3] > 0)
            g = gl.genotype(c[m], gl.fit_params(c[m])) if m.sum() > 10 else None
            het = (g["gt"] == "0/1") if g is not None else np.zeros(0, bool)
            out += ["sr>=%.2f n=%d het=%d med=%.3f" % (t, m.sum(), het.sum(), np.median(pe[m][het]) if het.any() else np.nan)]
        print(" | ".join(map(str, out)))
