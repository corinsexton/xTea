"""Per run and family: the germline AF_PE statistics a mosaic run can compute from its own
high-confidence (postF_hc) call list, B15 counts (fields 35/36/39/40 of the B15-recounted
inputs, joined by chrom/pos). Writes ceiling/germ_stats.tsv."""
import numpy as np, pandas as pd
R = "/n/scratch/users/c/cos689/xtea2_testdata/vafmix/runs"
I = "/n/scratch/users/c/cos689/xtea2_review/counting_v2/mosaic/inputs"
O = "/n/scratch/users/c/cos689/xtea2_review/mosaic_v2/ceiling/germ_stats.tsv"
rows = []
for n in ("BCM_300X", "UW_100X", "WashU_100X"):
    for f in ("Alu", "L1", "SVA"):
        inp = {}
        for ln in open(f"{I}/{n}_{f}.input_cv2.txt"):
            x = ln.split()
            inp[(x[0], x[1])] = [float(x[i]) for i in (35, 36, 39, 40)]
        hc = [ln.split() for ln in open(f"{R}/{n}/{n}/{f}/candidate_disc_filtered_cns.txt.high_confident.post_filtering.txt")]
        own = np.array([[float(x[i]) for i in (35, 36, 39, 40)] for x in hc]).reshape(-1, 4)
        c = np.array([inp[(x[0], x[1])] for x in hc if (x[0], x[1]) in inp]).reshape(-1, 4)
        with np.errstate(divide="ignore", invalid="ignore"):
            sr, pe = c[:, 0] / c[:, :2].sum(1), c[:, 2] / c[:, 2:].sum(1)
        gl = (sr >= 0.2) & (pe >= 0.2)
        het = gl & (pe < 0.8)
        rows.append(dict(run=n, fam=f, n_hc=len(hc), n_joined=len(c), n_gl=int(gl.sum()), n_het=int(het.sum()),
                         med_all=np.nanmedian(pe), med_gl=np.median(pe[gl]) if gl.any() else np.nan,
                         med_het=np.median(pe[het]) if het.any() else np.nan,
                         med_sr_gl=np.median(sr[gl]) if gl.any() else np.nan))
d = pd.DataFrame(rows)
d.to_csv(O, sep="\t", index=False, float_format="%.4f")
print(d.round(3).to_string(index=False))
