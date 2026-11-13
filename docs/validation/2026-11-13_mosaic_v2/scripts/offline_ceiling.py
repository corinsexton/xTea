"""Offline ceiling evaluation on postF_mosaic output with disc_upper opened to 1.0: apply a
per-record AF_PE ceiling (or any per-record rule) and score one-to-one (+-200 bp, same family)
against vafsite/sites/mix.tsv somatic sites; calls at listed HG005 germline sites counted too.
Module: records(name) -> DataFrame; score(df, keep_mask) -> per-family and pooled dict."""
import sys, numpy as np, pandas as pd
sys.path.insert(0, "/n/scratch/users/c/cos689/xtea2_testdata/gtval")
import score_gt as S
P = "/n/scratch/users/c/cos689/xtea2_review/mosaic_v2/pf_open"
MIX = "/n/scratch/users/c/cos689/xtea2_testdata/vafsite/sites/mix.tsv"
FAMS = {"L1": "L1", "Alu": "Alu", "SVA": "SVA"}
truth, hg005 = [], []
for ln in open(MIX):
    x = ln.rstrip("\n").split("\t")
    if x[0] == "chrom" or x[0] in S.SEX:
        continue
    s = S.Site(x[0], int(x[1]), x[2], None, None, float(x[4]) if x[4] != "NA" else -1.0)
    (truth if x[6] == "somatic" else hg005).append(s)


def records(name):
    rows = []
    for fam in FAMS:
        for ln in open(f"{P}/{name}/{fam}/mosaic_open.txt"):
            f = ln.split()
            a, b = int(f[2]), int(f[3])
            lo, hi = min(a, b), max(a, b)
            pos = hi if lo <= 0 else lo
            c = [float(f[i]) for i in (35, 36, 39, 40)]
            rows.append(dict(fam=fam, chrom=f[0], pos=pos, clip=c[0], fmap=c[1], disc=c[2], conc=c[3],
                             tsd=f[4], svlen=f[47] if len(f) > 47 else "-1", raw=ln))
    d = pd.DataFrame(rows)
    d["af_pe"] = d.disc / (d.disc + d.conc).where(lambda x: x > 0, np.nan)
    d["af_pe"] = d.af_pe.fillna(0.0)
    return d


def score(d, keep):
    d = d[np.asarray(keep, bool)]
    out, pool = {}, dict(n_truth=0, calls=0, TP=0, hg005=0)
    for fam, tfam in FAMS.items():
        t = [s for s in truth if s.family == tfam]
        sub = d[d.fam == fam]
        c = [S.Site(ch, int(p), tfam, None, None, -1.0) for ch, p in zip(sub.chrom, sub.pos)]
        m = S.match_one_to_one(t, c, 200)
        used = {ci for ci, _ in m.values()}
        rest = [c[i] for i in range(len(c)) if i not in used]
        hg = len(S.match_one_to_one([s for s in hg005 if s.family == tfam], rest, 200))
        out[fam] = dict(n_truth=len(t), calls=len(c), TP=len(m), hg005=hg)
        for k in pool:
            pool[k] += out[fam][k]
    for r in list(out.values()) + [pool]:
        r["FP"] = r["calls"] - r["TP"]
        p = r["TP"] / r["calls"] if r["calls"] else 0.0
        s = r["TP"] / r["n_truth"] if r["n_truth"] else 0.0
        r["F1"] = 2 * p * s / (p + s) if p + s else 0.0
    out["pooled"] = pool
    return out


if __name__ == "__main__":
    for n in ("BCM_300X", "UW_100X", "WashU_100X"):
        d = records(n)
        for v in (0.30, 0.35, 0.52):
            r = score(d, d.af_pe <= v)["pooled"]
            print(n, v, r)
