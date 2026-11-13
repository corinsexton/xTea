"""Per run and family: AF_PE of the germline hc calls (B15 counts), split by the pipeline
genotyper's GT; also AF_PE at calls matching listed HG005 germline sites."""
import sys, numpy as np, pandas as pd
sys.path.insert(0, "/n/scratch/users/c/cos689/xtea2_wt/vaf_site_v2/xtea")
import x_genotype_lik as gl
R = "/n/scratch/users/c/cos689/xtea2_testdata/vafmix/runs"
I = "/n/scratch/users/c/cos689/xtea2_review/counting_v2/mosaic/inputs"
MIX = "/n/scratch/users/c/cos689/xtea2_testdata/vafsite/sites/mix.tsv"
mix = pd.read_csv(MIX, sep="\t")
FAMMAP = {"Alu": "Alu", "L1": "L1", "SVA": "SVA", "LINE1": "L1", "ALU": "Alu"}
mix["fam"] = mix.family.map(FAMMAP)
rows = []
for n in ("BCM_300X", "UW_100X", "WashU_100X"):
    allc = []
    for f in ("Alu", "L1", "SVA"):
        inp = {}
        for ln in open(f"{I}/{n}_{f}.input_cv2.txt"):
            x = ln.split()
            inp[(x[0], x[1])] = [float(x[i]) for i in (35, 36, 39, 40)]
        hc = [ln.split() for ln in open(f"{R}/{n}/{n}/{f}/candidate_disc_filtered_cns.txt.high_confident.post_filtering.txt")]
        keys = [(x[0], x[1]) for x in hc if (x[0], x[1]) in inp]
        c = np.array([inp[k] for k in keys])
        allc.append((f, keys, c))
    # genotyper fitted on all families pooled (as the shared tail does)
    call = np.vstack([c for _, _, c in allc])
    P = gl.fit_params(call)
    for f, keys, c in allc + [("all", sum([k for _, k, _ in allc], []), call)]:
        g = gl.genotype(c, P)
        pe = g["af_pe"]
        het = (g["gt"] == "0/1") & np.isfinite(pe)
        h5 = mix[(mix.kind == "hg005") & ((mix.fam == f) | (f == "all"))]
        pos = {}
        for ch, p in zip(h5.chrom, h5.pos):
            pos.setdefault(ch, []).append(p)
        at5 = np.array([any(abs(int(k[1]) - q) <= 200 for q in pos.get(k[0], [])) for k in keys])
        rows.append(dict(run=n, fam=f, n_hc=len(keys), n_het=int(het.sum()), n_hom=int((g["gt"] == "1/1").sum()),
                         med_het_afpe=np.median(pe[het]) if het.any() else np.nan,
                         q25_het=np.nanpercentile(pe[het], 25) if het.any() else np.nan,
                         med_all_afpe=np.nanmedian(pe),
                         n_hg005=int(at5.sum()), med_hg005_afpe=np.nanmedian(pe[at5]) if at5.any() else np.nan,
                         p05_hg005_afpe=np.nanpercentile(pe[at5], 5) if at5.any() else np.nan,
                         pe_het=P["pe_het"], src=P["source"]))
d = pd.DataFrame(rows)
pd.set_option("display.width", 250)
print(d.round(3).to_string(index=False))
d.to_csv(sys.argv[1] if len(sys.argv) > 1 else "/dev/null", sep="\t", index=False)
