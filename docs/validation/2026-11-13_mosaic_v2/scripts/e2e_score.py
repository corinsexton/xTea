"""Score mosaic VCFs as sweep.py does (one-to-one, same family, +-200 bp, against the kind=somatic
sites of vafsite/sites/mix.tsv; calls at listed HG005 germline sites counted apart).
usage: e2e_score.py LABEL=VCF [LABEL=VCF ...]  -> TSV on stdout (per family and pooled)."""
import sys
sys.path.insert(0, "/n/scratch/users/c/cos689/xtea2_testdata/gtval")
import score_gt as S  # noqa: E402

MIX = "/n/scratch/users/c/cos689/xtea2_testdata/vafsite/sites/mix.tsv"
truth, hg005 = [], []
for ln in open(MIX):
    x = ln.rstrip("\n").split("\t")
    if x[0] == "chrom" or x[0] in S.SEX:
        continue
    s = S.Site(x[0], int(x[1]), x[2], None, None, float(x[4]) if x[4] != "NA" else -1.0)
    (truth if x[6] == "somatic" else hg005).append(s)

print("label\tfamily\tn_truth\tcalls\tTP\tFP\thg005\tF1")
for arg in sys.argv[1:]:
    lab, vcf = arg.split("=", 1)
    calls = S.load_calls([vcf])
    pool = dict(n_truth=0, calls=0, TP=0, hg005=0)
    rows = []
    for fam in S.FAMILIES:
        t = [s for s in truth if s.family == fam]
        c = [s for s in calls if s.family == fam]
        m = S.match_one_to_one(t, c, 200)
        used = {ci for ci, _ in m.values()}
        rest = [c[i] for i in range(len(c)) if i not in used]
        hg = len(S.match_one_to_one([s for s in hg005 if s.family == fam], rest, 200))
        d = dict(family=fam, n_truth=len(t), calls=len(c), TP=len(m), hg005=hg)
        rows.append(d)
        for k in pool:
            pool[k] += d[k]
    pool["family"] = "pooled"
    for d in rows + [pool]:
        fp = d["calls"] - d["TP"]
        f1 = 2 * d["TP"] / (d["n_truth"] + d["calls"]) if d["n_truth"] + d["calls"] else 0.0
        print("\t".join(map(str, [lab, d["family"], d["n_truth"], d["calls"], d["TP"], fp, d["hg005"],
                                  "%.4f" % f1])))
