#!/usr/bin/env python
"""Independent HG002 GT label from the GIAB Q100 v5.0q GRCh38 stvar benchmark, for every HG002
record of records.tsv.gz. Per haplotype within +-200 bp: 'ins' if it carries an insertion
>= max(50, 0.5*SVLEN); 'del' if a deletion >= 50 bp on that hap overlaps the site; else 'ref'.
q100 = hom (ins|ins), het (ins|ref), hemi (ins|del), absent, other. inbed = site within the
benchmark regions.
usage: q100_labels.py records.tsv.gz q100.vcf.gz bench.bed out.tsv.gz
"""
import gzip, sys, bisect
from collections import defaultdict
rec, vcf, bed, out = sys.argv[1:5]
W = 200
ins = defaultdict(list); dels = defaultdict(list)
for line in gzip.open(vcf, "rt"):
    if line[0] == "#":
        continue
    x = line.split("\t", 10)
    gt = x[9].split(":")[0].replace("/", "|").split("|")
    if len(gt) != 2:
        continue
    pos, ref = int(x[1]), x[3]
    for ai, alt in enumerate(x[4].split(","), 1):
        d = len(alt) - len(ref)
        haps = tuple(int(g == str(ai)) for g in gt)
        if d >= 50:
            ins[x[0]].append((pos, d, haps))
        elif d <= -50:
            dels[x[0]].append((pos, pos + len(ref), haps))
for c in ins:
    ins[c].sort()
B = defaultdict(list)
for line in open(bed):
    c, s, e = line.split()[:3]
    B[c].append((int(s), int(e)))
for c in B:
    B[c].sort()

def inbed(c, p):
    a = B.get(c, [])
    i = bisect.bisect_right(a, (p, 1 << 62)) - 1
    return i >= 0 and a[i][0] <= p < a[i][1]

fo = gzip.open(out, "wt")
fo.write("dataset\tidx\tq100\tq100_hapA\tq100_hapB\tq100_inslen\tinbed\n")
hdr = None
for line in gzip.open(rec, "rt"):
    x = line.rstrip("\n").split("\t")
    if hdr is None:
        hdr = {k: i for i, k in enumerate(x)}; continue
    if not x[0].startswith("HG002"):
        continue
    c, p = x[hdr["chrom"]], int(x[hdr["pos"]])
    try:
        sv = int(x[hdr["SVLEN"]])
    except ValueError:
        sv = 0
    mn = max(50, 0.5 * sv)
    st = ["ref", "ref"]; L = []
    a = ins.get(c, [])
    i = bisect.bisect_left(a, (p - W, -1, ()))
    while i < len(a) and a[i][0] <= p + W:
        if a[i][1] >= mn:
            for h in (0, 1):
                if a[i][2][h]:
                    st[h] = "ins"
            L.append(str(a[i][1]))
        i += 1
    for s, e, haps in dels.get(c, []):
        if s <= p <= e:
            for h in (0, 1):
                if haps[h] and st[h] == "ref":
                    st[h] = "del"
    k = sorted(st)
    lab = {("ins", "ins"): "hom", ("ins", "ref"): "het", ("del", "ins"): "hemi",
           ("ref", "ref"): "absent"}.get(tuple(k), "other")
    fo.write("\t".join([x[0], x[1], lab, st[0], st[1], ",".join(L) or ".",
                        str(int(inbed(c, p)))]) + "\n")
fo.close()
