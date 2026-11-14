#!/usr/bin/env python
"""Per variant: of the candidates passing the clip check with the filter off (sim/<name>.tsv),
how many near a truth site (same family, +-200 bp) vs not does the filter drop.
This is the clip-check stage, before all later filters -- a proxy; the arms give final calls."""
import bisect, sys
from collections import defaultdict
G = "/n/scratch/users/c/cos689/xtea2_testdata/gtval"
sys.path.insert(0, G)
import score_gt as S
H = "/n/data1/hms/dbmi/park/corinne/claude_sandbox/test_data/HG002/truth"
for n in sys.argv[1:]:
    T = H if n.startswith("HG002") else f"{G}/truth/truth_{n}.tsv"
    tr = defaultdict(list)
    for s in S.load_truth(T):
        tr[(s.chrom, s.family)].append(s.pos)
    for k in tr: tr[k].sort()
    def near(c, f, p):
        v = tr.get((c, f), []); i = bisect.bisect_left(v, p - 200)
        return i < len(v) and v[i] <= p + 200
    rows = [l.rstrip("\n").split("\t") for l in open(f"sim/{n}.tsv")]
    h = rows[0]; rows = rows[1:]
    nt = sum(near(r[1], r[0], int(r[2])) for r in rows); nf = len(rows) - nt
    out = [f"{n}: pass-off truth-near {nt}, other {nf}"]
    for i, col in enumerate(h):
        if not col.startswith("drop_"): continue
        dt = sum(1 for r in rows if r[i] == "1" and near(r[1], r[0], int(r[2])))
        df = sum(1 for r in rows if r[i] == "1") - dt
        out.append(f"   {col[5:]:11s} drops truth {dt:4d} ({100*dt/max(nt,1):4.1f}%)  other {df:4d} ({100*df/max(nf,1):4.1f}%)")
    print("\n".join(out))
