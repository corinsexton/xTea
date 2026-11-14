#!/usr/bin/env python
"""Breakdown of the low-MAPQ clipped reads counted near candidates (census/<name>.tsv.gz)."""
import gzip, sys
from collections import Counter
for path in sys.argv[1:]:
    n = path.split("/")[-1].replace(".tsv.gz", "")
    c = Counter(); tot = 0; diff = Counter()
    for i, ln in enumerate(gzip.open(path, "rt")):
        if i == 0: continue
        x = ln.rstrip("\n").split("\t")
        flag, mapq, side = int(x[3]), int(x[4]), x[5]
        ie, ix = x[10] == "1", x[11] == "1"
        kind = "supp" if flag & 0x800 else ("sec" if flag & 0x100 else "prim")
        xa = x[13]
        if ie != ix: diff[(side, "est_only" if ie else "exact_only")] += 1
        if not ie: continue
        tot += 1
        c[("kind", kind)] += 1; c[("side", side)] += 1
        c[("mapq", min(mapq, 4))] += 1
        xac = "noXA" if xa == "-" else ("XA_primary_only" if xa == "primary" else ("XA_nonprimary_only" if "primary" not in xa else "XA_mixed"))
        c[("XA", xac)] += 1
        c[("kind_XA", kind + ":" + xac)] += 1
    print(f"== {n}: {tot} (read,candidate) pairs counted (est)")
    for k in sorted(c): print(f"   {k[0]:8s} {str(k[1]):24s} {c[k]:7d} {100*c[k]/tot:5.1f}%")
    print("   est vs exact disagreements:", dict(diff))
