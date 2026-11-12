"""Leave-one-run-out for the mosaic disc_upper, criterion = mean of per-run pooled-family F1
(so UW_100X's ~6000-call FP mass does not dominate a summed pool). Plateau = values within
0.0005 of the max, contiguous with it; middle rounded to the grid."""
import csv, glob
from collections import defaultdict
R = ["BCM_300X", "UW_100X", "WashU_100X"]
d = defaultdict(dict)
for f in glob.glob("sweeps/*_disc_new*.tsv"):
    for r in csv.DictReader(open(f), delimiter="\t"):
        if r["family"] == "pooled":
            d[r["dataset"]][float(r["value"])] = r
grid = sorted(set.intersection(*[set(d[x]) for x in R]))
def choose(runs):
    sc = [(v, sum(float(d[x][v]["F1"]) for x in runs) / len(runs)) for v in grid]
    i = max(range(len(sc)), key=lambda k: sc[k][1]); b = sc[i][1]; lo = hi = i
    while lo > 0 and sc[lo-1][1] >= b - 0.0005: lo -= 1
    while hi < len(sc)-1 and sc[hi+1][1] >= b - 0.0005: hi += 1
    m = min(grid, key=lambda g: abs(g - (sc[lo][0] + sc[hi][0]) / 2))
    return m, b, sc[lo][0], sc[hi][0]
print("all three:", choose(R))
for x in R: print("  own optimum", x, choose([x]))
for h in R:
    m, b, lo, hi = choose([x for x in R if x != h])
    a, c = d[h][0.52], d[h][m]
    bins = [k for k in a if k.startswith("sens_")]
    print("hold %-10s chosen %.2f (plateau %.2f-%.2f) | held-out at 0.52: TP %s FP %s F1 %s | at %.2f: TP %s FP %s F1 %s | dTP %+d dFP %+d"
          % (h, m, lo, hi, a["TP"], a["FP"], a["F1"], m, c["TP"], c["FP"], c["F1"], int(c["TP"]) - int(a["TP"]), int(c["FP"]) - int(a["FP"])))
    print("     sens by VAF bin 0.52 -> %.2f: %s" % (m, "  ".join("%s %s->%s" % (k[5:], a[k], c[k]) for k in bins)))
