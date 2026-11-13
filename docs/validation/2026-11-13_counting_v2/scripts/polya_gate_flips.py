#!/usr/bin/env python
"""is_polyA_dominant_one_side (x_post_filter.py, nclip_half_cutoff 5, f_side_polyA_cutoff) evaluated on
the Alu/L1/SVA postF input (cns.txt.high_confident) of two arms; prints per dataset how many
one_side/one_half records flip, and whether each lost call in OUT_PREFIX.diff.tsv is one of them.
usage: polya_gate_flips.py OLD_RUNS NEW_RUNS DIFF_TSV"""
import os, sys
sys.path.insert(0, "/n/scratch/users/c/cos689/xtea2_wt/counting_v2/xtea")
import x_post_filter as XP
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from compare_arms import ORDER
old, new, diff = sys.argv[1:4]



import global_values; fcut = global_values.ONE_SIDE_POLYA_CUTOFF
def gate(r):
    lp, rp, el, er, lr = int(r[9]), int(r[10]), int(r[5]), int(r[6]), int(r[35])
    if el <= 0 and er <= 0: return True
    if el <= 0 and rp / er > fcut and lr - er >= 5: return True
    if er <= 0 and lp / el > fcut and lr - el >= 5: return True
    return False
def load(run, n, fam):
    p = os.path.join(run, n, n, fam, "candidate_disc_filtered_cns.txt.high_confident")
    return {(x[0], int(x[1])): x for x in (l.split("\t") for l in open(p)) if len(x) > 40} if os.path.exists(p) else {}
lost = {}
for l in open(diff):
    x = l.rstrip("\n").split("\t")
    if x[1] == "lost": lost[(x[0], x[2], x[3], int(x[4]))] = x[5]
print("f_side_polyA_cutoff", fcut)
print("dataset\tfamily\tone_side/half records\tgate off->on\ton->off\tlost calls\tlost explained (TP/FP)")
for n in ORDER:
    for fam in ("Alu", "L1", "SVA"):
        a, b = load(old, n, fam), load(new, n, fam)
        k = [s for s in set(a) & set(b) if "one_half" in a[s][32] or "one_side" in a[s][32]]
        on = [s for s in k if not gate(a[s]) and gate(b[s])]; off = [s for s in k if gate(a[s]) and not gate(b[s])]
        L = [s for (d, fm, c, p), st in lost.items() if d == n and fm == fam for s in [(c, p)]]
        ex = [s for s in L if s in on]
        tp = sum(lost[(n, fam, c, p)] == "TP" for c, p in ex)
        print("%s\t%s\t%d\t%d\t%d\t%d\t%d (%d/%d)" % (n, fam, len(k), len(on), len(off), len(L), len(ex), tp, len(ex) - tp))
