#!/usr/bin/env python
"""GQ calibration, old vs new arm: per dataset (and pooled over datasets), GT error rate within
GQ bins on truth-matched, GT-labeled calls, next to the error rate GQ predicts (mean 10^-GQ/10).

usage: gq_calib.py OLD_RUNS NEW_RUNS OUT.tsv
"""
import os
import sys
from collections import defaultdict

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from compare_arms import ORDER, G, H, S  # noqa: E402

BINS = [(0, 10), (10, 20), (20, 30), (30, 50), (50, 60), (60, 1e9)]
old_runs, new_runs, out = sys.argv[1:4]
acc = defaultdict(lambda: [0, 0, 0.0])  # (dataset|ALL, arm, bin) -> n, errors, sum predicted err
for n in ORDER:
    if not all(os.path.exists(os.path.join(r, n, n, n + ".vcf")) for r in (old_runs, new_runs)):
        continue
    truth = S.load_truth(H) if n.startswith("HG002") else \
        S.load_truth(os.path.join(G, "truth", "truth_%s.tsv" % n), True, True)
    for arm, r in (("old", old_runs), ("new", new_runs)):
        calls = S.load_calls([os.path.join(r, n, n, n + ".vcf")])
        res = S.score(truth, calls, 200)
        for fam in S.FAMILIES:
            for t, c in res[fam]["pairs"]:
                if t.gt not in ("het", "hom") or c.gq is None:
                    continue
                for lo, hi in BINS:
                    if lo <= c.gq < hi:
                        b = "%g-%s" % (lo, "inf" if hi > 1e8 else "%g" % hi)
                        for ds in (n, "ALL"):
                            v = acc[(ds, arm, b)]
                            v[0] += 1; v[1] += int(t.gt != c.gt); v[2] += 10 ** (-c.gq / 10.0)
with open(out, "w") as fo:
    fo.write("dataset\tarm\tgq_bin\tn\terrors\tobs_err\tpred_err\n")
    for (ds, arm, b), (k, e, p) in sorted(acc.items(), key=lambda x: (x[0][0], x[0][2], x[0][1])):
        fo.write("%s\t%s\t%s\t%d\t%d\t%.4f\t%.4f\n" % (ds, arm, b, k, e, e / k, p / k))
