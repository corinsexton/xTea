#!/usr/bin/env python
"""Per-record gate flips between two arms, on the postF input (cns2.txt.high_confident).

usage: gate_flips.py OLD_RUNS NEW_RUNS OUT.tsv
For every record present in both arms (chrom, pos, family): old/new AF_disc, AF_concod,
AF_PE and f_ef_disc, and whether each germline gate (0.075 on f_ef_clip, f_ef_disc,
f_clip_fmap, AF_PE; transduction AF_PE >= 0.08) flips. Prints per-dataset flip counts, the
old/new concordant ratio, and the cutoff on the new scale that keeps the same pass/fail set
(rank-matched: the new AF_PE value at the rank where old AF_PE crosses 0.075).
"""
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from compare_arms import ORDER, load_postf_input  # noqa: E402

old_runs, new_runs, out = sys.argv[1:4]
fo = open(out, "w")
fo.write("dataset\tfamily\tchrom\tpos\tsubtype\told_disc\told_conc\tnew_disc\tnew_conc\told_afpe\tnew_afpe"
         "\told_efdisc\tnew_efdisc\tflips\n")
print("dataset\trecords\tconc_old/new(sum)\tmed_afpe_old\tmed_afpe_new\tflip_afpe\tflip_efdisc\tflip_td08"
      "\tflip_any_two_side\tnew_cut_keep_0.075\tnew_cut_keep_0.08\tlike4like_0.075")
for n in ORDER:
    if not (os.path.isdir(os.path.join(old_runs, n)) and os.path.isdir(os.path.join(new_runs, n))):
        continue
    po, pn = load_postf_input(old_runs, n), load_postf_input(new_runs, n)
    keys = sorted(set(po) & set(pn))
    c = dict(afpe=0, efdisc=0, td=0, two=0)
    oa, na = [], []
    for k in keys:
        a, b = po[k], pn[k]
        fl = []
        if (a["af_pe"] > 0.075) != (b["af_pe"] > 0.075): fl.append("afpe"); c["afpe"] += 1
        if (a["ef_disc"] > 0.075) != (b["ef_disc"] > 0.075): fl.append("efdisc"); c["efdisc"] += 1
        if (a["af_pe"] >= 0.08) != (b["af_pe"] >= 0.08): fl.append("td08"); c["td"] += 1
        if fl and "two_side" in a["sub"]: c["two"] += 1
        oa.append(a["af_pe"]); na.append(b["af_pe"])
        fo.write("\t".join(map(str, [n, k[2], k[0], k[1], a["sub"], a["disc"], a["conc"], b["disc"], b["conc"],
                                     "%.4f" % a["af_pe"], "%.4f" % b["af_pe"], "%.4f" % a["ef_disc"],
                                     "%.4f" % b["ef_disc"], ",".join(fl) or "-"])) + "\n")
    oa, na = np.array(oa), np.array(na)
    so = sum(po[k]["conc"] for k in keys); sn = sum(pn[k]["conc"] for k in keys)

    def keep(cut):
        # the smallest new-scale cutoff t such that {new > t} has as many records as {old > cut}
        m = int((oa > cut).sum())
        s = np.sort(na)[::-1]
        return s[m] if m < len(s) else 0.0
    print("%s\t%d\t%.3f\t%.3f\t%.3f\t%d\t%d\t%d\t%d\t%.4f\t%.4f\t%.4f" % (
        n, len(keys), so / max(sn, 1), np.median(oa), np.median(na), c["afpe"], c["efdisc"], c["td"],
        c["two"], keep(0.075), keep(0.08), 2 * 0.075 / 1.075))
