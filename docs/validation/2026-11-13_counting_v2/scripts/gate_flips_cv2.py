#!/usr/bin/env python
"""Per-record germline gate flips between two arms on the postF input (cns2.txt.high_confident).

usage: gate_flips_cv2.py OLD_RUNS NEW_RUNS OUT.tsv
For every record present in both arms (chrom, pos, family) evaluates each is_qualified_rcd gate at
the current cutoffs (f_ef_clip, f_ef_disc, f_clip_fmap > 0.075; AF_PE > 0.14) and the transduction
gate (AF_PE >= 0.15), and records which flip. Prints per-dataset flip counts (all records, and
two_side records, the only branch where is_qualified_rcd is applied), plus how many records pass
all four gates in each arm.
"""
import os, sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from compare_arms import ORDER, load_postf_input, CUTS, TD_CUT  # noqa: E402

G = ("ef_clip", "ef_disc", "clip_fmap", "af_pe")
old_runs, new_runs, out = sys.argv[1:4]
fo = open(out, "w")
fo.write("dataset\tfamily\tchrom\tpos\tsubtype\t" + "\t".join("old_%s\tnew_%s" % (g, g) for g in G)
         + "\tflips\n")
print("dataset\trecords\t" + "\t".join("flip_" + g for g in G) + "\tflip_td\ttwo_side_flip_any\t"
      "two_side_pass_old\ttwo_side_pass_new\tafpe_up_only_flips\tafpe_down_only_flips")
for n in ORDER:
    if not (os.path.isdir(os.path.join(old_runs, n)) and os.path.isdir(os.path.join(new_runs, n))):
        continue
    po, pn = load_postf_input(old_runs, n), load_postf_input(new_runs, n)
    keys = sorted(set(po) & set(pn))
    c = {g: 0 for g in G}; c.update(td=0, two=0, po=0, pnw=0, up=0, dn=0)
    for k in keys:
        a, b = po[k], pn[k]
        fl = [g for g in G if (a[g] > CUTS[g]) != (b[g] > CUTS[g])]
        for g in fl: c[g] += 1
        if (a["af_pe"] >= TD_CUT) != (b["af_pe"] >= TD_CUT): fl.append("td"); c["td"] += 1
        if "af_pe" in fl:
            c["up" if b["af_pe"] > a["af_pe"] else "dn"] += 1
        if "two_side" in a["sub"]:
            if fl: c["two"] += 1
            c["po"] += all(a[g] > CUTS[g] for g in G); c["pnw"] += all(b[g] > CUTS[g] for g in G)
        fo.write("\t".join(map(str, [n, k[2], k[0], k[1], a["sub"]] + sum(
            (["%.4f" % a[g], "%.4f" % b[g]] for g in G), []) + [",".join(fl) or "-"])) + "\n")
    print("\t".join(map(str, [n, len(keys)] + [c[g] for g in G] + [c["td"], c["two"], c["po"], c["pnw"], c["up"], c["dn"]])))
