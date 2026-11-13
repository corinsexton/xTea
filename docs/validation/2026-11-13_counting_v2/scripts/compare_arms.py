#!/usr/bin/env python
"""Control vs new arm comparison (counting-v2; copied from the bug2 review, gates updated to the
2026-11-12 cutoffs: 0.075 on f_ef_clip/f_ef_disc/f_clip_fmap, AF_PE > 0.14, transduction AF_PE >= 0.15).

usage: compare_arms.py OLD_RUNS NEW_RUNS OUT_PREFIX [--postf-cut 0.075]

For each dataset present in both runs dirs (gtval layout <runs>/<name>/<name>/<name>.vcf):
  * scores both merged VCFs with score_gt (same truth/flags as score_final.sh);
  * lists calls present in only one arm (chrom, pos, family exact) with TP/FP status in that
    arm's own one-to-one matching;
  * attaches, per site, both arms' postF-input record fields (cns2.txt.high_confident): the
    subtype (field 32), AF_disc/AF_concod (39/40), and the four is_qualified_rcd ratios, and
    names the ratio(s) whose pass/fail differs between the arms.
Writes OUT_PREFIX.metrics.tsv and OUT_PREFIX.diff.tsv.
"""
import os
import sys
from collections import defaultdict

G = "/n/scratch/users/c/cos689/xtea2_testdata/gtval"
sys.path.insert(0, G)
import score_gt as S  # noqa: E402

H = "/n/data1/hms/dbmi/park/corinne/claude_sandbox/test_data/HG002/truth"
CUT = 0.075
CUTS = dict(ef_clip=0.075, ef_disc=0.075, clip_fmap=0.075, af_pe=0.14)
TD_CUT = 0.15
ORDER = ["HG00438", "HG02257", "HG02622", "HG002_HiSeq30x", "HG002_ds10x", "HG002_ds15x",
         "HG002_ds20x", "HG002_NovaSeq30x"]


def ratios(f):
    n_ef_clip = int(f[5]) + int(f[6]); n_ef_disc = int(f[7]) + int(f[8])
    n_clip, n_fmap, n_disc, n_conc = int(f[35]), int(f[36]), int(f[39]), int(f[40])
    return dict(ef_clip=n_ef_clip / n_clip if n_clip else 0.0,
                ef_disc=n_ef_disc / n_disc if n_disc else 0.0,
                clip_fmap=n_clip / (n_clip + n_fmap) if n_clip + n_fmap else 0.0,
                af_pe=n_disc / (n_disc + n_conc) if n_disc + n_conc else 0.0,
                disc=n_disc, conc=n_conc, sub=f[32])


def load_postf_input(run, name):
    m = {}
    for fam in S.FAMILIES:
        p = os.path.join(run, name, name, fam, "candidate_disc_filtered_cns2.txt.high_confident")
        if not os.path.exists(p):
            continue
        for ln in open(p):
            f = ln.split()
            if len(f) > 40:
                m[(f[0], int(f[1]), fam)] = ratios(f)
    return m


def load_features(run, name):
    """chrom,pos,fam -> (disc, concd) from every candidate's gntp.features (pre-filter)."""
    m = {}
    for fam in S.FAMILIES:
        p = os.path.join(run, name, name, fam, "candidate_disc_filtered_cns.txt.gntp.features")
        if not os.path.exists(p):
            continue
        for ln in open(p):
            x = ln.split()
            m[(x[0], int(x[1]), fam)] = (int(x[6]), int(x[7]))
    return m


def arm(run, name, truth):
    v = os.path.join(run, name, name, name + ".vcf")
    calls = S.load_calls([v])
    res = S.score(truth, calls, 200)
    rows = []
    S.report(res, name, rows)
    tp_keys = set()
    for fam in S.FAMILIES:
        for t, c in res[fam]["pairs"]:
            tp_keys.add((c.chrom, c.pos, c.family))
    keys = {(c.chrom, c.pos, c.family): c for c in calls}
    return rows, keys, tp_keys


def main():
    old_runs, new_runs, out = sys.argv[1:4]
    names = [n for n in ORDER if os.path.exists(os.path.join(old_runs, n, n, n + ".vcf"))
             and os.path.exists(os.path.join(new_runs, n, n, n + ".vcf"))]
    fm = open(out + ".metrics.tsv", "w")
    fm.write("dataset\tarm\tfamily\tn_truth\tcalls\tTP\tFP\tFN\tsens\tprec\tgt_lab\tgt_ok\tgt_acc"
             "\thet_het\thet_hom\thom_het\thom_hom\n")
    fd = open(out + ".diff.tsv", "w")
    fd.write("dataset\tchange\tfamily\tchrom\tpos\tstatus\tgt\tsubtype\told_disc\told_conc\told_afpe"
             "\tnew_disc\tnew_conc\tnew_afpe\told_efdisc\tnew_efdisc\tflipped_gates\n")
    for n in names:
        if n.startswith("HG002"):
            truth = S.load_truth(H)
        else:
            truth = S.load_truth(os.path.join(G, "truth", "truth_%s.tsv" % n), True, True)
        ro, ko, to = arm(old_runs, n, truth)
        rn, kn, tn = arm(new_runs, n, truth)
        for lab, rows in (("old", ro), ("new", rn)):
            for r in rows:
                _, fam, nt, nc, tp, nl, ok, hh, th, oo, tho, h_o, o_h = r
                fm.write("\t".join(map(str, [n, lab, fam, nt, nc, tp, nc - tp, nt - tp,
                                             "%.4f" % (tp / nt if nt else 0), "%.4f" % (tp / nc if nc else 0),
                                             nl, ok, "%.4f" % (ok / nl if nl else 0),
                                             hh, h_o, o_h, oo])) + "\n")
        po, pn = load_postf_input(old_runs, n), load_postf_input(new_runs, n)
        for change, keys, other, tps, src in (("gained", kn, ko, tn, kn), ("lost", ko, kn, to, ko)):
            for k in sorted(set(keys) - set(other)):
                a, b = po.get(k), pn.get(k)
                flipped = []
                if a and b:
                    for g in ("ef_clip", "ef_disc", "clip_fmap", "af_pe"):
                        if (a[g] > CUTS[g]) != (b[g] > CUTS[g]):
                            flipped.append(g)
                    if (a["af_pe"] >= TD_CUT) != (b["af_pe"] >= TD_CUT):
                        flipped.append("td_afpe")
                fd.write("\t".join(map(str, [
                    n, change, k[2], k[0], k[1], "TP" if k in tps else "FP", src[k].gt,
                    (a or b or {}).get("sub", "NA"),
                    a["disc"] if a else "NA", a["conc"] if a else "NA", "%.3f" % a["af_pe"] if a else "NA",
                    b["disc"] if b else "NA", b["conc"] if b else "NA", "%.3f" % b["af_pe"] if b else "NA",
                    "%.3f" % a["ef_disc"] if a else "NA", "%.3f" % b["ef_disc"] if b else "NA",
                    ",".join(flipped) or "-"])) + "\n")
    fm.close(); fd.close()


if __name__ == "__main__":
    main()
