#!/usr/bin/env python
"""Per-record table from the runs_b15 merged VCFs: INFO counts, current GT/GQ, truth GT label
(matched as score_final.sh / gq_calib.py do). All records are kept (the EM fit uses all of them);
truth is NA for unmatched, unlabeled or chrX/Y records.

usage: extract.py RUNS_DIR OUT.tsv.gz
"""
import gzip, os, re, sys
G = "/n/scratch/users/c/cos689/xtea2_testdata/gtval"
sys.path.insert(0, G)
import score_gt as S  # noqa
H = "/n/data1/hms/dbmi/park/corinne/claude_sandbox/test_data/HG002/truth"
ORDER = ["HG002_NovaSeq30x", "HG002_HiSeq30x", "HG002_ds20x", "HG002_ds15x", "HG002_ds10x",
         "HG00438", "HG02257", "HG02622"]
KEYS = ["AF_CLIP", "AF_FMAP", "AF_DISC", "AF_CONCORDNT", "SUBTYPE", "SVLEN", "LCOV", "RCOV",
        "REF_REP", "TD_SRC", "LCLIP", "RCLIP", "LDISC", "RDISC", "LPOLYA", "RPOLYA"]
runs, out = sys.argv[1:3]
with gzip.open(out, "wt") as fo:
    fo.write("\t".join(["dataset", "idx", "chrom", "pos", "family"] + KEYS +
                       ["gt", "gq", "truth"]) + "\n")
    for n in ORDER:
        v = os.path.join(runs, n, n, n + ".vcf")
        truth = S.load_truth(H) if n.startswith("HG002") else \
            S.load_truth(os.path.join(G, "truth", "truth_%s.tsv" % n), True, True)
        calls = S.load_calls([v])
        res = S.score(truth, calls, 200)
        lab = {}
        for fam in S.FAMILIES:
            for t, c in res[fam]["pairs"]:
                lab[(c.chrom, c.pos, c.family)] = t.gt or "NA"
        i = 0
        for line in open(v):
            if line.startswith("#"):
                continue
            x = line.rstrip("\n").split("\t")
            info = dict(S._INFO_RE.findall(x[7]))
            fmt = dict(zip(x[8].split(":"), x[9].split(":")))
            fam = S.VCF_ALT_TO_FAMILY.get(x[4].strip("<>"), "NA")
            fo.write("\t".join([n, str(i), x[0], x[1], fam] + [info.get(k, "NA") for k in KEYS] +
                               [S.norm_call_gt(fmt["GT"]), fmt["GQ"],
                                lab.get((x[0], int(x[1]), fam), "NA")]) + "\n")
            i += 1
