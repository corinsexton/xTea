#!/usr/bin/env python
"""Score the lowq-filter arms against truth and list calls each arm gains/loses vs a reference arm.

usage: compare.py OUT_PREFIX REF_ARM ARM=DIR [ARM=DIR ...]

Each DIR has the gtval layout DIR/<name>/<name>/<name>.vcf. A dataset missing from an arm's DIR
falls back to REF_ARM's call set (used where the offline replay showed the arm's clip-check drop
set equals the reference's, so its calls are identical by construction); such cells are marked
'=ref' in the metrics table. Truth/flags as gtval/score_final.sh; one-to-one, +-200 bp, same family.
Writes OUT_PREFIX.metrics.tsv (pooled TP/FP/recall/PPV per dataset x arm) and OUT_PREFIX.diff.tsv.
"""
import os
import sys

G = "/n/scratch/users/c/cos689/xtea2_testdata/gtval"
sys.path.insert(0, G)
import score_gt as S  # noqa: E402

H = "/n/data1/hms/dbmi/park/corinne/claude_sandbox/test_data/HG002/truth"
ORDER = ["HG002_NovaSeq30x", "HG002_HiSeq30x", "HG002_ds20x", "HG002_ds15x", "HG002_ds10x",
         "HG00438", "HG02257", "HG02622"]


def vcf(d, n):
    return os.path.join(d, n, n, n + ".vcf")


def main():
    out, ref = sys.argv[1], sys.argv[2]
    arms = dict(a.split("=", 1) for a in sys.argv[3:])
    fm = open(out + ".metrics.tsv", "w")
    fm.write("dataset\tarm\tsource\ttruth\tcalls\tTP\tFP\trecall\tPPV\n")
    fd = open(out + ".diff.tsv", "w")
    fd.write("dataset\tarm\tchange_vs_" + ref + "\tfamily\tchrom\tpos\tstatus\ttruth_pos\ttruth_gt\n")
    for n in ORDER:
        if not os.path.exists(vcf(arms[ref], n)):
            continue
        truth = S.load_truth(H) if n.startswith("HG002") else \
            S.load_truth(os.path.join(G, "truth", "truth_%s.tsv" % n), True, True)
        per = {}
        for a, d in arms.items():
            src = "run"
            p = vcf(d, n)
            if not os.path.exists(p):
                p, src = vcf(arms[ref], n), "=ref"
            calls = S.load_calls([p])
            res = S.score(truth, calls, 200)
            tp = {}
            for fam in S.FAMILIES:
                for t, c in res[fam]["pairs"]:
                    tp[(c.chrom, c.pos, c.family)] = t
            nt = sum(r["n_truth"] for r in res.values())
            nc = len(calls)
            k = len(tp)
            fm.write(f"{n}\t{a}\t{src}\t{nt}\t{nc}\t{k}\t{nc - k}\t{k / nt:.4f}\t{k / nc:.4f}\n")
            per[a] = ({(c.chrom, c.pos, c.family) for c in calls}, tp)
        rk, rtp = per[ref]
        for a in arms:
            if a == ref:
                continue
            ak, atp = per[a]
            for change, ks, tps in (("gained", ak - rk, atp), ("lost", rk - ak, rtp)):
                for key in sorted(ks):
                    t = tps.get(key)
                    fd.write(f"{n}\t{a}\t{change}\t{key[2]}\t{key[0]}\t{key[1]}\t"
                             f"{'TP' if t else 'FP'}\t{t.pos if t else '-'}\t{(t.gt or '-') if t else '-'}\n")
    fm.close()
    fd.close()


if __name__ == "__main__":
    main()
