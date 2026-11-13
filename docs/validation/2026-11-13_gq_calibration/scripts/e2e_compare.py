"""Compare the end-to-end VCFs (worktree code, e2e.sh) with the B15 VCFs and the offline result."""
import sys, pandas as pd
S = "/n/scratch/users/c/cos689/xtea2_review/gq_calib"
B15 = "/n/scratch/users/c/cos689/xtea2_testdata/gtval/runs_b15"
a = pd.read_csv(f"{S}/final/records_after.tsv.gz", sep="\t", keep_default_na=False)
recs = lambda f: [l.rstrip("\n").split("\t") for l in open(f) if not l.startswith("#")]
print("dataset\tn_before\tn_e2e\tcols1_8_same\tother_format_same\tgt_same_as_before\tgt_eq_offline\tgq_eq_offline")
for n in sys.argv[1:] or ["HG002_NovaSeq30x", "HG02257"]:
    o, w = recs(f"{B15}/{n}/{n}/{n}.vcf"), recs(f"{S}/e2e/{n}/{n}.vcf")
    x = a[a.dataset == n].sort_values("idx")
    gt = [r[9].split(":")[0] for r in w]; gq = [int(r[9].split(":")[1]) for r in w]
    print(n, len(o), len(w), all(p[:8] == q[:8] for p, q in zip(o, w)),
          all(p[9].split(":")[2:] == q[9].split(":")[2:] for p, q in zip(o, w)),
          [r[9].split(":")[0] for r in o] == gt,
          all((g == "1/1") == (e == "hom") for g, e in zip(gt, x.gt_after)),
          list(x.gq_after) == gq, sep="\t")
