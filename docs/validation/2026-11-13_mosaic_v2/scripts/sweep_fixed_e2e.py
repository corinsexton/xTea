#!/usr/bin/env python
"""(mosaic_v2 copy: drops --germline-hc from the postF_mosaic argv so a vaf-site-v2 run can be
rescored at a fixed ceiling.) Real end-to-end sweep of one post-filter cutoff on an existing run.

For each cutoff value: re-run the per-family post-filter step (postF_hc, or postF_mosaic) with the
tune code copy (tune/code/xtea, identical to the bug2_fix worktree except that the three cutoffs
are read from XTEA2_TUNE_AFPE / XTEA2_TUNE_TD / XTEA2_TUNE_DISC_UPPER), then the real shared tail
(run_xtea.run_shared_tail: locus merge -> gene -> genotype -> gVCF), then score the merged VCF.
Germline scoring is score_gt exactly as score_final.sh does it (1KG: --tier1 --require-agree).
Mosaic scoring: one-to-one match (same family, +-200 bp, score_gt.match_one_to_one) to the
kind=somatic sites of vafsite/sites/mix.tsv; every other call is an FP (calls at listed HG005
germline sites are counted separately too).
If the three per-family post-filter outputs are byte-identical to an earlier value's, the earlier
value's VCF is reused (the tail is deterministic given its input).

usage: sweep.py --run RUN_DIR --name NAME --mode germline|mosaic --param afpe|td|disc
                --values v1,v2,... --out OUT_DIR --tsv OUT.tsv [--afpe X --td X --disc X]
                [--mosaic-inputs FAM=PATH,...] [-n CORES]
"""
import argparse, bisect, hashlib, os, shlex, shutil, subprocess, sys, types
from collections import defaultdict

TUNE = "/n/scratch/users/c/cos689/xtea2_review/counting_v2/mosaic"
CODE = os.path.join(TUNE, "code", "xtea")
G = "/n/scratch/users/c/cos689/xtea2_testdata/gtval"
H = "/n/data1/hms/dbmi/park/corinne/claude_sandbox/test_data/HG002/truth"
MIX = "/n/scratch/users/c/cos689/xtea2_testdata/vafsite/sites/mix.tsv"
sys.path.insert(0, G)
import score_gt as S  # noqa: E402

ap = argparse.ArgumentParser()
ap.add_argument("--run", required=True); ap.add_argument("--name", required=True)
ap.add_argument("--mode", choices=("germline", "mosaic"), required=True)
ap.add_argument("--param", choices=("afpe", "td", "disc"), required=True)
ap.add_argument("--values", required=True)
ap.add_argument("--out", required=True); ap.add_argument("--tsv", required=True)
ap.add_argument("--afpe", default="0.14"); ap.add_argument("--td", default="0.148")
ap.add_argument("--disc", default="0.52")
ap.add_argument("--mosaic-inputs", default="")
ap.add_argument("-n", type=int, default=8)
ap.add_argument("--keep", action="store_true", help="keep every point's dir (default: only VCFs)")
a = ap.parse_args()

root = os.path.join(a.run, a.name)                       # sample_root of the source run
step = "postF_mosaic" if a.mode == "mosaic" else "postF_hc"
order = []
for ln in open(os.path.join(a.run, "driver.log")):
    x = ln.split()
    if len(x) >= 5 and x[-2] == step + ":" and x[-3] in ("L1", "Alu", "SVA"):
        order.append(x[-3])
assert order, "no %s lines in driver.log" % step
override = dict(kv.split("=", 1) for kv in a.mosaic_inputs.split(",") if kv)


def step_argv(fam):
    lines = open(os.path.join(root, "%s_%s.steps.log" % (a.name, fam))).read().split("\n")
    cmd = shlex.split(lines[lines.index("=== STEP %s:%s ===" % (fam, step)) + 1])
    cmd = cmd[2:]                                        # drop python and x_TEA_main.py
    if "--germline-hc" in cmd:                           # vaf-site-v2 option; the tune copy
        i = cmd.index("--germline-hc"); del cmd[i:i + 2]  # takes the ceiling from the env
    return cmd


# shared tail from the tune code copy (same code as the worktree)
sys.path.insert(0, CODE)
os.chdir(CODE)
import run_xtea as RX  # noqa: E402

cfg = {"reference": None, "gene_gff3": None}
for ln in open(os.path.join(root, "%s_xtea2.log" % a.name)):
    if "x_TEA_main.py --gene " in ln:
        c = shlex.split(ln); cfg["gene_gff3"] = c[c.index("-a") + 1]
    if "x_TEA_main.py --gVCF " in ln:
        c = shlex.split(ln); cfg["reference"] = c[c.index("--ref") + 1]
        bam_list = c[c.index("-b") + 1]
assert cfg["reference"] and cfg["gene_gff3"]
opt = types.SimpleNamespace(sample=a.name, cores=a.n, case_control=False, mosaic=(a.mode == "mosaic"))

# truth
if a.mode == "germline":
    if a.name.startswith("HG002"):
        truth = S.load_truth(H)
    else:
        truth = S.load_truth(os.path.join(G, "truth", "truth_%s.tsv" % a.name), True, True)
else:
    truth, hg005 = [], []
    for ln in open(MIX):
        x = ln.rstrip("\n").split("\t")
        if x[0] == "chrom" or x[0] in S.SEX:
            continue
        s = S.Site(x[0], int(x[1]), x[2], None, None, float(x[4]) if x[4] != "NA" else -1.0)
        (truth if x[6] == "somatic" else hg005).append(s)
VAF_BINS = [(0, 0.004, "0.25%"), (0.004, 0.015, "0.5-1.25%"), (0.015, 0.035, "1.5-3.25%"),
            (0.035, 0.065, "4-6.25%"), (0.065, 1, ">=7%")]


def score(vcf):
    calls = S.load_calls([vcf])
    rows = []
    if a.mode == "germline":
        res = S.score(truth, calls, 200)
        S.report(res, a.name, rows)
        out = []
        for r in rows:
            _, fam, nt, nc, tp, nl, ok, hh, th, oo, tho, h_o, o_h = r
            out.append(dict(family=fam, n_truth=nt, calls=nc, TP=tp, gt_lab=nl, gt_ok=ok))
        return out
    out, pool = [], dict(n_truth=0, calls=0, TP=0, hg005=0)
    pool_bins = defaultdict(lambda: [0, 0])
    for fam in S.FAMILIES:
        t = [s for s in truth if s.family == fam]
        c = [s for s in calls if s.family == fam]
        m = S.match_one_to_one(t, c, 200)
        used = {ci for ci, _ in m.values()}
        rest = [c[i] for i in range(len(c)) if i not in used]
        hg = len(S.match_one_to_one([s for s in hg005 if s.family == fam], rest, 200))
        d = dict(family=fam, n_truth=len(t), calls=len(c), TP=len(m), hg005=hg)
        for lo, hi, lab in VAF_BINS:
            idx = [i for i, s in enumerate(t) if lo <= s.extra < hi]
            k = sum(1 for i in idx if i in m)
            d["sens_" + lab] = "%d/%d" % (k, len(idx))
            pool_bins[lab][0] += k; pool_bins[lab][1] += len(idx)
        out.append(d)
        for k in pool:
            pool[k] += d[k]
    pool["family"] = "pooled"
    for lo, hi, lab in VAF_BINS:
        pool["sens_" + lab] = "%d/%d" % tuple(pool_bins[lab])
    out.append(pool)
    return out


os.makedirs(a.out, exist_ok=True)
env = dict(os.environ, XTEA2_TUNE_AFPE=a.afpe, XTEA2_TUNE_TD=a.td, XTEA2_TUNE_DISC_UPPER=a.disc)
var = {"afpe": "XTEA2_TUNE_AFPE", "td": "XTEA2_TUNE_TD", "disc": "XTEA2_TUNE_DISC_UPPER"}[a.param]
seen = {}
fout = open(a.tsv, "w")
hdr = None
for v in a.values.split(","):
    env[var] = v
    pdir = os.path.join(a.out, "v" + v)
    sroot = os.path.join(pdir, a.name)
    if os.path.exists(pdir):
        shutil.rmtree(pdir, ignore_errors=True)
    procs, outs, lfs = [], {}, []
    for fam in order:
        argv = step_argv(fam)
        tdir = os.path.join(sroot, fam); os.makedirs(tdir, exist_ok=True)
        o = os.path.join(tdir, os.path.basename(argv[argv.index("-o") + 1]))
        argv[argv.index("-o") + 1] = o
        if fam in override:
            argv[argv.index("-i") + 1] = override[fam]
        outs[fam] = o
        lf = open(o + ".stdout", "w"); lfs.append(lf)
        procs.append((fam, subprocess.Popen([sys.executable, RX.XTEA_MAIN] + argv, cwd=CODE, env=env,
                                            stdout=lf, stderr=subprocess.STDOUT)))
    for fam, p in procs:
        if p.wait() != 0:
            sys.exit("postF failed for %s %s=%s" % (fam, a.param, v))
    for lf in lfs:
        lf.close()
    h = hashlib.sha256()
    for fam in order:
        h.update(fam.encode()); h.update(open(outs[fam], "rb").read())
    key = h.hexdigest()[:16]
    vcf = os.path.join(sroot, a.name + ".vcf")
    if key in seen:
        reused = seen[key][0]
        res = seen[key][1]
        shutil.rmtree(pdir, ignore_errors=True)
    else:
        reused = "-"
        os.environ.update({k: env[k] for k in ("XTEA2_TUNE_AFPE", "XTEA2_TUNE_TD", "XTEA2_TUNE_DISC_UPPER")})
        with open(os.path.join(sroot, "tail.log"), "w") as lfh:
            RX.run_shared_tail(cfg, opt, RX.Path(sroot), {f: outs[f] for f in order}, bam_list, lfh)
        res = score(vcf)
        seen[key] = (v, res)
        if not a.keep:
            for fam in order:
                shutil.rmtree(os.path.join(sroot, fam), ignore_errors=True)
    for r in res:
        row = dict(dataset=a.name, param=a.param, value=v, postf_hash=key, reused_from=reused)
        row.update(r)
        row["FP"] = row["calls"] - row["TP"]; row["FN"] = row["n_truth"] - row["TP"]
        p = row["TP"] / row["calls"] if row["calls"] else 0.0
        s = row["TP"] / row["n_truth"] if row["n_truth"] else 0.0
        row["sens"] = "%.4f" % s; row["prec"] = "%.4f" % p
        row["F1"] = "%.4f" % (2 * p * s / (p + s) if p + s else 0.0)
        if hdr is None:
            hdr = list(row.keys()); fout.write("\t".join(hdr) + "\n")
        fout.write("\t".join(str(row.get(k, "")) for k in hdr) + "\n")
    fout.flush()
    print("%s %s=%s hash %s reused %s pooled TP %s calls %s" % (
        a.name, a.param, v, key, reused, res[-1]["TP"], res[-1]["calls"]), flush=True)
