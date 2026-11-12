#!/usr/bin/env python
"""Put the vafmix mosaic post-filter inputs on the new (pairs-counted-once) scale.

The vafmix runs were made before the bug-2 fix, so fields 36 (n_full_map), 39 (n_disc_pairs) and 40
(n_concd_pairs) of every postF_mosaic input record are on the old scale. At every record, the
features are recounted with the OLD code (ca5b843) and the NEW code (bug2_fix worktree), with the
same insert-size window and read length, and the fix's per-site effect is applied to the
pipeline's own value: new_field = max(0, pipeline_field + (new_recount - old_recount)).
(Applying the delta, not the raw recount, keeps every other difference between the recount and the
pipeline's genotyping pass out of the comparison.) Writes <out>/<name>_<fam>.input_new.txt and a
per-record table <out>/<name>_<fam>.recount.tsv.
"""
import argparse, importlib.util, os, shutil, sys, tempfile
from multiprocessing import Pool

ap = argparse.ArgumentParser()
ap.add_argument("--run"); ap.add_argument("--bam"); ap.add_argument("--ref")
ap.add_argument("--is", dest="i_is", type=int); ap.add_argument("--rlen", type=int)
ap.add_argument("-n", type=int, default=8); ap.add_argument("--out")
a = ap.parse_args()

NEW = "/n/scratch/users/c/cos689/xtea2_wt/bug2_fix/xtea/x_genotype_feature.py"
OLD = "/n/scratch/users/c/cos689/xtea2_wt/rlen_fix/xtea/x_genotype_feature.py"   # ca5b843


def load(name, path):
    sys.path.insert(0, os.path.dirname(path))
    spec = importlib.util.spec_from_file_location(name, path)
    m = importlib.util.module_from_spec(spec); spec.loader.exec_module(m)
    sys.path.pop(0)
    return m


sys.path.insert(0, os.path.dirname(NEW))
import global_values  # noqa: E402
new_mod = load("gf_new", NEW); old_mod = load("gf_old", OLD)
global_values.set_insert_size(a.i_is); global_values.set_read_length(a.rlen)
name = os.path.basename(a.run.rstrip("/"))
os.makedirs(a.out, exist_ok=True)


def count(mod, chrm, pos, wd):
    g = mod.XGenotyper(a.ref, wd, 1)
    g.collect_features_one_site(((chrm, pos, a.i_is), a.bam, wd))
    p = wd + "{0}_{1}".format(chrm, pos) + global_values.GNTP_FEATURE_SUFFIX
    x = open(p).read().split(); os.remove(p)
    return int(x[3]), int(x[6]), int(x[7])          # n_full_map, n_disc_pairs, n_concd_pairs


def work(site):
    wd = tempfile.mkdtemp(dir=a.out) + "/"
    o = count(old_mod, site[0], site[1], wd)
    n = count(new_mod, site[0], site[1], wd)
    shutil.rmtree(wd, ignore_errors=True)
    return site, o, n


import shlex  # noqa: E402
for fam in ("L1", "Alu", "SVA"):
    log = os.path.join(a.run, name, "%s_%s.steps.log" % (name, fam))
    lines = open(log).read().split("\n")
    cmd = shlex.split(lines[lines.index("=== STEP %s:postF_mosaic ===" % fam) + 1])
    src = cmd[cmd.index("-i") + 1]
    recs = [ln.rstrip("\n").split("\t") for ln in open(src) if ln.strip()]
    sites = sorted({(r[0], int(r[1])) for r in recs})
    with Pool(a.n) as p:
        res = {s: (o, n) for s, o, n in p.imap_unordered(work, sites, chunksize=8)}
    with open(os.path.join(a.out, "%s_%s.input_new.txt" % (name, fam)), "w") as fo, \
            open(os.path.join(a.out, "%s_%s.recount.tsv" % (name, fam)), "w") as ft:
        ft.write("chrom\tpos\tpipe_fmap\tpipe_disc\tpipe_conc\told_fmap\told_disc\told_conc"
                 "\tnew_fmap\tnew_disc\tnew_conc\tadj_fmap\tadj_disc\tadj_conc\n")
        for r in recs:
            o, n = res[(r[0], int(r[1]))]
            pipe = [int(r[36]), int(r[39]), int(r[40])]
            adj = [max(0, pipe[i] + n[i] - o[i]) for i in range(3)]
            ft.write("\t".join(map(str, [r[0], r[1]] + pipe + list(o) + list(n) + adj)) + "\n")
            r[36], r[39], r[40] = (str(v) for v in adj)
            fo.write("\t".join(r) + "\n")
    print(fam, "records", len(recs), "sites", len(sites), "->", src, file=sys.stderr, flush=True)
