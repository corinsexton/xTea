#!/usr/bin/env python
"""Move the HapMap-mixture mosaic post-filter inputs from the xTea2 0553edb counts onto the
counting-v2 counts (B15 right-clip map_pos fix + B16 supplementary AF clips).

Base inputs are the B14-scale records (bug2 tune/mosaic_inputs/<name>_<fam>.input_new.txt, i.e.
what 0553edb would feed postF_mosaic). At every record the features are recounted with the
CONTROL code (0553edb worktree) and the NEW code (counting_v2 worktree), with the same insert-size
window and read length, and the per-site change is applied to the base value:
new_field = max(0, base_field + (new_recount - ctrl_recount)) for fields 35 (n_af_clip),
36 (n_full_map), 39 (n_disc_pairs) and 40 (n_concd_pairs). Writes
<out>/<name>_<fam>.input_cv2.txt and a per-record table <out>/<name>_<fam>.recount.tsv.
"""
import argparse, importlib.util, os, shutil, sys, tempfile
from multiprocessing import Pool

ap = argparse.ArgumentParser()
ap.add_argument("--name"); ap.add_argument("--bam"); ap.add_argument("--ref")
ap.add_argument("--is", dest="i_is", type=int); ap.add_argument("--rlen", type=int)
ap.add_argument("-n", type=int, default=8); ap.add_argument("--out")
a = ap.parse_args()

NEW = "/n/scratch/users/c/cos689/xtea2_wt/counting_v2/xtea/x_genotype_feature.py"
OLD = "/n/scratch/users/c/cos689/xtea2_wt/cv2_ctrl/xtea/x_genotype_feature.py"   # 0553edb
BASE = "/n/scratch/users/c/cos689/xtea2_review/bug2/validation/tune/mosaic_inputs"
FIELDS = (35, 36, 39, 40)
FEAT = (2, 3, 6, 7)


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
os.makedirs(a.out, exist_ok=True)


def count(mod, chrm, pos, wd):
    g = mod.XGenotyper(a.ref, wd, 1)
    g.collect_features_one_site(((chrm, pos, a.i_is), a.bam, wd))
    p = wd + "{0}_{1}".format(chrm, pos) + global_values.GNTP_FEATURE_SUFFIX
    x = open(p).read().split(); os.remove(p)
    return tuple(int(x[i]) for i in FEAT)


def work(site):
    wd = tempfile.mkdtemp(dir=a.out) + "/"
    o = count(old_mod, site[0], site[1], wd)
    n = count(new_mod, site[0], site[1], wd)
    shutil.rmtree(wd, ignore_errors=True)
    return site, o, n


for fam in ("L1", "Alu", "SVA"):
    src = os.path.join(BASE, "%s_%s.input_new.txt" % (a.name, fam))
    recs = [ln.rstrip("\n").split("\t") for ln in open(src) if ln.strip()]
    sites = sorted({(r[0], int(r[1])) for r in recs})
    with Pool(a.n) as p:
        res = {s: (o, n) for s, o, n in p.imap_unordered(work, sites, chunksize=8)}
    with open(os.path.join(a.out, "%s_%s.input_cv2.txt" % (a.name, fam)), "w") as fo, \
            open(os.path.join(a.out, "%s_%s.recount.tsv" % (a.name, fam)), "w") as ft:
        ft.write("chrom\tpos\t" + "\t".join("%s_%d" % (k, f) for k in ("base", "ctrl", "new", "adj")
                                            for f in FIELDS) + "\n")
        for r in recs:
            o, n = res[(r[0], int(r[1]))]
            base = [int(r[f]) for f in FIELDS]
            adj = [max(0, base[i] + n[i] - o[i]) for i in range(len(FIELDS))]
            ft.write("\t".join(map(str, [r[0], r[1]] + base + list(o) + list(n) + adj)) + "\n")
            for f, v in zip(FIELDS, adj):
                r[f] = str(v)
            fo.write("\t".join(r) + "\n")
    print(fam, "records", len(recs), "sites", len(sites), "<-", src, file=sys.stderr, flush=True)
