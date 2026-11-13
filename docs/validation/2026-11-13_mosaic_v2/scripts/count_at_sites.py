#!/usr/bin/env python
"""Genotype-feature counts (AF_CLIP, AF_FMAP, AF_DISC, AF_CONCORDNT) at given sites, computed by
the pipeline's own XGenotyper.collect_features_one_site with the same globals the cns step sets
(READ_LENGTH = measured read length, DFT_IS = int(mean + 3 sd)).

The insertion position used for counting is the clip peak within +-WIN bp of the given position
(left-clip starts + right-clip ends, primary non-dup reads), as xTea's candidate positions are clip
peaks. No clipped read -> the given position.

Also writes <out>.isize: the measured read length, mean/sd as X_BasicInfo computes them, and an
empirical insert-size histogram (read1, proper pair, primary, TLEN>0) from random windows.

usage: count_at_sites.py --bam B --ref FA --sites sites.tsv --out out.tsv [-n 8]
sites.tsv: header line, first 3 columns chrom pos family; other columns are passed through.
"""
import argparse
import os
import random
import shutil
import sys
import tempfile
from multiprocessing import Pool

XTEA = os.environ.get("XTEA_DIR", "/n/scratch/users/c/cos689/xtea2_wt/vaf_site_v2/xtea")
sys.path.insert(0, XTEA)
import pysam  # noqa: E402
import global_values  # noqa: E402
from x_genotype_feature import XGenotyper  # noqa: E402
from x_basic_info import X_BasicInfo  # noqa: E402

WIN = 25
G = {}


def _end(a):
    e = a.reference_start
    for op, ln in a.cigartuples:
        if op in (0, 2, 3, 7, 8):
            e += ln
    return e


def peak(sam, chrom, pos):
    cnt = {}
    lc, rc = {}, {}
    for a in sam.fetch(chrom, max(0, pos - WIN - 200), pos + WIN + 200):
        if a.is_duplicate or a.is_unmapped or a.is_secondary or a.is_supplementary:
            continue
        ct = a.cigartuples
        if not ct:
            continue
        if ct[0][0] == 4:
            p = a.reference_start
            if abs(p - pos) <= WIN:
                cnt[p] = cnt.get(p, 0) + 1
                lc[p] = lc.get(p, 0) + 1
        if ct[-1][0] == 4:
            p = _end(a)
            if abs(p - pos) <= WIN:
                cnt[p] = cnt.get(p, 0) + 1
                rc[p] = rc.get(p, 0) + 1
    if not cnt:
        return pos, 0, "NA"
    best = max(cnt.items(), key=lambda kv: (kv[1], -abs(kv[0] - pos)))
    # TSD estimate: right-clip peak (left junction) minus left-clip peak (right junction)
    tsd = "NA"
    if lc and rc:
        lp = max(lc.items(), key=lambda kv: (kv[1], -abs(kv[0] - pos)))
        rp = max(rc.items(), key=lambda kv: (kv[1], -abs(kv[0] - pos)))
        if lp[1] >= 2 and rp[1] >= 2:
            tsd = str(rp[0] - lp[0])
    return best[0], best[1], tsd


def init(rlth, dft_is):
    global_values.set_read_length(rlth)
    global_values.set_insert_size(dft_is)


def work(args):
    i, chrom, pos = args
    sam = pysam.AlignmentFile(G["bam"], reference_filename=G["ref"])
    ppos, npk, tsd = peak(sam, chrom, pos)
    sam.close()
    wf = tempfile.mkdtemp(dir=G["tmp"]) + "/"
    xg = XGenotyper(G["ref"], wf, 1)
    xg.collect_features_one_site(((chrom, ppos, G["dft_is"]), G["bam"], wf))
    fn = wf + "%s_%d%s" % (chrom, ppos, global_values.GNTP_FEATURE_SUFFIX)
    with open(fn) as fh:
        f = fh.readline().split()
    shutil.rmtree(wf, ignore_errors=True)
    # chrm pos n_af_clip n_full_map n_l_raw n_r_raw n_disc n_concd n_disc_indel clip_lens polyA nchr
    return i, ppos, npk, tsd, int(f[2]), int(f[3]), int(f[6]), int(f[7]), f[9]


def isize_hist(bam, ref, n_win=400, seed=1):
    sam = pysam.AlignmentFile(bam, reference_filename=ref)
    rng = random.Random(seed)
    chroms = [(c, l) for c, l in zip(sam.references, sam.lengths)
              if c.replace("chr", "").isdigit()]
    tot = sum(l for _, l in chroms)
    vals = []
    for _ in range(n_win):
        x = rng.randrange(tot)
        for c, l in chroms:
            if x < l:
                break
            x -= l
        for a in sam.fetch(c, x, x + 20000):
            if (a.is_read1 and a.is_proper_pair and not a.is_secondary and not a.is_supplementary
                    and not a.is_duplicate and a.template_length > 0):
                vals.append(a.template_length)
            if len(vals) > 400000:
                break
    sam.close()
    return vals


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--bam", required=True)
    ap.add_argument("--ref", required=True)
    ap.add_argument("--sites", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("-n", type=int, default=8)
    ap.add_argument("--basic", help="existing basic_cov_is_rlth_info.txt (skip re-estimation)")
    ap.add_argument("--isize-from", help="existing <out>.isize: reuse its header and histogram")
    a = ap.parse_args()
    tmp = tempfile.mkdtemp(prefix="cnt_", dir=os.environ.get("TMPDIR", "/tmp"))
    old_isize = None
    if a.isize_from:
        with open(a.isize_from) as fh:
            old_isize = fh.read()
        hd = dict(kv.split("=") for kv in old_isize.splitlines()[0][1:].split())
        cov, rlth, mean_is, sd = float(hd["cov"]), float(hd["rlth"]), float(hd["mean_is"]), float(hd["sd"])
        if float(rlth).is_integer():
            rlth = int(float(rlth))
    elif a.basic:
        with open(a.basic) as fh:
            f = fh.readline().split("\t")
        cov, rlth, mean_is, sd = float(f[1]), float(f[2]), float(f[3]), float(f[4])
    else:
        bl = tmp + "/bam_list.txt"
        with open(bl, "w") as fh:
            fh.write(a.bam + "\tillumina\n")
        xb = X_BasicInfo(tmp + "/", a.n, a.ref)
        cov, rlth, mean_is, sd = xb.get_cov_is_rlth(bl, a.ref, 500, True)
    dft_is = int(mean_is + 3 * sd)
    sys.stderr.write("cov %.1f rlth %s mean_is %.1f sd %.1f DFT_IS %d\n" % (cov, rlth, mean_is, sd, dft_is))
    G.update(bam=a.bam, ref=a.ref, dft_is=dft_is, tmp=tmp)
    init(rlth, dft_is)
    rows = []
    with open(a.sites) as fh:
        hdr = fh.readline().rstrip("\n").split("\t")
        for ln in fh:
            if ln.strip():
                rows.append(ln.rstrip("\n").split("\t"))
    jobs = [(i, r[0], int(r[1])) for i, r in enumerate(rows)]
    with Pool(a.n, initializer=init, initargs=(rlth, dft_is)) as pool:
        res = sorted(pool.map(work, jobs, chunksize=4))
    with open(a.out, "w") as fo:
        fo.write("\t".join(hdr + ["cnt_pos", "n_peak", "tsd_est", "clip", "fmap", "disc", "conc", "clip_lens"]) + "\n")
        for (i, ppos, npk, tsd, c, fm, d, co, cl) in res:
            fo.write("\t".join(rows[i] + [str(ppos), str(npk), tsd, str(c), str(fm), str(d), str(co), cl]) + "\n")
    if old_isize is not None:
        with open(a.out + ".isize", "w") as fo:
            fo.write(old_isize)
        shutil.rmtree(tmp, ignore_errors=True)
        return
    vals = isize_hist(a.bam, a.ref)
    with open(a.out + ".isize", "w") as fo:
        fo.write("#cov=%.3f rlth=%s mean_is=%.3f sd=%.3f dft_is=%d\n" % (cov, rlth, mean_is, sd, dft_is))
        h = {}
        for v in vals:
            h[v] = h.get(v, 0) + 1
        for k in sorted(h):
            fo.write("%d\t%d\n" % (k, h[k]))
    shutil.rmtree(tmp, ignore_errors=True)


if __name__ == "__main__":
    main()
