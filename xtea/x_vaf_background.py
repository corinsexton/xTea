"""Per-sample background for mosaic VAF, measured at random control positions.

Background (b_SR, b_PE in x_vaf_site.p_alt) is the alt-supporting read rate at a site that has no
insertion: clipped reads that are not from the insertion and discordant pairs that are not from it.
It is library-specific. On the SMaHT HapMap mixture it is about 2.5x higher in the UW library than in
BCM, so one constant misestimates low-VAF calls on some libraries by about 2x.

This module counts the genotyping features (x_genotype_feature.XGenotyper, as the pipeline does
at its calls) at N random autosomal positions of the sample's own BAM. Each position is first moved
to its strongest clip peak within +-25 bp, as xTea's candidate positions are clip peaks. It then
fits b by maximum likelihood at v = 0 (x_vaf_site.fit_background) and multiplies by
x_vaf_site.CONTEXT_FACTOR: MEI call sites (in repeats, at clip peaks) carry more background than
random positions do. That factor was fitted once on the 300x BCM mixture BAM. Held-out check
(truth sites at 0.25-0.5% VAF):

                predicted (random x factor)   measured at MEI sites
    UW_100X     b_SR 0.057  b_PE 0.079        b_SR 0.045  b_PE 0.085
    WashU_100X  b_SR 0.019  b_PE 0.035        b_SR 0.023  b_PE 0.036

Cost: about 100 s for 2000 positions on 2 cores (100x BAM); the default is 1000 positions.

usage: x_vaf_background.py --bam B --ref FA --basic-info F [-n N] [--cores C] [--out F]
prints (and writes to --out) "b_SR,b_PE", the format --vaf-background takes.
"""
import argparse
import os
import random
import shutil
import sys
import tempfile
from multiprocessing import Pool

import numpy as np

try:
    from . import x_vaf_site as xvs
except ImportError:
    import x_vaf_site as xvs

N_CONTROL = 1000
PEAK_WIN = 25
_G = {}


def _ref_end(a):
    e = a.reference_start
    for op, ln in a.cigartuples:
        if op in (0, 2, 3, 7, 8):
            e += ln
    return e


def clip_peak(sam, chrom, pos, win=PEAK_WIN):
    """The position within +-win of pos with the most clipped read ends (left-clip starts and
    right-clip ends of primary, non-duplicate reads); pos itself when none is clipped."""
    cnt = {}
    for a in sam.fetch(chrom, max(0, pos - win - 200), pos + win + 200):
        if a.is_duplicate or a.is_unmapped or a.is_secondary or a.is_supplementary:
            continue
        ct = a.cigartuples
        if not ct:
            continue
        for clipped, p in ((ct[0][0] == 4, a.reference_start), (ct[-1][0] == 4, None)):
            if not clipped:
                continue
            p = _ref_end(a) if p is None else p
            if abs(p - pos) <= win:
                cnt[p] = cnt.get(p, 0) + 1
    if not cnt:
        return pos
    return max(cnt.items(), key=lambda kv: (kv[1], -abs(kv[0] - pos)))[0]


def random_positions(sam, n, seed=1):
    """n uniform positions on the autosomes (chr1-22 or 1-22), sorted."""
    chroms = [(c, l) for c, l in zip(sam.references, sam.lengths)
              if c.replace("chr", "").isdigit() and 1 <= int(c.replace("chr", "")) <= 22]
    tot = sum(l for _, l in chroms)
    rng = random.Random(seed)
    out = []
    for _ in range(n):
        x = rng.randrange(tot)
        for c, l in chroms:
            if x < l:
                break
            x -= l
        out.append((c, x))
    return sorted(out)


def _init(rlth, dft_is):
    import global_values
    global_values.set_read_length(rlth)
    global_values.set_insert_size(dft_is)


def _count_one(args):
    import pysam
    import global_values
    from x_genotype_feature import XGenotyper
    chrom, pos = args
    sam = pysam.AlignmentFile(_G["bam"], reference_filename=_G["ref"])
    ppos = clip_peak(sam, chrom, pos)
    sam.close()
    wf = tempfile.mkdtemp(dir=_G["tmp"]) + "/"
    try:
        XGenotyper(_G["ref"], wf, 1).collect_features_one_site(((chrom, ppos, _G["dft_is"]),
                                                                 _G["bam"], wf))
        with open(wf + "%s_%d%s" % (chrom, ppos, global_values.GNTP_FEATURE_SUFFIX)) as fh:
            f = fh.readline().split()
        # chrm pos n_af_clip n_full_map n_l_raw n_r_raw n_disc n_concd ...
        return [int(f[2]), int(f[3]), int(f[6]), int(f[7])]
    except (OSError, IndexError, ValueError):
        return [0, 0, 0, 0]
    finally:
        shutil.rmtree(wf, ignore_errors=True)


def control_counts(bam, ref, lib, n=N_CONTROL, cores=1, seed=1, tmp_dir=None):
    """(n, 4) AF_CLIP, AF_FMAP, AF_DISC, AF_CONCORDNT at n random autosomal clip peaks."""
    import pysam
    sam = pysam.AlignmentFile(bam, reference_filename=ref)
    pos = random_positions(sam, n, seed)
    sam.close()
    tmp = tempfile.mkdtemp(prefix="vafbg_", dir=tmp_dir)
    rlth, dft_is = int(round(lib.R)), int(lib.dft_is)
    _G.update(bam=bam, ref=ref, dft_is=dft_is, tmp=tmp)
    _init(rlth, dft_is)
    try:
        if cores > 1:
            with Pool(cores, initializer=_init, initargs=(rlth, dft_is)) as pool:
                rows = pool.map(_count_one, pos, chunksize=8)
        else:
            rows = [_count_one(p) for p in pos]
    finally:
        shutil.rmtree(tmp, ignore_errors=True)
    return np.array(rows, float).reshape(-1, 4)


def first_bam(bam_or_list):
    """A BAM/CRAM path, or the first entry of an xTea bam list (path [tab tech])."""
    if bam_or_list.endswith((".bam", ".cram")):
        return bam_or_list
    with open(bam_or_list) as fh:
        for ln in fh:
            if ln.strip():
                return ln.split()[0]
    raise ValueError("empty bam list: %s" % bam_or_list)


def sample_background(bam_or_list, ref, lib, n=N_CONTROL, cores=1, seed=1, tmp_dir=None):
    """Background odds for this sample's MEI calls: dict b_SR, b_PE (and n_control)."""
    c = control_counts(first_bam(bam_or_list), ref, lib, n=n, cores=cores, seed=seed,
                       tmp_dir=tmp_dir)
    b = xvs.control_background(c, lib)
    b["n_control"] = int(len(c))
    return b


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--bam", required=True, help="BAM/CRAM or an xTea bam list (first BAM used)")
    ap.add_argument("--ref", required=True)
    ap.add_argument("--basic-info", required=True, help="basic_cov_is_rlth_info.txt")
    ap.add_argument("-n", type=int, default=N_CONTROL, help="control positions (default %(default)s)")
    ap.add_argument("--cores", type=int, default=1)
    ap.add_argument("--seed", type=int, default=1)
    ap.add_argument("--out", default=None)
    a = ap.parse_args(argv)
    lib = xvs.read_basic_info(a.basic_info)
    if lib is None:
        sys.stderr.write("x_vaf_background: unusable --basic-info %s\n" % a.basic_info)
        return 1
    b = sample_background(a.bam, a.ref, lib, n=a.n, cores=a.cores, seed=a.seed,
                          tmp_dir=os.path.dirname(os.path.abspath(a.out)) if a.out else None)
    s = "%.6g,%.6g" % (b["b_SR"], b["b_PE"])
    if a.out:
        with open(a.out, "w") as fo:
            fo.write(s + "\n")
    print(s)
    return 0


if __name__ == "__main__":
    sys.exit(main())
