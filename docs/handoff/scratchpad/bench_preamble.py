"""Old vs new per-site preamble cost in collect_clipped_disc_reads_by_region.

OLD, per candidate site: BamInfo(...).is_chrm_contain_chr()  [1 BAM open + header walk]
                       + pysam.AlignmentFile(...)            [2nd BAM open]
                       + _get_chrm_id_name(samfile)          [get_tid over all 3,366 contigs]
NEW, per candidate site: two dict lookups against a per-process cache.

Single-threaded on purpose: this is a per-call constant, not a scaling measurement.
"""
import sys
import time

sys.path.insert(0, "/n/data1/hms/dbmi/park/corinne/claude_sandbox/xTea/xtea")

import pysam
from x_alignments import BamInfo, get_cached_bam, get_cached_bam_meta

BAM = "/n/scratch/users/c/cos689/xtea2_testdata/NA12878.wgs.bam"
REF = None
N = 20
N_SITES_L1 = 53953

print("python %s   pysam %s" % (sys.version.split()[0], pysam.__version__))


def old_preamble():
    bam_info = BamInfo(BAM, REF)
    b_with_chr = bam_info.is_chrm_contain_chr()
    samfile = pysam.AlignmentFile(BAM, "rb", reference_filename=REF)
    m_chrm = {}
    for schrm in samfile.references:
        m_chrm[samfile.get_tid(schrm)] = schrm
    m_chrm[-1] = "*"
    samfile.close()
    return b_with_chr, m_chrm


def new_preamble():
    b_with_chr, m_chrm = get_cached_bam_meta(BAM, REF)
    samfile = get_cached_bam(BAM, REF)
    return b_with_chr, m_chrm


def timeit(label, fn, n=N):
    fn()  # warm
    ts = []
    for _ in range(n):
        t0 = time.time()
        fn()
        ts.append(time.time() - t0)
    ts.sort()
    mean = sum(ts) / len(ts)
    print("%-14s mean %.5fs  median %.5fs   x %d L1 sites = %.0fs of worker CPU"
          % (label, mean, ts[len(ts) // 2], N_SITES_L1, mean * N_SITES_L1))
    return mean


t_old = timeit("OLD preamble", old_preamble)
t_new = timeit("NEW preamble", new_preamble)
print("saved per site: %.5fs   -> %.0fs of worker CPU over %d L1 sites (%.1f min at -n 8)"
      % (t_old - t_new, (t_old - t_new) * N_SITES_L1, N_SITES_L1,
         (t_old - t_new) * N_SITES_L1 / 8 / 60))
