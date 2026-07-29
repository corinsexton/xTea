"""Equivalence + fork-safety check for the new per-process BAM handle cache.

Three claims to prove, on the real 39 GB NA12878 WGS BAM:
  1. get_cached_bam_meta returns EXACTLY what the old per-site code computed
     (BamInfo.is_chrm_contain_chr, and XClipDisc._get_chrm_id_name).
  2. A reused handle returns the same reads as a fresh open, across repeated interleaved
     fetches -- i.e. reuse does not corrupt the file offset.
  3. Forked children do not share the parent's handle: each child's reads match a fresh open.
"""
import os
import sys
from multiprocessing import Pool

sys.path.insert(0, "/n/data1/hms/dbmi/park/corinne/claude_sandbox/xTea/xtea")

import pysam
from x_alignments import BamInfo, get_cached_bam, get_cached_bam_meta

BAM = "/n/scratch/users/c/cos689/xtea2_testdata/NA12878.wgs.bam"
REF = None
REGIONS = [("chr1", 1000000, 1002000), ("chr7", 5000000, 5002000), ("chr1", 1000000, 1002000)]


def fresh_fetch(region):
    chrm, s, e = region
    f = pysam.AlignmentFile(BAM, "rb", reference_filename=REF)
    out = [(a.query_name, a.reference_start, a.is_read1) for a in f.fetch(chrm, s, e)]
    f.close()
    return out


def cached_fetch(region):
    chrm, s, e = region
    f = get_cached_bam(BAM, REF)
    return [(a.query_name, a.reference_start, a.is_read1) for a in f.fetch(chrm, s, e)]


def child_task(region):
    return (os.getpid(), cached_fetch(region) == fresh_fetch(region), len(cached_fetch(region)))


# ---- claim 1: meta equivalence vs the old code paths -------------------------------------
b_with_chr_new, m_tid_new = get_cached_bam_meta(BAM, REF)

bam_info = BamInfo(BAM, REF)
b_with_chr_old = bam_info.is_chrm_contain_chr()

f = pysam.AlignmentFile(BAM, "rb", reference_filename=REF)
m_tid_old = {}
for schrm in f.references:
    m_tid_old[f.get_tid(schrm)] = schrm
m_tid_old[-1] = "*"
f.close()

print("claim 1  b_with_chr  old=%s new=%s  MATCH=%s"
      % (b_with_chr_old, b_with_chr_new, b_with_chr_old == b_with_chr_new))
print("claim 1  tid map     n=%d  MATCH=%s" % (len(m_tid_new), m_tid_old == m_tid_new))

# ---- claim 2: reuse does not corrupt the offset ------------------------------------------
ok2 = True
for region in REGIONS:
    a, b = cached_fetch(region), fresh_fetch(region)
    same = (a == b)
    ok2 = ok2 and same
    print("claim 2  %s:%d-%d  n=%d  MATCH=%s" % (region[0], region[1], region[2], len(b), same))

# ---- claim 3: fork safety ----------------------------------------------------------------
# The parent's cache is ALREADY warm at this point (claims 1-2 populated it), which is the
# dangerous case: if the pid were not in the key, children would inherit this handle.
pool = Pool(4)
results = pool.map(child_task, REGIONS * 2, 1)
pool.close()
pool.join()
pids = set(r[0] for r in results)
ok3 = all(r[1] for r in results)
print("claim 3  %d child tasks across %d pids, all match fresh open: %s"
      % (len(results), len(pids), ok3))
print("parent pid was %d (not among children: %s)" % (os.getpid(), os.getpid() not in pids))

print("\nALL PASS" if (b_with_chr_old == b_with_chr_new and m_tid_old == m_tid_new
                       and ok2 and ok3) else "\nFAILURE")
