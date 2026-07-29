"""Single-threaded pysam open/fetch microbenchmark.

Purpose: the transduction profiling found port's `open` phase in parse_td_sibling_by_site is
~1.8x master's (0.0928s vs 0.0509s per call) with ZERO code difference between the trees. If
that is environmental (pysam 0.17/py3.8 vs pysam 0.24/py3.10) it must reproduce here, with no
xTea code involved at all.

Measures exactly what the `open` phase does: an AlignmentFile open, a header-references touch
(what BamInfo.is_chrm_contain_chr costs), and a second open. Single-threaded on purpose --
this is a per-call constant, not a scaling measurement.
"""
import sys
import time

import pysam

BAM = "/n/scratch/users/c/cos689/xtea2_testdata/NA12878.wgs.bam"
REF = sys.argv[1] if len(sys.argv) > 1 else None
N = 30

print("python  %s" % sys.version.split()[0])
print("pysam   %s" % pysam.__version__)
print("htslib  %s" % getattr(pysam, "__htslib_version__", "?"))

# warm the page cache for the header so every timed iteration is fair
f = pysam.AlignmentFile(BAM, "rb", reference_filename=REF)
chrms = set(f.references)
f.close()
print("n_refs  %d" % len(chrms))


def timeit(label, fn, n=N):
    ts = []
    for _ in range(n):
        t0 = time.time()
        fn()
        ts.append(time.time() - t0)
    ts.sort()
    print("%-22s mean %.4fs  median %.4fs  min %.4fs  max %.4fs"
          % (label, sum(ts) / len(ts), ts[len(ts) // 2], ts[0], ts[-1]))


def open_close():
    f = pysam.AlignmentFile(BAM, "rb", reference_filename=REF)
    f.close()


def open_refs_close():
    f = pysam.AlignmentFile(BAM, "rb", reference_filename=REF)
    _ = set(f.references)
    f.close()


def open_fetch_close():
    f = pysam.AlignmentFile(BAM, "rb", reference_filename=REF)
    n = 0
    for _a in f.fetch("chr1", 1000000, 1000600):
        n += 1
    f.close()


timeit("open+close", open_close)
timeit("open+references+close", open_refs_close)
timeit("open+fetch600bp+close", open_fetch_close)
