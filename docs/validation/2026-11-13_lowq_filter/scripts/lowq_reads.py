#!/usr/bin/env python
"""Offline census of the low-MAPQ clipped reads the cns lowq filter counts.

usage: lowq_reads.py RUN_SAMPLE_DIR BAM REF READ_LENGTH OUT.tsv.gz

RUN_SAMPLE_DIR is <runs>/<name>/<name> (holds L1/Alu/SVA/candidate_list_from_disc.txt and
tmp/cns/low_mapq_clip_cnt.txt). For every candidate (union over families) it lists each
MAPQ<5 read with a soft clip at the first (left) or, failing that, last (right) CIGAR op --
the reads cns_emit.hpp considers -- with both position rules:
  est   = start (left) | start + READ_LENGTH - clip_len (right)   [master / current xTea2]
  exact = start (left) | reference_end (bam_endpos)               [suspect 1]
and the read's flag, MAPQ, XA/SA contigs, so the alt-contig question (suspect 2) can be asked.
It also writes OUT.counts.tsv: per family candidate, the scanner's count and our recount under
est, exact, and the alt-aware variants.
"""
import bisect
import gzip
import os
import sys
from collections import defaultdict

import pysam

NEAR = 35
LOWQ = 5
GATE = 12


def contig_class(c):
    if c.endswith("_alt"):
        return "alt"
    if c.startswith("HLA"):
        return "hla"
    if c.endswith("_decoy") or c == "chrEBV" or c == "hs37d5":
        return "decoy"
    if "_" in c:          # chrUn / random
        return "unplaced"
    return "primary"


def tag_contigs(r, tag):
    if not r.has_tag(tag):
        return []
    out = []
    for h in r.get_tag(tag).rstrip(";").split(";"):
        if h:
            out.append(h.split(",")[0])
    return out


def main():
    sdir, bam, ref, rl, out = sys.argv[1:6]
    rl = int(rl)
    fams = [f for f in ("L1", "Alu", "SVA") if os.path.exists(f"{sdir}/{f}/candidate_list_from_disc.txt")]
    cands = defaultdict(set)
    fam_cands = {}
    for f in fams:
        s = []
        for ln in open(f"{sdir}/{f}/candidate_list_from_disc.txt"):
            x = ln.split()
            cands[x[0]].add(int(x[1]))
            s.append((x[0], int(x[1])))
        fam_cands[f] = s
    cands = {c: sorted(v) for c, v in cands.items()}

    def near(c, p):
        v = cands.get(c)
        if not v:
            return []
        lo = bisect.bisect_left(v, p - NEAR + 1)
        hi = bisect.bisect_left(v, p + NEAR)
        return v[lo:hi]

    af = pysam.AlignmentFile(bam, reference_filename=ref)
    alt_ref = any(contig_class(c) == "alt" for c in af.references)
    fo = gzip.open(out, "wt")
    fo.write("chrom\tcand\tqname\tflag\tmapq\tside\tclip\tcigar\test_pos\texact_pos\tin_est\tin_exact"
             "\tXA_n\tXA_classes\tSA_classes\n")
    cnt = defaultdict(lambda: defaultdict(int))  # (c,C) -> variant -> n
    margin = NEAR + 2000
    for c, v in cands.items():
        # merge windows
        wins = []
        for C in v:
            b, e = max(0, C - margin), C + NEAR + 1
            if wins and b <= wins[-1][1]:
                wins[-1][1] = max(wins[-1][1], e)
            else:
                wins.append([b, e])
        for b, e in wins:
            for r in af.fetch(c, b, e):
                if r.flag & (0x400 | 0x4):
                    continue
                if r.mapping_quality >= LOWQ:
                    continue
                if r.reference_start < b:      # owned by the previous window
                    continue
                cig = r.cigartuples
                if not cig:
                    continue
                if cig[0][0] == 4:
                    side, clip = "L", cig[0][1]
                    est = exact = r.reference_start
                elif cig[-1][0] == 4:
                    side, clip = "R", cig[-1][1]
                    est = r.reference_start + rl - clip
                    exact = r.reference_end
                else:
                    continue
                ne, nx = near(c, est), near(c, exact)
                if not ne and not nx:
                    continue
                xa = tag_contigs(r, "XA")
                sa = tag_contigs(r, "SA")
                xac = sorted({contig_class(z) for z in xa})
                sac = sorted({contig_class(z) for z in sa})
                # alt-aware: the read's alternative hits are all on non-primary contigs
                xa_nonprim_only = bool(xa) and all(contig_class(z) != "primary" for z in xa)
                for C in sorted(set(ne) | set(nx)):
                    ie, ix = C in ne, C in nx
                    fo.write(f"{c}\t{C}\t{r.query_name}\t{r.flag}\t{r.mapping_quality}\t{side}\t{clip}\t"
                             f"{r.cigarstring}\t{est}\t{exact}\t{int(ie)}\t{int(ix)}\t{len(xa)}\t"
                             f"{','.join(xac) or '-'}\t{','.join(sac) or '-'}\n")
                    k = cnt[(c, C)]
                    if ie:
                        k["est"] += 1
                        if not xa_nonprim_only:
                            k["est_alt"] += 1
                        if not (r.flag & 0x900):
                            k["est_prim"] += 1
                    if ix:
                        k["exact"] += 1
                        if not xa_nonprim_only:
                            k["exact_alt"] += 1
                        if not (r.flag & 0x900):
                            k["exact_prim"] += 1
    fo.close()
    with open(out.replace(".tsv.gz", "") + ".counts.tsv", "w") as fc:
        cols = ["est", "exact", "est_alt", "exact_alt", "est_prim", "exact_prim"]
        fc.write("family\tchrom\tpos\tscanner\t" + "\t".join(cols) + "\n")
        for f in fams:
            sc = {}
            p = f"{sdir}/{f}/tmp/cns/low_mapq_clip_cnt.txt"
            if os.path.exists(p):
                for ln in open(p):
                    x = ln.split()
                    sc[(x[0], int(x[1]))] = int(x[2])
            for (c, C) in fam_cands[f]:
                k = cnt.get((c, C), {})
                fc.write(f"{f}\t{c}\t{C}\t{sc.get((c, C), 'NA')}\t" +
                         "\t".join(str(k.get(z, 0)) for z in cols) + "\n")
    print(f"alt_ref={alt_ref}", file=sys.stderr)


if __name__ == "__main__":
    main()
