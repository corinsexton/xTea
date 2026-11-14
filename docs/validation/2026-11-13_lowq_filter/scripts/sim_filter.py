#!/usr/bin/env python
"""Offline replay of the cns clip-consistency check with different low-MAPQ clip counts.

usage: sim_filter.py WT RUN_SAMPLE_DIR COUNTS_TSV OUT_TSV

Re-parses each family's tmp/cns/temp_clip.sam exactly as call_MEIs_consensus does and runs
check_clip_alone_consistency2 once with m_lowq_clip=None (filter off) and once per variant
column of COUNTS_TSV (from lowq_reads.py). Writes one row per candidate that passes the check
with the filter off: family, chrom, pos, cnt_clip (realigned clips in the check), and for each
variant its low-q count and whether the filter drops the site.
"""
import os
import re
import sys
from collections import defaultdict


def main():
    wt, sdir, counts, out = sys.argv[1:5]
    sys.path.insert(0, os.path.join(wt, "xtea"))
    import global_values
    from x_clip_disc_filter import XClipDiscFilter

    cols = None
    cnts = defaultdict(dict)
    for ln in open(counts):
        x = ln.rstrip("\n").split("\t")
        if cols is None:
            cols = x[3:]
            continue
        cnts[x[0]][(x[1], int(x[2]))] = dict(zip(cols, [int(v) if v != "NA" else 0 for v in x[3:]]))

    fo = open(out, "w")
    fo.write("family\tchrom\tpos\tcnt_clip\t" + "\t".join(f"n_{c}\tdrop_{c}" for c in cols) + "\n")
    for fam in ("L1", "Alu", "SVA"):
        wf = f"{sdir}/{fam}/tmp/cns/"
        if not os.path.exists(wf + "temp_clip.sam"):
            continue
        log = open(f"{sdir}/{os.path.basename(sdir)}_{fam}.steps.log").read()
        m = re.search(r"Filter \(on cns\) cutoff: (\d+) and (\d+)", log)
        nclip = int(m.group(1))
        import family_params
        f = XClipDiscFilter(f"{sdir}/bam_list.txt", wf, 1, None,
                            family=family_params.from_sva_flag(fam == "SVA"))
        lpos, rpos, _, _, _, _ = f.parse_clip_realignment_consensus(
            wf + "temp_clip.sam", global_values.MIN_CLIP_MAPPED_RATIO)
        idist = global_values.CLIP_SEARCH_WINDOW
        ratio = global_values.CLIP_CONSISTENT_RATIO
        base = f.check_clip_alone_consistency2(lpos, rpos, idist, ratio, nclip, idist, None)
        # cnt_clip as the check sees it: left+right peak clip counts
        lchk = f.check_one_side_clip_consistency(lpos, idist, ratio, idist)
        rchk = f.check_one_side_clip_consistency(rpos, idist, ratio, idist)
        res = {}
        for c in cols:
            m_lowq = defaultdict(dict)
            for (ch, p), d in cnts[fam].items():
                m_lowq[ch][p] = d[c]
            res[c] = f.check_clip_alone_consistency2(lpos, rpos, idist, ratio, nclip, idist, m_lowq)
        for ch in base:
            for p in base[ch]:
                n = 0
                if ch in lchk and p in lchk[ch]:
                    n += lchk[ch][p][-1]
                if ch in rchk and p in rchk[ch]:
                    n += rchk[ch][p][-1]
                row = [fam, ch, str(p), str(n)]
                for c in cols:
                    row += [str(cnts[fam].get((ch, p), {}).get(c, 0)),
                            str(int(not (ch in res[c] and p in res[c][ch])))]
                fo.write("\t".join(row) + "\n")
    fo.close()


if __name__ == "__main__":
    main()
