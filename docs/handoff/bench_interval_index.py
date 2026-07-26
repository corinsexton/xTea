"""Characterize the repeat-annotation query problem and benchmark candidate structures.

Question being answered: the vendored AVL interval tree is 96% of postF. What does the
workload ACTUALLY require, and what is the right structure for it?

HISTORICAL — this script measures the OLD vendored AVL tree against its replacement, so it
requires `xtea/intervaltree.py` + `interval.py` + `node.py`, which were deleted in step 2.
To re-run it, check the tree out first:

    git checkout 4a9e824 -- xtea/intervaltree.py xtea/interval.py xtea/node.py

Results are recorded in docs/clipdisc_contract_spec.md §4 P1-EVALUATED.
"""
import sys, os, time, random
from bisect import bisect_right

sys.path.insert(0, "/n/data1/hms/dbmi/park/corinne/claude_sandbox/xTea/xtea")

ALU = "/home/cos689/xTea/rep_lib_annotation/Alu/hg38/hg38_Alu.out"
L1 = "/home/cos689/xTea/rep_lib_annotation/LINE/hg38/hg38_L1.fa.out"


def t(label, fn):
    s = time.time()
    r = fn()
    e = time.time() - s
    print(f"  {label:<46s} {e:8.2f} s")
    return r, e


# ---------------------------------------------------------------- parsing
def parse_python(path, extnd, min_len):
    """Exactly what load_rmsk_annotation_with_extnd_div_with_lenth_cutoff does."""
    m = {}
    with open(path) as fin:
        for line in fin:
            if len(line.rstrip()) < 1:
                continue
            fields = line.split()
            if len(fields) < 15:
                continue
            if fields[0].rstrip() == "SW" or fields[0].rstrip() == "score":
                continue
            mismatch_rate = float(fields[1])
            chrm = fields[4]
            start_pos = int(fields[5])
            end_pos = int(fields[6])
            if abs(end_pos - start_pos) < min_len:
                continue
            sub_type = fields[9]
            family = fields[10]
            b_rc = fields[8] == "C"
            if b_rc:
                csn_start, csn_end = int(fields[13]), int(fields[12])
            else:
                csn_start, csn_end = int(fields[11]), int(fields[12])
            e_s = start_pos - extnd
            e_e = end_pos + extnd
            m.setdefault(chrm, {}).setdefault(e_s, []).append(
                (e_e, mismatch_rate, b_rc, sub_type, family, csn_start, csn_end))
    return m


def parse_numpy(path, extnd, min_len):
    """Same fields, but column-parsed with pandas."""
    import pandas as pd
    import numpy as np
    df = pd.read_csv(path, sep=r"\s+", header=None, engine="c",
                     usecols=[1, 4, 5, 6, 8, 9],
                     names=["div", "chrm", "start", "end", "strand", "subtype"],
                     dtype={1: str, 4: str, 5: str, 6: str, 8: str, 9: str},
                     comment=None, skip_blank_lines=True, on_bad_lines="skip")
    df = df[pd.to_numeric(df["start"], errors="coerce").notna()]
    start = df["start"].astype(np.int64).to_numpy()
    end = df["end"].astype(np.int64).to_numpy()
    div = df["div"].astype(np.float64).to_numpy()
    chrm = df["chrm"].to_numpy()
    keep = np.abs(end - start) >= min_len
    return chrm[keep], start[keep] - extnd, end[keep] + extnd, div[keep]


# ---------------------------------------------------------------- structures
class SortedArrayIndex:
    """Plan's proposal: per-chrom sorted starts + prefix-max-end backward scan.
    Returns ALL overlapping intervals (so it can apply a selection rule)."""
    __slots__ = ("starts", "ends", "divs", "pmax")

    def __init__(self, starts, ends, divs):
        order = sorted(range(len(starts)), key=lambda i: starts[i])
        self.starts = [starts[i] for i in order]
        self.ends = [ends[i] for i in order]
        self.divs = [divs[i] for i in order]
        pm, run = [], -1
        for e in self.ends:
            if e > run:
                run = e
            pm.append(run)
        self.pmax = pm

    def hits(self, p):
        hi = bisect_right(self.starts, p)
        out = []
        ends, pmax = self.ends, self.pmax
        i = hi - 1
        while i >= 0 and pmax[i] > p:
            if ends[i] > p:
                out.append(i)
            i -= 1
        out.reverse()
        return out


class MergedBooleanIndex:
    """Boolean-only: merge overlaps into disjoint runs, store as one flat
    alternating boundary array. Query = ONE bisect, parity test."""
    __slots__ = ("bounds",)

    def __init__(self, starts, ends):
        iv = sorted(zip(starts, ends))
        bounds = []
        cs, ce = None, None
        for s, e in iv:
            if cs is None:
                cs, ce = s, e
            elif s <= ce:
                if e > ce:
                    ce = e
            else:
                bounds.append(cs); bounds.append(ce)
                cs, ce = s, e
        if cs is not None:
            bounds.append(cs); bounds.append(ce)
        self.bounds = bounds

    def contains(self, p):
        # inside iff index is odd (begin <= p < end, half-open)
        return (bisect_right(self.bounds, p) & 1) == 1


def main():
    which = sys.argv[1] if len(sys.argv) > 1 else "Alu"
    if which == "Alu":
        path, extnd, min_len = ALU, 100, 200
    else:
        path, extnd, min_len = L1, 400, 125

    print(f"\n{'='*70}\n{which}  extnd={extnd}  min_len={min_len}\n{'='*70}")

    print("\n[PARSE]")
    m, t_py = t("current Python line loop", lambda: parse_python(path, extnd, min_len))
    (chrm_a, st_a, en_a, dv_a), t_np = t("pandas/numpy column parse",
                                         lambda: parse_numpy(path, extnd, min_len))
    n_total = sum(len(v) for v in m.values())
    print(f"    -> {n_total:,} distinct start keys across {len(m)} chroms"
          f" | numpy rows {len(st_a):,} | parse speedup {t_py/t_np:.1f}x")

    # build per-chrom python lists from the dict (the real input to indexing)
    per_chrom = {}
    for c, mp in m.items():
        ss = sorted(mp.keys())
        per_chrom[c] = (ss, [mp[s][0][0] for s in ss], [mp[s][0][1] for s in ss])

    print("\n[BUILD]")
    try:
        from intervaltree import IntervalTree

        def build_avl():
            trees = {}
            for c, (ss, ee, _) in per_chrom.items():
                tr = IntervalTree()
                for s, e in zip(ss, ee):
                    tr.addi(s, e)
                trees[c] = tr
            return trees
        avl, t_avl = t("vendored AVL IntervalTree (current)", build_avl)
    except Exception as ex:
        print("    AVL build failed:", ex)
        avl, t_avl = None, float("nan")

    sai, t_sai = t("sorted arrays + prefix-max-end", lambda: {
        c: SortedArrayIndex(ss, ee, dd) for c, (ss, ee, dd) in per_chrom.items()})
    mbi, t_mbi = t("merged disjoint boundary array (bool only)", lambda: {
        c: MergedBooleanIndex(ss, ee) for c, (ss, ee, _) in per_chrom.items()})

    nb = sum(len(v.bounds) for v in mbi.values())
    print(f"    -> merged index is {nb//2:,} disjoint runs"
          f" (vs {n_total:,} raw intervals, {100*(1-nb/2/n_total):.1f}% fewer)")
    if t_avl == t_avl:
        print(f"    -> build speedup: sorted {t_avl/t_sai:.0f}x, merged {t_avl/t_mbi:.0f}x")

    # ------------------------------------------------------------ data shape
    print("\n[DATA SHAPE]")
    lens = [e - s for c, (ss, ee, _) in per_chrom.items() for s, e in zip(ss, ee)]
    lens.sort()
    print(f"    interval length: median {lens[len(lens)//2]:,}  "
          f"p99 {lens[int(len(lens)*0.99)]:,}  max {lens[-1]:,}")

    # multi-overlap rate over every interval start (B8's measurement)
    n_multi = n_q = 0
    for c, (ss, _, _) in per_chrom.items():
        idx = sai[c]
        for p in ss:
            h = idx.hits(p)
            n_q += 1
            if len(h) > 1:
                n_multi += 1
    print(f"    multi-overlap: {n_multi:,}/{n_q:,} = {100*n_multi/n_q:.1f}% of queries")

    # ------------------------------------------------------------ query
    print("\n[QUERY]  (100k random positions, hot-path pattern)")
    random.seed(7)
    big = max(per_chrom, key=lambda c: len(per_chrom[c][0]))
    lo, hi = per_chrom[big][0][0], per_chrom[big][0][-1]
    qs = [random.randint(lo, hi) for _ in range(100_000)]

    if avl:
        def q_avl():
            tr = avl[big]
            n = 0
            for p in qs:
                if len(tr[p]):
                    n += 1
            return n
        n1, t1 = t("AVL tree point query (current)", q_avl)
    else:
        n1, t1 = None, float("nan")

    def q_sai():
        idx = sai[big]
        n = 0
        for p in qs:
            if idx.hits(p):
                n += 1
        return n
    n2, t2 = t("sorted arrays + scan (returns all hits)", q_sai)

    def q_mbi():
        idx = mbi[big]
        n = 0
        for p in qs:
            if idx.contains(p):
                n += 1
        return n
    n3, t3 = t("merged boundary, 1 bisect (bool only)", q_mbi)

    def q_np():
        import numpy as np
        b = np.array(mbi[big].bounds, dtype=np.int64)
        a = np.array(qs, dtype=np.int64)
        return int(((np.searchsorted(b, a, side="right") & 1) == 1).sum())
    n4, t4 = t("numpy vectorized batch (bool only)", q_np)

    print(f"    -> agreement: AVL {n1} | sorted {n2} | merged {n3} | numpy {n4}")
    if t1 == t1:
        print(f"    -> query speedup vs AVL: sorted {t1/t2:.1f}x, "
              f"merged {t1/t3:.1f}x, numpy {t1/t4:.0f}x")


main()
