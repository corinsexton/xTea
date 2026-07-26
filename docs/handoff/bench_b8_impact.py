"""Re-measure B8's impact at the REAL postF parameters.

The spec's B8 table used L1 min_len=300, but x_post_filter.py:415 overrides it with
L1_min_ref_copy_len=125 for the L1 branch. Re-measure, and also quantify how often the
min-divergence rule (the chosen fix) differs from what the set happened to return.

HISTORICAL — this script measures the OLD vendored AVL tree against its replacement, so it
requires `xtea/intervaltree.py` + `interval.py` + `node.py`, which were deleted in step 2.
To re-run it, check the tree out first:

    git checkout 4a9e824 -- xtea/intervaltree.py xtea/interval.py xtea/node.py

Results are recorded in docs/clipdisc_contract_spec.md §4 P1-EVALUATED.
"""
import sys, time
from bisect import bisect_right

sys.path.insert(0, "/n/data1/hms/dbmi/park/corinne/claude_sandbox/xTea/xtea")
from intervaltree import IntervalTree

ALU = "/home/cos689/xTea/rep_lib_annotation/Alu/hg38/hg38_Alu.out"
L1 = "/home/cos689/xTea/rep_lib_annotation/LINE/hg38/hg38_L1.fa.out"
CUTOFF = 15.0  # global_values.REP_DIVERGENT_CUTOFF


def parse(path, extnd, min_len):
    m = {}
    with open(path) as fin:
        for line in fin:
            if len(line.rstrip()) < 1:
                continue
            f = line.split()
            if len(f) < 15 or f[0].rstrip() in ("SW", "score"):
                continue
            s, e = int(f[5]), int(f[6])
            if abs(e - s) < min_len:
                continue
            m.setdefault(f[4], {}).setdefault(s - extnd, []).append((e + extnd, float(f[1])))
    return m


def run(name, path, extnd, min_len):
    print(f"\n=== {name}  extnd={extnd}  min_len={min_len} ===")
    m = parse(path, extnd, min_len)
    n_iv = sum(len(v) for v in m.values())

    n_q = n_multi = 0
    n_setfirst_ne_min = 0     # set's pick != min-divergence pick
    n_verdict_flip = 0        # the two picks land on opposite sides of the 15% cutoff
    n_setfirst_ne_left = 0    # set's pick != leftmost (the spec's measured column)

    for c, mp in m.items():
        starts = sorted(mp.keys())
        ends = [mp[s][0][0] for s in starts]
        divs = [mp[s][0][1] for s in starts]
        pm, run_max = [], -1
        for e in ends:
            if e > run_max:
                run_max = e
            pm.append(run_max)

        tree = IntervalTree()
        for s, e in zip(starts, ends):
            tree.addi(s, e)

        for p in starts:
            n_q += 1
            hi = bisect_right(starts, p)
            idx = []
            i = hi - 1
            while i >= 0 and pm[i] > p:
                if ends[i] > p:
                    idx.append(i)
                i -= 1
            idx.reverse()
            if len(idx) <= 1:
                continue
            n_multi += 1

            # what the CURRENT code returns: first element of a set iteration
            st = next(iter(tree[p]))[0]
            div_set = divs[starts.index(st)] if False else mp[st][0][1]

            # the chosen fix: minimum divergence, tie-broken leftmost
            best = idx[0]
            for j in idx[1:]:
                if divs[j] < divs[best]:
                    best = j
            div_min = divs[best]
            st_min = starts[best]
            st_left = starts[idx[0]]

            if st != st_min:
                n_setfirst_ne_min += 1
            if st != st_left:
                n_setfirst_ne_left += 1
            # verdict = "in a young copy" -> call dropped
            if (div_set < CUTOFF) != (div_min < CUTOFF):
                n_verdict_flip += 1

    print(f"  intervals                     {n_iv:,}")
    print(f"  queries                       {n_q:,}")
    print(f"  multi-overlap                 {n_multi:,} ({100*n_multi/n_q:.1f}%)")
    print(f"  set-first != leftmost         {n_setfirst_ne_left:,} ({100*n_setfirst_ne_left/n_q:.1f}%)")
    print(f"  set-first != min-divergence   {n_setfirst_ne_min:,} ({100*n_setfirst_ne_min/n_q:.1f}%)")
    print(f"  FILTER VERDICT FLIPS          {n_verdict_flip:,} ({100*n_verdict_flip/n_q:.2f}%)")


run("Alu", ALU, 100, 200)
run("L1 (real postF param)", L1, 400, 125)
run("L1 (spec's recorded param)", L1, 400, 300)
