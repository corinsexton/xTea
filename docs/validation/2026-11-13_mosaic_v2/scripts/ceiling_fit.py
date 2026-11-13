"""Per-sample mosaic AF_PE ceiling vs the fixed 0.35, from the real disc_upper sweeps
(0.20-0.70 in 0.01 steps; each point a postF_mosaic -> merge -> gVCF rerun scored one-to-one,
+-200 bp, same family, against vafsite mix.tsv somatic sites; HG005-site calls counted apart).
A family's result depends only on that family's ceiling, so a per-family ceiling c_f is scored
by taking each family's row at round(c_f, 2) and summing (the cross-family locus merge is not
re-run; checked against a real run separately).
Rules: fixed c; k x m_f with m_f a per-family germline AF_PE statistic from the run's own
postF_hc calls (ceiling/germ_stats.tsv); k x m_pooled. Fit on BCM_300X + WashU_100X (max mean
per-run F1), test on UW_100X, and leave one run out. Writes ceiling/*.tsv."""
import glob, numpy as np, pandas as pd
B = "/n/scratch/users/c/cos689/xtea2_review/mosaic_v2"
RUNS = ["BCM_300X", "UW_100X", "WashU_100X"]
FAMS = ["Alu", "L1", "SVA"]
sw = pd.concat([pd.read_csv(f, sep="\t") for f in glob.glob(f"{B}/sweeps/*_fine*.tsv")])
sw = sw[sw.family.isin(FAMS)].drop_duplicates(["dataset", "family", "value"])
sw["v"] = (sw.value * 100).round().astype(int)
T = {(r.dataset, r.family, r.v): r for r in sw.itertuples()}
gs = pd.read_csv(f"{B}/ceiling/germ_stats.tsv", sep="\t")


def score_run(run, ceil):
    """ceil: {fam: ceiling}. Pooled TP/FP/hg005/F1 for one run."""
    tp = fp = hg = nt = 0
    for f in FAMS:
        v = int(round(min(max(ceil[f], 0.20), 0.70) * 100))
        r = T[(run, f, v)]
        tp += r.TP; fp += r.FP; hg += r.hg005; nt += r.n_truth
    return dict(TP=tp, FP=fp, hg005=hg, F1=2 * tp / (nt + tp + fp))


def m_of(run, stat):
    g = gs[gs.run == run].set_index("fam")
    if stat.endswith("_pooled"):
        s = stat[:-7]
        # pooled median approximated by the n-weighted mean of family medians is not a median;
        # use the het/gl count weights only to pick the dominant family -> report both rules
        w = g["n_het" if s == "med_het" else "n_gl"]
        val = float((g[s] * w).sum() / w.sum())
        return {f: val for f in FAMS}
    return {f: float(g.loc[f, stat]) for f in FAMS}


def rule(name, p, run):
    if name == "fixed":
        return {f: p for f in FAMS}
    m = m_of(run, name)
    return {f: p * m[f] for f in FAMS}


GRID = {"fixed": np.round(np.arange(0.20, 0.701, 0.01), 2)}
for s in ("med_het", "med_gl", "med_het_pooled"):
    GRID[s] = np.round(np.arange(0.40, 1.001, 0.01), 2)


def fit(name, train):
    best = None
    for p in GRID[name]:
        f1 = np.mean([score_run(r, rule(name, p, r))["F1"] for r in train])
        if best is None or f1 > best[1] + 1e-12:
            best = (p, f1)
    return best[0]


rows = []
for name in GRID:
    for test in RUNS + ["UW_100X*"]:
        if test == "UW_100X*":
            continue
        train = [r for r in RUNS if r != test]
        p = fit(name, train)
        c = rule(name, p, test)
        s = score_run(test, c)
        rows.append(dict(rule=name, held_out=test, param=p, **{f"ceil_{f}": round(c[f], 3) for f in FAMS}, **s))
loro = pd.DataFrame(rows)
loro.to_csv(f"{B}/ceiling/loro.tsv", sep="\t", index=False, float_format="%.4f")
print("== leave one run out (fit on the other two, max mean per-run F1); UW_100X row = fit on BCM+WashU")
print(loro.round(4).to_string(index=False))
# all runs at the pooled best and at 0.35, per-run table
rows = []
for name in GRID:
    p_all = fit(name, RUNS)
    for r in RUNS:
        for lab, p in (("fit_all", p_all),) + ((("0.35", 0.35),) if name == "fixed" else ()):
            c = rule(name, p, r)
            rows.append(dict(rule=name, setting=lab, param=p, run=r, **{f"ceil_{f}": round(c[f], 3) for f in FAMS},
                             **score_run(r, c)))
tab = pd.DataFrame(rows)
tab.to_csv(f"{B}/ceiling/per_run.tsv", sep="\t", index=False, float_format="%.4f")
print("== per run"); print(tab.round(4).to_string(index=False))
# per-run optimum of a fixed ceiling, and its ratio to med_het
rows = []
for r in RUNS:
    f1 = {p: score_run(r, {f: p for f in FAMS})["F1"] for p in GRID["fixed"]}
    pb = max(f1, key=f1.get)
    within = [p for p in f1 if f1[p] >= f1[pb] - 0.002]
    rows.append(dict(run=r, best=pb, F1_best=f1[pb], F1_035=f1[0.35], plateau_lo=min(within), plateau_hi=max(within)))
print("== per-run fixed optimum (plateau: F1 within 0.002 of best)"); print(pd.DataFrame(rows).round(4).to_string(index=False))
pd.DataFrame(rows).to_csv(f"{B}/ceiling/per_run_optimum.tsv", sep="\t", index=False, float_format="%.4f")
