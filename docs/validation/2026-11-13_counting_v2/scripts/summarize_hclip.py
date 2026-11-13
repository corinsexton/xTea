#!/usr/bin/env python
"""Summarize measure_hclip.py output: SR AF shift if hard-clipped supplementary junction reads
(per read, not already a soft clip) were added to the clip count.
usage: summarize_hclip.py SITES.tsv [...]"""
import sys
import numpy as np
import pandas as pd

def sr(clip, fmap):
    n = clip + fmap
    return np.where(n > 0, clip / np.maximum(n, 1), np.nan)

rows = []
for p in sys.argv[1:]:
    d = pd.read_csv(p, sep="\t")
    name = p.split("/")[-1].replace(".sites.tsv", "")
    # current: dominant side by soft counts
    L = d.n_l + d.h_l + d.s_supp_l; R = d.n_r + d.h_r + d.s_supp_r
    left = L >= R
    new_clip = np.where(left, L, R)
    new_fmap = np.where(left, d.fmap - d.lfm, d.fmap - d.rfm)
    d["sr_old"] = sr(d.af_clip, d.af_fmap)
    d["sr_new"] = sr(new_clip, new_fmap)
    d["added"] = new_clip - d.af_clip
    for (fam, gt), g in list(d.groupby(["family", "truth_gt"])) + [(("ALL", "all"), d)]:
        rows.append(dict(dataset=name, family=fam, truth_gt=gt, n=len(g),
                         sites_with_hclip=int(((g.h_l + g.h_r) > 0).sum()),
                         sites_with_supp=int(((g.h_l + g.h_r + g.s_supp_l + g.s_supp_r) > 0).sum()), hclip_reads=int((g.h_l + g.h_r).sum()), dup_with_soft=int((g.h_l_dup + g.h_r_dup).sum()),
                         soft_supp=int((g.s_supp_l + g.s_supp_r).sum()),
                         soft_clip_reads=int(g.af_clip.sum()), mean_added=round(g.added.mean(), 3),
                         med_sr_old=round(np.nanmedian(g.sr_old), 3), med_sr_new=round(np.nanmedian(g.sr_new), 3),
                         mean_sr_old=round(np.nanmean(g.sr_old), 3), mean_sr_new=round(np.nanmean(g.sr_new), 3)))
out = pd.DataFrame(rows)
out["d_med"] = (out.med_sr_new - out.med_sr_old).round(3)
out["d_mean"] = (out.mean_sr_new - out.mean_sr_old).round(3)
pd.set_option("display.width", 250); pd.set_option("display.max_rows", 500)
print(out.to_string(index=False))
