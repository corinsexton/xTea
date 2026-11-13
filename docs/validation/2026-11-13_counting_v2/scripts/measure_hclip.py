#!/usr/bin/env python
"""Issue B measurement: hard-clipped supplementary junction reads at final calls.

usage: measure_hclip.py RUN_DIR NAME BAM REF OUT_PREFIX

For every final call in RUN_DIR/NAME/NAME/NAME.vcf that matches truth (score_gt, 200 bp,
one-to-one), re-reads the BAM at the genotype-feature position (the gntp.features row of that
family) exactly as collect_features_one_site does for the split-read AF, and additionally counts
supplementary records clipped at the breakpoint (|clip pos - ins_pos| < CLIP_EXACT_CLIP_SLACK):

  h_l / h_r      hard-clipped (op 5) supplementaries, left / right clip at the breakpoint, whose
                 READ (query name + R1/R2) is not already counted as a soft clip on either side
  h_l_dup/h_r_dup  same, but the read is already counted through a soft-clipped primary
  hq_l / hq_r    as h_l / h_r but keyed by query name only (template), for comparison
  s_supp_l/_r    soft-clipped (op 4) supplementaries at the breakpoint, deduplicated like h_l/h_r
                 (bwa -Y style BAMs; none expected with bwa defaults)
  hraw_l/hraw_r  hard-clipped supplementaries within TSD_CUTOFF (the raw-clip-count window)

Checks that the soft recount equals the pipeline's n_af_clip / n_full_map (features fields 2/3).
Writes OUT_PREFIX.sites.tsv (one row per TP call).
"""
import os
import sys

import pysam

G = "/n/scratch/users/c/cos689/xtea2_testdata/gtval"
sys.path.insert(0, G)
import score_gt as S  # noqa: E402

H = "/n/data1/hms/dbmi/park/corinne/claude_sandbox/test_data/HG002/truth"
SLACK = 3        # CLIP_EXACT_CLIP_SLACK
TSD = 100        # TSD_CUTOFF
HALF = 45        # BWA_HALF_READ_MIN_SCORE
EXT = 550        # DFT_IS (fetch half-window); value used by the pipeline is logged per run


def ref_end(cig, start):
    p = start
    for op, ln in cig[:-1]:
        if op in (4, 5, 1):
            continue
        p += ln
    return p


def count_site(bam, chrm, ins):
    r = dict(n_l=0, n_r=0, fmap=0, lfm=0, rfm=0, h_l=0, h_r=0, h_l_dup=0, h_r_dup=0,
             hq_l=0, hq_r=0, s_supp_l=0, s_supp_r=0, hraw_l=0, hraw_r=0, h_l_len=[], h_r_len=[])
    soft = set()      # reads (qname, is_read1) with a soft clip counted at the breakpoint
    soft_q = set()
    hard = []         # (side, key, qname, cliplen)
    st = max(1, ins - EXT)
    for a in bam.fetch(chrm, st, ins + EXT):
        if a.is_duplicate or a.is_secondary or a.is_unmapped:
            continue
        cig = a.cigartuples
        if not cig:
            continue
        mp = a.reference_start
        key = (a.query_name, a.is_read1)
        if len(cig) == 1 and cig[0][0] == 0:
            e = mp + cig[0][1]
            if mp <= ins <= e:
                r["fmap"] += 1
                if abs(mp - ins) < HALF:
                    r["rfm"] += 1
                elif abs(e - ins) < HALF:
                    r["lfm"] += 1
        if a.is_supplementary:
            if cig[0][0] in (4, 5) and abs(mp - ins) < SLACK:
                hard.append(("l", key, a.query_name, cig[0][1], cig[0][0]))
            if cig[0][0] == 5 and abs(mp - ins) < TSD:
                r["hraw_l"] += 1
            if cig[-1][0] in (4, 5):
                cp = ref_end(cig, mp)
                if abs(cp - ins) < SLACK:
                    hard.append(("r", key, a.query_name, cig[-1][1], cig[-1][0]))
                if cig[-1][0] == 5 and abs(cp - ins) < TSD:
                    r["hraw_r"] += 1
            continue
        if cig[0][0] == 4 and abs(mp - ins) < SLACK:
            r["n_l"] += 1; soft.add(key); soft_q.add(a.query_name)
        if cig[-1][0] == 4:
            cp = ref_end(cig, mp)
            if abs(cp - ins) < SLACK:
                r["n_r"] += 1; soft.add(key); soft_q.add(a.query_name)
    seen = set(); seen_q = set()
    for side, key, q, ln, op in hard:
        if key in soft:
            r["h_%s_dup" % side] += 1
        elif key not in seen:
            seen.add(key)
            if op == 5:
                r["h_%s" % side] += 1; r["h_%s_len" % side].append(ln)
            else:
                r["s_supp_%s" % side] += 1
        if q not in soft_q and q not in seen_q:
            seen_q.add(q); r["hq_%s" % side] += 1
    return r


def load_feat(run, name):
    m = {}
    for fam in S.FAMILIES:
        p = os.path.join(run, name, name, fam, "candidate_disc_filtered_cns.txt.gntp.features")
        for ln in open(p):
            x = ln.split()
            m[(fam, x[0], int(x[1]))] = (int(x[2]), int(x[3]))
    return m


def main():
    run, name, bam_p, ref, out = sys.argv[1:6]
    truth = S.load_truth(H) if name.startswith("HG002") else \
        S.load_truth(os.path.join(G, "truth", "truth_%s.tsv" % name), True, True)
    calls = S.load_calls([os.path.join(run, name, name, name + ".vcf")])
    res = S.score(truth, calls, 200)
    feat = load_feat(run, name)
    vend = {}
    for ln in open(os.path.join(run, name, name, name + ".vcf")):
        if ln.startswith("#"):
            continue
        x = ln.split("\t")
        info = dict(kv.split("=", 1) for kv in x[7].split(";") if "=" in kv)
        fam_v = S.VCF_ALT_TO_FAMILY.get(x[4].strip("<>"))
        vend[(fam_v, x[0], int(x[1]))] = (int(info.get("END", x[1])), int(info["AF_CLIP"]), int(info["AF_FMAP"]))
    bam = pysam.AlignmentFile(bam_p, reference_filename=ref)
    cols = ["family", "chrom", "pos", "feat_pos", "truth_gt", "call_gt", "af_clip", "af_fmap",
            "re_clip", "re_fmap", "n_l", "n_r", "fmap", "lfm", "rfm", "h_l", "h_r", "h_l_dup",
            "h_r_dup", "hq_l", "hq_r", "s_supp_l", "s_supp_r", "hraw_l", "hraw_r", "h_lens"]
    nmiss = nbad = 0
    with open(out + ".sites.tsv", "w") as fo:
        fo.write("\t".join(cols) + "\n")
        for fam in S.FAMILIES:
            for t, c in res[fam]["pairs"]:
                fp = None
                end, v_clip, v_fmap = vend[(fam, c.chrom, c.pos)]
                for q in (end, c.pos):
                    if feat.get((fam, c.chrom, q)) == (v_clip, v_fmap):
                        fp = q; break
                if fp is None:
                    nmiss += 1; continue
                af_clip, af_fmap = feat[(fam, c.chrom, fp)]
                r = count_site(bam, c.chrom, fp)
                if r["n_l"] < r["n_r"]:
                    re_clip, re_fmap = r["n_r"], r["fmap"] - r["rfm"]
                else:
                    re_clip, re_fmap = r["n_l"], r["fmap"] - r["lfm"]
                if (re_clip, re_fmap) != (af_clip, af_fmap):
                    nbad += 1
                fo.write("\t".join(map(str, [fam, c.chrom, c.pos, fp, t.gt, c.gt, af_clip, af_fmap,
                                             re_clip, re_fmap] + [r[k] for k in cols[10:-1]]
                                       + [":".join(map(str, r["h_l_len"] + r["h_r_len"])) or "none"])) + "\n")
    print("no features row: %d; recount mismatches: %d" % (nmiss, nbad))


if __name__ == "__main__":
    main()
