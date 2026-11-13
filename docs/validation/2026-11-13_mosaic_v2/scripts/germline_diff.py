"""Diff the germline VCFs re-made from the gtval runs_b15 inputs by xTea2 (base) and
vaf-site-v2 (v2): call set, GT, GQ, every other field; VAF / VAF_CI and the ##xtea2_vaf_params
header reported separately. Writes germline_check/diff.tsv."""
import re, sys
B = "/n/scratch/users/c/cos689/xtea2_review/mosaic_v2/germline_check"
DS = ["HG002_HiSeq30x", "HG002_NovaSeq30x", "HG002_ds10x", "HG002_ds15x", "HG002_ds20x", "HG00438", "HG02257", "HG02622"]
VAFKEYS = {"VAF", "VAF_CI", "VAF_LO", "VAF_HI"}


def load(p):
    hdr, recs = [], {}
    for ln in open(p):
        if ln.startswith("##"):
            hdr.append(ln.rstrip()); continue
        if ln.startswith("#"):
            continue
        f = ln.rstrip("\n").split("\t")
        info = dict(kv.split("=", 1) if "=" in kv else (kv, True) for kv in f[7].split(";"))
        fmt = dict(zip(f[8].split(":"), f[9].split(":"))) if len(f) > 9 else {}
        recs[(f[0], f[1], f[3], f[4])] = (f[:7], info, fmt)
    return hdr, recs


print("dataset\tn_base\tn_v2\tcalls_only_one\tGT_diff\tGQ_diff\tother_diff\tVAF_diff\tVAF_CI_diff\thdr_diff_lines")
for n in DS:
    hb, rb = load(f"{B}/base/{n}/{n}/{n}.vcf")
    hv, rv = load(f"{B}/v2/{n}/{n}/{n}.vcf")
    only = len(set(rb) ^ set(rv))
    gt = gq = oth = vaf = ci = 0
    otherkeys = set()
    for k in set(rb) & set(rv):
        (cb, ib, fb), (cv, iv, fv) = rb[k], rv[k]
        gt += fb.get("GT") != fv.get("GT")
        gq += fb.get("GQ") != fv.get("GQ")
        vaf += (ib.get("VAF"), fb.get("VAF")) != (iv.get("VAF"), fv.get("VAF"))
        ci += (ib.get("VAF_CI"), fb.get("VAF_CI")) != (iv.get("VAF_CI"), fv.get("VAF_CI"))
        d = cb != cv
        for key in (set(ib) | set(iv)) - VAFKEYS:
            if ib.get(key) != iv.get(key):
                d = True; otherkeys.add("INFO/" + key)
        for key in (set(fb) | set(fv)) - VAFKEYS - {"GT", "GQ"}:
            if fb.get(key) != fv.get(key):
                d = True; otherkeys.add("FORMAT/" + key)
        oth += d
    hd = [l for l in set(hb) ^ set(hv)]
    hdk = sorted({re.match(r"##([^=]+)", l).group(1) for l in hd})
    print(f"{n}\t{len(rb)}\t{len(rv)}\t{only}\t{gt}\t{gq}\t{oth}{'(' + ','.join(sorted(otherkeys)) + ')' if otherkeys else ''}\t{vaf}\t{ci}\t{len(hd)}:{','.join(hdk)}")
