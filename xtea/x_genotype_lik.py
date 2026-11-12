"""Closed-form likelihood genotyper for germline MEI calls, with per-sample EM calibration.

Replaces the ONNX random forest (x_genotype_classify_sklearn / genotype_model_6_7_2024.onnx),
which reached 75% GT accuracy on HG002 and called only 10% of true homozygotes 1/1.

Model
-----
Each call carries two pairs of read counts, both written by x_genotype_feature:

    SR (split reads):  alt = AF_CLIP (clipped at the breakpoint), ref = AF_FMAP (fully mapped)
    PE (read pairs):   alt = AF_DISC (discordant),                ref = AF_CONCORDNT

Neither alt fraction is 0.5 / 1.0 for a het / hom insertion: alt reads are under-counted
(mapping and clipping bias), and a hom call still has "ref" reads (reads spanning the TSD,
pairs whose insert straddles a short insertion). So each read type t and genotype g get their
own expected alt fraction b_t[g]. The two genotypes considered are 0/1 and 1/1; a called
insertion is never 0/0 here.

    log L(g) = log BetaBin(clip | clip+fmap, mean b_SR[g], overdispersion rho)
             + log BetaBin(disc | disc+conc, mean b_PE[g], overdispersion rho)

The four biases and the hom fraction w are estimated per sample by EM over all of that
sample's calls (no labels needed); that is what makes the genotyper portable across
libraries, read lengths and insert-size distributions. The E-step uses the beta-binomial
above; the M-step sets each b_t[g] to the responsibility-weighted pooled alt fraction (the
binomial MLE, i.e. a method-of-moments step for the beta-binomial mean). When there are too few
informative calls, or EM lands on a degenerate solution, the HG002 defaults are used instead.

Per-call outputs (FORMAT fields)
--------------------------------
    GT      0/1 or 1/1 (./. when the call has no SR and no PE reads)
    GQ      phred of the posterior error after an error floor: perr' = eps + (1 - 2 eps) perr,
            eps = GQ_ERR_FLOOR = 0.001, so GQ tops out at 30. The floor stands for sites whose
            counts do not reflect the genotype (complex or nested insertions, duplicate calls
            at one site); it does not change GT. Conditional on the site carrying the
            insertion: 0/1 vs 1/1 only, 0/0 is not modeled
    GTL     log10 likelihoods of (0/1, 1/1). Named GTL, not GL, because VCF reserves GL
            for Number=G (all three diploid genotypes)
    AF_SR   raw clip/(clip+fmap)
    AF_PE   raw disc/(disc+conc)
    VAF     maximum-likelihood allele fraction from x_vaf_site: the expected SR and PE alt
            fractions of each call are worked out from the read length, the insert-size
            distribution, the insertion length and the TSD, so a short Alu and a full-length
            L1 at the same VAF are no longer read off the same line. Needs the sample's read
            length and insert size (basic_info / --basic-info); without them VAF falls back to
            the bias line below.
    VAF_CI  95% profile-likelihood interval for VAF (x_vaf_site).

    Fallback (no read length / insert size): each raw fraction is mapped through its type's
    fitted bias line -- piecewise linear through (0, 0), (b_het, 0.5), (b_hom, 1.0) -- and
    the two are averaged weighted by depth; VAF_CI is Beta(alt+0.5, ref+0.5) on the pooled
    counts mapped through the depth-weighted bias line.

The fitted parameters are written to the VCF header as ##xtea2_genotype_params=<...> and
##xtea2_vaf_params=<...>.

Mosaic mode (mosaic=True / --mosaic) writes no genotype: GT is ./., GQ and GTL are missing,
and VAF / VAF_CI come from x_vaf_site with a background term and calibration on the sample's
germline calls (see x_vaf_site.estimate).

Post-hoc use
------------
    python -m xtea.x_genotype_lik in.vcf -o out.vcf      (or: xtea2-genotype in.vcf -o out.vcf)

re-genotypes an existing xTea2 VCF from its INFO counts (AF_CLIP, AF_FMAP, AF_DISC,
AF_CONCORDNT) and replaces the FORMAT/sample columns. All records in the file are fitted
together, so pass the merged (all-family) VCF to get a per-sample fit.

Dependencies: numpy and scipy only. This module imports nothing else from xTea, so it works
both as `xtea.x_genotype_lik` and as the flat `x_genotype_lik` the pipeline scripts use.
"""
import argparse
import gzip
import math
import os
import sys

import numpy as np
from scipy.special import betaincinv, betaln, gammaln

try:                                    # xtea.x_genotype_lik (package / CLI entry point)
    from . import x_vaf_site as xvs
except ImportError:                     # flat x_genotype_lik (pipeline scripts)
    import x_vaf_site as xvs

# HG002 (30x NovaSeq PCR-free) EM fit, used when a sample cannot support its own fit.
# Fit on the full v2 rerun with the AF_FMAP fix, read pairs counted once each on both sides of
# AF_PE, and the pair checks on right-clipped reads measured from the alignment start (1494
# calls, all informative); the raw fit was sr 0.4411/0.9577, pe 0.5703/0.9137, hom_prior 0.3083.
DEFAULT_PARAMS = {"sr_het": 0.441, "sr_hom": 0.958, "pe_het": 0.570, "pe_hom": 0.914,
                  "hom_prior": 0.308}
DEFAULT_RHO = 0.05
GQ_CAP = 60
# Per-call error floor on GQ (see gq_from_error). Fit on HG002 NovaSeq 30x (GIAB Q100 v5.0q
# labels) plus HG00438 (HPRC tier1): MLE 0.00084, 95% CI 0.00013-0.0025; rounded to 0.001.
# The beta-binomial alone cannot reach this floor without raising GT errors.
GQ_ERR_FLOOR = 0.001
MIN_INFORMATIVE = 50        # fewer informative calls than this -> defaults
# Parameter ranges. Bias estimates are clamped into these inside EM; a het bias or a hom
# prior that ends on a bound, or a het/hom gap below MIN_SEP, marks the fit as degenerate.
# The hom-bias upper bound is a cap, not a failure: SR hom fractions near 1 are expected.
HET_RANGE = (0.10, 0.80)
HOM_RANGE = (0.40, 0.995)
PRIOR_RANGE = (0.02, 0.95)
MIN_SEP = 0.10
EM_MAX_ITER = 500
EM_TOL = 1e-7
_KEYS = ("sr_het", "sr_hom", "pe_het", "pe_hom")


# ------------------------------------------------------------------------------------------
# likelihood
# ------------------------------------------------------------------------------------------
def betabinom_logpmf(k, n, p, rho=DEFAULT_RHO):
    """log P(k | n) for a beta-binomial with mean p and intra-class correlation rho.

    rho <= 0 gives the binomial. n = 0 gives 0 (the read type contributes nothing)."""
    k = np.asarray(k, float)
    n = np.asarray(n, float)
    lc = gammaln(n + 1) - gammaln(k + 1) - gammaln(n - k + 1)
    if rho <= 0:
        with np.errstate(divide="ignore", invalid="ignore"):
            out = lc + np.where(k > 0, k * np.log(p), 0.0) + \
                np.where(n - k > 0, (n - k) * np.log1p(-p), 0.0)
        return out
    a = p * (1.0 - rho) / rho
    b = (1.0 - p) * (1.0 - rho) / rho
    return lc + betaln(k + a, n - k + b) - betaln(a, b)


def site_loglik(counts, params, rho=DEFAULT_RHO):
    """Natural-log likelihoods, shape (N, 2): column 0 = 0/1, column 1 = 1/1.

    counts: array (N, 4) of clip, fmap, disc, conc."""
    c = np.asarray(counts, float).reshape(-1, 4)
    n_sr = c[:, 0] + c[:, 1]
    n_pe = c[:, 2] + c[:, 3]
    out = np.empty((len(c), 2))
    for j, g in enumerate(("het", "hom")):
        out[:, j] = (betabinom_logpmf(c[:, 0], n_sr, params["sr_" + g], rho)
                     + betabinom_logpmf(c[:, 2], n_pe, params["pe_" + g], rho))
    return out


def _posterior_hom(ll, hom_prior):
    lp_het = ll[:, 0] + math.log(1.0 - hom_prior)
    lp_hom = ll[:, 1] + math.log(hom_prior)
    return np.exp(lp_hom - np.logaddexp(lp_het, lp_hom))


# ------------------------------------------------------------------------------------------
# EM
# ------------------------------------------------------------------------------------------
def _clamp(x, lo_hi):
    return float(min(max(x, lo_hi[0]), lo_hi[1]))


def run_em(counts, rho=DEFAULT_RHO, init=None, max_iter=EM_MAX_ITER, tol=EM_TOL):
    """Unconstrained-start EM for the two-genotype mixture. Returns (params, n_iter, loglik).

    Bias estimates are clamped into HET_RANGE / HOM_RANGE each M-step and the prior into
    PRIOR_RANGE; whether the result is usable is decided by fit_params."""
    c = np.asarray(counts, float).reshape(-1, 4)
    P = dict(init or DEFAULT_PARAMS)
    prev = -np.inf
    ll_tot = prev
    it = 0
    for it in range(1, max_iter + 1):
        ll = site_loglik(c, P, rho)
        lp = ll + np.log([1.0 - P["hom_prior"], P["hom_prior"]])
        z = np.logaddexp(lp[:, 0], lp[:, 1])
        ll_tot = float(z.sum())
        r = np.exp(lp[:, 1] - z)                     # responsibility for 1/1
        newP = {"hom_prior": _clamp(r.mean(), PRIOR_RANGE)}
        for t, (ia, ib) in (("sr", (0, 1)), ("pe", (2, 3))):
            for g, w in (("het", 1.0 - r), ("hom", r)):
                den = (w * (c[:, ia] + c[:, ib])).sum()
                est = (w * c[:, ia]).sum() / den if den > 0 else P[t + "_" + g]
                newP[t + "_" + g] = _clamp(est, HET_RANGE if g == "het" else HOM_RANGE)
        delta = max(abs(newP[k] - P[k]) for k in P)
        P = newP
        if delta < tol and abs(ll_tot - prev) < 1e-6 * max(1.0, abs(ll_tot)):
            break
        prev = ll_tot
    return P, it, ll_tot


def _degenerate_reason(P):
    for t in ("sr", "pe"):
        h, o = P[t + "_het"], P[t + "_hom"]
        if h >= o - MIN_SEP:
            return "%s_het>=%s_hom-%.2f" % (t, t, MIN_SEP)
        if h <= HET_RANGE[0] or h >= HET_RANGE[1]:
            return "%s_het_at_bound" % t
        if o <= HOM_RANGE[0]:
            return "%s_hom_at_bound" % t
    if P["hom_prior"] <= PRIOR_RANGE[0] or P["hom_prior"] >= PRIOR_RANGE[1]:
        return "hom_prior_at_bound"
    return None


def fit_params(counts, rho=DEFAULT_RHO, min_informative=MIN_INFORMATIVE):
    """Per-sample fit. Returns a dict with the four biases, hom_prior, rho, source
    ("em" or "default"), reason (why defaults were used, or "."), n_informative, n_iter.

    A call is informative when it has at least one SR and one PE read. EM is started from
    the HG002 defaults and from two other starts; the highest-likelihood non-degenerate
    solution wins."""
    c = np.asarray(counts, float).reshape(-1, 4)
    inf = (c[:, 0] + c[:, 1] > 0) & (c[:, 2] + c[:, 3] > 0)
    n_inf = int(inf.sum())
    out = dict(DEFAULT_PARAMS, rho=rho, gq_err_floor=GQ_ERR_FLOOR, source="default", reason=".",
               n_informative=n_inf, n_iter=0)
    if n_inf < min_informative:
        out["reason"] = "n_informative<%d" % min_informative
        return out
    starts = [DEFAULT_PARAMS,
              {"sr_het": 0.45, "sr_hom": 0.85, "pe_het": 0.30, "pe_hom": 0.65, "hom_prior": 0.2},
              {"sr_het": 0.55, "sr_hom": 0.99, "pe_het": 0.45, "pe_hom": 0.90, "hom_prior": 0.4}]
    best, reason = None, None
    for s in starts:
        P, it, ll = run_em(c[inf], rho, init=s)
        why = _degenerate_reason(P)
        if why is None and (best is None or ll > best[2]):
            best = (P, it, ll)
        elif why is not None and reason is None:
            reason = why
    if best is None:
        out["reason"] = "degenerate_em:" + reason
        return out
    out.update(best[0])
    out.update(source="em", n_iter=best[1])
    return out


# ------------------------------------------------------------------------------------------
# per-site outputs
# ------------------------------------------------------------------------------------------
def _bias_map(f, b_het, b_hom):
    """Piecewise-linear map of a raw alt fraction through (0,0), (b_het,0.5), (b_hom,1)."""
    f = np.asarray(f, float)
    b_het = np.asarray(b_het, float)
    b_hom = np.asarray(b_hom, float)
    lo = 0.5 * f / b_het
    hi = 0.5 + 0.5 * (f - b_het) / (b_hom - b_het)
    return np.clip(np.where(f <= b_het, lo, hi), 0.0, 1.0)


def gq_from_error(perr, err_floor=GQ_ERR_FLOOR, gq_cap=GQ_CAP):
    """GQ (int array) from the model's posterior probability that the called GT is wrong.

    perr' = err_floor + (1 - 2 err_floor) * perr: with probability 2 err_floor the counts are
    taken to say nothing about the genotype (a 50:50 call), so GQ <= -10 log10(err_floor). The
    map is monotone and keeps perr' < 0.5 whenever perr < 0.5, so the GT call is unchanged."""
    p = err_floor + (1.0 - 2.0 * err_floor) * np.asarray(perr, float)
    with np.errstate(divide="ignore"):
        gq = np.minimum(gq_cap, -10.0 * np.log10(np.maximum(p, 1e-300)))
    return np.floor(gq + 1e-9).astype(int)


def genotype(counts, params, gq_cap=GQ_CAP):
    """Genotype every call. Returns a dict of arrays: gt (str), gq (int), gl (N,2 log10; the GTL field),
    af_sr, af_pe, vaf, vaf_lo, vaf_hi (NaN where undefined), post_hom."""
    c = np.asarray(counts, float).reshape(-1, 4)
    rho = params.get("rho", DEFAULT_RHO)
    ll = site_loglik(c, params, rho)
    post = _posterior_hom(ll, params["hom_prior"])
    n_sr = c[:, 0] + c[:, 1]
    n_pe = c[:, 2] + c[:, 3]
    none = (n_sr + n_pe) == 0
    is_hom = post >= 0.5
    perr = np.where(is_hom, 1.0 - post, post)
    gq = gq_from_error(perr, params.get("gq_err_floor", GQ_ERR_FLOOR), gq_cap)
    gq = np.where(none, 0, gq).astype(int)
    gt = np.where(none, "./.", np.where(is_hom, "1/1", "0/1"))

    with np.errstate(divide="ignore", invalid="ignore"):
        af_sr = np.where(n_sr > 0, c[:, 0] / n_sr, np.nan)
        af_pe = np.where(n_pe > 0, c[:, 2] / n_pe, np.nan)
    v_sr = _bias_map(np.nan_to_num(af_sr), params["sr_het"], params["sr_hom"])
    v_pe = _bias_map(np.nan_to_num(af_pe), params["pe_het"], params["pe_hom"])
    n_all = n_sr + n_pe
    with np.errstate(divide="ignore", invalid="ignore"):
        vaf = np.where(none, np.nan, (n_sr * v_sr + n_pe * v_pe) / n_all)
        w_sr = np.where(none, 0.5, n_sr / n_all)
    b_het = w_sr * params["sr_het"] + (1 - w_sr) * params["pe_het"]
    b_hom = w_sr * params["sr_hom"] + (1 - w_sr) * params["pe_hom"]
    alt = c[:, 0] + c[:, 2]
    ref = c[:, 1] + c[:, 3]
    q_lo = betaincinv(alt + 0.5, ref + 0.5, 0.025)
    q_hi = betaincinv(alt + 0.5, ref + 0.5, 0.975)
    vaf_lo = np.where(none, np.nan, _bias_map(q_lo, b_het, b_hom))
    vaf_hi = np.where(none, np.nan, _bias_map(q_hi, b_het, b_hom))
    return {"gt": gt, "gq": gq, "gl": ll / math.log(10.0), "af_sr": af_sr, "af_pe": af_pe,
            "vaf": vaf, "vaf_lo": vaf_lo, "vaf_hi": vaf_hi, "post_hom": post}


# ------------------------------------------------------------------------------------------
# formatting
# ------------------------------------------------------------------------------------------
FORMAT_KEYS = ("GT", "GQ", "GTL", "AF_SR", "AF_PE", "VAF", "VAF_CI")

FORMAT_HEADER_LINES = (
    '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype (0/1 or 1/1; ./. when the call has no split-read or read-pair counts)">',
    '##FORMAT=<ID=GQ,Number=1,Type=Integer,Description="Genotype quality: phred-scaled posterior probability that GT is wrong, with a per-call error floor (gq_err_floor in ##xtea2_genotype_params; 0.001 gives GQ <= 30). Conditional on the site carrying the insertion: it compares 0/1 with 1/1 only; 0/0 is not modeled">',
    '##FORMAT=<ID=GTL,Number=2,Type=Float,Description="log10 likelihoods of the two modeled genotypes, in the order 0/1,1/1, under the beta-binomial genotype model (see ##xtea2_genotype_params). Not the VCF GL tag (Number=G): 0/0 is not modeled">',
    '##FORMAT=<ID=AF_SR,Number=1,Type=Float,Description="Raw split-read alt fraction AF_CLIP/(AF_CLIP+AF_FMAP)">',
    '##FORMAT=<ID=AF_PE,Number=1,Type=Float,Description="Raw read-pair alt fraction AF_DISC/(AF_DISC+AF_CONCORDNT)">',
    '##FORMAT=<ID=VAF,Number=1,Type=Float,Description="Allele fraction: maximum-likelihood estimate from AF_CLIP/AF_FMAP and AF_DISC/AF_CONCORDNT under a per-site read-yield model (read length, insert size, SVLEN, TSD; see ##xtea2_vaf_params). Without ##xtea2_vaf_params: AF_SR and AF_PE mapped through the fitted het/hom bias lines and averaged by depth">',
    '##FORMAT=<ID=VAF_CI,Number=2,Type=Float,Description="95% interval for VAF: profile-likelihood interval under the read-yield model (lower bound 0 when the data cannot exclude 0). Without ##xtea2_vaf_params: Beta(alt+0.5,ref+0.5) quantiles mapped through the bias line">',
)
PARAMS_HEADER_KEY = "##xtea2_genotype_params="
VAF_HEADER_KEY = xvs.VAF_HEADER_KEY


def params_header_line(P):
    """One ##xtea2_genotype_params=<...> line (no trailing newline)."""
    return (PARAMS_HEADER_KEY + "<model=betabinom,source={source},reason={reason},"
            "n_informative={n_informative},rho={rho:g},sr_het={sr_het:.4f},sr_hom={sr_hom:.4f},"
            "pe_het={pe_het:.4f},pe_hom={pe_hom:.4f},hom_prior={hom_prior:.4f},"
            "gq_cap={cap},gq_err_floor={floor:g}>").format(
                cap=GQ_CAP, floor=P.get("gq_err_floor", GQ_ERR_FLOOR), **P)


def parse_params_header_line(line):
    """Inverse of params_header_line (numbers back to float/int)."""
    body = line.strip()[len(PARAMS_HEADER_KEY):].strip("<>")
    out = {}
    for kv in body.split(","):
        k, v = kv.split("=", 1)
        try:
            out[k] = int(v) if k in ("n_informative", "gq_cap") else float(v)
        except ValueError:
            out[k] = v
    return out


def _fmt(x, nd=3):
    return "." if x is None or (isinstance(x, float) and math.isnan(x)) else "%.*f" % (nd, x)


def format_values(res, i):
    """FORMAT values for call i, in FORMAT_KEYS order, as strings. A mosaic result (gq None)
    has GT ./. and missing GQ/GTL but still carries AF_SR, AF_PE, VAF and VAF_CI."""
    tail = [_fmt(float(res["af_sr"][i])), _fmt(float(res["af_pe"][i])),
            _fmt(float(res["vaf"][i])),
            "%s,%s" % (_fmt(float(res["vaf_lo"][i])), _fmt(float(res["vaf_hi"][i])))]
    if res.get("gq") is None:
        return ["./.", ".", ".,."] + tail
    if res["gt"][i] == "./.":
        return ["./.", "0", ".,.", ".", ".", ".", ".,."]
    gl = res["gl"][i]
    return [str(res["gt"][i]), str(int(res["gq"][i])), "%.2f,%.2f" % (gl[0], gl[1])] + tail


def mosaic_result(counts):
    """A genotype()-shaped result with no genotype: GT ./., GQ/GTL missing, raw AFs, VAF unset
    (filled by with_site_vaf)."""
    c = np.asarray(counts, float).reshape(-1, 4)
    n_sr = c[:, 0] + c[:, 1]
    n_pe = c[:, 2] + c[:, 3]
    nan = np.full(len(c), np.nan)
    with np.errstate(divide="ignore", invalid="ignore"):
        af_sr = np.where(n_sr > 0, c[:, 0] / n_sr, np.nan)
        af_pe = np.where(n_pe > 0, c[:, 2] / n_pe, np.nan)
    return {"gt": np.full(len(c), "./."), "gq": None, "gl": None, "af_sr": af_sr,
            "af_pe": af_pe, "vaf": nan, "vaf_lo": nan.copy(), "vaf_hi": nan.copy()}


def with_site_vaf(res, counts, lib, S, t, mosaic=False, calib=None, background=None,
                  rho=None, min_informative=MIN_INFORMATIVE):
    """Replace res's VAF / VAF_CI with x_vaf_site's per-site estimate. Returns the
    ##xtea2_vaf_params dict, or None (and leaves res unchanged) when lib is None."""
    if lib is None:
        return None
    vaf, lo, hi, VP = xvs.estimate(counts, lib, S, t, mosaic=mosaic, calib=calib,
                                   background=background, rho=rho,
                                   min_informative=min_informative)
    res["vaf"], res["vaf_lo"], res["vaf_hi"] = vaf, lo, hi
    return VP


# ------------------------------------------------------------------------------------------
# pipeline step: raw xTea result rows (fields 35/36/39/40) -> rows + 7 genotype columns
# ------------------------------------------------------------------------------------------
RAW_IDX = (35, 36, 39, 40)          # AF_CLIP, AF_FMAP, AF_DISC, AF_CONCORDNT in a raw row
PARAMS_SIDECAR_SUFFIX = ".genotype_params"


def _read_raw(sf):
    with open(sf) as fin:
        return [ln.rstrip("\n") for ln in fin if ln.strip()]


def genotype_raw_file(sf_in, sf_out, rho=DEFAULT_RHO, basic_info=None, type_map=None,
                      rtype=None, mosaic=False, calib_files=None, background=None):
    """The --gntp_classify step. Fits once over every row of sf_in (one sample, all
    families), appends the FORMAT_KEYS values as 7 tab-separated columns (fields 50-56, read
    back by x_gvcf), and writes the params header lines to sf_out + PARAMS_SIDECAR_SUFFIX.

    basic_info   X_BasicInfo's basic_cov_is_rlth_info.txt (read length, insert size); without
                 it VAF / VAF_CI use the old bias line
    type_map     locus_table type map (merged path), or rtype (per-type path): the family,
                 for the default insertion length when SVLEN is missing
    mosaic       no genotype (GT ./.); VAF with background, calibrated on calib_files
    calib_files  [(raw file, rtype)] of the sample's germline calls (mosaic calibration)
    background   {"b_SR", "b_PE"} for mosaic; default x_vaf_site.MOSAIC_BACKGROUND"""
    lines = _read_raw(sf_in)
    rows = [ln.split("\t") for ln in lines]
    counts = np.array([[int(f[i]) for i in RAW_IDX] for f in rows], float).reshape(-1, 4)
    lib = xvs.read_basic_info(basic_info) if basic_info else None
    if basic_info and lib is None:
        print("[x_genotype_lik] WARNING: no read length / insert size in %s; VAF uses the "
              "bias line" % basic_info)
    tmap = xvs.load_type_map(type_map) if type_map else None
    VP = None
    if mosaic:
        P, res = None, mosaic_result(counts)
    else:
        P = fit_params(counts, rho)
        res = genotype(counts, P)
    if lib is not None:
        _, S, t = xvs.raw_site_inputs(rows, rtype, tmap)
        calib = None
        if calib_files:
            parts = []
            for sf, rt in calib_files:
                if os.path.isfile(sf):
                    parts.append(xvs.raw_site_inputs([ln.split("\t") for ln in _read_raw(sf)],
                                                     rt))
            if parts:
                calib = tuple(np.concatenate([p[k] for p in parts]) for k in range(3))
        VP = with_site_vaf(res, counts, lib, S, t, mosaic=mosaic, calib=calib,
                           background=background)
    elif mosaic:
        print("[x_genotype_lik] WARNING: mosaic VAF needs the read length and insert size "
              "(basic_info); VAF left empty")
    with open(sf_out, "w") as fout:
        for i, ln in enumerate(lines):
            fout.write(ln + "\t" + "\t".join(format_values(res, i)) + "\n")
    with open(sf_out + PARAMS_SIDECAR_SUFFIX, "w") as fp:
        if P is not None:
            fp.write(params_header_line(P) + "\n")
        if VP is not None:
            fp.write(xvs.vaf_params_header_line(VP) + "\n")
    if P is not None:
        print("[x_genotype_lik] %d calls, %s fit (%s): sr %.3f/%.3f pe %.3f/%.3f hom_prior %.3f"
              % (len(lines), P["source"], P["reason"], P["sr_het"], P["sr_hom"], P["pe_het"],
                 P["pe_hom"], P["hom_prior"]))
    if VP is not None:
        print("[x_genotype_lik] VAF (%s, %s calibration, %d calls): %s"
              % (VP["mode"], VP["source"], VP["n_calib"], xvs.vaf_params_header_line(VP)))
    return P


# ------------------------------------------------------------------------------------------
# post-hoc CLI: re-genotype an existing xTea2 VCF
# ------------------------------------------------------------------------------------------
INFO_KEYS = ("AF_CLIP", "AF_FMAP", "AF_DISC", "AF_CONCORDNT")


def _info_dict(s):
    return dict(kv.split("=", 1) for kv in s.split(";") if "=" in kv)


def _open_text(sf):
    """Plain or gzip/bgzip-compressed text, by magic bytes."""
    with open(sf, "rb") as fb:
        gz = fb.read(2) == b"\x1f\x8b"
    return gzip.open(sf, "rt") if gz else open(sf)


DEFAULT_SAMPLE = "SAMPLE"
# FORMAT header lines dropped on re-genotyping: the keys written here plus GL, which an early
# build of this genotyper wrote instead of GTL.
_STALE_FORMAT_IDS = FORMAT_KEYS + ("GL",)


ALT_FAMILY = {"INS:ME:ALU": "ALU", "INS:ME:LINE1": "LINE1", "INS:ME:SVA": "SVA",
              "INS:ME:HERV-K": "HERV"}


def library_from_header(head):
    """x_vaf_site.Library from an input VCF's ##xtea2_vaf_params line, or None."""
    for h in head:
        if h.startswith(VAF_HEADER_KEY):
            V = xvs.parse_vaf_params_header_line(h)
            try:
                return xvs.Library(V["read_len"], V["mean_is"], V["sd_is"])
            except (KeyError, TypeError, ValueError):
                return None
    return None


def regenotype_vcf(sf_in, sf_out, rho=DEFAULT_RHO, min_informative=MIN_INFORMATIVE,
                   sample=None, lib=None, mosaic=False, background=None, vaf_rho=None):
    """Re-genotype every record of a single-sample xTea2 VCF (plain or .gz) from its INFO
    counts. Output: the 8 fixed columns, FORMAT and one sample column. The sample keeps its
    name unless `sample` is given; a sites-only input gets `sample` or DEFAULT_SAMPLE. Old
    FORMAT lines for the keys written here and any old params lines are replaced; the old
    FORMAT/sample values are discarded. Multi-sample input raises ValueError.

    VAF uses x_vaf_site when the read length and insert size are known: `lib`, else the
    input's ##xtea2_vaf_params line; otherwise the bias line. mosaic=True writes no genotype.
    Returns the genotype params (None in mosaic mode)."""
    head, recs, chrom_line = [], [], None
    with _open_text(sf_in) as fin:
        for line in fin:
            line = line.rstrip("\n")
            if line.startswith("##"):
                head.append(line)
            elif line.startswith("#"):
                chrom_line = line.split("\t")
            elif line.strip():
                recs.append(line.split("\t"))
    if chrom_line is None:
        raise ValueError("%s has no #CHROM header line; not a VCF?" % sf_in)
    samples = chrom_line[9:]
    if len(samples) > 1:
        raise ValueError("%s has %d samples (%s); xtea2-genotype fits one sample at a time. "
                         "Split it first, e.g. bcftools view -s NAME."
                         % (sf_in, len(samples), ",".join(samples)))
    s_name = sample or (samples[0] if samples else DEFAULT_SAMPLE)

    counts = np.zeros((len(recs), 4))
    for j, f in enumerate(recs):
        if len(f) < 8:
            raise ValueError("%s: record %d has %d columns, expected at least 8"
                             % (sf_in, j + 1, len(f)))
        info = _info_dict(f[7])
        missing = [k for k in INFO_KEYS if k not in info]
        if missing:
            raise ValueError("%s:%s lacks INFO %s; not an xTea2 VCF?" % (f[0], f[1], missing))
        counts[j] = [int(info[k]) for k in INFO_KEYS]
    if mosaic:
        P, res = None, mosaic_result(counts)
    else:
        P = fit_params(counts, rho, min_informative)
        res = genotype(counts, P)
    if lib is None:
        lib = library_from_header(head)
    VP = None
    if lib is not None:
        S, t = np.zeros(len(recs)), np.zeros(len(recs))
        for j, f in enumerate(recs):
            info = _info_dict(f[7])
            fam = ALT_FAMILY.get(f[4].strip("<>"))
            S[j] = xvs.parse_svlen(info.get("SVLEN"), fam)
            t[j] = xvs.parse_tsdlen(info.get("TSDLEN"))
        VP = with_site_vaf(res, counts, lib, S, t, mosaic=mosaic, background=background,
                           rho=vaf_rho, min_informative=min_informative)
    elif mosaic:
        sys.stderr.write("[xtea2-genotype] WARNING: no read length / insert size (pass "
                         "--basic-info or --read-length/--insert-mean/--insert-sd); VAF left "
                         "empty\n")
    drop = (tuple("##FORMAT=<ID=%s," % k for k in _STALE_FORMAT_IDS)
            + (PARAMS_HEADER_KEY, VAF_HEADER_KEY))
    new_head = [h for h in head if not h.startswith(drop)]
    new_head += list(FORMAT_HEADER_LINES)
    if P is not None:
        new_head.append(params_header_line(P))
    if VP is not None:
        new_head.append(xvs.vaf_params_header_line(VP))
    new_head.append("\t".join(["#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO",
                               "FORMAT", s_name]))
    with open(sf_out, "w") as fout:
        for h in new_head:
            fout.write(h + "\n")
        for j, f in enumerate(recs):
            fout.write("\t".join(f[:8] + [":".join(FORMAT_KEYS),
                                          ":".join(format_values(res, j))]) + "\n")
    return P


def main(argv=None):
    ap = argparse.ArgumentParser(
        prog="xtea2-genotype",
        description="Re-genotype an xTea2 VCF from its INFO read counts (AF_CLIP, AF_FMAP, "
                    "AF_DISC, AF_CONCORDNT) with the per-sample EM-calibrated beta-binomial "
                    "genotyper. Fit over all records in the file: pass the merged VCF.")
    ap.add_argument("vcf", help="input xTea2 VCF, plain or .gz (one sample, or sites-only)")
    ap.add_argument("-o", "--output", required=True, help="output VCF (uncompressed)")
    ap.add_argument("--sample", default=None,
                    help="sample column name (default: the input's sample, or %s for a "
                         "sites-only input)" % DEFAULT_SAMPLE)
    ap.add_argument("--rho", type=float, default=DEFAULT_RHO,
                    help="beta-binomial overdispersion (default %(default)s)")
    ap.add_argument("--min-informative", type=int, default=MIN_INFORMATIVE,
                    help="fewer informative calls than this -> HG002 default parameters "
                         "(default %(default)s)")
    ap.add_argument("--basic-info", default=None,
                    help="xTea2 tmp_basic/basic_cov_is_rlth_info.txt of the sample (read length "
                         "and insert size, for VAF). Default: the input's ##xtea2_vaf_params")
    ap.add_argument("--read-length", type=float, default=None, help="read length (with "
                    "--insert-mean and --insert-sd, instead of --basic-info)")
    ap.add_argument("--insert-mean", type=float, default=None)
    ap.add_argument("--insert-sd", type=float, default=None)
    ap.add_argument("--mosaic", action="store_true",
                    help="no genotype (GT ./.): VAF and VAF_CI only, with a background term")
    ap.add_argument("--vaf-background", default=None, metavar="B_SR,B_PE|FILE",
                    help="background alt odds for --mosaic, or a file from x_vaf_background.py "
                         "(default %s,%s)"
                         % (xvs.MOSAIC_BACKGROUND["b_SR"], xvs.MOSAIC_BACKGROUND["b_PE"]))
    a = ap.parse_args(argv)
    lib = None
    if a.basic_info:
        lib = xvs.read_basic_info(a.basic_info)
        if lib is None:
            ap.error("no read length / insert size in %s" % a.basic_info)
    elif a.read_length is not None:
        if a.insert_mean is None or a.insert_sd is None:
            ap.error("--read-length needs --insert-mean and --insert-sd")
        lib = xvs.Library(a.read_length, a.insert_mean, a.insert_sd)
    P = regenotype_vcf(a.vcf, a.output, a.rho, a.min_informative, a.sample, lib=lib,
                       mosaic=a.mosaic, background=parse_background(a.vaf_background))
    if P is not None:
        sys.stderr.write(params_header_line(P) + "\n")
    return 0


def parse_background(s):
    """"b_SR,b_PE", or a file whose first line is that (x_vaf_background --out) -> dict with a
    "source" label, or None."""
    if not s:
        return None
    src = "user"
    if os.path.isfile(s):
        with open(s) as fh:
            s, src = fh.readline().strip(), "file"
    b = [float(x) for x in s.split(",")]
    if len(b) != 2 or min(b) < 0:
        raise ValueError("--vaf-background takes two non-negative numbers, b_SR,b_PE")
    return {"b_SR": b[0], "b_PE": b[1], "source": src}


if __name__ == "__main__":
    sys.exit(main())
