"""Per-site VAF for MEI calls from short-read counts, from the expected read yield of each allele.

The four counts xTea writes for each call (x_genotype_feature) are

    SR  alt = AF_CLIP       reads soft-clipped at the dominant breakpoint
        ref = AF_FMAP       reads fully mapped across the breakpoint (one M), minus the 44-bp
                            window x_genotype_feature leaves out on the non-dominant side
    PE  alt = AF_DISC       reads within DFT_IS whose mate maps > 2 kb away or to another chrom
        ref = AF_CONCORDNT  reads of near pairs whose span covers the breakpoint

A het insertion does not give 50% alt reads, and how far from 50% depends on the site. A short
insertion gives fewer discordant pairs, and pairs that span it are counted as concordant; a
short insert size relative to the read length turns discordant pairs into clipped reads. So each
read type gets an expected number of alt-allele and ref-allele reads per unit fragment-start
density. These come from the read length R, the insert-size distribution P(F), the insertion
length S and the TSD length t.

Geometry
--------
For one junction (the other is its mirror), take the alt haplotype with the left flank ending at
0, the insertion at [0, S) and the right flank starting at S. A read whose end reaches more than c
bases into the insertion is soft-clipped, so long as at least a = R - K of its bases are in the
flank. K is the longest clip a flank-primary read can carry; with fewer flank bases the read maps
to a TE copy and is the "TE-primary" mate of a discordant pair. For a fragment of length F, per
unit fragment-start density:

    disc   n_D(F) = max(0, min(F - R, S + 2a - R))     anchor in the flank, mate TE-primary
    span   n_S(F) = max(0, F - 2R - S + 2c)            one read anchored on each side
    ref    n_C(F) = F, folded at DFT_IS                fragments covering the breakpoint

Expected reads per haplotype (both junctions for PE, the dominant one for SR):

    D_alt = 2 E[n_D]        C_alt = 2 E[n_S] + l_PE + 2t        C_ref = 2 E[n_C]
    A_SR  = 2 (K - c)       B_SR  = 2 (R - 44)                  L_SR  = l_SR

l_PE counts concordant reads near an alt allele that the geometry above misses (pairs whose read
ends just past the breakpoint, hard-clipped supplementary records of TE-primary reads). 2t counts
pairs that start inside the duplicated TSD, and l_SR counts fully mapped reads at an alt allele.
For an allele fraction v

    p_SR(v) = (v e_SR A_SR + b_SR B_SR) / (v e_SR A_SR + b_SR B_SR + (1 - v) B_SR + v L_SR)
    p_PE(v) = (v e_PE D_alt + b_PE C_ref) / (v e_PE D_alt + b_PE C_ref + v C_alt + (1 - v) C_ref)

e_SR and e_PE are the sample's detection efficiencies relative to the geometry. b_SR and b_PE are
background alt reads (clips and chimeric or far-mate pairs found near any position) relative to
the ref yield. VAF is the maximum-likelihood v under

    log L(v) = log BetaBin(clip | clip + fmap, p_SR(v), rho)
             + log BetaBin(disc | disc + conc, p_PE(v), rho)

with a 95% profile-likelihood interval (log L within 1.92 of the maximum). No genotype is needed,
so the same code serves germline and mosaic calls.

Calibration (fit_calibration) fits e_SR, e_PE, l_SR, l_PE and the hom fraction w by maximum
likelihood over a sample's calls, as an unlabeled mixture of v = 0.5 and v = 1. Background is not
part of that fit: a germline het or hom carries far more alt reads than background, so b = 0 there.
For mosaic calls b matters: at 1% VAF it is about as large as the signal. It is fitted at v = 0 on
counts at 1000 random autosomal clip peaks of the sample's BAM (fit_background, counted by
x_vaf_background) and scaled by CONTEXT_FACTOR, the excess background of MEI call sites over random
positions. Without a BAM, MOSAIC_BACKGROUND (the 300x BCM HapMap mixture) is used.

Dependencies: numpy and scipy only.
"""
import math

import numpy as np
from scipy.optimize import minimize
from scipy.special import betaln, gammaln

BWA_HALF_READ_MIN_SCORE = 45      # x_genotype_feature fmap window (|start - pos| < 45)
DEFAULT_RHO = 0.05                # germline; coverage of the 95% interval 0.91-0.98 on 6 genomes
MOSAIC_RHO = 0.02                 # mosaic; 0.05 biases low-VAF estimates upward (see README_NOTES)
DEFAULT_TSD = 14                  # median TSD length of xTea2 calls (all families)
MAX_SVLEN = 20000                 # longer "insertion lengths" (orphan transductions) are not lengths
MIN_INFORMATIVE = 50
# Shortest clip that counts (c) and longest clip a flank-primary read can carry (K). K is set by the
# aligner, not the read length: the clip-length plateau gave K = 63-81 at R = 151 and R = 250.
CLIP_C = 2.0
CLIP_K = 70.0
# HG002 30x NovaSeq PCR-free, mixture fit over its 1568 truth-matched calls with the pipeline's
# insert-size estimate (normal, mean 419, sd 178), R = 151, K = 70, rho = 0.05. Used when a sample cannot support its own
# fit. On five other genomes the fitted values range e_SR 1.0-1.4, e_PE 0.9-1.5, l_PE 23-145.
DEFAULT_CAL = {"e_SR": 1.21, "e_PE": 1.33, "l_SR": 8.5, "l_PE": 137.0, "b_SR": 0.0, "b_PE": 0.0}
DEFAULT_HOM_FRAC = 0.27
FIT_KEYS = ("e_SR", "e_PE", "l_SR", "l_PE")
FAMILY_DEFAULT_LEN = {"ALU": 300, "LINE1": 6000, "L1": 6000, "SVA": 2000, "HERV": 7000}
GENERIC_DEFAULT_LEN = 300

# VAF grid: log-spaced below 0.1, linear above.
V_GRID = np.unique(np.concatenate([np.logspace(-4, -1, 301), np.linspace(0.1, 1.0, 901)]))
CHI2_1_95_HALF = 1.920729


# ------------------------------------------------------------------------------------------
# sample geometry
# ------------------------------------------------------------------------------------------
class Library:
    """Read length and insert-size distribution of one sample."""

    def __init__(self, read_len, mean_is, sd_is):
        self.R = float(read_len)
        self.mean_is = float(mean_is)
        self.sd_is = max(float(sd_is), 1.0)
        self.dft_is = self.mean_is + 3.0 * self.sd_is          # x_genotype_feature window
        self.F, self.w = normal_isize(self.mean_is, self.sd_is)

    def __repr__(self):
        return "Library(R=%g, mean_is=%.1f, sd_is=%.1f)" % (self.R, self.mean_is, self.sd_is)


def read_basic_info(sf):
    """Library from X_BasicInfo's basic_cov_is_rlth_info.txt (bam, cov, rlth, mean_is, sd); the
    read-count-weighted mean when the file has one line per BAM. None if unusable."""
    rows = []
    try:
        with open(sf) as fin:
            for ln in fin:
                f = ln.split()
                if len(f) >= 5:
                    rows.append([float(x) for x in f[1:5]])
    except (OSError, ValueError):
        return None
    rows = [r for r in rows if r[1] > 0 and r[2] > 0 and r[3] > 0]
    if not rows:
        return None
    a = np.array(rows)
    wt = a[:, 0] / a[:, 0].sum() if a[:, 0].sum() > 0 else np.full(len(a), 1.0 / len(a))
    return Library(*(wt @ a[:, 1:4]))


def normal_isize(mean, sd, step=2.0):
    """Discretized normal insert-size pmf (support, probs), truncated at 4 sd and at 1 bp."""
    sd = max(float(sd), 1.0)
    f = np.arange(max(1.0, mean - 4 * sd), mean + 4 * sd + step, step)
    w = np.exp(-0.5 * ((f - mean) / sd) ** 2)
    return f, w / w.sum()


# ------------------------------------------------------------------------------------------
# yields and expected alt fractions
# ------------------------------------------------------------------------------------------
def yields(lib, S, t, cal):
    """Expected reads per unit fragment-start density, per haplotype, for arrays S (insertion
    length) and t (TSD length). Returns a dict of arrays A_SR, B_SR, L_SR, D_alt, C_alt, C_ref."""
    S = np.atleast_1d(np.asarray(S, float))
    t = np.clip(np.atleast_1d(np.asarray(t, float)), 0, None)
    R, F, w = lib.R, lib.F[None, :], lib.w
    K = min(CLIP_K, R - 1.0)
    a = R - K
    n_d = np.maximum(0.0, np.minimum(F - R, S[:, None] + 2 * a - R))
    n_s = np.maximum(0.0, F - 2 * R - S[:, None] + 2 * CLIP_C)
    n_c = np.where(lib.F <= lib.dft_is, lib.F, np.maximum(0.0, 2 * lib.dft_is - lib.F))
    d_alt = 2.0 * (n_d @ w)
    n = len(S)
    return {"A_SR": np.full(n, 2.0 * max(K - CLIP_C, 1.0)),
            "B_SR": np.full(n, 2.0 * max(R - (BWA_HALF_READ_MIN_SCORE - 1), 1.0)),
            "L_SR": np.full(n, float(cal["l_SR"])),
            "D_alt": d_alt,
            "C_alt": 2.0 * (n_s @ w) + cal["l_PE"] + 2.0 * t,
            "C_ref": np.full(n, 2.0 * float(n_c @ w))}


def p_alt(v, y, cal):
    """Expected alt fractions (p_SR, p_PE) at VAF v (scalar or array broadcastable with y)."""
    v = np.asarray(v, float)
    a_sr = v * cal["e_SR"] * y["A_SR"] + cal["b_SR"] * y["B_SR"]
    p_sr = a_sr / (a_sr + (1 - v) * y["B_SR"] + v * y["L_SR"])
    a_pe = v * cal["e_PE"] * y["D_alt"] + cal["b_PE"] * y["C_ref"]
    p_pe = a_pe / (a_pe + v * y["C_alt"] + (1 - v) * y["C_ref"])
    return np.clip(p_sr, 1e-9, 1 - 1e-9), np.clip(p_pe, 1e-9, 1 - 1e-9)


# ------------------------------------------------------------------------------------------
# likelihood
# ------------------------------------------------------------------------------------------
def betabinom_logpmf(k, n, p, rho=DEFAULT_RHO):
    """log P(k | n) for a beta-binomial with mean p and intra-class correlation rho (rho <= 0:
    binomial). Broadcasts; n = 0 gives 0."""
    k = np.asarray(k, float)
    n = np.asarray(n, float)
    p = np.asarray(p, float)
    lc = gammaln(n + 1) - gammaln(k + 1) - gammaln(n - k + 1)
    if rho <= 0:
        return lc + k * np.log(p) + (n - k) * np.log1p(-p)
    s = (1.0 - rho) / rho
    return lc + betaln(k + p * s, n - k + (1 - p) * s) - betaln(p * s, (1 - p) * s)


def loglik(counts, v, y, cal, rho):
    """log L at VAF v for every site; counts (N, 4) clip, fmap, disc, conc."""
    c = np.asarray(counts, float).reshape(-1, 4)
    p_sr, p_pe = p_alt(v, y, cal)
    return (betabinom_logpmf(c[:, 0], c[:, 0] + c[:, 1], p_sr, rho)
            + betabinom_logpmf(c[:, 2], c[:, 2] + c[:, 3], p_pe, rho))


def site_vaf(counts, y, cal, rho=DEFAULT_RHO, grid=V_GRID):
    """MLE VAF and 95% profile interval per site. Returns vaf, lo, hi (NaN with no reads; lo = 0
    when the interval reaches the bottom of the grid)."""
    c = np.asarray(counts, float).reshape(-1, 4)
    N = len(c)
    vaf, lo, hi = np.full(N, np.nan), np.full(N, np.nan), np.full(N, np.nan)
    if N == 0:
        return vaf, lo, hi
    yy = {k: v[:, None] for k, v in y.items()}
    p_sr, p_pe = p_alt(grid[None, :], yy, cal)
    ll = (betabinom_logpmf(c[:, :1], c[:, :1] + c[:, 1:2], p_sr, rho)
          + betabinom_logpmf(c[:, 2:3], c[:, 2:3] + c[:, 3:4], p_pe, rho))
    has = c.sum(1) > 0
    j = np.argmax(ll, 1)
    inside = ll >= (ll[np.arange(N), j] - CHI2_1_95_HALF)[:, None]
    first = np.argmax(inside, 1)
    last = len(grid) - 1 - np.argmax(inside[:, ::-1], 1)
    vaf[has], lo[has], hi[has] = grid[j][has], grid[first][has], grid[last][has]
    lo[has & (first == 0)] = 0.0
    return vaf, lo, hi


# ------------------------------------------------------------------------------------------
# calibration
# ------------------------------------------------------------------------------------------
def fit_calibration(counts, lib, S, t, rho=DEFAULT_RHO, min_informative=MIN_INFORMATIVE,
                    use=None):
    """Per-sample fit of FIT_KEYS and the hom fraction w, maximizing
    sum log[(1 - w) L(0.5) + w L(1)] over the calls in `use` (default: all) that have both SR and
    PE reads. Returns a dict: the calibration, w, source ("fit" or "default"), reason, n_used."""
    c = np.asarray(counts, float).reshape(-1, 4)
    S = np.asarray(S, float)
    t = np.asarray(t, float)
    m = (c[:, 0] + c[:, 1] > 0) & (c[:, 2] + c[:, 3] > 0)
    if use is not None:
        m &= np.asarray(use, bool)
    out = dict(DEFAULT_CAL, w=DEFAULT_HOM_FRAC, source="default", reason=".", n_used=int(m.sum()))
    if m.sum() < min_informative:
        out["reason"] = "n_informative<%d" % min_informative
        return out
    c, S, t = c[m], S[m], t[m]
    y0 = yields(lib, S, t, dict(DEFAULT_CAL, l_PE=0.0))

    def unpack(x):
        cal = dict(DEFAULT_CAL, **{k: float(np.exp(v)) for k, v in zip(FIT_KEYS, x[:-1])})
        y = dict(y0, C_alt=y0["C_alt"] + cal["l_PE"], L_SR=np.full(len(S), cal["l_SR"]))
        return cal, 1.0 / (1.0 + math.exp(-x[-1])), y

    def nll(x):
        cal, w, y = unpack(x)
        return -np.logaddexp(loglik(c, 0.5, y, cal, rho) + math.log(1 - w),
                             loglik(c, 1.0, y, cal, rho) + math.log(w)).sum()

    x0 = np.r_[np.log([DEFAULT_CAL[k] for k in FIT_KEYS]),
               math.log(DEFAULT_HOM_FRAC / (1 - DEFAULT_HOM_FRAC))]
    best = None
    for d in (0.0, 0.5, -0.5):
        start = x0 + np.r_[np.full(len(FIT_KEYS), d), 0.0]
        r = minimize(nll, start, method="Nelder-Mead",
                     options={"maxiter": 4000, "xatol": 1e-3, "fatol": 1e-3})
        if best is None or r.fun < best.fun:
            best = r
    cal, w, _ = unpack(best.x)
    why = _implausible(cal, w)
    if why:
        out["reason"] = "implausible_fit:" + why
        return out
    out.update(cal, w=w, source="fit")
    return out


# A fit outside these ranges means the mixture did not separate (or the call set is not
# germline-like); the defaults are safer. The upper ends are ~3x above what 7 germline genomes
# gave. Low efficiencies are real: the GIAB HG002 2x250 BAM has few discordant pairs (mates of
# repeat reads left unmapped) and fits e_PE 0.14, which recovers het 0.51 / hom 1.0, while the
# default e_PE 1.33 put its hets at 0.25.
_PLAUSIBLE = {"e_SR": (0.05, 4.0), "e_PE": (0.05, 4.0), "l_SR": (0.0, 60.0), "l_PE": (0.0, 600.0)}


def _implausible(cal, w):
    for k, (lo, hi) in _PLAUSIBLE.items():
        if not lo <= cal[k] <= hi:
            return "%s=%.3g" % (k, cal[k])
    if not 0.02 <= w <= 0.95:
        return "w=%.3g" % w
    return None


def background_from_counts(counts):
    """Background (b_SR, b_PE) from counts at control positions without an insertion: the
    beta-binomial means of clip/(clip+fmap) and disc/(disc+conc), as odds (p / (1 - p)), which is
    what p_alt(0) = b / (1 + b) inverts. Positions with no reads are skipped. The pooled mean is
    used, not the median: background is heavy-tailed, and the likelihood needs its mean."""
    c = np.asarray(counts, float).reshape(-1, 4)
    out = {}
    for key, (i, j) in (("b_SR", (0, 1)), ("b_PE", (2, 3))):
        n = c[:, i] + c[:, j]
        k = c[:, i][n > 0]
        n = n[n > 0]
        p = float(k.sum() / n.sum()) if n.sum() > 0 else 0.0
        out[key] = p / (1.0 - p) if p < 1 else 0.0
    return out


# ------------------------------------------------------------------------------------------
# inputs
# ------------------------------------------------------------------------------------------
def parse_svlen(s, family=None):
    """Insertion length from a raw or VCF value; missing, non-positive (xTea writes -1 when it
    has none) or > MAX_SVLEN -> the family default (300 when the family is unknown)."""
    try:
        x = int(float(s))
    except (TypeError, ValueError):
        x = 0
    if x <= 0 or x > MAX_SVLEN:
        x = FAMILY_DEFAULT_LEN.get(str(family).upper(), GENERIC_DEFAULT_LEN)
    return x


def parse_tsdlen(s):
    """TSD length; missing, negative or > 50 -> DEFAULT_TSD."""
    try:
        x = int(float(s))
    except (TypeError, ValueError):
        return DEFAULT_TSD
    return x if 0 <= x <= 50 else DEFAULT_TSD


# ------------------------------------------------------------------------------------------
# raw xTea result rows
# ------------------------------------------------------------------------------------------
RAW_COUNT_IDX = (35, 36, 39, 40)    # AF_CLIP, AF_FMAP, AF_DISC, AF_CONCORDNT
RAW_TSD_IDX = 4
RAW_SVLEN_IDX = 47
RTYPE_FAMILY = {1: "LINE1", 2: "ALU", 4: "SVA", 8: "HERV"}
BASIC_INFO_NAME = "basic_cov_is_rlth_info.txt"


def family_of_rtype(rtype):
    """Family name for an xTea rtype bitmask (lowest set bit), or None."""
    try:
        r = int(rtype)
    except (TypeError, ValueError):
        return None
    for bit in sorted(RTYPE_FAMILY):
        if r & bit:
            return RTYPE_FAMILY[bit]
    return None


def raw_call_position(f):
    """min(lclip_pos, rclip_pos), or the other one when that is <= 0. The same key x_gvcf and
    locus_table.call_position use for the type map."""
    a, b = int(f[2]), int(f[3])
    lo, hi = min(a, b), max(a, b)
    return hi if lo <= 0 else lo


def load_type_map(sf):
    """locus_table type map -> {(chrom, pos): rtype}."""
    m = {}
    with open(sf) as fin:
        for ln in fin:
            f = ln.split()
            if len(f) >= 3:
                m[(f[0], int(f[1]))] = int(f[2])
    return m


def raw_site_inputs(rows, rtype=None, type_map=None):
    """counts (N, 4), insertion length S and TSD length t from split raw rows. The family (for
    the default length) comes from type_map at the call position, else from rtype."""
    counts = np.array([[float(f[i]) for i in RAW_COUNT_IDX] for f in rows], float).reshape(-1, 4)
    S, t = np.zeros(len(rows)), np.zeros(len(rows))
    fam0 = family_of_rtype(rtype)
    for j, f in enumerate(rows):
        fam = fam0
        if type_map:
            r = type_map.get((f[0], raw_call_position(f)))
            fam = family_of_rtype(r) or fam
        S[j] = parse_svlen(f[RAW_SVLEN_IDX] if len(f) > RAW_SVLEN_IDX else None, fam)
        t[j] = parse_tsdlen(f[RAW_TSD_IDX] if len(f) > RAW_TSD_IDX else None)
    return counts, S, t


def germline_like(counts, min_af=0.2):
    """Calls whose raw SR and PE alt fractions are both >= min_af: what a mosaic run can
    calibrate on. A germline het gives ~0.3-0.45 on both."""
    c = np.asarray(counts, float).reshape(-1, 4)
    with np.errstate(divide="ignore", invalid="ignore"):
        sr = c[:, 0] / (c[:, 0] + c[:, 1])
        pe = c[:, 2] / (c[:, 2] + c[:, 3])
    return (sr >= min_af) & (pe >= min_af)


# ------------------------------------------------------------------------------------------
# one entry point for the pipeline and the CLI
# ------------------------------------------------------------------------------------------
# Background for mosaic runs, as odds per ref read (see p_alt). ML fit on the 276 truth sites at
# 0.25-0.5% VAF of the 300x BCM HapMap mixture (rho 0.02, default calibration). Background is
# sample-specific (about 2.5x higher in the UW library), so the pipeline measures it per sample
# at random control positions (x_vaf_background, --vaf-background auto); this constant is the
# fallback when that is not available.
MOSAIC_BACKGROUND = {"b_SR": 0.017, "b_PE": 0.037}

# MEI call sites (in repeats, at clip peaks) carry more background than random autosomal clip
# peaks do. Ratio of the two ML fits on the 300x BCM mixture (MEI truth sites at <= 0.5% VAF vs
# 2000 random positions): b_SR 0.0172 / 0.0092, b_PE 0.0365 / 0.0281.
CONTEXT_FACTOR = {"b_SR": 1.87, "b_PE": 1.30}
CONTROL_SVLEN = 300.0
CONTROL_TSD = float(DEFAULT_TSD)


def fit_background(counts, lib, cal=None, rho=MOSAIC_RHO, v=0.0):
    """ML background odds (b_SR, b_PE) from counts at sites of known VAF v (default 0: control
    positions). SR and PE are fitted separately (each depends only on its own b). The yields use
    an ALU-length insertion and the default TSD; at v = 0 they cancel except through B and C_ref."""
    from scipy.optimize import minimize_scalar
    c = np.asarray(counts, float).reshape(-1, 4)
    cal = dict(DEFAULT_CAL if cal is None else cal)
    n = len(c)
    y = yields(lib, np.full(n, CONTROL_SVLEN), np.full(n, CONTROL_TSD), cal)
    out = {}
    for j, key in enumerate(("b_SR", "b_PE")):
        k, m = c[:, 2 * j], c[:, 2 * j] + c[:, 2 * j + 1]
        if m.sum() == 0:
            out[key] = 0.0
            continue

        def f(lb):
            p = p_alt(v, y, dict(cal, **{key: float(np.exp(lb))}))[j]
            return -betabinom_logpmf(k, m, p, rho).sum()
        out[key] = float(np.exp(minimize_scalar(f, bounds=(-12, 0), method="bounded").x))
    return out


def control_background(counts, lib, cal=None, rho=MOSAIC_RHO):
    """Background for MEI calls from counts at random control positions: fit_background x
    CONTEXT_FACTOR."""
    b = fit_background(counts, lib, cal=cal, rho=rho)
    return {k: b[k] * CONTEXT_FACTOR[k] for k in ("b_SR", "b_PE")}


def estimate(counts, lib, S, t, mosaic=False, calib=None, background=None, rho=None,
             min_informative=MIN_INFORMATIVE, bg_source=None):
    """VAF and 95% interval for every call.

    counts, S, t   the calls to estimate (arrays as in yields / site_vaf)
    calib          (counts, S, t) to calibrate on; default: the calls themselves. Germline mode
                   fits on all of them, mosaic mode only on germline_like ones (a mosaic call
                   set has none, so the pipeline passes the sample's germline calls here)
    background     dict b_SR, b_PE (optional "source" label); default 0 for germline,
                   MOSAIC_BACKGROUND for mosaic
    rho            beta-binomial overdispersion; default DEFAULT_RHO / MOSAIC_RHO
    bg_source      label for the header (e.g. "control:n=1000", "user"); default "none" for
                   germline, "default" for mosaic without a background

    Returns (vaf, lo, hi, params) where params is the dict for vaf_params_header_line."""
    if rho is None:
        rho = MOSAIC_RHO if mosaic else DEFAULT_RHO
    cc, cs, ct = calib if calib is not None else (counts, S, t)
    use = germline_like(cc) if mosaic else None
    fit = fit_calibration(cc, lib, cs, ct, rho=DEFAULT_RHO, min_informative=min_informative,
                          use=use)
    cal = {k: fit[k] for k in FIT_KEYS}
    bg = background if background is not None else (MOSAIC_BACKGROUND if mosaic
                                                     else {"b_SR": 0.0, "b_PE": 0.0})
    if bg_source is None:
        bg_source = (background.get("source", "user") if background is not None
                     else ("default" if mosaic else "none"))
    cal.update(b_SR=float(bg["b_SR"]), b_PE=float(bg["b_PE"]))
    y = yields(lib, S, t, cal)
    vaf, lo, hi = site_vaf(counts, y, cal, rho=rho)
    params = dict(cal, mode="mosaic" if mosaic else "germline", source=fit["source"],
                  reason=fit["reason"], n_calib=fit["n_used"], hom_frac=fit["w"], rho=rho,
                  read_len=lib.R, mean_is=lib.mean_is, sd_is=lib.sd_is, K=min(CLIP_K, lib.R - 1),
                  c=CLIP_C, bg_source=bg_source)
    return vaf, lo, hi, params


VAF_HEADER_KEY = "##xtea2_vaf_params="
_HEADER_FIELDS = ("mode", "source", "reason", "n_calib", "read_len", "mean_is", "sd_is", "K", "c",
                  "rho", "e_SR", "e_PE", "l_SR", "l_PE", "b_SR", "b_PE", "bg_source", "hom_frac")


def vaf_params_header_line(P):
    """##xtea2_vaf_params=<...> (no newline); the read-length and insert-size fields let a later
    xtea2-genotype run rebuild the Library from the VCF alone."""
    def fmt(k):
        v = P.get(k, "NA")
        return "%s=%s" % (k, ("%.6g" % v) if isinstance(v, (float, np.floating)) else v)
    return VAF_HEADER_KEY + "<model=site_yield_betabinom," + ",".join(fmt(k) for k in _HEADER_FIELDS) + ">"


def parse_vaf_params_header_line(line):
    body = line.strip()[len(VAF_HEADER_KEY):].strip("<>")
    out = {}
    for kv in body.split(","):
        k, v = kv.split("=", 1)
        try:
            out[k] = float(v)
        except ValueError:
            out[k] = v
    return out
