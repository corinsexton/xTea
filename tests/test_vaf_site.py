"""Per-site VAF model (x_vaf_site) and its use in the genotype step (x_genotype_lik)."""
import os
import sys

import numpy as np
import pytest
from scipy.stats import betabinom

sys.path.insert(0, os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "xtea"))
import x_genotype_lik as gl  # noqa: E402
import x_vaf_site as xv  # noqa: E402

from .test_gvcf_multi import _row, bam_file, ref_fasta  # noqa: E402,F401  (fixtures)

LIB = xv.Library(151, 420, 100)
CAL = dict(xv.DEFAULT_CAL)


def _simulate(v, S, t, cal=CAL, lib=LIB, rho=xv.DEFAULT_RHO, depth=(30, 60), seed=0):
    """Counts drawn from the model itself: SR/PE depth Poisson, alt reads beta-binomial."""
    rng = np.random.default_rng(seed)
    v = np.broadcast_to(np.asarray(v, float), np.shape(S))
    y = xv.yields(lib, S, t, cal)
    p_sr, p_pe = xv.p_alt(v, y, cal)
    s = (1 - rho) / rho
    n_sr = rng.poisson(depth[0], len(S))
    n_pe = rng.poisson(depth[1], len(S))
    k_sr = betabinom.rvs(n_sr, p_sr * s, (1 - p_sr) * s, random_state=rng)
    k_pe = betabinom.rvs(n_pe, p_pe * s, (1 - p_pe) * s, random_state=rng)
    return np.c_[k_sr, n_sr - k_sr, k_pe, n_pe - k_pe].astype(float)


# ---- geometry ------------------------------------------------------------------------------
def test_yields_follow_insertion_length():
    S = np.array([50, 150, 300, 1000, 6000])
    y = xv.yields(LIB, S, np.full(5, 14), CAL)
    assert np.all(np.diff(y["D_alt"]) >= 0)            # longer insertion -> more discordant
    assert y["D_alt"][-1] == pytest.approx(y["D_alt"][-2], rel=0.05)   # saturates past F
    assert np.all(np.diff(y["C_alt"]) <= 0)            # fewer pairs span it
    assert np.all(y["C_ref"] == y["C_ref"][0])         # ref yield does not depend on S
    # a short Alu reads lower than a full-length L1 at the same VAF
    p_pe = xv.p_alt(0.5, y, CAL)[1]
    assert p_pe[1] < p_pe[-1]


def test_short_inserts_give_fewer_discordant_pairs():
    short, normal = xv.Library(151, 250, 60), xv.Library(151, 450, 100)
    S, t = np.array([300.0]), np.array([14.0])
    assert xv.yields(short, S, t, CAL)["D_alt"][0] < xv.yields(normal, S, t, CAL)["D_alt"][0]


def test_p_alt_background_and_monotone():
    y = xv.yields(LIB, np.array([300.0]), np.array([14.0]), CAL)
    cal = dict(CAL, b_SR=0.02, b_PE=0.04)
    p0 = xv.p_alt(0.0, y, cal)
    assert p0[0][0] == pytest.approx(0.02 / 1.02) and p0[1][0] == pytest.approx(0.04 / 1.04)
    v = np.linspace(0, 1, 50)
    p_sr, p_pe = xv.p_alt(v, y, cal)
    assert np.all(np.diff(p_sr) > 0) and np.all(np.diff(p_pe) > 0)


# ---- estimator -----------------------------------------------------------------------------
@pytest.mark.parametrize("v", [0.05, 0.25, 0.5, 1.0])
def test_site_vaf_recovers_truth(v):
    rng = np.random.default_rng(1)
    S = rng.choice([280.0, 2000.0, 6000.0], 400)
    t = np.full(400, 14.0)
    c = _simulate(v, S, t, depth=(60, 120), seed=int(v * 100))
    vaf, lo, hi = xv.site_vaf(c, xv.yields(LIB, S, t, CAL), CAL)
    assert np.median(vaf) == pytest.approx(v, rel=0.15)
    assert 0.85 <= np.mean((lo <= v) & (v <= hi)) <= 1.0
    assert np.all((lo <= vaf) & (vaf <= hi))


def test_site_vaf_without_reads_is_missing():
    y = xv.yields(LIB, np.array([300.0]), np.array([14.0]), CAL)
    vaf, lo, hi = xv.site_vaf(np.zeros((1, 4)), y, CAL)
    assert np.isnan(vaf[0]) and np.isnan(lo[0]) and np.isnan(hi[0])


def test_low_vaf_interval_reaches_zero():
    y = xv.yields(LIB, np.array([300.0]), np.array([14.0]), CAL)
    vaf, lo, hi = xv.site_vaf(np.array([[0, 200, 0, 600.0]]), y, CAL, rho=xv.MOSAIC_RHO)
    assert lo[0] == 0.0 and hi[0] < 0.05


# ---- calibration ---------------------------------------------------------------------------
def test_fit_calibration_recovers_germline_sample():
    truth = {"e_SR": 0.9, "e_PE": 1.6, "l_SR": 6.0, "l_PE": 60.0, "b_SR": 0.0, "b_PE": 0.0}
    rng = np.random.default_rng(2)
    n = 1500
    S = rng.choice([280.0, 2000.0, 6000.0], n, p=[0.8, 0.1, 0.1])
    t = rng.integers(5, 20, n).astype(float)
    hom = rng.random(n) < 0.3
    c = _simulate(np.where(hom, 1.0, 0.5), S, t, cal=truth, seed=3)
    fit = xv.fit_calibration(c, LIB, S, t)
    assert fit["source"] == "fit"
    assert fit["e_SR"] == pytest.approx(0.9, rel=0.15)
    assert fit["e_PE"] == pytest.approx(1.6, rel=0.15)
    assert fit["w"] == pytest.approx(0.3, abs=0.06)
    # and with that fit the VAFs land at 0.5 / 1
    cal = {k: fit[k] for k in xv.FIT_KEYS}
    cal.update(b_SR=0.0, b_PE=0.0)
    vaf = xv.site_vaf(c, xv.yields(LIB, S, t, cal), cal)[0]
    assert np.median(vaf[~hom]) == pytest.approx(0.5, abs=0.05)
    assert np.median(vaf[hom]) > 0.9


def test_fit_calibration_falls_back_on_few_calls():
    S, t = np.full(20, 300.0), np.full(20, 14.0)
    fit = xv.fit_calibration(_simulate(0.5, S, t), LIB, S, t)
    assert fit["source"] == "default" and fit["reason"].startswith("n_informative")
    assert fit["e_SR"] == xv.DEFAULT_CAL["e_SR"]


def test_background_from_counts():
    b = xv.background_from_counts(np.array([[1, 99, 2, 98], [0, 0, 0, 0], [1, 99, 2, 98.0]]))
    assert b["b_SR"] == pytest.approx(0.01 / 0.99) and b["b_PE"] == pytest.approx(0.02 / 0.98)


def test_estimate_mosaic_calibrates_on_germline_like_only():
    rng = np.random.default_rng(4)
    S = np.full(300, 300.0)
    t = np.full(300, 14.0)
    germ = _simulate(np.where(rng.random(300) < 0.3, 1.0, 0.5), S, t, seed=5)
    mos = _simulate(0.05, S[:100], t[:100], cal=dict(CAL, **xv.MOSAIC_BACKGROUND), seed=6,
                    rho=xv.MOSAIC_RHO, depth=(200, 400))
    vaf, lo, hi, P = xv.estimate(mos, LIB, S[:100], t[:100], mosaic=True, calib=(germ, S, t))
    assert P["mode"] == "mosaic" and P["source"] == "fit" and P["rho"] == xv.MOSAIC_RHO
    assert P["b_SR"] == xv.MOSAIC_BACKGROUND["b_SR"]
    assert np.median(vaf) == pytest.approx(0.05, rel=0.25)
    # the mosaic calls themselves are not germline-like: nothing to calibrate on
    P2 = xv.estimate(mos, LIB, S[:100], t[:100], mosaic=True)[3]
    assert P2["source"] == "default" and P2["n_calib"] == 0


# ---- inputs --------------------------------------------------------------------------------
def test_parse_svlen_and_tsd():
    assert xv.parse_svlen("312", "ALU") == 312
    assert xv.parse_svlen("-1") == xv.GENERIC_DEFAULT_LEN
    assert xv.parse_svlen("77142134", "LINE1") == 6000           # orphan transduction junk
    assert xv.parse_svlen("0", "SVA") == 2000
    assert xv.parse_svlen(None) == xv.GENERIC_DEFAULT_LEN
    assert xv.parse_tsdlen("-1") == xv.DEFAULT_TSD
    assert xv.parse_tsdlen("16") == 16
    assert xv.parse_tsdlen("NULL") == xv.DEFAULT_TSD


def test_read_basic_info(tmp_path):
    sf = tmp_path / "basic_cov_is_rlth_info.txt"
    sf.write_text("a.bam\t30\t150\t400\t100\nb.bam\t10\t150\t480\t140\n")
    lib = xv.read_basic_info(str(sf))
    assert lib.R == 150 and lib.mean_is == pytest.approx(420) and lib.sd_is == pytest.approx(110)
    sf.write_text("a.bam\t30\t0\t400\t100\n")
    assert xv.read_basic_info(str(sf)) is None
    assert xv.read_basic_info(str(tmp_path / "missing.txt")) is None


def _raw(chrm, pos, counts, svlen, tsd=14):
    f = ["0"] * 50
    f[0], f[2], f[3], f[4] = chrm, str(pos), str(pos + 15), str(tsd)
    for i, x in zip(xv.RAW_COUNT_IDX, counts):
        f[i] = str(int(x))
    f[47] = str(svlen)
    return f


def test_raw_site_inputs_family_from_type_map(tmp_path):
    rows = [_raw("chr1", 100, (5, 5, 5, 20), -1), _raw("chr1", 500, (5, 5, 5, 20), 300, tsd=-1)]
    sf = tmp_path / "map.tsv"
    sf.write_text("chr1\t100\t4\nchr1\t500\t2\n")
    c, S, t = xv.raw_site_inputs(rows, rtype=1, type_map=xv.load_type_map(str(sf)))
    assert list(S) == [2000, 300] and list(t) == [14, xv.DEFAULT_TSD]
    assert c[0].tolist() == [5, 5, 5, 20]
    assert xv.raw_site_inputs(rows, rtype=1)[1][0] == 6000        # no map: rtype's family


def test_vaf_params_header_roundtrip():
    P = xv.estimate(_simulate(0.5, np.full(80, 300.0), np.full(80, 14.0)), LIB,
                    np.full(80, 300.0), np.full(80, 14.0))[3]
    back = xv.parse_vaf_params_header_line(xv.vaf_params_header_line(P))
    assert back["mode"] == "germline" and back["read_len"] == 151
    for k in ("mean_is", "sd_is", "e_SR", "e_PE", "l_PE", "rho"):
        assert back[k] == pytest.approx(P[k], rel=1e-5)


# ---- genotype step -------------------------------------------------------------------------
def _germline_raw(tmp_path, n=200, seed=7):
    rng = np.random.default_rng(seed)
    S = rng.choice([280.0, 6000.0], n)
    hom = rng.random(n) < 0.3
    c = _simulate(np.where(hom, 1.0, 0.5), S, np.full(n, 14.0), seed=seed)
    sf = tmp_path / "gene.txt"
    sf.write_text("".join("\t".join(_raw("chr1", 1000 + 100 * i, c[i], int(S[i]))) + "\n"
                          for i in range(n)))
    info = tmp_path / "basic_cov_is_rlth_info.txt"
    info.write_text("x.bam\t30\t%g\t%g\t%g\n" % (LIB.R, LIB.mean_is, LIB.sd_is))
    return str(sf), str(info), hom


def _vaf_cols(sf):
    rows = [ln.rstrip("\n").split("\t") for ln in open(sf)]
    k = gl.FORMAT_KEYS
    return [dict(zip(k, r[50:50 + len(k)])) for r in rows]


def test_genotype_step_uses_site_vaf(tmp_path):
    sf, info, hom = _germline_raw(tmp_path)
    out = str(tmp_path / "gntp.txt")
    gl.genotype_raw_file(sf, out, basic_info=info, rtype=2)
    side = open(out + gl.PARAMS_SIDECAR_SUFFIX).read()
    assert gl.PARAMS_HEADER_KEY in side and gl.VAF_HEADER_KEY in side
    cols = _vaf_cols(out)
    vaf = np.array([float(c["VAF"]) for c in cols])
    assert np.median(vaf[~hom]) == pytest.approx(0.5, abs=0.07)
    assert np.median(vaf[hom]) > 0.9
    for c in cols:
        lo, hi = map(float, c["VAF_CI"].split(","))
        assert lo <= float(c["VAF"]) <= hi and c["GT"] in ("0/1", "1/1")


def test_genotype_step_without_basic_info_keeps_bias_line(tmp_path):
    sf, _, _ = _germline_raw(tmp_path)
    out = str(tmp_path / "gntp.txt")
    gl.genotype_raw_file(sf, out)
    assert gl.VAF_HEADER_KEY not in open(out + gl.PARAMS_SIDECAR_SUFFIX).read()
    rows = [ln.split("\t") for ln in open(sf)]
    counts = np.array([[int(f[i]) for i in gl.RAW_IDX] for f in rows], float)
    res = gl.genotype(counts, gl.fit_params(counts))
    assert [c["VAF"] for c in _vaf_cols(out)] == [gl._fmt(float(x)) for x in res["vaf"]]


def test_genotype_step_mosaic(tmp_path):
    germ, info, _ = _germline_raw(tmp_path, n=300)
    S = np.full(60, 300.0)
    c = _simulate(0.05, S, np.full(60, 14.0), cal=dict(CAL, **xv.MOSAIC_BACKGROUND),
                  rho=xv.MOSAIC_RHO, depth=(200, 400), seed=8)
    sf = tmp_path / "mosaic.txt"
    sf.write_text("".join("\t".join(_raw("chr2", 1000 + 100 * i, c[i], 300)) + "\n"
                          for i in range(60)))
    out = str(tmp_path / "gntp.txt")
    assert gl.genotype_raw_file(str(sf), out, basic_info=info, rtype=2, mosaic=True,
                                calib_files=[(germ, 2)]) is None
    side = open(out + gl.PARAMS_SIDECAR_SUFFIX).read()
    assert gl.PARAMS_HEADER_KEY not in side
    P = xv.parse_vaf_params_header_line(side.strip())
    assert P["mode"] == "mosaic" and P["source"] == "fit"
    cols = _vaf_cols(out)
    assert all(x["GT"] == "./." and x["GQ"] == "." and x["GTL"] == ".,." for x in cols)
    assert np.median([float(x["VAF"]) for x in cols]) == pytest.approx(0.05, rel=0.3)


def test_gvcf_header_carries_vaf_params(tmp_path, ref_fasta, bam_file):
    from x_gvcf import gVCF
    sf_raw = tmp_path / "gene.txt"
    sf_raw.write_text(_row("chr1", 100, 110) + "\tgene_info\n")
    info = tmp_path / "basic.txt"
    info.write_text("x.bam\t30\t150\t420\t100\n")
    sf_gntp = str(tmp_path / "gntp.txt")
    gl.genotype_raw_file(str(sf_raw), sf_gntp, basic_info=str(info), rtype=1)
    sf_vcf = str(tmp_path / "out.vcf")
    gVCF().cvt_raw_rslt_to_gvcf("SAMPLE", bam_file, sf_gntp, 1, ref_fasta, sf_vcf)
    head = [ln for ln in open(sf_vcf) if ln.startswith("##")]
    assert sum(h.startswith(gl.VAF_HEADER_KEY) for h in head) == 1
    assert sum(h.startswith(gl.PARAMS_HEADER_KEY) for h in head) == 1


def test_cli_reads_library_from_header(tmp_path):
    rng = np.random.default_rng(9)
    S = rng.choice([280.0, 6000.0], 200)
    hom = rng.random(200) < 0.3
    c = _simulate(np.where(hom, 1.0, 0.5), S, np.full(200, 14.0), seed=9)
    sf_in, sf_a, sf_b = (str(tmp_path / x) for x in ("in.vcf", "a.vcf", "b.vcf"))
    with open(sf_in, "w") as fo:
        fo.write("##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n")
        for i in range(200):
            fo.write("chr1\t%d\t.\tA\t<INS:ME:ALU>\t.\tPASS\tSVLEN=%d;TSDLEN=14;AF_CLIP=%d;"
                     "AF_FMAP=%d;AF_DISC=%d;AF_CONCORDNT=%d\n" % ((1000 + i, S[i]) + tuple(c[i])))
    # no library anywhere: bias-line VAF, no vaf header
    gl.main([sf_in, "-o", sf_a])
    assert not any(ln.startswith(gl.VAF_HEADER_KEY) for ln in open(sf_a))
    # library on the command line -> header -> picked up again from the header
    gl.main([sf_in, "-o", sf_a, "--read-length", "151", "--insert-mean", "420",
             "--insert-sd", "100"])
    gl.main([sf_a, "-o", sf_b])
    assert open(sf_a).read() == open(sf_b).read()
    vaf = np.array([float(ln.split("\t")[9].split(":")[5]) for ln in open(sf_b)
                    if not ln.startswith("#")])
    assert np.median(vaf[~hom]) == pytest.approx(0.5, abs=0.07)
    # mosaic: no GT, VAF still written
    gl.main([sf_a, "-o", sf_b, "--mosaic", "--vaf-background", "0.01,0.02"])
    rec = [ln.split("\t") for ln in open(sf_b) if not ln.startswith("#")]
    assert all(r[9].startswith("./.:.:.,.:") for r in rec)
    vh = [ln for ln in open(sf_b) if ln.startswith(gl.VAF_HEADER_KEY)][0]
    assert xv.parse_vaf_params_header_line(vh)["b_PE"] == pytest.approx(0.02)


# ---- per-sample background at control positions ------------------------------------------------
def test_fit_background_recovers_simulated_background():
    b = {"b_SR": 0.02, "b_PE": 0.05}
    n = 1500
    c = _simulate(0.0, np.full(n, xv.CONTROL_SVLEN), np.full(n, xv.CONTROL_TSD),
                  cal=dict(CAL, **b), rho=xv.MOSAIC_RHO, depth=(40, 80), seed=11)
    fit = xv.fit_background(c, LIB)
    assert fit["b_SR"] == pytest.approx(0.02, rel=0.15)
    assert fit["b_PE"] == pytest.approx(0.05, rel=0.15)
    mei = xv.control_background(c, LIB)
    for k in ("b_SR", "b_PE"):
        assert mei[k] == pytest.approx(fit[k] * xv.CONTEXT_FACTOR[k])
    # no reads at all -> no background
    assert xv.fit_background(np.zeros((5, 4)), LIB) == {"b_SR": 0.0, "b_PE": 0.0}


def test_estimate_records_background_source():
    S, t = np.full(5, 300.0), np.full(5, 14.0)
    c = np.tile([5, 95, 5, 95.0], (5, 1))
    assert xv.estimate(c, LIB, S, t)[3]["bg_source"] == "none"
    assert xv.estimate(c, LIB, S, t, mosaic=True)[3]["bg_source"] == "default"
    P = xv.estimate(c, LIB, S, t, mosaic=True,
                    background={"b_SR": 0.03, "b_PE": 0.06, "source": "control:n=1000"})[3]
    assert P["bg_source"] == "control:n=1000" and P["b_PE"] == 0.06
    assert "bg_source=control:n=1000" in xv.vaf_params_header_line(P)
    assert xv.parse_vaf_params_header_line(xv.vaf_params_header_line(P))["bg_source"] == "control:n=1000"


def test_parse_background_value_or_file(tmp_path):
    assert gl.parse_background(None) is None
    assert gl.parse_background("0.01,0.02") == {"b_SR": 0.01, "b_PE": 0.02, "source": "user"}
    f = tmp_path / "bg.txt"
    f.write_text("0.03,0.04\n")
    assert gl.parse_background(str(f)) == {"b_SR": 0.03, "b_PE": 0.04, "source": "file"}
    with pytest.raises(ValueError):
        gl.parse_background("0.01")


@pytest.fixture
def control_bam(tmp_path):
    """chr1 (4 kb) with 100 bp reads every 10 bp, proper pairs (insert 300), every 5th read
    soft-clipped by 30 bp: a uniform background of clipped reads and no discordant pairs."""
    import pysam
    ref = tmp_path / "ref.fa"
    seq = ("ACGTTGCAAG" * 400)
    ref.write_text(">chr1\n" + "\n".join(seq[i:i + 60] for i in range(0, len(seq), 60)) + "\n")
    pysam.faidx(str(ref))
    sf = str(tmp_path / "ctl.bam")
    header = {"HD": {"VN": "1.6", "SO": "coordinate"}, "SQ": [{"SN": "chr1", "LN": len(seq)}]}
    reads = []
    for i, start in enumerate(range(0, len(seq) - 400, 10)):
        for mate in (0, 1):
            a = pysam.AlignedSegment()
            a.query_name = "r%d" % i
            pos = start if mate == 0 else start + 200
            a.reference_id, a.reference_start = 0, pos
            a.next_reference_id, a.next_reference_start = 0, start + 200 if mate == 0 else start
            a.query_sequence = seq[pos:pos + 100]
            a.query_qualities = pysam.qualitystring_to_array("I" * 100)
            a.flag = (1 | 2 | (32 if mate == 0 else 16) | (64 if mate == 0 else 128))
            a.mapping_quality = 60
            a.template_length = 300 if mate == 0 else -300
            a.cigartuples = [(4, 30), (0, 70)] if (i % 5 == 0 and mate == 0) else [(0, 100)]
            reads.append(a)
    reads.sort(key=lambda r: r.reference_start)
    with pysam.AlignmentFile(sf, "wb", header=header) as out:
        for a in reads:
            out.write(a)
    pysam.index(sf)
    return sf, str(ref)


def test_control_counts_on_a_bam(control_bam, tmp_path):
    import x_vaf_background as xb
    bam, ref = control_bam
    lib = xv.Library(100, 300, 30)
    c = xb.control_counts(bam, ref, lib, n=12, cores=1, seed=3, tmp_dir=str(tmp_path))
    assert c.shape == (12, 4)
    assert c[:, 1].sum() > 0 and c[:, 3].sum() > 0     # full-map reads, concordant pairs
    assert c[:, 2].sum() == 0                          # no discordant pairs in this BAM
    bl = tmp_path / "bams.txt"
    bl.write_text(bam + "\tillumina\n")
    assert xb.first_bam(str(bl)) == bam
    b = xb.sample_background(str(bl), ref, lib, n=12, tmp_dir=str(tmp_path))
    assert b["n_control"] == 12 and b["b_PE"] == pytest.approx(0.0, abs=1e-4)
