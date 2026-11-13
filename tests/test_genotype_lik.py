"""Unit tests for xtea/x_genotype_lik.py, the beta-binomial genotyper.

Covers: EM parameter recovery on synthetic counts, the low-n and degenerate fallbacks to the
HG002 defaults, the GQ cap, the post-hoc CLI round trip on a small VCF, and the pipeline
step (genotype_raw_file) feeding x_gvcf's FORMAT columns.
"""
import os
import sys

import numpy as np
import pytest

sys.path.insert(0, os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "xtea"))

import x_genotype_lik as gl  # noqa: E402

from .test_gvcf_multi import _row, bam_file, ref_fasta  # noqa: E402,F401  (fixtures)


TRUE = {"sr_het": 0.45, "sr_hom": 0.93, "pe_het": 0.33, "pe_hom": 0.72, "hom_prior": 0.35}


def _simulate(n, P, rho=gl.DEFAULT_RHO, seed=0, sr_depth=20, pe_depth=25):
    """counts (n, 4) and the true hom indicator, drawn from the model itself."""
    rng = np.random.default_rng(seed)
    hom = rng.random(n) < P["hom_prior"]
    out = np.zeros((n, 4))
    for t, (ia, ib), depth in (("sr", (0, 1), sr_depth), ("pe", (2, 3), pe_depth)):
        p = np.where(hom, P[t + "_hom"], P[t + "_het"])
        a = p * (1 - rho) / rho
        b = (1 - p) * (1 - rho) / rho
        q = rng.beta(a, b)
        d = rng.poisson(depth, n) + 1
        k = rng.binomial(d, q)
        out[:, ia], out[:, ib] = k, d - k
    return out, hom


def test_em_recovers_parameters():
    counts, hom = _simulate(3000, TRUE)
    P = gl.fit_params(counts)
    assert P["source"] == "em", P["reason"]
    for k, v in TRUE.items():
        assert abs(P[k] - v) < 0.03, (k, P[k], v)
    res = gl.genotype(counts, P)
    acc = np.mean((res["gt"] == "1/1") == hom)
    assert acc > 0.95


def test_betabinom_pmf_normalizes():
    from scipy.special import logsumexp
    for n in (0, 1, 7, 40):
        k = np.arange(n + 1)
        tot = logsumexp(gl.betabinom_logpmf(k, np.full_like(k, n), 0.3))
        assert abs(tot) < 1e-9


def test_low_n_falls_back_to_defaults():
    counts, _ = _simulate(30, TRUE)
    P = gl.fit_params(counts)
    assert P["source"] == "default"
    assert P["reason"].startswith("n_informative<")
    for k, v in gl.DEFAULT_PARAMS.items():
        assert P[k] == v


def test_uninformative_calls_do_not_count():
    counts, _ = _simulate(200, TRUE)
    counts[:, 2:] = 0                       # no PE reads anywhere -> nothing informative
    P = gl.fit_params(counts)
    assert P["n_informative"] == 0 and P["source"] == "default"


def test_degenerate_fit_falls_back_to_defaults():
    # Every call has the same alt fraction: the two components cannot separate.
    counts = np.tile([10, 10, 10, 10], (500, 1)).astype(float)
    P = gl.fit_params(counts)
    assert P["source"] == "default"
    assert P["reason"].startswith("degenerate_em:")
    for k, v in gl.DEFAULT_PARAMS.items():
        assert P[k] == v


def test_gq_cap_and_no_read_call():
    P = dict(gl.DEFAULT_PARAMS, rho=gl.DEFAULT_RHO)
    counts = np.array([[400, 0, 400, 0],     # overwhelming hom evidence
                       [200, 200, 150, 250],  # overwhelming het evidence
                       [0, 0, 0, 0]], float)  # nothing at all
    res = gl.genotype(counts, P)
    assert list(res["gt"]) == ["1/1", "0/1", "./."]
    gq_max = int(-10 * np.log10(gl.GQ_ERR_FLOOR) + 1e-9)    # the error floor's ceiling
    assert res["gq"][0] >= gq_max - 1 and res["gq"][1] >= gq_max - 1   # 29 or 30
    assert res["gq"].max() <= min(gl.GQ_CAP, gq_max)
    assert res["gq"][2] == 0
    assert gl.format_values(res, 2)[0] == "./."
    lo, hi = res["vaf_lo"][1], res["vaf_hi"][1]
    assert 0 <= lo < hi <= 1


def test_params_header_roundtrip():
    P = gl.fit_params(_simulate(1000, TRUE)[0])
    line = gl.params_header_line(P)
    back = gl.parse_params_header_line(line)
    assert back["source"] == "em" and back["gq_cap"] == gl.GQ_CAP
    assert back["gq_err_floor"] == gl.GQ_ERR_FLOOR
    for k in ("sr_het", "sr_hom", "pe_het", "pe_hom", "hom_prior"):
        assert abs(back[k] - P[k]) < 1e-4


def _vcf_text(counts, samples=("S1",)):
    out = ["##fileformat=VCFv4.2\n",
           '##FORMAT=<ID=GT,Number=1,Type=String,Description="old">\n',
           '##FORMAT=<ID=GL,Number=2,Type=Float,Description="stale, from an early build">\n',
           "##xtea2_genotype_params=<stale>\n",
           "\t".join(["#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO"]
                     + (["FORMAT"] + list(samples) if samples else [])) + "\n"]
    for i, (a, b, c, d) in enumerate(counts.astype(int)):
        info = "SVTYPE=INS;AF_CLIP=%d;AF_FMAP=%d;AF_DISC=%d;AF_CONCORDNT=%d" % (a, b, c, d)
        f = ["chr1", str(1000 + 10 * i), ".", "A", "<INS:ME:ALU>", ".", "PASS", info]
        if samples:
            f += ["GT"] + ["0/1"] * len(samples)
        out.append("\t".join(f) + "\n")
    return "".join(out)


def _write_vcf(sf, counts, samples=("S1",)):
    with open(sf, "w") as fout:
        fout.write(_vcf_text(counts, samples))


def _read_out(sf):
    head, recs = [], []
    with open(sf) as fin:
        for ln in fin:
            (head if ln.startswith("#") else recs).append(ln.rstrip("\n"))
    return head, recs


def test_cli_roundtrip(tmp_path, capsys):
    counts, hom = _simulate(400, TRUE, seed=3)
    sf_in, sf_out = str(tmp_path / "in.vcf"), str(tmp_path / "out.vcf")
    _write_vcf(sf_in, counts)
    assert gl.main([sf_in, "-o", sf_out]) == 0
    assert gl.PARAMS_HEADER_KEY in capsys.readouterr().err

    head, recs = _read_out(sf_out)
    assert head[-1].split("\t") == ["#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER",
                                     "INFO", "FORMAT", "S1"]
    assert not any(h.startswith("##FORMAT=<ID=GL,") for h in head)
    params = [h for h in head if h.startswith(gl.PARAMS_HEADER_KEY)]
    assert len(params) == 1 and "<stale>" not in params[0]
    assert gl.parse_params_header_line(params[0])["source"] == "em"
    for k in gl.FORMAT_KEYS:
        assert sum(h.startswith("##FORMAT=<ID=%s," % k) for h in head) == 1
    assert len(recs) == 400
    gts = []
    for r in recs:
        f = r.split("\t")
        assert len(f) == 10 and f[8] == ":".join(gl.FORMAT_KEYS)
        v = f[9].split(":")
        assert len(v) == len(gl.FORMAT_KEYS)
        assert 0 <= int(v[1]) <= gl.GQ_CAP
        gts.append(v[0])
    assert np.mean((np.array(gts) == "1/1") == hom) > 0.9

    # Re-running on its own output is stable.
    sf_out2 = str(tmp_path / "out2.vcf")
    gl.main([sf_out, "-o", sf_out2])
    assert open(sf_out).read() == open(sf_out2).read()


def test_cli_rejects_vcf_without_counts(tmp_path):
    sf_in = str(tmp_path / "bad.vcf")
    with open(sf_in, "w") as fout:
        fout.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n")
        fout.write("chr1\t5\t.\tA\t<INS>\t.\tPASS\tSVTYPE=INS\n")
    with pytest.raises(ValueError):
        gl.main([sf_in, "-o", str(tmp_path / "o.vcf")])


def test_raw_step_feeds_gvcf(tmp_path, ref_fasta, bam_file):
    from x_gvcf import gVCF
    sf_raw = tmp_path / "gene.txt"
    sf_raw.write_text(_row("chr1", 100, 110) + "\tgene_info\n")
    sf_gntp = str(tmp_path / "gntp.txt")
    P = gl.genotype_raw_file(str(sf_raw), sf_gntp)
    assert P["source"] == "default"                       # one call -> defaults
    row = open(sf_gntp).read().rstrip("\n").split("\t")
    assert len(row) == 50 + len(gl.FORMAT_KEYS)
    assert os.path.exists(sf_gntp + gl.PARAMS_SIDECAR_SUFFIX)

    sf_vcf = str(tmp_path / "out.vcf")
    gVCF().cvt_raw_rslt_to_gvcf("SAMPLE", bam_file, sf_gntp, 1, ref_fasta, sf_vcf)
    lines = open(sf_vcf).read().splitlines()
    assert any(h.startswith(gl.PARAMS_HEADER_KEY) for h in lines)
    for k in gl.FORMAT_KEYS:
        assert any(h.startswith("##FORMAT=<ID=%s," % k) for h in lines)
    rec = [ln for ln in lines if not ln.startswith("#")]
    assert len(rec) == 1
    f = rec[0].split("\t")
    assert f[8] == ":".join(gl.FORMAT_KEYS)
    assert f[9].split(":") == row[50:]


def test_cli_sites_only_input_gets_one_sample(tmp_path):
    counts, _ = _simulate(100, TRUE, seed=4)
    sf_in = str(tmp_path / "sites.vcf")
    _write_vcf(sf_in, counts, samples=())
    for argv, name in (([], gl.DEFAULT_SAMPLE), (["--sample", "HG002"], "HG002")):
        sf_out = str(tmp_path / ("o_%s.vcf" % name))
        gl.main([sf_in, "-o", sf_out] + argv)
        head, recs = _read_out(sf_out)
        assert head[-1].split("\t")[8:] == ["FORMAT", name]
        assert all(len(r.split("\t")) == 10 for r in recs)


def test_cli_rejects_multi_sample(tmp_path):
    counts, _ = _simulate(60, TRUE, seed=5)
    sf_in = str(tmp_path / "two.vcf")
    _write_vcf(sf_in, counts, samples=("A", "B"))
    with pytest.raises(ValueError, match="2 samples"):
        gl.main([sf_in, "-o", str(tmp_path / "o.vcf")])


def test_cli_reads_gzip(tmp_path):
    import gzip
    counts, _ = _simulate(200, TRUE, seed=6)
    sf_plain, sf_gz = str(tmp_path / "in.vcf"), str(tmp_path / "in.vcf.gz")
    _write_vcf(sf_plain, counts)
    with gzip.open(sf_gz, "wt") as fout:
        fout.write(_vcf_text(counts))
    gl.main([sf_plain, "-o", str(tmp_path / "a.vcf")])
    gl.main([sf_gz, "-o", str(tmp_path / "b.vcf")])
    assert open(str(tmp_path / "a.vcf")).read() == open(str(tmp_path / "b.vcf")).read()
