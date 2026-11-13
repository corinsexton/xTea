"""GQ error floor and the per-count beta-binomial in xtea/x_genotype_lik.py.

The beta-binomial reduces to the binomial as rho -> 0 and matches scipy's parameterization;
gq_from_error is monotone, tops out at -10 log10(GQ_ERR_FLOOR), reduces to the plain capped
phred when the floor is 0, and never changes the GT call.
"""
import os
import sys

import numpy as np
from scipy import stats

sys.path.insert(0, os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "xtea"))

import x_genotype_lik as gl  # noqa: E402


def test_betabinom_tends_to_binomial():
    k = np.arange(0, 31)
    n = np.full_like(k, 30)
    for p in (0.05, 0.44, 0.96):
        ref = stats.binom.logpmf(k, 30, p)
        assert np.allclose(gl.betabinom_logpmf(k, n, p, rho=0), ref, atol=1e-10)
        assert np.allclose(gl.betabinom_logpmf(k, n, p, rho=1e-9), ref, atol=1e-5)


def test_betabinom_matches_scipy():
    k = np.arange(0, 21)
    n = np.full_like(k, 20)
    p, rho = 0.57, 0.05
    a, b = p * (1 - rho) / rho, (1 - p) * (1 - rho) / rho
    assert np.allclose(gl.betabinom_logpmf(k, n, p, rho), stats.betabinom.logpmf(k, 20, a, b))


def test_betabinom_known_value():
    # n = 1: P(k=1) is the mean whatever rho is.
    for rho in (0.0, 0.05, 0.3):
        assert abs(np.exp(gl.betabinom_logpmf(1, 1, 0.3, rho)) - 0.3) < 1e-12
    # n = 2, k = 2: E[q^2] = p^2 + rho p (1 - p).
    p, rho = 0.4, 0.1
    assert abs(np.exp(gl.betabinom_logpmf(2, 2, p, rho)) - (p * p + rho * p * (1 - p))) < 1e-12


def test_gq_floor_ceiling_and_zero_floor():
    gq_max = int(-10 * np.log10(gl.GQ_ERR_FLOOR) + 1e-9)
    assert gl.gq_from_error(0.0)[()] == gq_max == 30
    assert gl.gq_from_error(1e-30)[()] == gq_max
    perr = np.array([0.5, 0.1, 0.01, 1e-3, 1e-5, 1e-9])
    plain = np.floor(np.minimum(gl.GQ_CAP, -10 * np.log10(perr)) + 1e-9).astype(int)
    assert list(gl.gq_from_error(perr, err_floor=0.0)) == list(plain)
    # Low GQ is barely touched by the floor.
    assert list(gl.gq_from_error([0.5, 0.1, 0.01])) == [3, 9, 19]


def test_gq_monotone_in_error():
    perr = np.sort(np.concatenate([np.logspace(-15, np.log10(0.5), 400), [0.0]]))
    for floor in (0.0, 1e-4, gl.GQ_ERR_FLOOR, 0.01):
        gq = gl.gq_from_error(perr, err_floor=floor)
        assert np.all(np.diff(gq) <= 0)


def test_floor_does_not_change_gt():
    rng = np.random.default_rng(3)
    n = 2000
    counts = np.column_stack([rng.integers(0, 30, n), rng.integers(0, 30, n),
                              rng.integers(0, 60, n), rng.integers(0, 60, n)]).astype(float)
    P0 = dict(gl.DEFAULT_PARAMS, rho=gl.DEFAULT_RHO, gq_err_floor=0.0)
    P1 = dict(P0, gq_err_floor=gl.GQ_ERR_FLOOR)
    r0, r1 = gl.genotype(counts, P0), gl.genotype(counts, P1)
    assert list(r0["gt"]) == list(r1["gt"])
    assert np.array_equal(r0["gl"], r1["gl"])          # GTL is the count likelihood only
    assert np.all(r1["gq"] <= r0["gq"])
    assert r1["gq"].max() <= 30


def test_fit_params_records_floor():
    P = gl.fit_params(np.zeros((0, 4)))
    assert P["gq_err_floor"] == gl.GQ_ERR_FLOOR
    assert "gq_err_floor=0.001" in gl.params_header_line(P)
