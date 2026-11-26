"""Shield breach: the breakpoint is found where F turns down, the line wins without a bend, the null p and the
identity gap behave."""

import numpy as np
import pytest

from entropy_age.experiments.shield_breach import identity_gap, interval, line_null_p, piecewise_fit


def _bent(seed: int, n: int = 300):
    rng = np.random.default_rng(seed)
    b = np.sort(rng.uniform(0, 1, n))
    f = 0.5 - 0.8 * np.maximum(0, b - 0.6) + rng.normal(0, 0.01, n)
    return b, f, rng


def test_breakpoint_found_where_f_starts_falling():
    b, f, _ = _bent(1)
    fit = piecewise_fit(b, f)
    assert fit["tau"] == pytest.approx(0.6, abs=0.05)
    assert abs(fit["slope_before"]) < 0.05
    assert fit["slope_change"] == pytest.approx(-0.8, abs=0.1)
    assert fit["delta_bic"] > 10


def test_breakpoint_with_a_covariate_held():
    b, f, rng = _bent(2)
    age = rng.normal(0, 1, b.size)
    fit = piecewise_fit(b, f + 0.3 * age, age[:, None])
    assert fit["tau"] == pytest.approx(0.6, abs=0.05)
    assert fit["slope_change"] == pytest.approx(-0.8, abs=0.1)


def test_a_straight_line_gives_no_supported_bend():
    rng = np.random.default_rng(3)
    b = rng.uniform(0, 1, 300)
    f = 0.2 - 0.3 * b + rng.normal(0, 0.05, b.size)
    assert piecewise_fit(b, f)["delta_bic"] < 2


def test_line_null_p_is_calibrated_under_a_line():
    ps = []
    for seed in range(12):
        rng = np.random.default_rng(100 + seed)
        b = rng.uniform(0, 1, 200)
        ps.append(line_null_p(b, 0.2 - 0.3 * b + rng.normal(0, 0.05, b.size), None, 39, rng))
    assert min(ps) > 0.1 and np.median(ps) > 0.25  # no false bend under a line, at these seeds


def test_line_null_rejects_a_real_bend():
    b, f, rng = _bent(4)
    assert line_null_p(b, f, None, 49, rng) <= 0.02


def test_interval_needs_enough_draws():
    assert np.isnan(interval(np.arange(5.0))[0])
    lo, hi = interval(np.linspace(0, 1, 1001))
    assert lo == pytest.approx(0.025) and hi == pytest.approx(0.975)


def test_identity_gap_sign_follows_faster_drift_at_members():
    rng = np.random.default_rng(5)
    n, p = 120, 600
    ages = rng.uniform(20, 80, n)
    base = rng.uniform(0.2, 0.8, p)
    member = np.zeros(p, bool)
    member[:40] = True
    spread = np.where(member, 0.004, 0.0005)  # members spread away from their level faster with age
    x = np.clip(base + rng.normal(0, 1, (n, p)) * spread * (ages[:, None] - 15), 0.01, 0.99)
    island = rng.integers(0, 2, p).astype(bool)
    gap = identity_gap(x, ages, np.zeros((n, 0)), member, np.ones(p), island)
    assert gap > 0.5
