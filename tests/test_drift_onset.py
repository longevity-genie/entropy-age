"""RM29: development removed by the trend, the change point found, and the line preferred when there is no bend."""

import numpy as np
import pytest

from entropy_age.experiments.drift_onset import hinge_fit, noise_per_animal, trend_residuals


def test_hinge_finds_a_flat_then_rising_noise():
    rng = np.random.default_rng(10)
    age = np.sort(rng.uniform(0, 20, 300))
    d = 1.0 + 0.05 * np.maximum(0, age - 5) + rng.normal(0, 0.02, age.size)
    fit = hinge_fit(age, d)
    assert fit["tau"] == pytest.approx(5, abs=1)
    assert abs(fit["slope_before"]) < 0.01 and fit["slope_after"] == pytest.approx(0.05, abs=0.01)
    assert fit["delta_bic"] > 10 and fit["onset"]


def test_a_straight_line_gives_no_supported_bend():
    rng = np.random.default_rng(11)
    age = np.sort(rng.uniform(0, 20, 300))
    d = 1.0 + 0.03 * age + rng.normal(0, 0.05, age.size)
    assert not hinge_fit(age, d)["onset"]


def test_trend_removes_development_and_leaves_noise_growth():
    rng = np.random.default_rng(12)
    n, p = 400, 200
    age = rng.uniform(0.05, 10, n)
    level = rng.uniform(0.2, 0.8, p)
    development = 0.15 * np.tanh(3 * np.log(age + 0.1))[:, None] * rng.choice([-1, 1], p)  # directed, early
    noise_sd = 0.01 + 0.004 * np.maximum(0, age - 2)  # noise rises after age 2
    x = np.clip(level + development + rng.normal(size=(n, p)) * noise_sd[:, None], 0, 1)
    res = trend_residuals(x, age, maturity=2.0)
    assert np.abs(res).mean() < np.abs(x - x.mean(axis=0)).mean() / 2  # development absorbed
    d = noise_per_animal(x, res, np.arange(p))
    fit = hinge_fit(age, d)
    assert fit["onset"] and fit["tau"] == pytest.approx(2, abs=1.2)


def test_development_alone_makes_no_bend_at_maturity():
    """The null that matters: steep early directed change and steady noise must not look like a late onset."""
    rng = np.random.default_rng(13)
    n, p = 400, 200
    age = rng.uniform(0.05, 10, n)
    level = rng.uniform(0.2, 0.8, p)
    development = 0.15 * np.tanh(3 * np.log(age + 0.1))[:, None] * rng.choice([-1, 1], p)
    x = np.clip(level + development + rng.normal(0, 0.02, (n, p)), 0, 1)
    res = trend_residuals(x, age, maturity=2.0)
    fit = hinge_fit(age, noise_per_animal(x, res, np.arange(p)))
    assert not fit["onset"] and fit["delta_bic"] < 6
