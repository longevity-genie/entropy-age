"""RM32: drift toward each CpG's own equilibrium."""

import numpy as np

from entropy_age.experiments.equilibrium_drift import adjust, bounded_relaxation, pooled_rates


def test_pooled_rates_recover_gain_and_loss():
    rng = np.random.default_rng(0)
    level = rng.uniform(0, 1, 500)
    slope = 0.02 * (1 - level) - 0.05 * level + rng.normal(0, 1e-4, level.size)
    r = pooled_rates(slope, level)
    assert abs(r["k_gain"] - 0.02) < 1e-3 and abs(r["k_loss"] - 0.05) < 1e-3
    assert abs(r["z"] - 0.02 / 0.07) < 0.01


def test_bounded_relaxation_recovers_equilibrium_and_holds_bounds():
    ages = np.linspace(20, 90, 300)
    t = ages - ages.min()
    bending = 0.6 + (0.1 - 0.6) * np.exp(-0.05 * t)
    beyond = np.clip(0.95 + 0.004 * t, 0, 1.2)  # a straight rise whose free fit would put z above 1
    fit = bounded_relaxation(np.column_stack([bending, beyond]), ages)
    assert abs(fit["equilibrium"][0] - 0.6) < 0.02 and abs(fit["rate"][0] - 0.05) < 0.01
    assert fit["identified"][0]
    assert 0 <= fit["equilibrium"][1] <= 1
    assert abs(fit["k_gain"][0] - fit["equilibrium"][0] * fit["rate"][0] * 10) < 1e-9


def test_adjust_removes_covariate_keeps_age():
    rng = np.random.default_rng(1)
    ages = rng.uniform(20, 80, 400)
    batch = (rng.random(400) < 0.5).astype(float)
    x = (0.3 + 0.002 * ages + 0.1 * batch + rng.normal(0, 0.005, 400))[:, None]
    y = adjust(x, ages, batch[:, None])
    slope = np.polyfit(ages, y[:, 0], 1)[0]
    assert abs(slope - 0.002) < 2e-4
    assert (
        abs(y[batch == 1].mean() - y[batch == 0].mean() - 0.002 * (ages[batch == 1].mean() - ages[batch == 0].mean()))
        < 0.003
    )
