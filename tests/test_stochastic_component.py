"""RM23: Tong-style simulation modes and the blood-cell deconvolution."""

import numpy as np
import pytest

from entropy_age.experiments.stochastic_component import deconvolve, simulate, step_probability


def test_rates_follow_the_effect_only_in_tong_mode_and_the_shield_slows_set_cpgs():
    effect = np.array([0.0, 0.1, -0.3])
    shielded = np.array([False, True, False])
    tong = step_probability("tong", effect, 1.0, shielded)
    assert tong[0] == 0 and tong[2] > tong[1]
    common = step_probability("direction", effect, 1.0, shielded)
    assert np.allclose(common, common[0])
    sh = step_probability("shield", effect, 1.0, shielded)
    assert sh[1] < sh[0] == sh[2]
    with pytest.raises(ValueError):
        step_probability("other", effect, 1.0, shielded)


def test_directed_modes_move_with_the_effect_and_noise_does_not():
    rng = np.random.default_rng(0)
    ground = np.full((5, 2), 0.3)
    effect = np.array([0.2, -0.2])
    ages = np.full(200, 60.0)
    shielded = np.zeros(2, dtype=bool)
    tong = simulate(ground, ages, 30.0, effect, shielded, "tong", 20.0, 0.02, rng).mean(axis=0)
    assert tong[0] > 0.33 and tong[1] < 0.27
    noise = simulate(ground, ages, 30.0, effect, shielded, "noise", 20.0, 0.02, rng).mean(axis=0)
    assert abs(noise[0] - noise[1]) < 0.02  # no direction from the effect


def test_deconvolution_recovers_mixture_fractions():
    rng = np.random.default_rng(1)
    reference = rng.uniform(0, 1, (200, 3))
    truth = np.array([[0.6, 0.3, 0.1], [0.1, 0.1, 0.8]])
    x = truth @ reference.T
    assert np.allclose(deconvolve(x, reference), truth, atol=1e-6)
