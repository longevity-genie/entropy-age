"""RM27: a treatment that undoes directed aging, or the noise, is told apart; an inert one is not credited."""

import numpy as np

from entropy_age.experiments.intervention_split import split_test


def _cohort(rng, directed_back: float, noise_back: float, n=(12, 14, 12), p=2000):
    """Young, aged and treated mice; the first 100 CpGs shift with age, the next 100 get noisier."""
    base = rng.uniform(0.2, 0.8, p)
    shift = np.zeros(p)
    shift[:100] = 0.15
    sd_young, sd_aged = np.full(p, 0.02), np.full(p, 0.02)
    sd_aged[100:200] = 0.08

    def draw(k, mean_shift, sd):
        return base + mean_shift + rng.normal(size=(k, p)) * sd

    young = draw(n[0], 0, sd_young)
    aged = draw(n[1], shift, sd_aged)
    sd_treated = sd_aged - noise_back * (sd_aged - sd_young)
    treated = draw(n[2], shift * (1 - directed_back), sd_treated)
    x = np.vstack([young, aged, treated])
    groups = np.array(["Young_ctrl"] * n[0] + ["Aged_ctrl"] * n[1] + ["Aged_x"] * n[2])
    return x, groups


def test_directed_reversal_is_found_and_noise_left():
    rng = np.random.default_rng(16)
    x, g = _cohort(rng, directed_back=0.5, noise_back=0.0)
    r = split_test(x, g, "Aged_x", rng, splits=6, permutations=60)
    assert r["directed_reversal"] - r["directed_null"] > 0.3 and r["directed_p"] < 0.05
    assert r["noise_p"] > 0.05


def test_noise_reversal_is_found_and_directed_left():
    rng = np.random.default_rng(17)
    x, g = _cohort(rng, directed_back=0.0, noise_back=0.8)
    r = split_test(x, g, "Aged_x", rng, splits=6, permutations=60)
    assert r["noise_p"] < 0.05 and r["noise_reversal"] > r["noise_null"]
    assert r["directed_p"] > 0.05


def test_an_inert_treatment_is_not_credited():
    rng = np.random.default_rng(18)
    x, g = _cohort(rng, directed_back=0.0, noise_back=0.0)
    r = split_test(x, g, "Aged_x", rng, splits=6, permutations=60)
    assert r["directed_p"] > 0.05 and r["noise_p"] > 0.05
