"""RM16 step 2: the marrow anomaly checks."""

import numpy as np

from entropy_age.experiments.marrow_anomaly import donor_key, gaps, mix_correlation


def _cohort(rng, set_rate: float, n: int = 40, p: int = 400):
    age = np.sort(rng.uniform(30, 80, n))
    level = rng.uniform(0.05, 0.95, p)
    rate = np.full(p, 0.004)
    in_set = np.zeros(p, bool)
    in_set[:60] = True
    rate[in_set] = set_rate
    noise = rng.normal(size=(n, p)) * (rate * (age[:, None] - 30))
    x = np.clip(level + noise, 0.001, 0.999)
    return x, age, in_set


def test_gap_is_negative_when_the_set_drifts_less():
    rng = np.random.default_rng(1)
    x, age, in_set = _cohort(rng, 0.001)
    out = gaps(x, age, np.zeros(x.shape[1], bool), {"housekeeping": (in_set, [])})
    assert out[("dispersion_slope", "housekeeping")] < 0
    x, age, in_set = _cohort(rng, 0.004)
    flat = gaps(x, age, np.zeros(x.shape[1], bool), {"housekeeping": (in_set, [])})
    assert flat[("dispersion_slope", "housekeeping")] > out[("dispersion_slope", "housekeeping")]


def test_donors_match_by_age_and_sex():
    keys = donor_key(np.array([38.3, 38.3, 41.5]), np.array(["0", "1", "0"]))
    assert keys.tolist() == ["38.3|0", "38.3|1", "41.5|0"]


def test_mix_correlation_ignores_missing():
    s = np.array([1.0, 2.0, np.nan, 4.0])
    d = np.array([2.0, 4.0, 1.0, 8.0])
    assert abs(mix_correlation(s, d) - 1.0) < 1e-12
