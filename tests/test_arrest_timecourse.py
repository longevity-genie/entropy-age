import numpy as np

from entropy_age.experiments.arrest_timecourse import drift_test, slopes


def test_slopes_and_drift_test_separate_a_trend_from_noise():
    rng = np.random.default_rng(3)
    t = np.arange(6.0)
    noise = rng.normal(0, 0.01, (6, 4000))
    assert np.allclose(slopes(np.outer(t, [0.5, -1.0]), t), [0.5, -1.0])
    flat = drift_test(0.5 + noise, t)
    trend = drift_test(0.5 + noise + np.outer(t, np.full(4000, 0.01)), t)
    assert flat["reorderings"] == 720 and flat["p_reorder"] > 0.05
    assert trend["p_reorder"] <= 2 / 720 and trend["ratio_to_null"] > 5  # only the true order and its mirror reach it
