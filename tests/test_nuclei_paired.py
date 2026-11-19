import numpy as np

from entropy_age.experiments.nuclei_paired import paired_difference


def test_paired_difference_scales_by_the_full_sample_sd_and_tests_the_difference():
    key = ("kl_slope", "housekeeping")
    glia, neurons = {key: (-0.2, 2.0)}, {key: (-0.05, 1.0)}
    rng = np.random.default_rng(1)
    boots = [({key: (-0.2 + e, 9.9)}, {key: (-0.05 + e / 2, 9.9)}) for e in rng.normal(0, 0.01, 200)]
    (row,) = paired_difference(glia, neurons, boots)
    assert np.isclose(row["gap_glia"], -0.1) and np.isclose(row["gap_neurons"], -0.05)  # each over its own full SD
    assert np.isclose(row["difference"], -0.05)
    assert row["ci_low"] < -0.05 < row["ci_high"] < 0  # shared noise cancels in the paired difference
    assert row["p_bootstrap"] <= 0.01 and row["resamples"] == 200
