import numpy as np
import polars as pl

from entropy_age.experiments.senescence_modes import paired, sign_flip


def test_sign_flip_is_exact_and_paired_takes_within_chip_differences():
    assert sign_flip(np.array([1.0, 1.0, 1.0, 1.0])) == 2 / 16  # only all-plus and all-minus reach the observed mean
    values = pl.DataFrame({"chip": ["a", "a", "b", "b", "c", "c", "d", "d"], "arm": ["R", "X"] * 4,
                           "F_h": [-0.3, -0.1, -0.4, -0.2, -0.2, 0.0, -0.5, -0.3]})  # fmt: skip
    row = paired(values, "R", "X", "F_h")
    assert row["chips"] == 4 and np.isclose(row["difference"], -0.2) and row["chips_negative"] == 4


def test_noise_floor_removes_pure_noise_from_the_squared_scores():
    from entropy_age.experiments.senescence_modes import noise_floor, sample_scores

    rng = np.random.default_rng(0)
    n = 40_000
    truth = np.full(n, 0.5)
    controls = truth + rng.normal(0, 0.02, (2, n))
    treated = truth + rng.normal(0, 0.02, (1, n))  # no real drift at all
    reference = controls.mean(axis=0)
    in_set = (np.arange(n) % 4 == 0).astype(float)
    sets = {"h": (in_set, None)}
    strata = np.zeros(n, dtype=int).astype(str)
    scores = sample_scores(treated, reference, strata, sets, noise_floor(controls, reference))
    assert scores["B_h"][0] > 0.03  # the uncorrected deviation is all noise
    assert abs(scores["B2_h"][0]) < 2e-4  # and the floor removes it
