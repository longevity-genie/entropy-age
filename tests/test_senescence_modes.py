import numpy as np
import polars as pl

from entropy_age.experiments.senescence_modes import paired, sign_flip


def test_sign_flip_is_exact_and_paired_takes_within_chip_differences():
    assert sign_flip(np.array([1.0, 1.0, 1.0, 1.0])) == 2 / 16  # only all-plus and all-minus reach the observed mean
    values = pl.DataFrame({"chip": ["a", "a", "b", "b", "c", "c", "d", "d"], "arm": ["R", "X"] * 4,
                           "F_h": [-0.3, -0.1, -0.4, -0.2, -0.2, 0.0, -0.5, -0.3]})  # fmt: skip
    row = paired(values, "R", "X", "F_h")
    assert row["chips"] == 4 and np.isclose(row["difference"], -0.2) and row["chips_negative"] == 4
