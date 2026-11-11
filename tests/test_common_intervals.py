"""RM9: resamples rescaled by the full sample's spread, not their own."""

import polars as pl
import pytest

from entropy_age.experiments.common import rescale


def test_rescale_divides_every_resample_by_the_full_sample_sd():
    point = pl.DataFrame(
        {"tissue": ["a", "a", "b"], "metric": ["m", "m", "m"], "group": ["x", "y", "x"], "sd": [2.0, 2.0, 4.0]}
    )
    boot = pl.DataFrame(
        {"tissue": ["a", "b", "a"], "metric": ["m"] * 3, "group": ["x", "x", "y"], "matched_diff": [1.0, 1.0, -2.0]}
    )
    out = rescale(boot, point, ["tissue", "metric"])
    assert out["matched_diff_sd"].to_list() == [0.5, 0.25, -1.0]
    assert "sd" not in out.columns


def test_rescale_refuses_a_resample_without_a_full_sample_sd():
    point = pl.DataFrame({"tissue": ["a"], "metric": ["m"], "sd": [2.0]})
    boot = pl.DataFrame({"tissue": ["b"], "metric": ["m"], "matched_diff": [1.0]})
    with pytest.raises(ValueError):
        rescale(boot, point, ["tissue", "metric"])
