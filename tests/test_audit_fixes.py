"""Fixes from the outside audit of 2026-11-18."""

import numpy as np
import polars as pl
import pytest

from entropy_age.experiments.common import rescale
from entropy_age.experiments.global_fdr import p_columns
from entropy_age.stochasticity import baseline_strata, drift_metrics, relative_entropy


@pytest.mark.parametrize("ref", [0.0, 0.001, 0.01, 0.5, 0.99, 0.999, 1.0])
def test_relative_entropy_is_zero_at_the_reference_and_rises_away(ref):
    assert relative_entropy(np.array([ref]), np.array([ref]))[0] == pytest.approx(0.0, abs=1e-12)
    for direction in (-1, 1):
        steps = np.clip(ref + direction * np.array([0.0, 0.002, 0.005, 0.02, 0.1]), 0, 1)
        d = relative_entropy(steps, np.full(steps.size, ref))
        assert np.all(np.diff(d) >= -1e-12)


def test_kl_slope_is_positive_for_a_site_leaving_the_floor():
    rng = np.random.default_rng(0)
    ages = np.linspace(20, 80, 200)
    x = np.clip(0.001 + 0.004 * (ages - 20) / 60 + rng.normal(0, 1e-4, ages.size), 0, 1)[:, None]
    assert drift_metrics(x, ages, ["c"], spearman=False)["kl_slope"][0] > 0


def test_baseline_strata_keep_matching_with_a_missing_value():
    v = np.array([*np.linspace(0, 1, 10), np.nan])
    s = baseline_strata(v, 5)
    assert s[-1] == 5 and len(set(s[:-1])) == 5


def test_rescale_refuses_conflicting_denominators():
    point = pl.DataFrame({"t": ["a", "a"], "sd": [1.0, 2.0]})
    boot = pl.DataFrame({"t": ["a"], "matched_diff": [0.5]})
    with pytest.raises(ValueError):
        rescale(boot, point, ["t"])


def test_p_columns_read_every_raw_p_and_no_adjusted_one():
    cols = ["p", "p_permutation", "directed_p", "B_p_welch", "spearman_p", "p_bh", "q_global", "adj_p", "sd", "P_core_essential"]
    assert p_columns(cols) == ["p", "p_permutation", "directed_p", "B_p_welch", "spearman_p"]
