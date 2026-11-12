"""RM11: a clock's linear score and the partial correlation of shield and acceleration."""

import numpy as np
import polars as pl
import pytest

from entropy_age.experiments.shield_clocks import linear_score, partial_r


def test_linear_score_uses_present_cpgs_and_fills_gaps_with_the_mean():
    x = np.array([[0.1, 0.5, np.nan], [0.3, 0.7, 0.4], [0.2, 0.6, 0.6]])
    w = pl.DataFrame({"cpg": ["a", "c", "zz"], "weight": [2.0, 1.0, 9.0]})
    score, used = linear_score(x, ["a", "b", "c"], w)
    assert used == 2
    assert score.tolist() == pytest.approx([0.2 + 0.5, 0.6 + 0.4, 0.4 + 0.6])


def test_partial_r_finds_a_link_beyond_age_and_none_where_only_age_links():
    rng = np.random.default_rng(0)
    n = 500
    age = rng.uniform(20, 80, n)
    shared = rng.normal(size=n)
    a = 0.05 * age + shared + rng.normal(size=n)
    b = 0.03 * age + shared + rng.normal(size=n)
    design = np.column_stack([np.ones(n), age])
    r, p, _ = partial_r(a, b, design)
    assert r > 0.4 and p < 1e-10
    c = 0.05 * age + rng.normal(size=n)
    d = 0.03 * age + rng.normal(size=n)
    r0, p0, _ = partial_r(c, d, design)
    assert abs(r0) < 0.12


def test_linear_score_skips_a_cpg_no_sample_has():
    """A clock CpG absent from the array (all NaN) must not turn every score into NaN (R42)."""
    x = np.array([[0.1, np.nan], [0.3, np.nan]])
    w = pl.DataFrame({"cpg": ["a", "b"], "weight": [2.0, 5.0]})
    score, used = linear_score(x, ["a", "b"], w)
    assert used == 1 and np.allclose(score, [0.2, 0.6])
