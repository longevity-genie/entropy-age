"""RM26: coupling of two tissues' accelerations over shared donors."""

import numpy as np
import polars as pl
import pytest

from entropy_age.experiments.organ_coupling import MIN_SHARED, coupling


def _tissue(donors, values):
    return pl.DataFrame({"donor": donors, "acc": values})


def test_coupling_pairs_by_donor_not_by_row_order():
    rng = np.random.default_rng(3)
    pace = rng.normal(size=60)
    donors = [f"d{i}" for i in range(60)]
    a = _tissue(donors, pace + rng.normal(0, 0.5, 60))
    order = rng.permutation(60)
    b = _tissue([donors[i] for i in order], (pace + rng.normal(0, 0.5, 60))[order])
    r = coupling(a, b, "acc")
    assert r["n"] == 60 and r["r"] > 0.6 and r["low"] < r["r"] < r["high"]


def test_independent_tissues_give_an_interval_around_zero_and_few_donors_none():
    rng = np.random.default_rng(4)
    donors = [f"d{i}" for i in range(200)]
    r = coupling(_tissue(donors, rng.normal(size=200)), _tissue(donors, rng.normal(size=200)), "acc")
    assert r["low"] < 0 < r["high"]
    few = donors[: MIN_SHARED - 1]
    assert coupling(_tissue(few, np.arange(len(few)) * 1.0), _tissue(few, np.arange(len(few)) * 1.0), "acc") is None


def test_missing_values_are_dropped_per_pair():
    donors = [f"d{i}" for i in range(30)]
    a = _tissue(donors, [float(i) for i in range(30)])
    b = _tissue(donors, [None] * 5 + [float(i) for i in range(5, 30)])
    assert coupling(a, b, "acc")["n"] == 25
    assert coupling(a, b, "acc")["r"] == pytest.approx(1.0)
