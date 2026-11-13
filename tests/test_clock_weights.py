"""R64: published clocks' weights on protected CpGs."""

import numpy as np
import polars as pl

from entropy_age.experiments.clock_weights import clock_summary, cpg_stats


def test_cpg_stats_slope_and_spread():
    rng = np.random.default_rng(0)
    ages = np.linspace(20, 80, 200)
    rising = 0.01 * ages + rng.normal(0, 0.001, ages.size)
    spreading = 0.5 + rng.normal(0, 1, ages.size) * ages * 1e-3
    st = cpg_stats(np.column_stack([rising, spreading]), ages)
    assert abs(st["slope"][0] - 0.1) < 0.005
    assert st["spread_slope"][1] > 5 * abs(st["spread_slope"][0])
    assert st["young"][0] < 0.4


def test_weight_beyond_slope_detects_extra_weight():
    rng = np.random.default_rng(1)
    n = 300
    prot = np.arange(n) < 60
    slope = np.exp(rng.normal(-3, 0.5, n))
    weight = 0.01 / slope * np.where(prot, 2.0, 1.0) * np.exp(rng.normal(0, 0.1, n))
    g = pl.DataFrame(
        {
            "weight": weight,
            "slope": slope,
            "spread_slope": rng.normal(0, 1, n),
            "young": rng.uniform(0, 1, n),
            "protected": prot,
        }
    )
    r = clock_summary(g, 50, np.random.default_rng(2))
    assert 1.8 < r["weight_beyond_slope"] < 2.2 and r["ci_low"] > 1.5 and r["p"] < 1e-6
    assert r["protected_cpgs"] == 60


def test_too_few_protected_is_skipped():
    g = pl.DataFrame(
        {
            "weight": [1.0] * 8,
            "slope": [0.1] * 8,
            "spread_slope": [0.0] * 8,
            "young": [0.2] * 8,
            "protected": [True] * 2 + [False] * 6,
        }
    )
    assert clock_summary(g, 10, np.random.default_rng(0)) is None
