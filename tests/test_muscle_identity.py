import numpy as np
import polars as pl

from entropy_age.experiments.muscle_identity import (
    axis_projection,
    band_gaps,
    differential,
    identity_gene_cpgs,
    matched_partners,
    null_projection,
    slope_interval,
)


def test_axis_projection_places_the_ends_at_zero_and_one():
    origin = np.array([0.1, 0.9, 0.5, 0.2])
    target = np.array([0.8, 0.2, 0.5, 0.6])
    x = np.vstack([origin, target, (origin + target) / 2, origin + 0.25 * (target - origin)])
    assert np.allclose(axis_projection(x, origin, target), [0.0, 1.0, 0.5, 0.25])


def test_axis_projection_ignores_missing_values():
    origin, target = np.array([0.1, 0.9, 0.3]), np.array([0.8, 0.2, np.nan])
    x = np.array([[0.8, 0.2, 0.0], [np.nan, 0.2, 0.3]])
    assert np.allclose(axis_projection(x, origin, target), [1.0, 1.0])


def test_differential_needs_threshold_and_confirmation_of_sign():
    a = np.array([0.9, 0.9, 0.9, 0.5, np.nan])
    b = np.array([0.1, 0.1, 0.1, 0.4, 0.1])
    confirm = np.array([0.5, -0.5, 0.05, 0.5, 0.5])
    assert differential(a, b, 0.3).tolist() == [True, True, True, False, False]
    assert differential(a, b, 0.3, confirm, 0.1).tolist() == [True, False, False, False, False]


def test_identity_gene_cpgs_keeps_promoters_and_enhancers_of_listed_genes():
    links = pl.DataFrame(
        {
            "cpg": ["c1", "c2", "c3", "c4", "c1"],
            "symbol": ["DES", "DES", "DES", "GAPDH", "DES"],
            "region": ["TSS200", "Body", "Body", "TSS200", "Body"],
        }
    )
    out = identity_gene_cpgs(links, {"c2", "c4"})
    assert out["cpg"].to_list() == ["c1", "c2"]
    assert dict(out.select("cpg", "kind").iter_rows()) == {"c1": "promoter", "c2": "enhancer"}


def test_slope_interval_recovers_a_known_slope():
    rng = np.random.default_rng(1)
    age = rng.uniform(20, 80, 200)
    sex = rng.integers(0, 2, 200).astype(float)
    y = 0.003 * age + 0.1 * sex + rng.normal(0, 0.02, 200)
    r = slope_interval(y, age, sex[:, None], 200, np.random.default_rng(2))
    assert abs(r["slope"] - 0.003) < 0.0005
    assert r["ci_low"] < 0.003 < r["ci_high"] and r["p"] < 1e-10


def test_null_projection_follows_the_matched_partners():
    origin = np.array([0.1, 0.9, 0.1, 0.9])
    target = np.array([0.6, 0.4, 0.1, 0.9])
    x = np.array([[0.2, 0.9, 0.6, 0.9]])  # the first background CpG moved fully toward the target direction
    strata = np.array(["a", "b", "a", "b"])
    in_axis = np.array([True, True, False, False])
    partners = matched_partners(strata, in_axis, ~in_axis, 5, np.random.default_rng(0))
    assert (partners == [[2, 3]] * 5).all()
    assert np.allclose(null_projection(x, origin, target, partners, np.flatnonzero(in_axis)), [0.5])


def test_band_gaps_sees_a_set_that_drifts_faster_than_matched_cpgs():
    rng = np.random.default_rng(3)
    n = 4000
    baseline = rng.uniform(0.05, 0.95, n)
    in_set = np.zeros(n, dtype=bool)
    in_set[:400] = True
    metrics = pl.DataFrame(
        {m: rng.normal(0, 1, n) + 2.0 * in_set for m in ("kl_slope", "dispersion_slope", "entropy_slope")}
    )
    metrics = metrics.with_columns(baseline=pl.Series(baseline))
    rows = band_gaps(metrics, rng.integers(0, 2, n).astype(bool), {"set": in_set})
    gap = next(r for r in rows if r["band"] == "all" and r["metric"] == "kl_slope")["gap_sd"]
    assert gap > 1.0
