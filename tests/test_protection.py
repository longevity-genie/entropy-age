"""The per-person protection score: contrast weights, the reference's composition fit, the gene resampling."""

import numpy as np
import polars as pl
import pytest

from entropy_age.experiments.per_person import age_fits, lifespan_scores, paired_turnover
from entropy_age.experiments.tissue_contrast import Prepared
from entropy_age.protection import (
    contrast_weights,
    deviations,
    fit_reference,
    matched_scores,
    resampled_gaps,
)


def test_weights_give_the_set_minus_matched_mean_per_stratum():
    d = np.array([[1.0, 5.0, 3.0, 2.0, 6.0, 9.0]])
    in_set = np.array([1.0, 1.0, 0.0, 1.0, 0.0, 0.0])
    strata = np.array(["a", "a", "a", "b", "b", "c"])  # c holds no set feature and must not count
    p, b = contrast_weights(in_set, strata)
    assert (d @ p)[0] == pytest.approx(2 / 3 * (3 - 3) + 1 / 3 * (2 - 6))  # weights: the set's share per stratum
    assert (d @ b)[0] == pytest.approx(2 / 3 * 3 + 1 / 3 * 6)


def test_a_reference_component_is_removed_from_a_sample():
    rng = np.random.default_rng(1)
    n, p = 60, 400
    loading = rng.normal(0, 0.05, p)
    level = rng.uniform(0.2, 0.8, p)
    factor = rng.normal(0, 1, n)
    x = level + np.outer(factor, loading) + rng.normal(0, 0.005, (n, p))
    background = np.ones(p, dtype=bool)
    island = np.zeros(p, dtype=bool)
    sample = level + 3.0 * loading + rng.normal(0, 0.005, (1, p))
    raw = deviations(sample, fit_reference(x, island, background).level).mean()
    reference = fit_reference(x, island, background, components=1)
    adjusted = deviations(reference.adjust(sample), reference.level).mean()
    assert adjusted < 0.3 * raw


def test_gene_resamples_centre_on_the_point_gap_and_vanish_when_genes_are_alike():
    rng = np.random.default_rng(2)
    d = rng.gamma(2.0, 1.0, (4, 300))
    strata = np.repeat(["s0", "s1", "s2"], 100)
    genes = [np.arange(j, j + 3) for j in range(0, 60, 3)] + [np.arange(j, j + 3) for j in range(100, 160, 3)]
    in_set = np.zeros(300, dtype=bool)
    for g in genes:
        in_set[g] = True
    point, _ = matched_scores(d, in_set, strata)
    draws = resampled_gaps(d, genes, strata, 400, np.random.default_rng(3))
    assert np.abs(draws.mean(axis=0) - point).max() < 4 * draws.std(axis=0).max() / np.sqrt(400) + 1e-9
    assert (draws.std(axis=0) > 1e-3).all()
    flat = np.where(in_set, 1.0, 2.0) * np.ones_like(d)  # every gene alike: only the weights could move P
    assert resampled_gaps(flat, genes, strata, 50, np.random.default_rng(4)).std(axis=0).max() < 1e-12


def _cohort(n_young: int, protect: float) -> tuple[Prepared, pl.DataFrame, pl.DataFrame]:
    """Set CpGs drift ``protect`` times as fast as the others; ages 20-90."""
    rng = np.random.default_rng(5)
    p = 600
    features = [f"cg{j}" for j in range(p)]
    ages = np.r_[rng.uniform(20, 35, n_young), rng.uniform(40, 90, 120)]
    level = rng.uniform(0.1, 0.9, p)
    in_set = np.arange(p) < 150
    rate = np.where(in_set, protect, 1.0) * 0.002
    x = np.clip(level + rng.normal(0, 1, (ages.size, p)) * rate * (ages[:, None] - 15), 0, 1)
    probes = pl.DataFrame({"cpg": features, "cpg_island": [j % 2 == 0 for j in range(p)]})
    membership = pl.DataFrame(
        {"feature": features[:150], "group": "housekeeping", "gene": [f"g{j // 3}" for j in range(150)]}
    )
    prepared = Prepared("blood", x, ages, None, features, np.arange(ages.size), None)
    return prepared, probes, membership


def test_the_score_falls_with_age_where_the_set_drifts_less():
    p, probes, membership = _cohort(30, protect=0.4)
    scores = lifespan_scores(p, probes, membership, 0, 30, np.random.default_rng(6))
    assert scores.height == 120 and (scores["age"] > 35).all()
    fit = {r["measure"]: r for r in age_fits(scores, None)}
    assert fit["P"]["slope_per_decade"] < 0 and fit["P"]["p"] < 0.01
    assert 0 < fit["P"]["reliability"] <= 1


def test_too_few_young_samples_is_refused():
    p, probes, membership = _cohort(10, protect=0.4)
    with pytest.raises(ValueError, match="fewer than 20"):
        lifespan_scores(p, probes, membership, 0, 0, None)


def test_turnover_counts_donors_with_blood_below_each_region():
    rows = []
    for donor, (blood, cortex) in enumerate([(-0.3, -0.1), (-0.2, -0.25), (-0.4, 0.0)]):
        for tissue, f in (("whole blood", blood), ("frontal cortex", cortex)):
            rows.append({"donor": donor, "tissue": tissue, "components": 0, "F_housekeeping": f})
    table = paired_turnover(pl.DataFrame(rows))
    row = table.filter((pl.col("region") == "frontal cortex") & (pl.col("measure") == "F")).row(0, named=True)
    assert (row["donors"], row["blood_below"]) == (3, 2)
    assert row["median_diff"] == pytest.approx(-0.2)
