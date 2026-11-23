import polars as pl

from entropy_age.experiments.series import groups_of, sample_table, specs


def test_every_description_file_parses_and_names_its_columns():
    for name, spec in specs().items():
        assert {"accession", "age", "groups"} <= set(spec), name
        assert all(isinstance(filters, dict) for filters in spec["groups"].values()), name


def test_groups_keep_one_passing_array_per_donor():
    chars = pl.DataFrame(
        {
            "sample": ["a", "b", "c", "d", "e"],
            "age": ["30", "40", "50", "x", "60"],
            "Sex": ["F", "F", "M", "M", "M"],
            "dx": ["Control", "Control", "Case", "Control", "Control"],
        }
    )
    spec = {
        "age": "age",
        "sex": "Sex",
        "donor": "title:^(\\w+)_",
        "hold": ["dx"],
        "groups": {"controls": {"dx": "Control"}, "all": {}},
    }
    table = sample_table(chars, ["D1_x", "D1_y", "D2_x", "D3_x", "D4_x"], spec)
    assert table["donor"].to_list() == ["D1", "D1", "D2", "D4"]  # d has no age
    groups = groups_of(table, spec, usable=["b", "c", "e"], exclude={"D4"})
    assert groups["controls"]["sample"].to_list() == ["b"]  # a failed detection, D4 is excluded
    assert groups["all"]["sample"].to_list() == ["b", "c"]


def test_donor_interval_brackets_the_point_estimate():
    import numpy as np

    from entropy_age.experiments.series import donor_intervals
    from entropy_age.protection import matched_scores
    from entropy_age.stochasticity import drift_metrics

    rng = np.random.default_rng(3)
    donors, cpgs = 80, 600
    ages = rng.uniform(20, 80, donors)
    in_set = np.arange(cpgs) < 150
    drift = np.where(in_set, 0.0004, 0.0012)  # the set gains spread a third as fast
    x = np.clip(0.3 + rng.normal(0, 1, (donors, cpgs)) * (0.02 + drift * (ages[:, None] - 20)), 0.01, 0.99)
    feats = [f"c{i}" for i in range(cpgs)]
    island = np.zeros(cpgs, bool)
    cov = np.empty((donors, 0))
    out = donor_intervals(x, ages, cov, feats, island, {"set": (in_set, None)}, 60, rng)
    lo, hi, kept = out[("all", "dispersion_slope", "set")]
    assert kept == 60
    full = drift_metrics(x, ages, feats, covariates=None, spearman=False)
    values = full["dispersion_slope"].to_numpy()
    strata = np.char.add(
        "False|",
        __import__("entropy_age.stochasticity", fromlist=["x"])
        .baseline_strata(full["baseline"].to_numpy(), 10)
        .astype(str),
    )
    point = float(matched_scores(values[None, :], in_set, strata)[0][0]) / float(values.std())
    assert lo < point < hi and point < 0
    assert abs((lo + hi) / 2 - point) < 0.35 * abs(point)  # not pulled toward zero


def test_sample_table_refuses_titles_without_a_donor():
    import pytest

    chars = pl.DataFrame({"sample": ["a", "b"], "age": ["30", "40"], "Sex": ["F", "M"]})
    spec = {"age": "age", "sex": "Sex", "donor": "title:^(D\\d+)_", "groups": {"all": {}}}
    with pytest.raises(ValueError, match="do not match"):
        sample_table(chars, ["D1_x", "other"], spec)
