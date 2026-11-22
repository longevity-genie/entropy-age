import polars as pl

from entropy_age.experiments.astrocyte_passage import bulk_cultures, htert_cultures


def test_bulk_cultures_keeps_unsorted_cultures_with_doublings():
    chars = pl.DataFrame(
        {
            "sample": ["a", "b", "c", "d"],
            "FACS.code": ["Bulk", "Sen", "Bulk", "Bulk"],
            "donor": ["1", "1", "2", "1"],
            "cpd": ["3.5", "9", "NA", "0"],
        }
    )
    kept = bulk_cultures(chars)
    assert kept["sample"].to_list() == ["d", "a"]
    assert kept["doublings"].to_list() == [0.0, 3.5]


def test_htert_cultures_keeps_immortalised_cultures_by_replicate_line():
    chars = pl.DataFrame(
        {
            "sample": ["a", "b", "c", "d"],
            "htert": ["1", "0", "1", "1"],
            "replicate": ["2", "3", "1", "1"],
            "cpd": ["45.8", "6.0", "72.2", "29.5"],
        }
    )
    kept = htert_cultures(chars)
    assert kept["sample"].to_list() == ["d", "c", "a"]
    assert kept["donor"].to_list() == ["line 1", "line 1", "line 2"]


def test_culture_bootstrap_brackets_a_planted_gap():
    import numpy as np

    from entropy_age.experiments.astrocyte_passage import culture_bootstrap
    from entropy_age.experiments.common import dummies
    from entropy_age.experiments.tissue_contrast import Prepared

    rng = np.random.default_rng(1)
    n, m = 24, 400
    doublings = np.tile(np.linspace(0, 20, n // 2), 2)
    lines = np.repeat(["a", "b"], n // 2)
    slope = np.where(np.arange(m) < 100, 0.002, 0.006)  # the set drifts a third as fast
    x = np.clip(0.2 + doublings[:, None] * slope[None, :] + rng.normal(0, 0.01, (n, m)), 0.01, 0.99)
    cov = dummies(lines)
    prepared = Prepared("t", x, doublings, cov, [f"c{i}" for i in range(m)], lines, cov)
    in_set = np.arange(m) < 100
    out = culture_bootstrap(prepared, lines, {"set": (in_set, None)}, np.repeat("s", m), 20, rng)
    assert out[0]["resamples"] == 20 and out[0]["ci_high"] < 0
