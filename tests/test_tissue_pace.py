"""RM25: tissue paces from result tables and the agreement of their orders."""

import polars as pl
import pytest

from entropy_age.experiments.tissue_pace import agreement, paces


def _table(values: dict[str, float], kind: str = "abs_slope", **extra) -> pl.DataFrame:
    n = len(values)
    return pl.DataFrame(
        {
            "metric": [kind] * n + ["entropy_slope"] * n,
            "group": ["housekeeping"] * (2 * n),
            "tissue": list(values) * 2,
            "median_rest": list(values.values()) + [9.0] * n,
            **{k: [v] * (2 * n) for k, v in extra.items()},
        }
    )


def test_paces_keep_the_kind_and_filter():
    t = pl.concat([_table({"a": 1.0, "b": 2.0}, species="mouse"), _table({"a": 5.0, "b": 6.0}, species="rat")])
    p = paces(t, {"species": "rat"})
    assert p.sort("tissue")["pace"].to_list() == [5.0, 6.0]
    assert p["kind"].unique().to_list() == ["abs_slope"]


def test_agreement_is_spearman_over_shared_tissues_only():
    a = paces(_table({"a": 1.0, "b": 2.0, "c": 3.0, "d": 4.0, "x": 0.0}), {})
    b = paces(_table({"a": 10.0, "b": 20.0, "c": 30.0, "d": 40.0, "y": 0.0}), {})
    (r,) = agreement(a, b)
    assert r["tissues"] == 4 and r["rho"] == pytest.approx(1.0)
    reversed_ = paces(_table({"a": 4.0, "b": 3.0, "c": 2.0, "d": 1.0}), {})
    assert agreement(a, reversed_)[0]["rho"] == pytest.approx(-1.0)


def test_too_few_shared_tissues_give_no_row():
    a = paces(_table({"a": 1.0, "b": 2.0, "c": 3.0}), {})
    assert agreement(a, a) == []
