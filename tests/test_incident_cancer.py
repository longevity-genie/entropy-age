"""RM28: the logistic fit and the outcome definitions of EPIC-Italy."""

import numpy as np
import polars as pl
import pytest

from entropy_age.experiments.incident_cancer import logit, outcome_rows


def test_logit_recovers_a_simulated_odds_ratio():
    rng = np.random.default_rng(14)
    n = 5000
    z = rng.normal(size=n)
    y = (rng.uniform(size=n) < 1 / (1 + np.exp(-(-0.3 + 0.5 * z)))).astype(float)
    beta, se = logit(y, np.column_stack([np.ones(n), z]))
    assert beta[1] == pytest.approx(0.5, abs=4 * se[1])


def test_outcomes_split_by_site_and_lag():
    rng = np.random.default_rng(15)
    n = 300
    icd = np.array([""] * 150 + ["C50"] * 100 + ["C18"] * 50, dtype=object)
    ttd = np.where(icd == "", np.nan, rng.uniform(0, 10, n))
    people = pl.DataFrame(
        {
            "age": rng.uniform(40, 70, n),
            "female": np.ones(n, dtype=bool),
            "icd": [c or None for c in icd],
            "years_to_diagnosis": ttd,
        }
    )
    score = (icd == "C50") * 1.0 + rng.normal(0, 1, n)  # higher in breast cancer cases only
    rows = {(r["outcome"], r["lag_over_2y"]): r for r in outcome_rows(people, {"m": score})}
    assert rows[("breast", False)]["cases"] == 100 and rows[("colorectal", False)]["cases"] == 50
    assert rows[("breast", True)]["cases"] == int(((icd == "C50") & (ttd > 2)).sum())
    assert rows[("breast", False)]["odds_ratio"] > 2 and rows[("breast", False)]["p"] < 1e-6
    assert rows[("colorectal", False)]["p"] > 0.01
