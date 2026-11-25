import numpy as np
import polars as pl

from entropy_age.experiments import mitotic_clocks


def _coef():
    return {
        "epiTOC2": pl.DataFrame({"cpg": ["a", "b"], "delta": [0.01, 0.02], "beta0": [0.0, 0.5]}),
        "RepliTali": pl.DataFrame({"term": ["(Intercept)", "a", "c"], "coefficient": [10.0, 2.0, -4.0]}),
        "epiCMIT": pl.DataFrame({"cpg": ["a", "c"], "class": ["hyper", "hypo"]}),
    }


def test_formulas():
    betas = pl.DataFrame({"cpg": ["a", "b", "c"], "s1": [0.1, 0.6, 0.8], "s2": [0.3, 0.5, 0.2]})
    out = mitotic_clocks.score(betas, ["s1", "s2"], _coef()).pivot(on="clock", index="sample", values="prediction")
    s1 = out.filter(pl.col("sample") == "s1")
    # epiTOC2: 2 * mean((0.1-0)/(0.01*1), (0.6-0.5)/(0.02*0.5)) = 2 * mean(10, 10) = 20
    assert np.isclose(s1["epiTOC2"][0], 20.0)
    assert np.isclose(s1["RepliTali"][0], 10 + 2 * 0.1 - 4 * 0.8)
    assert np.isclose(s1["epiCMIT"][0], max(0.1, 1 - 0.8))


def test_coverage_rule_drops_a_clock():
    betas = pl.DataFrame({"cpg": ["a", "b"], "s1": [0.1, 0.6]})  # RepliTali's heavier CpG c is absent
    out = mitotic_clocks.score(betas, ["s1"], _coef())
    assert "RepliTali" not in out["clock"].to_list()
    assert "epiTOC2" in out["clock"].to_list()
