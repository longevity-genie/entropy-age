import numpy as np
import polars as pl

from entropy_age.experiments import mitotic_clocks


def _coef():
    toc = pl.DataFrame({"cpg": ["a", "b"], "delta": [0.01, 0.02], "beta0": [0.0, 0.5]})
    return {
        "epiTOC2": toc,
        "epiTOC3": toc,
        "HypoClock": pl.DataFrame({"cpg": ["a", "c"]}),
        "RepliTali": pl.DataFrame({"term": ["(Intercept)", "a", "c"], "coefficient": [10.0, 2.0, -4.0]}),
        "epiCMIT": pl.DataFrame({"cpg": ["a", "c"], "class": ["hyper", "hypo"]}),
        "EpiTOC1": pl.DataFrame({"cpg": ["z"]}),
    }


def _wide(scores: pl.DataFrame) -> pl.DataFrame:
    return scores.pivot(on="clock", index="sample", values="prediction")


def test_formulas():
    betas = pl.DataFrame({"cpg": ["a", "b", "c"], "s1": [0.1, 0.6, 0.8], "s2": [0.3, 0.5, 0.2]})
    scores, cover = mitotic_clocks.score(betas, ["s1", "s2"], _coef())
    s1 = _wide(scores).filter(pl.col("sample") == "s1")
    # epiTOC2: 2 * mean((0.1-0)/(0.01*1), (0.6-0.5)/(0.02*0.5)) = 2 * mean(10, 10) = 20
    assert np.isclose(s1["epiTOC2"][0], 20.0)
    assert np.isclose(s1["HypoClock"][0], (0.1 + 0.8) / 2)
    assert np.isclose(s1["RepliTali"][0], 10 + 2 * 0.1 - 4 * 0.8)
    assert np.isclose(s1["epiCMIT"][0], max(0.1, 1 - 0.8))
    assert cover["passed"].all()


def test_coverage_rule_drops_a_clock():
    betas = pl.DataFrame({"cpg": ["a", "b"], "s1": [0.1, 0.6]})  # RepliTali's heavier CpG c is absent
    scores, cover = mitotic_clocks.score(betas, ["s1"], _coef())
    assert "RepliTali" not in scores["clock"].to_list()
    assert "epiTOC2" in scores["clock"].to_list()
    assert not cover.filter(pl.col("clock") == "RepliTali")["passed"][0]


def test_coverage_is_judged_per_array():
    # c is measured in one array of ten: pooled availability would pass, per-array shares do not
    n = 10
    cols = {f"s{i}": [0.1, 0.6, (0.5 if i == 0 else None)] for i in range(n)}
    betas = pl.DataFrame({"cpg": ["a", "b", "c"], **cols})
    scores, _ = mitotic_clocks.score(betas, list(cols), _coef())
    assert "RepliTali" not in scores["clock"].to_list()
    assert "epiCMIT" not in scores["clock"].to_list()  # its hypo arm (c) fails on its own


def test_cpgs_include_epitoc1():
    assert {"a", "b", "c", "z"} == mitotic_clocks.cpgs(_coef())


def test_array_below_its_share_gets_no_score():
    # eleven arrays measure everything; one misses b, half of epiTOC2's weight: the group passes, that array is NaN
    cols = {f"s{i}": [0.1, 0.6, 0.5] for i in range(11)}
    cols["bad"] = [0.1, None, 0.5]
    betas = pl.DataFrame({"cpg": ["a", "b", "c"], **cols})
    scores, cover = mitotic_clocks.score(betas, list(cols), _coef())
    toc = scores.filter(pl.col("clock") == "epiTOC2")
    assert np.isnan(toc.filter(pl.col("sample") == "bad")["prediction"][0])
    assert np.isfinite(toc.filter(pl.col("sample") == "s0")["prediction"][0])
    assert cover.filter(pl.col("clock") == "epiTOC2")["passed"][0]
