import numpy as np
import polars as pl

from entropy_age.experiments.mesa_expression import (
    NO_PROBE,
    coupling,
    cpg_expression,
    expression_bins,
    gene_level,
    gradient,
    stratum_residual,
)


def test_gene_level_keeps_the_best_expressed_probe():
    values = np.array([[1.0, 1.0], [5.0, 7.0], [2.0, 2.0], [9.0, 9.0]])
    genes, expr = gene_level(["p1", "p2", "p3", "p4"], values, {"p1": "A", "p2": "A", "p3": "B"})
    assert genes == ["A", "B"]
    assert expr.tolist() == [[5.0, 7.0], [2.0, 2.0]]


def test_cpg_expression_takes_promoter_links_and_the_best_gene():
    links = pl.DataFrame(
        {
            "cpg": ["c1", "c1", "c2", "c3"],
            "symbol": ["A", "B", "C", "A"],
            "promoter": [True, True, True, False],
        }
    )
    level, gene = cpg_expression(["c1", "c2", "c3"], links, {"A": 3.0, "B": 8.0})
    assert gene == ["B", None, None]
    assert level[0] == 8.0 and np.isnan(level[1]) and np.isnan(level[2])


def test_expression_bins_marks_missing():
    level = np.array([np.nan, *np.arange(20.0)])
    bins = expression_bins(level, 10)
    assert bins[0] == NO_PROBE
    assert set(bins[1:]) == {str(i) for i in range(1, 11)}
    assert bins[1] == "1" and bins[-1] == "10"


def test_stratum_residual_centres_each_stratum():
    r = stratum_residual(np.array([1.0, 3.0, 10.0, 20.0]), np.array(["a", "a", "b", "b"]))
    assert r.tolist() == [-1.0, 1.0, -5.0, 5.0]


def test_gradient_finds_less_drift_at_higher_expression():
    rng = np.random.default_rng(0)
    n = 2000
    level = rng.random(n)
    bins = expression_bins(level, 10)
    values = 1.0 - level + rng.normal(0, 0.1, n)  # expressed genes drift less
    genes = [f"g{i}" for i in range(n)]
    rows = gradient(values, np.zeros(n, dtype=str), bins, genes, np.ones(n, bool))
    assert len(rows) == 10
    assert rows[0]["relative_slope_sd"] > 0 > rows[-1]["relative_slope_sd"]
    assert rows[0]["spearman"] < -0.8


def test_coupling_recovers_a_negative_link_and_a_null_near_zero():
    rng = np.random.default_rng(1)
    people, genes = 200, 50
    expr = rng.normal(size=(genes, people))
    meth = -0.5 * expr + rng.normal(size=(genes, people))
    design = np.column_stack([np.ones(people), rng.normal(size=people)])
    row = coupling(meth, expr, design, 50, rng)
    assert row["mean_r"] < -0.3 and row["share_negative"] > 0.95
    assert abs(row["null_mean"]) < 0.05 and row["p_permutation"] < 0.05
