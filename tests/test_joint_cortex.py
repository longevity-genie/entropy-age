import numpy as np
import polars as pl

from entropy_age.experiments.joint_cortex import NO_EXPRESSION, expression_class, joint_variants, usable_share


def test_expression_class_takes_the_highest_promoter_gene():
    links = pl.DataFrame(
        {"cpg": ["a", "a", "b", "c"], "symbol": ["G1", "G2", "G3", "G9"], "promoter": [True, True, False, True]}
    )
    genes = pl.DataFrame({"symbol": ["G1", "G2", "G3"], "decile": [2, 7, 5]})
    out = expression_class(["a", "b", "c", "d"], links, genes)
    assert out.tolist() == ["7", NO_EXPRESSION, NO_EXPRESSION, NO_EXPRESSION]  # b: no promoter link; c: unlisted


def test_joint_variants_nest_and_usable_share_counts_shared_strata():
    island = np.array([True, True, True, True])
    base = np.array([0.05, 0.05, 0.05, 0.05])
    expr = np.array(["1", "1", "2", "2"])
    state = np.array(["x", "x", "x", "y"])
    cells = np.array(["0|0"] * 4)
    v = joint_variants(island, base, expr, state, cells)
    assert np.unique(v["island x level"]).size == 1
    assert np.unique(v["+ expression + state + sorted ATAC"]).size == 3
    in_set = np.array([True, False, True, False])
    assert usable_share(in_set, v["island x level"]) == 1.0
    assert usable_share(in_set, v["+ expression + state + sorted ATAC"]) == 0.5  # the third CpG has no control left
