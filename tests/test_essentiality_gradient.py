"""RM13 (1): essentiality beyond chromatin and expression."""

import numpy as np
import polars as pl

from entropy_age.experiments.essentiality_gradient import cpg_entrez, demean, depmap_effects, slope


def test_slope_absorbs_strata():
    rng = np.random.default_rng(2)
    strata = np.repeat(np.array(["a", "b", "c"]), 200)
    x = rng.normal(size=600)
    y = 0.5 * x + np.repeat([10.0, -3.0, 4.0], 200) + rng.normal(scale=0.1, size=600)
    assert abs(slope(y, x, np.empty((600, 0)), strata) - 0.5) < 0.02
    assert np.allclose(demean(np.array([1.0, 3.0, 5.0]), np.array(["a", "a", "b"])), [-1.0, 1.0, 0.0])


def test_depmap_columns_give_entrez(tmp_path):
    path = tmp_path / "d.csv"
    path.write_text(",A1BG (1),AARS1 (16)\nACH-1,0.1,-1.0\nACH-2,-0.1,-2.0\n")
    d = depmap_effects(path)
    assert d["entrez"].to_list() == [1, 16] and np.allclose(d["effect"].to_numpy(), [0.0, -1.5])


def test_cpg_entrez_keeps_single_gene_cpgs_and_old_symbols():
    links = pl.DataFrame(
        {"cpg": ["c1", "c2", "c2", "c3"], "symbol": ["AARS", "A1BG", "AARS", "A1BG"], "promoter": [True] * 4}
    )
    hgnc = pl.DataFrame({"symbol": ["AARS1", "A1BG"], "previous": [["AARS"], []], "entrez": [16, 1]})
    got = dict(cpg_entrez(links, hgnc).sort("cpg").iter_rows())
    assert got == {"c1": 16, "c3": 1}
