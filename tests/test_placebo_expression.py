import numpy as np
import polars as pl

from entropy_age.experiments.placebo_expression import matched_draw, with_deciles


def test_deciles_split_evenly():
    g = with_deciles(pl.DataFrame({"symbol": [f"g{i}" for i in range(100)], "log_tpm": np.arange(100.0)}))
    assert g.group_by("decile").len()["len"].to_list() == [10] * 10
    assert g.filter(pl.col("symbol") == "g99")["decile"][0] == 9


def test_draw_keeps_decile_counts_and_reports_shortfall():
    pool = pl.DataFrame({"symbol": ["a", "b", "c", "d"], "decile": [0, 0, 9, 1]})
    real = pl.DataFrame({"symbol": ["x", "y", "z"], "decile": [0, 9, 9]})
    genes, short = matched_draw(real, pool, np.random.default_rng(0))
    assert sorted(genes)[-1] == "c" and len(genes) == 2
    assert short == 1
    assert len(set(genes) & {"a", "b"}) == 1
