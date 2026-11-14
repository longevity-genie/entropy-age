"""RM34: placebo gene sets."""

import numpy as np
import polars as pl

from entropy_age.experiments.placebo_sets import placebo_membership, summarise


def test_placebo_draws_genes_outside_the_real_sets():
    links = pl.DataFrame(
        {
            "cpg": [f"c{i}" for i in range(12)],
            "symbol": ["A", "A", "B", "C", "C", "D", "E", "F", "G", "G", "H", "H"],
            "promoter": [True] * 11 + [False],
        }
    )
    real = pl.DataFrame({"feature": ["c0", "c1"], "group": ["housekeeping"] * 2, "gene": ["A", "A"]})
    feats = {f"c{i}" for i in range(12)}
    mem = placebo_membership(links, real, feats, 3, "placebo", np.random.default_rng(0))
    assert not mem["feature"].is_in(["c0", "c1", "c11"]).any()
    genes = links.filter(pl.col("cpg").is_in(mem["feature"].implode()))["symbol"].unique()
    assert genes.len() == 3 and set(mem["group"]) == {"placebo"}


def test_summarise_ranks_real_against_placebos():
    real = pl.DataFrame({"tissue": ["t"], "group": ["g"], "metric": ["m"], "matched_diff_sd": [-0.3]})
    placebo = pl.DataFrame(
        {
            "tissue": ["t"] * 100,
            "set": ["g"] * 100,
            "metric": ["m"] * 100,
            "matched_diff_sd": np.linspace(-0.1, 0.1, 100),
        }
    )
    s = summarise(real, placebo)
    assert s["share_placebo_below_real"][0] == 0 and abs(s["placebo_mean"][0]) < 1e-9
