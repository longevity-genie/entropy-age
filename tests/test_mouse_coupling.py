import numpy as np
import polars as pl

from entropy_age.experiments.mouse_coupling import pair_r, shuffled


def test_pair_r_over_shared_mice_and_shuffle_keeps_groups():
    rng = np.random.default_rng(0)
    v = rng.normal(size=12)
    a = pl.DataFrame({"donor": [str(i) for i in range(12)], "group": ["Y"] * 6 + ["A"] * 6, "B_x": v})
    b = a.with_columns(pl.col("B_x") * 2.0)
    r, n = pair_r(a, b, "B_x")
    assert n == 12 and abs(r - 1.0) < 1e-12
    s = shuffled(b, rng)
    assert sorted(s.filter(pl.col("group") == "Y")["donor"].to_list()) == [str(i) for i in range(6)]
