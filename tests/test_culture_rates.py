import numpy as np
import polars as pl

from entropy_age.experiments.culture_rates import joint_rates


def test_joint_rates_separate_divisions_from_days():
    rng = np.random.default_rng(0)
    rows = []
    for e in range(4):
        for i in range(10):
            pd_, days = float(i), float(i * (1 + e))  # arms differ in days per doubling
            rows.append(
                {
                    "subexperiment": f"x{e}",
                    "donor": "d",
                    "doublings": pd_,
                    "days": days,
                    "B_x": 0.1 * e + 0.02 * pd_ + 0.001 * days + rng.normal(0, 1e-4),
                }
            )
    r = joint_rates(pl.DataFrame(rows), "B_x")
    assert abs(r["per_10_doublings"] - 0.2) < 0.01 and abs(r["per_10_days"] - 0.01) < 0.005
