"""RM30: per-gene coupling of methylation and expression residuals, signed and in spread."""

import numpy as np
import pytest

from entropy_age.experiments.meth_expr_spread import columnwise_r, slope_coupling, spread_coupling


def test_columnwise_r_matches_numpy_per_column():
    rng = np.random.default_rng(7)
    a, b = rng.normal(size=(50, 3)), rng.normal(size=(50, 3))
    b[:, 0] = a[:, 0] * 2 + 1
    r = columnwise_r(a, b)
    assert r[0] == pytest.approx(1.0)
    assert r[1] == pytest.approx(np.corrcoef(a[:, 1], b[:, 1])[0, 1])


def test_spread_coupling_sees_shared_noise_and_its_null_does_not():
    rng = np.random.default_rng(8)
    donors, genes = 120, 300
    noisy = rng.gamma(1.0, 1.0, size=(donors, genes))  # how unusual each donor is at each gene
    meth = noisy * rng.choice([-1, 1], size=(donors, genes)) + rng.normal(0, 0.3, (donors, genes))
    expr = -0.5 * meth + noisy * rng.choice([-1, 1], size=(donors, genes))  # silencing plus shared spread
    c = spread_coupling(meth, expr, rng, permutations=50)
    assert c["signed_mean_r"] < -0.1 and c["signed_share_negative"] > 0.8
    assert c["abs_mean_r"] > c["abs_null_mean"] + 5 * c["abs_null_sd"] and c["abs_p"] < 0.05
    unrelated = spread_coupling(rng.normal(size=(donors, genes)), rng.normal(size=(donors, genes)), rng, 50)
    assert abs(unrelated["abs_mean_r"] - unrelated["abs_null_mean"]) < 4 * unrelated["abs_null_sd"]


def test_slope_coupling_holds_the_level_fixed():
    rng = np.random.default_rng(9)
    n = 2000
    level_m, level_e = rng.uniform(size=n), rng.uniform(size=n)
    shared_level = rng.uniform(size=n)
    # both slopes follow a shared level only: no coupling once each slope's own level is held
    ms, es = shared_level + rng.normal(0, 0.1, n), shared_level + rng.normal(0, 0.1, n)
    r_raw, _ = slope_coupling(ms, level_m, es, level_e)
    r_held, _ = slope_coupling(ms, shared_level, es, shared_level)
    assert r_raw > 0.5 and abs(r_held) < 0.1
