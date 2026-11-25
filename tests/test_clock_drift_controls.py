import numpy as np
import polars as pl

from entropy_age.experiments.clock_drift_controls import (
    corr,
    independent,
    matched_placebo,
    panel_indices,
    partial_r,
    placebo_pool,
)


def test_panel_indices_level_compression_and_spread():
    x = np.array([[0.5, 0.5, 0.5, 0.5], [0.0, 1.0, 0.0, 1.0], [0.1, 0.1, 0.9, 0.9]])
    ind = panel_indices(x)
    assert np.allclose(ind["level"], [0.5, 0.5, 0.5])
    assert np.allclose(ind["compression"], [0.0, 0.5, 0.4])  # mean |beta - 0.5|: lower is more compressed
    assert np.allclose(ind["sd"], [0.0, 0.5, 0.4])


def test_compression_falls_when_betas_shrink_toward_half():
    rng = np.random.default_rng(1)
    x = rng.uniform(0, 1, (1, 5000))
    shrunk = 0.5 + 0.7 * (x - 0.5)
    assert panel_indices(shrunk)["compression"][0] < panel_indices(x)["compression"][0]
    assert np.isclose(panel_indices(shrunk)["level"][0], 0.5 + 0.7 * (panel_indices(x)["level"][0] - 0.5))


def test_partial_r_removes_a_shared_driver():
    rng = np.random.default_rng(2)
    z = rng.normal(size=2000)
    x = z + rng.normal(0, 0.3, z.size)
    y = z + rng.normal(0, 0.3, z.size)
    assert corr(x, y) > 0.8
    assert abs(partial_r(x, y, z[:, None])) < 0.1


def test_partial_r_keeps_a_direct_link():
    rng = np.random.default_rng(3)
    z = rng.normal(size=2000)
    x = rng.normal(size=z.size)
    y = 0.8 * x + z + rng.normal(0, 0.3, z.size)
    assert partial_r(x, y, z[:, None]) > 0.9


def test_matched_placebo_reproduces_stratum_counts():
    rng = np.random.default_rng(4)
    strata = rng.choice(np.array(["a", "b", "c"]), 600)
    real = np.zeros(600, dtype=bool)
    real[:60] = True
    pool = [np.arange(j, min(j + 4, 600)) for j in range(60, 600, 4)]
    mask, genes, fill = matched_placebo(strata, real, pool, np.random.default_rng(0))
    assert fill == 1.0
    assert not (mask & real).any()
    for s in "abc":
        assert (mask & (strata == s)).sum() == (real & (strata == s)).sum()
    assert sum(g.size for g in genes) == mask.sum()


def test_placebo_pool_skips_set_genes_and_members():
    kept = ["c1", "c2", "c3", "c4"]
    links = pl.DataFrame(
        {
            "cpg": ["c1", "c2", "c3", "c4", "c4"],
            "symbol": ["HK", "G1", "G1", "G2", "HK"],
            "promoter": [True, True, True, True, True],
        }
    )
    member = np.array([True, False, True, False])
    pool = placebo_pool(kept, member, links, {"HK"})
    assert [g.tolist() for g in pool] == [[1], [3]]


def test_independent_takes_largest_non_brain_group_per_series():
    t = pl.DataFrame(
        {
            "series": ["s1", "s1", "s2", "s3"],
            "group": ["all", "controls", "cortex", "x"],
            "n": [100, 40, 80, 30],
            "brain": [False, False, True, False],
            "sensitivity": [False, False, False, True],
        }
    )
    assert independent(t)[["series", "group"]].rows() == [("s1", "all")]


def test_loo_drift_matches_loo_scores():
    from entropy_age.experiments.clock_drift_controls import loo_drift
    from entropy_age.experiments.health_states import loo_scores

    rng = np.random.default_rng(5)
    x = np.clip(rng.beta(2, 5, (12, 400)) + rng.normal(0, 0.02, (12, 400)), 0, 1)
    island = rng.random(400) < 0.5
    a = np.zeros(400, dtype=bool)
    a[:40] = True
    b = np.zeros(400, dtype=bool)
    b[40:90] = True
    sets = {"hk": (a, [np.arange(40)]), "other": (b, [np.arange(40, 90)])}
    ours, theirs = loo_drift(x, island, sets), loo_scores(x, island, sets)
    for key in ("P_hk", "B_hk", "F_hk", "F_other"):
        assert np.allclose(ours[key], theirs[key])


def test_loo_drift_drop_leaves_cpgs_out_of_the_contrast():
    from entropy_age.experiments.clock_drift_controls import loo_drift

    rng = np.random.default_rng(6)
    x = rng.beta(2, 5, (10, 300))
    island = np.zeros(300, dtype=bool)
    a = np.zeros(300, dtype=bool)
    a[:30] = True
    hk = np.zeros(300, dtype=bool)
    hk[30:80] = True
    x[:, 30:80] += rng.normal(0, 0.2, (10, 50))  # noisy CpGs: in the look-alikes they raise B
    sets = {"pl": (a, [np.arange(30)]), "hk": (hk, [np.arange(30, 80)])}
    kept = loo_drift(x, island, sets)
    dropped = loo_drift(x, island, sets, {"pl": hk})
    assert not np.allclose(kept["B_pl"], dropped["B_pl"])
    assert np.allclose(kept["F_hk"], dropped["F_hk"])
    assert np.allclose(loo_drift(x, island, sets, {"pl": np.zeros(300, dtype=bool)})["F_pl"], kept["F_pl"])
