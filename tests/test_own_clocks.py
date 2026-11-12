"""RM11: equal-size pools, the look-alike pool matched stratum by stratum."""

import numpy as np

from entropy_age.experiments.own_clocks import pools


def test_pools_have_equal_size_and_lookalikes_match_the_strata():
    features = [f"c{i}" for i in range(60)]
    members = {"c0", "c1", "c2", "c30"}
    strata = np.array(["a"] * 30 + ["b"] * 30)
    p = pools(features, members, strata, np.random.default_rng(0))
    assert p["shielded"].tolist() == [0, 1, 2, 30]
    assert p["look-alike"].size == p["random"].size == 4
    assert sorted(strata[p["look-alike"]].tolist()) == ["a", "a", "a", "b"]
    assert not set(p["look-alike"]) & set(p["shielded"]) and not set(p["random"]) & set(p["shielded"])


def test_pools_trim_both_sides_where_background_is_scarce():
    features = [f"c{i}" for i in range(5)]
    strata = np.array(["a", "a", "a", "b", "a"])  # stratum a: set c0, c1, c2; background c4 only
    p = pools(features, {"c0", "c1", "c2"}, strata, np.random.default_rng(0))
    assert p["shielded"].size == p["look-alike"].size == p["random"].size == 1
    assert p["look-alike"].tolist() == [4]
