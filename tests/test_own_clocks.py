"""RM11: equal-size pools, the look-alike pool matched stratum by stratum."""

import numpy as np
import pytest

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


def test_pools_refuse_a_stratum_without_enough_background():
    with pytest.raises(ValueError):
        pools(["c0", "c1", "c2"], {"c0", "c1"}, np.array(["a", "a", "b"]), np.random.default_rng(0))
