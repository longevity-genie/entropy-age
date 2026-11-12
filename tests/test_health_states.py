"""RM24: health-state effects on the shield scores, and leave-one-out scoring."""

import numpy as np
import pytest

from entropy_age.experiments.health_states import Contrast, _enough, loo_scores, outcome_effects


def test_effect_is_the_outcome_coefficient_net_of_age_and_standardised():
    rng = np.random.default_rng(1)
    n = 80
    case = np.arange(n) < 30
    age = rng.uniform(40, 90, n) + 10 * case
    f = 0.1 * case + 0.005 * age + rng.normal(0, 0.02, n)
    rows = outcome_effects({"F_housekeeping": f, "P_housekeeping": f, "B_housekeeping": age}, case, age[:, None])
    by = {r["measure"]: r for r in rows}
    assert set(by) == {"F", "B"}
    assert by["F"]["effect"] == pytest.approx(0.1, abs=0.02)
    assert by["F"]["effect_sd"] == pytest.approx(0.1 / f.std(ddof=1), rel=0.2)
    assert by["B"]["p"] > 0.01  # B is age itself: nothing left for the outcome


def test_loo_scores_leave_each_person_out_of_their_reference():
    rng = np.random.default_rng(2)
    x = np.clip(rng.uniform(0.05, 0.95, 400)[None, :] + rng.normal(0, 0.03, (12, 400)), 0, 1)
    in_set = np.zeros(400, dtype=bool)
    in_set[:40] = True
    sets = {"housekeeping": (in_set, [np.arange(i, i + 4) for i in range(0, 40, 4)])}
    island = np.zeros(400, dtype=bool)
    outlier = x.copy()
    outlier[0, :40] = np.clip(outlier[0, :40] + 0.3, 0, 1)  # person 0's set CpGs drift far
    plain, shifted = loo_scores(x, island, sets), loo_scores(outlier, island, sets)
    assert set(plain) >= {"P_housekeeping", "B_housekeeping", "F_housekeeping", "S_housekeeping"}
    assert shifted["S_housekeeping"][0] > plain["S_housekeeping"][0]
    assert shifted["P_housekeeping"][0] > shifted["P_housekeeping"][1:].max()


def test_small_groups_are_left_out():
    x = np.zeros((20, 3))
    few = Contrast("c", "t", "o", x, np.arange(20) < 3, x[:, :1], ["a", "b", "c"], True)
    many = Contrast("c", "t", "o", x, np.arange(20) < 10, x[:, :1], ["a", "b", "c"], True)
    stages = Contrast("c", "t", "o", x, np.arange(20) % 4 * 1.0, x[:, :1], ["a", "b", "c"], False)
    assert not _enough(few) and _enough(many) and _enough(stages)
