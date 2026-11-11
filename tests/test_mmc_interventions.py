"""RM12: scoring treated and control animals alike, against the controls."""

import numpy as np

from entropy_age.experiments.consortium_interventions import arm_scores, chip_design, same_chip


def test_controls_are_scored_leave_one_out_and_a_weaker_shield_shows_in_treated():
    rng = np.random.default_rng(0)
    n, p = 24, 400
    in_set = np.arange(p) < 100
    treated = np.arange(n) >= 12
    level = rng.uniform(0.2, 0.8, p)
    noise = np.where(in_set, 0.01, 0.03)
    x = level + rng.normal(0, 1, (n, p)) * noise
    x[np.ix_(treated, in_set)] += rng.normal(
        0, 0.03, (int(treated.sum()), int(in_set.sum()))
    )  # treated set CpGs drift as much as the rest
    genes = [np.arange(j, j + 4) for j in range(0, 100, 4)]
    scores = arm_scores(np.clip(x, 0, 1), treated, np.zeros(p, dtype=bool), {"housekeeping": (in_set, genes)})
    f = scores["F_housekeeping"]
    assert f.shape == (n,)
    assert f[~treated].mean() < -0.3  # controls: set CpGs drift far less
    assert f[treated].mean() > f[~treated].mean() + 0.2  # treated: the shield is weaker


def test_with_no_treatment_controls_and_treated_drift_alike():
    # a control left in its own reference would sit closer to it and look less drifted
    rng = np.random.default_rng(1)
    n, p = 24, 2000
    in_set = np.arange(p) < 400
    treated = np.arange(n) >= 12
    x = np.clip(rng.uniform(0.2, 0.8, p) + rng.normal(0, 0.03, (n, p)), 0, 1)
    genes = [np.arange(j, j + 4) for j in range(0, 400, 4)]
    b = arm_scores(x, treated, np.zeros(p, dtype=bool), {"housekeeping": (in_set, genes)})["B_housekeeping"]
    assert abs(b[~treated].mean() / b[treated].mean() - 1) < 0.03


def test_same_chip_keeps_treated_and_only_controls_sharing_their_chips():
    chips = np.array(["a", "a", "b", "b", "c", "c", "a"])
    treated = np.array([True, False, True, False, False, False, False])
    is_control = np.array([False, True, False, True, True, True, False])  # the last is another arm
    assert same_chip(chips, treated, is_control).tolist() == [True, True, True, True, False, False, False]


def test_chip_design_has_one_column_per_extra_chip_and_none_for_one_chip():
    d = chip_design(np.array(["b", "a", "b", "c"]))
    assert d.tolist() == [[1, 0], [0, 0], [1, 0], [0, 1]]
    assert chip_design(np.array(["a", "a"])).shape == (2, 0)


def test_signed_shift_says_which_way_set_cpgs_moved():
    # set CpGs pushed up in treated animals: S rises; pushed down: S falls; P rises either way
    rng = np.random.default_rng(2)
    n, p = 24, 400
    in_set = np.arange(p) < 100
    treated = np.arange(n) >= 12
    base = rng.uniform(0.2, 0.8, p) + rng.normal(0, 0.02, (n, p))
    genes = [np.arange(j, j + 4) for j in range(0, 100, 4)]
    sets = {"housekeeping": (in_set, genes)}
    out = {}
    for sign in (1, -1):
        x = base.copy()
        x[np.ix_(treated, in_set)] += sign * 0.05
        out[sign] = arm_scores(np.clip(x, 0, 1), treated, np.zeros(p, dtype=bool), sets)
    for sign in (1, -1):
        s = out[sign]["S_housekeeping"]
        assert sign * (s[treated].mean() - s[~treated].mean()) > 0.1
        pp = out[sign]["P_housekeeping"]
        assert pp[treated].mean() > pp[~treated].mean() + 0.05


def test_composition_components_remove_a_shift_shared_by_every_cpg():
    # treated animals differ only along one composition axis that moves all CpGs; removing it leaves no effect
    rng = np.random.default_rng(3)
    n, p = 40, 1200
    in_set = np.arange(p) < 200
    treated = np.arange(n) >= 20
    loading = rng.normal(0, 1, p)
    axis = rng.normal(0, 1, n) + 3.0 * treated
    x = np.clip(rng.uniform(0.3, 0.7, p) + 0.02 * np.outer(axis, loading) + rng.normal(0, 0.005, (n, p)), 0, 1)
    genes = [np.arange(j, j + 4) for j in range(0, 200, 4)]
    sets = {"housekeeping": (in_set, genes)}
    raw = arm_scores(x, treated, np.zeros(p, dtype=bool), sets)["B_housekeeping"]
    adj = arm_scores(x, treated, np.zeros(p, dtype=bool), sets, components=1)["B_housekeeping"]
    assert raw[treated].mean() > 2 * raw[~treated].mean()  # composition looks like drift
    assert adj[treated].mean() < 1.3 * adj[~treated].mean()  # removed
