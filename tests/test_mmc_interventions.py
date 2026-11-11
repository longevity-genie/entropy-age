"""RM12: scoring treated and control animals alike, against the controls."""

import numpy as np

from entropy_age.experiments.consortium_interventions import arm_scores


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
