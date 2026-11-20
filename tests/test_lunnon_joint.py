import numpy as np

from entropy_age.data.lunnon_joint import joint_betas, quantile_normalise


def test_quantile_normalise_gives_every_array_one_distribution_and_keeps_ranks():
    rng = np.random.default_rng(0)
    x = rng.gamma(2.0, 1000.0, (500, 4)) * np.array([1.0, 2.0, 0.5, 3.0])  # arrays differ in scale only
    q = quantile_normalise(x)
    assert np.allclose(np.sort(q, axis=0), np.sort(q, axis=0)[:, [0]])
    for j in range(4):
        assert (np.argsort(q[:, j]) == np.argsort(x[:, j])).all()


def test_quantile_normalise_keeps_missing_values_missing():
    x = np.arange(40.0).reshape(10, 4)
    x[3, 1] = np.nan
    q = quantile_normalise(x)
    assert np.isnan(q[3, 1]) and np.isfinite(q).sum() == 39


def test_joint_betas_removes_an_array_wide_scale_and_masks_failed_detection():
    rng = np.random.default_rng(1)
    truth = rng.random(400)
    total = rng.gamma(4.0, 2000.0, 400)
    m, u = truth * total, (1 - truth) * total
    scale = np.array([1.0, 3.0])  # the second array is three times brighter in the unmethylated channel only
    mm, uu = np.column_stack([m, m]), np.column_stack([u, u * scale[1]])
    design = np.array(["I"] * 200 + ["II"] * 200)
    detection = np.zeros((400, 2))
    detection[5, 0] = 0.5
    beta = joint_betas(mm, uu, design, detection)
    assert np.isnan(beta[5, 0])
    ok = np.isfinite(beta).all(axis=1)
    assert np.abs(beta[ok, 0] - beta[ok, 1]).max() < 1e-9  # the channel imbalance is gone
    raw = mm / (mm + uu + 100)
    assert np.abs(raw[:, 0] - raw[:, 1]).max() > 0.2
