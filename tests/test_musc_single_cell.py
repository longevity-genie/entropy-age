import numpy as np

from entropy_age.experiments.musc_single_cell import between_cell_variance, readings


def test_between_cell_variance_removes_read_sampling():
    rng = np.random.default_rng(0)
    cells, promoters = 400, 30
    n = rng.integers(4, 30, (cells, promoters)).astype(float)
    same = rng.binomial(n.astype(int), 0.3).astype(float)  # every cell at 30 %: no variance between cells
    levels = np.clip(rng.normal(0.3, 0.1, (cells, promoters)), 0.01, 0.99)  # cells differ, SD 0.1
    varied = rng.binomial(n.astype(int), levels).astype(float)
    _, tau_same, used = between_cell_variance(same, n)
    _, tau_varied, _ = between_cell_variance(varied, n)
    assert (used == cells).all()
    assert abs(tau_same.mean()) < 0.001
    assert abs(tau_varied.mean() - 0.01) < 0.002


def test_readings_report_divergence_gained_in_old_cells():
    rng = np.random.default_rng(1)
    cells, promoters = 300, 20
    old = np.arange(cells) >= cells // 2
    n = rng.integers(5, 25, (cells, promoters)).astype(float)
    sd = np.where(old[:, None], 0.12, 0.04)
    levels = np.clip(rng.normal(0.4, sd, (cells, promoters)), 0.01, 0.99)
    k = rng.binomial(n.astype(int), levels).astype(float)
    out = readings(k, n, old)
    assert np.nanmean(out["divergence"]) > 0.03 and np.nanmean(out["abs_change"]) < 0.1
