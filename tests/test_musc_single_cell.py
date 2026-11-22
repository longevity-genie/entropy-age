import numpy as np

from entropy_age.experiments.musc_single_cell import between_cell_variance, readings


def test_between_cell_variance_uses_cells_with_enough_reads():
    k = np.array([[2.0, 0.0], [4.0, 1.0], [0.0, 3.0], [8.0, 0.0]])
    n = np.array([[4.0, 2.0], [8.0, 4.0], [4.0, 4.0], [8.0, 1.0]])
    p, variance, cells, reads = between_cell_variance(k, n)
    assert cells.tolist() == [4, 2] and reads.tolist() == [6.0, 4.0]
    assert np.isclose(p[0], 14 / 24) and np.isclose(variance[1], np.var([0.25, 0.75], ddof=1))


def test_readings_report_divergence_gained_in_old_cells():
    rng = np.random.default_rng(1)
    cells, promoters = 300, 20
    old = np.arange(cells) >= cells // 2
    n = rng.integers(5, 25, (cells, promoters)).astype(float)
    sd = np.where(old[:, None], 0.12, 0.04)
    levels = np.clip(rng.normal(0.4, sd, (cells, promoters)), 0.01, 0.99)
    k = rng.binomial(n.astype(int), levels).astype(float)
    out = readings(k, n, old)
    assert np.nanmean(out["divergence"]) > 0.008 and np.nanmean(out["abs_change"]) < 0.1
