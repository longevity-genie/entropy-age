import numpy as np

from entropy_age.experiments.atac_cortex import excess_dispersion_slopes


def test_excess_dispersion_slope_sees_growing_variation_and_ignores_read_depth():
    rng = np.random.default_rng(0)
    n, m = 400, 60
    age = rng.uniform(20, 90, n)
    total = rng.uniform(0.8e6, 1.2e6, n)
    depth = np.where(np.arange(m) < 30, 2e-5, 4e-4)  # two classes of peaks: few reads and many reads
    grows = (np.arange(m) % 2) == 0  # in half of each class the biological variation grows with age
    cv2 = np.where(grows[None, :], 0.02 + 0.002 * (age[:, None] - 20), 0.05)
    rate = rng.gamma(1 / cv2, cv2) * depth[None, :] * total[:, None]
    x = rng.poisson(rate).astype(float)
    slopes = excess_dispersion_slopes(x, total, age, np.empty((n, 0)))
    assert slopes[grows].mean() > 5 * abs(slopes[~grows].mean())
    low, high = slopes[grows & (np.arange(m) < 30)].mean(), slopes[grows & (np.arange(m) >= 30)].mean()
    assert 0.5 < low / high < 2.0  # the same growth is read alike in shallow and deep peaks
    assert abs(slopes[~grows & (np.arange(m) < 30)].mean()) < 0.0005
