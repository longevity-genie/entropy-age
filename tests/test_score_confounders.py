"""RM13 (5): what the per-person score picks up."""

import numpy as np

from entropy_age.experiments.score_confounders import partial_r2


def test_partial_r2_counts_only_the_new_block():
    rng = np.random.default_rng(4)
    n = 2000
    age = rng.uniform(20, 80, n)
    batch = rng.integers(0, 2, n).astype(float)
    y = 0.05 * age + 1.0 * batch + rng.normal(size=n)
    base = np.column_stack([np.ones(n), age])
    r2 = partial_r2(y, base, batch[:, None])
    assert 0.15 < r2 < 0.25
    assert partial_r2(y, base, rng.normal(size=(n, 1))) < 0.01
