"""RM31: lineage-specific age slopes recovered from mixed samples."""

import numpy as np

from entropy_age.experiments.lineage_shield import lineage_slopes


def test_slopes_are_recovered_per_lineage():
    rng = np.random.default_rng(22)
    n = 800
    age = rng.uniform(20, 90, n)
    f = rng.dirichlet([6, 3, 1], n)  # three lineages per person
    true = np.array([[0.002, 0.0], [0.0, -0.003], [0.001, 0.001]])  # lineage x CpG
    base = np.array([[0.3, 0.6], [0.5, 0.4], [0.7, 0.5]])
    x = sum(f[:, [k]] * (base[k] + true[k] * (age - age.mean())[:, None]) for k in range(3))
    x = x + rng.normal(0, 0.005, x.shape)
    s = lineage_slopes(x, age, f)
    assert s.shape == (3, 2)
    assert np.allclose(s[:2], true[:2], atol=4e-4)
