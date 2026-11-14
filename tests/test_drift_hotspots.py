"""RM34 item 4: hotspot enrichment against a level-matched background."""

import numpy as np

from entropy_age.experiments.drift_hotspots import enrichment


def test_enrichment_is_one_at_random_and_high_when_concentrated():
    rng = np.random.default_rng(0)
    n = 20_000
    strata = rng.integers(0, 10, n)
    member = rng.random(n) < 0.1
    hot = rng.random(n) < 0.05
    ratio, z, _, _ = enrichment(hot, member, strata)
    assert abs(ratio - 1) < 0.2 and abs(z) < 4
    hot2 = hot | (member & (rng.random(n) < 0.3))
    ratio2, z2, _, _ = enrichment(hot2, member, strata)
    assert ratio2 > 3 and z2 > 10
