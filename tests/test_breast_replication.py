import numpy as np

from entropy_age.experiments.breast_replication import standardised_within_chip


def test_standardised_within_chip_removes_chip_offsets():
    v = np.array([1.0, 2.0, 11.0, 12.0])
    z = standardised_within_chip(v, np.array(["a", "a", "b", "b"]))
    assert np.allclose(z, [-1, 1, -1, 1])
