"""RM13 (4): the species shield off the floor."""

import numpy as np

from entropy_age.experiments.species_floor import off_floor


def test_off_floor_uses_the_youngest_quarter():
    ages = np.array([1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0])
    x = np.array([[0.02, 0.5, 0.97]] * 2 + [[0.5, 0.5, 0.5]] * 6)
    assert off_floor(x, ages).tolist() == [False, True, False]
