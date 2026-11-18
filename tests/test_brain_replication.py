import numpy as np
import polars as pl

from entropy_age.experiments.brain_replication import complete_matrix


def test_complete_matrix_keeps_only_cpgs_measured_in_every_array():
    frame = pl.DataFrame({"cpg": ["a", "b", "c"], "s1": [0.1, None, 0.3], "s2": [0.2, 0.5, 0.4]})
    x, feats = complete_matrix(frame, ["a", "b", "c", "d"], ["s1", "s2"])
    assert feats == ["a", "c"]  # b is missing in s1, d is not on the array
    assert np.allclose(x, [[0.1, 0.3], [0.2, 0.4]])
