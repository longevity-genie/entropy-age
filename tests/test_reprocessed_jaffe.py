"""RM18 (3): the Jaffe controls used for the reprocessing check."""

import polars as pl

from entropy_age.experiments.reprocessed_jaffe import adult_controls


def test_adult_controls_keep_one_sample_per_donor():
    s = pl.DataFrame(
        {
            "sample": ["a", "b", "c", "d"],
            "donor": ["x", "x", "y", "z"],
            "age": [30.0, 30.0, 10.0, 50.0],
            "group": ["Control", "Control", "Control", "Schizo"],
            "drop_sample": [False, False, False, False],
            "best_qc": [False, True, True, True],
            "neuron_fraction": [0.3, 0.3, 0.3, 0.3],
        }
    )
    assert adult_controls(s)["sample"].to_list() == ["b"]
