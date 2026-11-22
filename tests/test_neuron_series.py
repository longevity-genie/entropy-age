import polars as pl

from entropy_age.experiments.neuron_series import ALL, CONTROLS, neuron_groups


def test_neuron_groups_split_controls_and_drop_missing_ages():
    chars = pl.DataFrame(
        {
            "sample": ["a", "b", "c", "d", "e"],
            "tissuebank.id": ["1", "2", "3", "4", "1"],
            "age": ["40", "68", "unknown", "55", "40"],
            "Sex": ["Male", "Female", "Male", "Male", "Male"],
            "dist.dx": ["Bipolar", "Control", "Control", "Schizophrenia", "Bipolar"],
        }
    )
    groups = neuron_groups(chars)
    assert groups[ALL]["sample"].to_list() == ["a", "b", "d"]  # c has no age; e is a second array of donor 1
    assert groups[CONTROLS]["sample"].to_list() == ["b"]
    assert groups[ALL]["age"].to_list() == [40.0, 68.0, 55.0]
