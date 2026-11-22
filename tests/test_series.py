import polars as pl

from entropy_age.experiments.series import groups_of, sample_table, specs


def test_every_description_file_parses_and_names_its_columns():
    for name, spec in specs().items():
        assert {"accession", "age", "groups"} <= set(spec), name
        assert all(isinstance(filters, dict) for filters in spec["groups"].values()), name


def test_groups_keep_one_passing_array_per_donor():
    chars = pl.DataFrame(
        {
            "sample": ["a", "b", "c", "d", "e"],
            "age": ["30", "40", "50", "x", "60"],
            "Sex": ["F", "F", "M", "M", "M"],
            "dx": ["Control", "Control", "Case", "Control", "Control"],
        }
    )
    spec = {
        "age": "age",
        "sex": "Sex",
        "donor": "title:^(\\w+)_",
        "hold": ["dx"],
        "groups": {"controls": {"dx": "Control"}, "all": {}},
    }
    table = sample_table(chars, ["D1_x", "D1_y", "D2_x", "D3_x", "D4_x"], spec)
    assert table["donor"].to_list() == ["D1", "D1", "D2", "D4"]  # d has no age
    groups = groups_of(table, spec, usable=["b", "c", "e"], exclude={"D4"})
    assert groups["controls"]["sample"].to_list() == ["b"]  # a failed detection, D4 is excluded
    assert groups["all"]["sample"].to_list() == ["b", "c"]
