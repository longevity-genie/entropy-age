import polars as pl

from entropy_age.experiments.brain_regions import region_groups


def test_region_groups_split_by_region_and_drop_known_donors():
    chars = pl.DataFrame(
        {
            "sample": ["a", "b", "c", "d", "e"],
            "brain region": ["PFC", "PFC", "AMG", "PFC", "PFC"],
            "age at death": ["30.5", "61", "30.5", "45", "45"],
            "Sex": ["F", "M", "F", "M", "M"],
            "primary dx": ["Control", "MDD", "Control", "Control", "Control"],
        }
    )
    groups = region_groups(chars, ["Br1_PFC", "Br2_PFC", "Br1_AMG", "Br3_PFC", "Br3_PFC"], {"Br2"})
    assert groups["prefrontal cortex, controls"]["sample"].to_list() == ["a", "d"]  # e repeats Br3, Br2 is left out
    assert groups["prefrontal cortex, every donor (diagnosis held)"]["donor"].to_list() == ["Br1", "Br3"]
    assert groups["amygdala, controls"]["sample"].to_list() == ["c"]
