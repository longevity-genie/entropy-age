import pytest

from entropy_age.experiments.muscle_training import SERIES, design


def test_design_reads_condition_and_person_from_titles():
    gs = design(["a", "b"], ["Muscle_PRE_GS2_Rep1", "Muscle_4WP_GS2_Rep1"], SERIES["GSE171140"]["pattern"])
    assert gs["condition"].to_list() == ["PRE", "4WP"] and gs["person"].to_list() == ["GS2", "GS2"]
    sb = design(
        ["c", "d"],
        ["SkM_Epi_Mem_1: Muscle_Baseline_Participant_001_Rep1", "SkM_Epi_Mem_27: Muscle_Acute_P006_Rep1"],
        SERIES["GSE114763"]["pattern"],
    )
    assert sb["condition"].to_list() == ["Baseline", "Acute"] and sb["person"].to_list() == ["1", "6"]


def test_design_refuses_a_title_it_cannot_read():
    with pytest.raises(ValueError, match="does not match"):
        design(["a"], ["something else"], SERIES["GSE171140"]["pattern"])
