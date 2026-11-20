import numpy as np
import pytest

from entropy_age.experiments import common
from entropy_age.experiments.chromatin_match import KINDS, strata_variants


def test_blood_state_default_and_switch(monkeypatch):
    monkeypatch.delenv("ENTROPY_AGE_BLOOD_STATE", raising=False)
    monkeypatch.delenv("ENTROPY_AGE_BETAS", raising=False)
    assert common.blood_state() == "E062"
    assert common.own_epigenome("blood (Hannum)", "E062") == "E062"
    assert common.source_suffix() == ""
    monkeypatch.setenv("ENTROPY_AGE_BLOOD_STATE", "E030")
    assert common.own_epigenome("blood (Hannum)", "E062") == "E030"
    assert common.own_epigenome("cortex (Jaffe)", "E073") == "E073"
    assert common.source_suffix() == "_e030"


def test_blood_state_refuses_unknown(monkeypatch):
    monkeypatch.setenv("ENTROPY_AGE_BLOOD_STATE", "E999")
    with pytest.raises(ValueError):
        common.blood_state()


def test_joint_blood_matching_is_finer_than_either_map():
    island = np.zeros(6, bool)
    baseline = np.linspace(0.1, 0.2, 6)
    states = {
        "E062": np.array(["a", "a", "a", "b", "b", "b"]),
        "E030": np.array(["a", "b", "a", "b", "a", "b"]),
        "E073": np.array(["a"] * 6),
    }
    v = strata_variants(island, baseline, states, "blood (Hannum)")
    joint = v["+ both blood states (E062, E030)"]
    for other in ("+ own state (E062)", "+ own state (E030)"):
        for s in np.unique(joint):
            assert np.unique(v[other][joint == s]).size == 1
    # the figure colours a matching by its kind: every matching must have one
    assert {name.split(" (")[0] for name in v} <= set(KINDS)


def test_conversion_switch_and_group_table(monkeypatch):
    from entropy_age.data.reprocessed import conversion_by_group

    monkeypatch.delenv("ENTROPY_AGE_BLOOD_STATE", raising=False)
    monkeypatch.delenv("ENTROPY_AGE_BETAS", raising=False)
    monkeypatch.setenv("ENTROPY_AGE_GCT", "1")
    assert common.conversion_terms() and common.source_suffix() == "_gct"
    gct = np.array([1.0, 1.1, 1.2, 1.5, 1.6, 1.7])
    row = conversion_by_group(gct, np.array([False, False, False, True, True, True]))
    assert abs(row["gct_carriers"] - 1.6) < 1e-12 and abs(row["gct_controls"] - 1.1) < 1e-12
    assert abs(row["difference_sd"] - 5.0) < 1e-9 and row["p_welch"] < 0.01


def test_human_annotation_switch(monkeypatch):
    for k in ("ENTROPY_AGE_BLOOD_STATE", "ENTROPY_AGE_BETAS", "ENTROPY_AGE_GCT", "ENTROPY_AGE_LUNNON"):
        monkeypatch.delenv(k, raising=False)
    assert not common.human_annotation()
    monkeypatch.setenv("ENTROPY_AGE_SPECIES_ANNOTATION", "human")
    assert common.human_annotation() and common.source_suffix() == "_humanann"
