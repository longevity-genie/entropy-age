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
