"""RM8: candidate pools by chromatin, and the same-donor brain-minus-blood offset of a clock."""

import numpy as np
import pytest

from entropy_age.experiments.chromatin_clocks import lunnon_offsets, pools


def test_pools_need_the_state_in_both_tissues():
    blood = np.array(["1_TssA", "1_TssA", "13_ReprPC", "13_ReprPC", "15_Quies"])
    brain = np.array(["2_TssAFlnk", "13_ReprPC", "10_TssBiv", "1_TssA", "15_Quies"])
    p = pools([f"cg{j}" for j in range(5)], blood, brain)
    assert p["active in both"].tolist() == [0]
    assert p["Polycomb in both"].tolist() == [2]
    assert p["any"].tolist() == [0, 1, 2, 3, 4]


def test_offsets_compare_each_donor_with_their_own_blood():
    regions = ("frontal cortex", "superior temporal gyrus", "entorhinal cortex", "cerebellum")
    predictions = {"whole blood": np.array([50.0, 60.0])} | {r: np.array([72.0, 55.0]) for r in regions}
    ages = {"whole blood": np.array([50.0, 70.0])} | {r: np.array([70.0, 50.0]) for r in regions}
    donors = {"whole blood": np.array(["a", "b"])} | {r: np.array(["b", "a"]) for r in regions}
    row = lunnon_offsets(predictions, ages, donors)[0]
    # donor b: region +2, blood -10 -> +12; donor a: region +5, blood 0 -> +5
    assert row["donors"] == 2 and row["brain_minus_blood"] == pytest.approx(8.5)
