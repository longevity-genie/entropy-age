"""RM8: published clocks' CpGs, their drift against matched CpGs, and their tissue offsets."""

import numpy as np
import pytest

from entropy_age.data.clocks import git_blob_sha1, read_clock
from entropy_age.experiments.clock_cpgs import chromatin_shares, clock_gap, tissue_offsets


def test_a_clock_file_loses_its_intercept_and_byte_order_mark(tmp_path):
    path = tmp_path / "clock.csv"
    path.write_bytes("﻿CpGmarker,CoefficientTraining\nintercept,30.5\ncg0001,0.5\nch.1.2,1\ncg0002,-0.25\n".encode())
    frame = read_clock(path)
    assert frame["cpg"].to_list() == ["cg0001", "cg0002"]
    assert frame["weight"].to_list() == [0.5, -0.25]


def test_the_blob_hash_is_gits():
    assert git_blob_sha1(b"") == "e69de29bb2d1d6434b8b29ae775ad8c2e48c5391"


def test_clock_cpgs_that_drift_faster_give_a_positive_gap():
    rng = np.random.default_rng(0)
    n = 3000
    values = rng.normal(0, 1, n)
    clock = np.zeros(n, dtype=bool)
    clock[:200] = True
    values[clock] += 0.5
    other = np.zeros(n, dtype=bool)
    other[200:300] = True  # another clock's CpGs: in neither side
    values[other] += 10
    strata = rng.integers(0, 5, n).astype(str)
    gap = clock_gap(values, clock, ~clock & ~other, strata, 100, 1)
    sd = values[~other].std()
    assert gap["gap_sd"] == pytest.approx(0.5 / sd, abs=0.15)
    assert gap["ci_low"] < gap["gap_sd"] < gap["ci_high"] and gap["ci_low"] > 0


def test_chromatin_shares_count_active_and_polycomb_states():
    states = np.array(["1_TssA", "2_TssAFlnk", "13_ReprPC", "15_Quies"])
    assert chromatin_shares(states, np.ones(4, dtype=bool)) == {"active": 0.5, "polycomb": 0.25}


def test_tissue_offsets_pair_donors_and_scale_by_blood():
    blood = np.array([0.0, 2.0, 4.0])  # SD 2
    scores = {
        "whole blood": blood,
        **{
            r: np.array([3.0, 4.0])
            for r in ("frontal cortex", "superior temporal gyrus", "entorhinal cortex", "cerebellum")
        },
    }
    donors = {
        "whole blood": np.array(["a", "b", "c"]),
        **{r: np.array(["c", "a"]) for r in scores if r != "whole blood"},
    }
    row = tissue_offsets(scores, donors)[0]
    assert row["donors"] == 2
    assert row["offset_sd"] == pytest.approx(((3 - 4) + (4 - 0)) / 2 / 2)
