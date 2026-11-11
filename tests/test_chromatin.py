"""RM7: chromatin states of 450k CpGs, and the gap left once they are matched on."""

import gzip

import numpy as np
import polars as pl
import pytest

from entropy_age.data.chromatin import cpg_states, read_segmentation
from entropy_age.experiments.chromatin_match import matched_gaps, strata_variants


def test_a_cpg_takes_the_state_of_the_segment_holding_it():
    segments = pl.DataFrame(
        {"chr": ["1", "1", "2"], "start": [0, 100, 0], "end": [100, 200, 50], "state": ["1_TssA", "15_Quies", "7_Enh"]}
    )
    # 1-based positions: 100 is the last base of the first segment, 101 the first of the second
    probes = pl.DataFrame({"chr": ["1", "1", "1", "2", "X"], "pos": [100, 101, 250, 10, 5]})
    assert cpg_states(probes, segments).to_list() == ["1_TssA", "15_Quies", None, "7_Enh", None]


def test_a_segmentation_outside_the_model_is_refused(tmp_path):
    path = tmp_path / "E000_15_coreMarks_mnemonics.bed.gz"
    with gzip.open(path, "wt") as handle:
        handle.write("chr1\t0\t100\t1_TssA\nchr1\t100\t200\t16_Unknown\n")
    with pytest.raises(ValueError, match="16_Unknown"):
        read_segmentation(path)


def test_matching_on_the_state_that_explains_the_gap_removes_it():
    rng = np.random.default_rng(0)
    n = 4000
    active = rng.random(n) < 0.3
    values = np.where(active, 0.2, 1.0) + rng.normal(0, 0.1, n)  # active chromatin drifts less
    in_set = active & (rng.random(n) < 0.5)  # the set is a sample of active CpGs, nothing more
    genes = [np.flatnonzero(in_set)[i : i + 2] for i in range(0, int(in_set.sum()), 2)]
    states = {
        "E062": np.where(active, "1_TssA", "15_Quies"),
        "E073": rng.choice(["1_TssA", "15_Quies"], n),
        "E030": np.where(active, "1_TssA", "15_Quies"),
    }
    variants = strata_variants(np.zeros(n, dtype=bool), rng.random(n), states, "blood (Hannum)")
    assert set(variants) == {"island x level", "+ own state (E062)", "+ other state (E073)", "+ own state (E030)"}
    rows = {r["matching"]: r for r in matched_gaps(values, {"housekeeping": (in_set, genes)}, variants, 50, 1)}
    assert rows["island x level"]["gap_sd"] < -1.0
    assert abs(rows["+ own state (E062)"]["gap_sd"]) < 0.05
    assert rows["+ own state (E062)"]["share_removed"] > 0.95
    assert rows["+ other state (E073)"]["share_removed"] < 0.2
