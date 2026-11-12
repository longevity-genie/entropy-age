"""RM18: probe design in the strata, and Zhou et al.'s mask."""

import gzip

import polars as pl
import pytest

from entropy_age.data.annotation import read_probe_mask
from entropy_age.experiments.probe_robustness import keep_unmasked, strata

PROBES = pl.DataFrame(
    {
        "cpg": ["a", "b", "c", "d"],
        "cpg_island": [True, True, False, False],
        "probe_type": ["I", "II", "I", "II"],
    }
)
METRICS = pl.DataFrame({"feature": ["a", "b", "c", "d"], "baseline": [0.05, 0.05, 0.05, 0.05]})


def test_probe_type_splits_otherwise_equal_strata():
    plain = strata(METRICS, PROBES, by_type=False)["stratum"].to_list()
    typed = strata(METRICS, PROBES, by_type=True)["stratum"].to_list()
    assert plain[0] == plain[1] and plain[2] == plain[3]  # island x decile only
    assert typed[0] != typed[1] and typed[2] != typed[3]  # design separates them
    assert len(set(typed)) == 4


def test_strata_refuse_a_missing_probe_type():
    probes = PROBES.with_columns(probe_type=pl.Series(["I", None, "I", "II"]))
    with pytest.raises(ValueError):
        strata(METRICS, probes, by_type=True)


def test_mask_reads_m_general_and_drops_only_masked_probes(tmp_path):
    path = tmp_path / "mask.tsv.gz"
    with gzip.open(path, "wt") as handle:
        handle.write("Probe_ID\tmask\tmaskUniq\tM_general\na\tx\tx\tTRUE\nb\ty\ty\tFALSE\n")
    mask = read_probe_mask(path)
    assert dict(mask.iter_rows()) == {"a": True, "b": False}
    kept = keep_unmasked(METRICS, set(mask.filter(pl.col("masked"))["cpg"]))
    assert kept["feature"].to_list() == ["b", "c", "d"]
