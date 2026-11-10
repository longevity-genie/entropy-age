import gzip
from pathlib import Path

import numpy as np
import polars as pl
import pytest

from entropy_age.data.cohorts import _fibres, _gtex
from entropy_age.experiments.human_tissues import (
    ADJUSTMENTS,
    FIBRES,
    TURNOVER,
    FibreSet,
    _boot_fibres,
    _boot_gtex,
    _subsample,
    summarise_fibres,
    summarise_gtex,
)
from entropy_age.experiments.tissue_contrast import Prepared, _run_parallel

DONORS = '\t"collaborator_participant_id: GTEX-A"\t"collaborator_participant_id: GTEX-B"'
GTEX_MATRIX = """!Series_title\t"gtex toy"
!Sample_title\t"GTEX-A-1"\t"GTEX-B-1"\t"GTEX-A-2"
!Sample_geo_accession\t"GSM1"\t"GSM2"\t"GSM3"
!Sample_characteristics_ch1{donors}\t"collaborator_participant_id: GTEX-A"
!Sample_characteristics_ch1\t"Sex: 1"\t"Sex: 2"\t"Sex: 1"
!Sample_characteristics_ch1\t"autolysis_score: Mild"\t"autolysis_score: NA"\t"autolysis_score: Severe"
!Sample_characteristics_ch1\t"tissue: Muscle - Skeletal"\t"tissue: Whole Blood"\t"tissue: Colon - Transverse"
!Sample_characteristics_ch1\t"age: 20-29"\t"age: 60-69"\t"age: 20-29"
!Sample_characteristics_ch1\t"smoker_status: Never"\t"smoker_status: Current"\t"smoker_status: Never"
!Sample_characteristics_ch1\t"sample_plate: P1"\t"sample_plate: P2"\t"sample_plate: P1"
!series_matrix_table_begin
"ID_REF"\t"GSM1"\t"GSM2"\t"GSM3"
!series_matrix_table_end
""".replace("{donors}", DONORS)


def _gz(path: Path, text: str) -> Path:
    with gzip.open(path, "wt") as handle:
        handle.write(text)
    return path


def test_gtex_sheet_maps_titles_bins_and_tissues(tmp_path: Path):
    matrix = _gz(tmp_path / "m.txt.gz", GTEX_MATRIX)
    betas = _gz(tmp_path / "b.txt.gz", '"","GTEX-A-2","GTEX-A-1","GTEX-B-1"\n"cg1",0.1,0.2,0.3\n"cg2",0.4,NA,0.6\n')
    samples, frame = _gtex(matrix, betas)
    assert samples["tissue"].to_list() == ["skeletal muscle", "whole blood", "colon"]
    assert samples["age"].to_list() == [25.0, 65.0, 25.0]  # bin midpoints
    assert samples["sex"].to_list() == ["male", "female", "male"]
    assert samples["donor"].to_list() == ["GTEX-A", "GTEX-B", "GTEX-A"]
    assert frame.columns == ["cpg", "GSM3", "GSM1", "GSM2"]  # titles -> accessions
    assert frame["GSM1"].to_list()[1] is None and frame["GSM3"].to_list() == pytest.approx([0.1, 0.4])

    stray = _gz(tmp_path / "s.txt.gz", '"","GTEX-A-1","GTEX-Z-9"\n"cg1",0.1,0.2\n')
    with pytest.raises(ValueError, match="not GEO sample titles"):
        _gtex(matrix, stray)
    odd = _gz(tmp_path / "o.txt.gz", GTEX_MATRIX.replace("Whole Blood", "Brain - Cortex"))
    with pytest.raises(ValueError, match="unknown GTEx tissues"):
        _gtex(odd, betas)


FIBRE_MATRIX = """!Series_title\t"fibre toy"
!Sample_title\t"s1_pre_WM_rep1"\t"f2_pre_WM_rep1"\t"TIIx_pre_rep1"
!Sample_geo_accession\t"GSM1"\t"GSM2"\t"GSM3"
!Sample_characteristics_ch1\t"Sex: male"\t"Sex: female"\t"Sex: not provided"
!Sample_characteristics_ch1\t"sample_id: s1"\t"sample_id: f2"\t"sample_id: Mix"
!series_matrix_table_begin
"ID_REF"\t"GSM1"\t"GSM2"\t"GSM3"
!series_matrix_table_end
"""


def test_fibre_columns_parse_drop_pools_and_average_replicate_probes(tmp_path: Path):
    matrix = _gz(tmp_path / "m.txt.gz", FIBRE_MATRIX)
    header = "Probe_ID,T2X.T2X-Pre,WM.s1-preWM-2,T1.s1-preT1,T2.f2-12wT2,WM.s1-4wWM,T1.Pooled-T1\n"
    rows = "cg2_BC21,0.9,0.5,0.5,0.5,0.5,0.9\ncg1_BC21,0.9,0.2,0.4,0.6,0.8,0.9\ncg1_TC21,0.9,0.4,0.6,0.8,1.0,0.9\n"
    betas = _gz(tmp_path / "b.csv.gz", header + rows)
    samples, frame = _fibres(matrix, betas)
    assert samples["sample"].to_list() == ["WM.s1-preWM-2", "T1.s1-preT1", "T2.f2-12wT2", "WM.s1-4wWM"]
    assert samples["fibre"].to_list() == ["WM", "T1", "T2", "WM"]
    assert samples["timepoint"].to_list() == ["pre", "pre", "12 weeks", "4 weeks"]
    assert samples["replicate"].to_list() == [2, 1, 1, 1]
    assert samples["sex"].to_list() == ["male", "male", "female", "male"]
    assert frame["cpg"].to_list() == ["cg1", "cg2"]
    assert frame.filter(pl.col("cpg") == "cg1").row(0)[1:] == pytest.approx((0.3, 0.5, 0.7, 0.9))

    bad = _gz(tmp_path / "x.csv.gz", "Probe_ID,WM.s1-preWM,T3.s1-preT3\ncg1_BC21,0.1,0.2\n")
    with pytest.raises(ValueError, match="do not parse"):
        _fibres(matrix, bad)


def test_subsample_draws_distinct_rows_of_the_same_size():
    rng = np.random.default_rng(0)
    for n in (47, 224):
        rows = _subsample(rng, n, 32)
        assert rows.size == 32 and np.unique(rows).size == 32 and rows.max() < n
        assert np.all(np.diff(rows) > 0)
    assert _subsample(rng, 20, 32).tolist() == list(range(20))  # a smaller tissue keeps everyone


def _synthetic_gtex(rng: np.random.Generator, features: list[str], members: int, sizes: dict[str, int]):
    """Every CpG drifts toward 0.5 (entropy gain); set members drift less, except in muscle."""
    out = {}
    for tissue, n in sizes.items():
        ages = rng.uniform(25, 75, n)
        level = rng.choice([0.08, 0.92], len(features)) + rng.uniform(-0.05, 0.05, len(features))
        rate = rng.uniform(0.004, 0.008, len(features))
        if TURNOVER[tissue] != "post-mitotic":
            rate[:members] *= 0.3
        drift = (0.5 - level) * np.clip(rate * (ages[:, None] - 20), 0, 1)
        x = np.clip(level + drift + 0.02 * rng.normal(size=(n, len(features))), 0, 1)
        sex = (rng.random(n) < 0.5).astype(float)[:, None]
        p = Prepared(tissue, x, ages, sex, features, np.array([f"{tissue}{i}" for i in range(n)]), sex)
        for adjustment in ADJUSTMENTS:
            out[(tissue, adjustment)] = p
    return out


def test_gtex_bootstrap_finds_the_unprotected_tissue_and_is_reproducible():
    rng = np.random.default_rng(5)
    features = [f"cg{i}" for i in range(400)]
    probes = pl.DataFrame({"cpg": features, "cpg_island": [i % 2 == 0 for i in range(400)]})
    membership = pl.DataFrame({"feature": features[:80], "group": "housekeeping", "gene": "g"})
    sizes = {"colon": 120, "lung": 90, "skeletal muscle": 40}
    tissues = _synthetic_gtex(rng, features, 80, sizes)
    state = {"tissues": tissues, "probes": probes, "membership": membership, "size": 30}
    serial = pl.DataFrame(_run_parallel(_boot_gtex, list(range(12)), 1, state))
    parallel = pl.DataFrame(_run_parallel(_boot_gtex, list(range(12)), 2, state))
    keys = ["boot", "tissue", "adjustment", "metric", "group"]
    assert serial.sort(keys).equals(parallel.sort(keys))
    assert set(serial["tissue"].to_list()) == set(sizes)

    point = serial.filter(pl.col("boot") == 0).drop("boot")
    summary = summarise_gtex(serial, point)
    entropy = summary.filter((pl.col("metric") == "entropy_slope") & (pl.col("adjustment") == "covariates"))
    for contrast in ("renewing - post-mitotic", "slow - post-mitotic", "colon - skeletal muscle"):
        row = entropy.filter(pl.col("contrast") == contrast)
        assert row.height == 1 and row["ci_high"][0] < 0, contrast  # more protection than in muscle


def test_fibre_bootstrap_compares_the_same_people():
    rng = np.random.default_rng(9)
    features = [f"cg{i}" for i in range(300)]
    probes = pl.DataFrame({"cpg": features, "cpg_island": [i % 3 == 0 for i in range(300)]})
    membership = pl.DataFrame({"feature": features[:60], "group": "housekeeping", "gene": "g"})
    n = 14
    level = rng.uniform(0.2, 0.8, len(features))
    x = {}
    for fibre in FIBRES:
        sd = np.full(len(features), 0.05)
        if fibre == "WM":
            sd[:60] = 0.015  # protected only where some cells divide
        x[fibre] = np.clip(level + sd * rng.normal(size=(n, len(features))), 0, 1)
    fs = FibreSet("pre", x, np.array(["male", "female"] * 7), np.array([f"p{i}" for i in range(n)]), features, level)
    state = {"sets": {"pre": fs}, "probes": probes, "membership": membership}
    boot = pl.DataFrame(_run_parallel(_boot_fibres, list(range(30)), 1, state))
    spread = boot.group_by("fibre", "group").agg(pl.col("matched_diff_sd").std())["matched_diff_sd"]
    assert spread.min() > 0  # every resample draws other people
    summary = summarise_fibres(boot)
    assert set(summary["contrast"].to_list()) == {"WM - T1", "WM - T2"}
    assert (summary["ci_high"] < 0).all()
