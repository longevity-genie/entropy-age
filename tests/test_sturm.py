"""GSE179847 (Sturm et al. 2022): the sample sheet and the replicate flasks of the culture test."""

import gzip
from pathlib import Path

import polars as pl
import pytest

from entropy_age.data.cohorts import _characteristics, _sturm
from entropy_age.experiments.passaged_cells import _replicate_groups

# Three samples: a baseline control, a DEX flask split from a culture at the treatment's start, and a
# SURF1 control. The second carries an extra characteristic ("notes") ahead of the others, as on GEO,
# so its keys sit one row below the first sample's.
_KEYS = [
    ("cell_line", ["HC1", "HC1", "SURF1_1"]),
    ("cell_line_group", ["HC1_c", "HC1_dex", "SURF1_1_c"]),
    ("clinical_condition", ["Normal", "Normal", "SURF1_Mutation"]),
    ("Sex", ["male", "male", "female"]),
    ("donor_age", ["29", "29", "0"]),
    ("study_part", ["2", "2", "2"]),
    ("treatments", ["Control", "DEX", "Control"]),
    ("percent_oxygen", ["21", "21", "21"]),
    ("passage", ["9", "9", "5"]),
    ("days_grown_udays", ["24.5", "24.5", "10"]),
    ("days_after_treatment_udays", ["NA", "-2", "NA"]),
    ("population_doublings_utotal_divisions", ["17.9", "17.9", "6.1"]),
    ("sample_plate", ["P1", "P1", "P2"]),
    ("slide", ["201", "202", "203"]),
    ("array", ["R01C01", "R02C01", "R03C01"]),
]


def _series(path: Path) -> Path:
    rows = [[f"{k}: {v[i]}" for k, v in _KEYS] for i in range(3)]
    rows[1].insert(3, "notes: Treatment Start")
    rows[0].append("")
    rows[2].append("")
    lines = ['!Sample_geo_accession\t"GSM1"\t"GSM2"\t"GSM3"']
    for r in range(len(rows[1])):
        lines.append("!Sample_characteristics_ch1\t" + "\t".join(f'"{rows[i][r]}"' for i in range(3)))
    with gzip.open(path, "wt") as handle:
        handle.write("\n".join(lines) + "\n!series_matrix_table_begin\n!series_matrix_table_end\n")
    return path


def _betas(path: Path, columns: list[str]) -> Path:
    header = ['""'] + [f'"{c} {kind}"' for c in columns for kind in ("beta", "Detection Pval")]
    body = [f'"cg{j:08d}",' + ",".join(f"{0.1 * (j + 1):.2f},0.001" for _ in columns) for j in range(3)]
    with gzip.open(path, "wt") as handle:
        handle.write(",".join(header) + "\n" + "\n".join(body) + "\n")
    return path


def test_characteristics_follow_each_sample_s_keys(tmp_path: Path):
    meta = _characteristics(_series(tmp_path / "s.txt.gz"))
    assert meta["sample"].to_list() == ["GSM1", "GSM2", "GSM3"]
    assert meta["treatments"].to_list() == ["Control", "DEX", "Control"]  # shifted one row in GSM2
    assert meta["notes"].to_list() == [None, "Treatment Start", None]


def test_sturm_sheet_names_arms_cultures_and_columns(tmp_path: Path):
    series = _series(tmp_path / "s.txt.gz")
    samples, betas = _sturm(series, _betas(tmp_path / "b.csv.gz", ["201_R01C01", "202_R02C01", "203_R03C01"]))
    by = {r["sample"]: r for r in samples.iter_rows(named=True)}
    assert by["GSM1"]["subexperiment"] == "Baseline profiling"
    assert by["GSM2"]["subexperiment"] == "Study part 2"
    assert by["GSM3"]["subexperiment"] == "SURF1 mutation"
    assert (by["GSM1"]["agent"], by["GSM2"]["agent"], by["GSM1"]["oxygen"]) == ("control", "DEX", "Ambient")
    assert by["GSM2"]["culture"] == "HC1|9|24.5"  # before the treatment: the culture it was split from
    assert by["GSM1"]["culture"] == "HC1_c"
    assert by["GSM2"]["doublings"] == pytest.approx(17.9)
    assert betas.columns == ["cpg", "GSM1", "GSM2", "GSM3"]  # detection p-values dropped, barcodes renamed


def test_sturm_refuses_a_matrix_that_misses_a_sample(tmp_path: Path):
    series = _series(tmp_path / "s.txt.gz")
    with pytest.raises(ValueError, match="1 samples have no column"):
        _sturm(series, _betas(tmp_path / "b.csv.gz", ["201_R01C01", "202_R02C01"]))


def test_replicates_are_flasks_of_one_culture_when_the_cohort_names_it():
    unit = pl.DataFrame(
        {
            "row": [0, 1, 2, 3],
            "days": [5.0, 5.0, 9.0, 9.0],
            "arm": ["agent=DEX", "agent=MitoQ", "agent=DEX", "agent=MitoQ"],
            "culture": ["HC5|4|5", "HC5|4|5", "dex", "mitoq"],
        }
    )
    groups, day = _replicate_groups(unit)
    assert day == 5.0 and [g.tolist() for g in groups] == [[0, 1]]
    groups, day = _replicate_groups(unit.drop("culture"))  # by arm: no arm holds two samples of one day
    assert groups == [] and day is None
