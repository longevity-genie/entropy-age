import gzip
from pathlib import Path

import numpy as np
import polars as pl
import pytest

from entropy_age.data.geo import read_series_matrix, series_matrix_url
from entropy_age.data.matrix import FeatureMatrix
from entropy_age.published import LinearClock

SERIES = """!Series_title\t"toy"
!Sample_title\t"a"\t"b"\t"c"
!Sample_geo_accession\t"GSM1"\t"GSM2"\t"GSM3"
!Sample_source_name_ch1\t"blood_control"\t"blood_case"\t"blood_control"
!Sample_characteristics_ch1\t"gender: male"\t"gender: female"\t"gender: male"
!Sample_characteristics_ch1\t"age: 16"\t"age: 40"\t"age: 88"
!Sample_characteristics_ch1\t"male"\t"female"\t"male"
!series_matrix_table_begin
"ID_REF"\t"GSM1"\t"GSM2"\t"GSM3"
"cg1"\t0.1\t0.2\t0.3
"cg2"\t0.9\t\t0.7
!series_matrix_table_end
"""


def test_series_matrix_parsing(tmp_path: Path):
    path = tmp_path / "toy_series_matrix.txt.gz"
    with gzip.open(path, "wt") as handle:
        handle.write(SERIES)
    parsed = read_series_matrix(path)
    assert parsed.samples["sample"].to_list() == ["GSM1", "GSM2", "GSM3"]
    assert parsed.samples["age"].to_list() == ["16", "40", "88"]
    assert "characteristic_3" in parsed.samples.columns
    assert parsed.values.shape == (2, 4)
    assert parsed.values["GSM2"][1] is None


def test_series_url():
    assert series_matrix_url("gse41037").endswith("GSE41nnn/GSE41037/matrix/GSE41037_series_matrix.txt.gz")
    with pytest.raises(ValueError):
        series_matrix_url("GSM1")


def test_feature_matrix_imputes_sporadic_gaps():
    frame = pl.DataFrame(
        {"cpg": ["a", "b", "c"], "s1": [0.1, None, 0.5], "s2": [0.2, 0.8, None], "s3": [0.3, 0.6, None]}
    )
    strict = FeatureMatrix.from_frame(frame)
    assert strict.features == ["a"]
    lenient = FeatureMatrix.from_frame(frame, max_missing=0.4)
    assert lenient.features == ["a", "b"]
    assert lenient.values[0, 1] == pytest.approx(0.7)  # median of 0.8 and 0.6
    assert lenient.select_features(["b"]).values.shape == (3, 1)


def test_linear_clock(tmp_path: Path):
    path = tmp_path / "toy.parquet"
    pl.DataFrame({"feature": ["(Intercept)", "cg2", "cg1"], "coef": [1.0, 2.0, 3.0]}).write_parquet(path)
    clock = LinearClock.read(path)
    assert clock.features == ["cg1", "cg2"]
    matrix = FeatureMatrix(np.array([[1.0, 0.0, 9.0]]), ["s"], ["cg2", "x", "cg1"])
    assert clock.predict_matrix(matrix)[0] == pytest.approx(1 + 3 * 9 + 2 * 1)
