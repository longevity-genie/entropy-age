"""RM14: reading a series matrix's people and betas, and the carrier effect."""

import gzip

import numpy as np
import pytest

from entropy_age.experiments.syndromes import effects, read_betas, read_people, with_known_age

MATRIX = """!Series_title\t"x"
!Sample_geo_accession\t"GSM1"\t"GSM2"\t"GSM3"
!Sample_characteristics_ch1\t"gender: female"\t"gender: male"\t"gender: male"
!Sample_characteristics_ch1\t"age (years): 4"\t"age (years): 12"\t"age (years): NA"
!Sample_characteristics_ch1\t"disease state: Control"\t"disease state: Kabuki"\t"disease state: Control"
!Sample_characteristics_ch1\t"tissue: whole blood"\t"tissue: whole blood"\t"tissue: whole blood"
!series_matrix_table_begin
"ID_REF"\t"GSM1"\t"GSM2"\t"GSM3"
"cg01"\t0.1\t0.2\t0.3
"cg02"\t0.4\t\t0.6
"cg03"\t0.7\t0.8\t0.9
!series_matrix_table_end
"""


@pytest.fixture
def matrix(tmp_path):
    path = tmp_path / "m.txt.gz"
    with gzip.open(path, "wt") as handle:
        handle.write(MATRIX)
    return path


def test_people_carry_disease_sex_and_age_with_missing_age_null(matrix):
    people = read_people(matrix)
    assert people["sample"].to_list() == ["GSM1", "GSM2", "GSM3"]
    assert people["disease"].to_list() == ["Control", "Kabuki", "Control"]
    assert people["sex"].to_list() == ["female", "male", "male"]
    assert people["age"].to_list() == [4.0, 12.0, None]


def test_betas_follow_the_requested_samples_and_features_with_gaps_as_nan(matrix):
    x = read_betas(matrix, ["cg03", "cg02", "cg99"], ["GSM3", "GSM2"])
    assert x.shape == (2, 3)
    assert x[0, :2].tolist() == [0.9, 0.6]
    assert x[1, 0] == 0.8 and np.isnan(x[1, 1])
    assert np.isnan(x[:, 2]).all()


def test_effect_is_the_carrier_difference_net_of_covariates():
    rng = np.random.default_rng(0)
    n = 60
    carrier = np.arange(n) < 20
    age = rng.uniform(0, 30, n) + 20 * carrier  # carriers older: the covariate matters
    f = -0.3 + 0.2 * carrier + 0.01 * age + rng.normal(0, 0.01, n)
    rows = effects({"F_housekeeping": f, "P_housekeeping": f, "S_housekeeping": -f}, carrier, age[:, None])
    by = {r["measure"]: r for r in rows}
    assert set(by) == {"F", "S"}  # P is not an effect measure
    assert by["F"]["effect"] == pytest.approx(0.2, abs=0.01)
    assert by["S"]["effect"] == pytest.approx(-0.2, abs=0.01)
    assert by["F"]["control_mean"] == pytest.approx(f[~carrier].mean())


def test_people_without_an_age_are_left_out_only_where_the_series_gives_ages(matrix):
    people = read_people(matrix)
    assert with_known_age(people)["sample"].to_list() == ["GSM1", "GSM2"]
    no_ages = people.with_columns(age=None)
    assert with_known_age(no_ages).height == 3
