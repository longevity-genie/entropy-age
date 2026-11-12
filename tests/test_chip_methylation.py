"""RM16: patients' clonal-haematopoiesis flags and the beta table keyed by sample title."""

import gzip

import numpy as np
import pytest

from entropy_age.experiments.chip_methylation import read_beta_table, read_patients

MATRIX = """!Sample_title\t"Patient_1"\t"Patient_2"\t"Patient_3"
!Sample_geo_accession\t"GSM1"\t"GSM2"\t"GSM3"
!Sample_characteristics_ch1\t"age: 63"\t"age: 61"\t"age: 40"
!Sample_characteristics_ch1\t"Sex: male"\t"Sex: female"\t"Sex: male"
!Sample_characteristics_ch1\t"dnmt3a mutation: FALSE"\t"dnmt3a mutation: TRUE"\t"dnmt3a mutation: FALSE"
!Sample_characteristics_ch1\t"tet2 mutation: FALSE"\t"tet2 mutation: TRUE"\t"tet2 mutation: FALSE"
!Sample_characteristics_ch1\t"clonal hematopoiesis: FALSE"\t"clonal hematopoiesis: TRUE"\t"clonal hematopoiesis: FALSE"
"""


@pytest.fixture
def files(tmp_path):
    m, b = tmp_path / "m.txt.gz", tmp_path / "b.txt.gz"
    with gzip.open(m, "wt") as h:
        h.write(MATRIX)
    with gzip.open(b, "wt") as h:
        h.write("ID_REF\tPatient_1\tPatient_2\tPatient_3\ncg1\t0.1\t0.2\t0.3\ncg2\t0.4\tNA\t0.6\n")
    return m, b


def test_patients_carry_flags_by_title(files):
    p = read_patients(files[0]).sort("sample")
    assert p["sample"].to_list() == ["Patient_1", "Patient_2", "Patient_3"]
    assert p["ch"].to_list() == [False, True, False]
    assert p["dnmt3a"].to_list() == [False, True, False] and p["tet2"].to_list() == [False, True, False]
    assert p["disease"].to_list() == ["Control", "CH", "Control"]
    assert p["age"].to_list() == [63.0, 61.0, 40.0]


def test_beta_table_follows_requested_samples(files):
    x = read_beta_table(files[1], ["cg2", "cg1", "cg9"], ["Patient_3", "Patient_2"])
    assert x[0, :2].tolist() == [0.6, 0.3]
    assert np.isnan(x[1, 0]) and x[1, 1] == 0.2 and np.isnan(x[:, 2]).all()
