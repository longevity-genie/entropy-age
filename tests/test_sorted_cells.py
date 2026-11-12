"""RM31: reading MESA's normalized M-value tables and people."""

import gzip

import numpy as np
import pytest

from entropy_age.experiments.sorted_cells import read_betas, read_people


def test_betas_come_from_m_values_for_the_named_participants(tmp_path):
    path = tmp_path / "n.txt.gz"
    with gzip.open(path, "wt") as h:
        h.write('"ID_REF"\t"100001.Mvalue"\t"100001.detectionPval"\t"100004.Mvalue"\t"100004.detectionPval"\n')
        h.write('"cg01"\t0\t0\t1\t0\n"cg02"\t-1\t0\tNA\t0\n"cg03"\t2\t0\t3\t0\n')
    x = read_betas(path, ["cg03", "cg01", "cg99"], ["100004", "100001"])
    assert x[0, 0] == pytest.approx(8 / 9) and x[1, 0] == pytest.approx(0.8)
    assert x[0, 1] == pytest.approx(2 / 3) and x[1, 1] == pytest.approx(0.5)
    assert np.isnan(x[:, 2]).all()


def test_people_carry_the_participant_id_age_and_health(tmp_path):
    path = tmp_path / "m.txt.gz"
    with gzip.open(path, "wt") as h:
        h.write(
            '!Sample_title\t"100001_peripheral_CD14 [methylation]"\t"100004_peripheral_CD14 [methylation]"\n'
            '!Sample_geo_accession\t"GSM1"\t"GSM2"\n'
            '!Sample_characteristics_ch1\t"age: 59"\t"age: 63"\n'
            '!Sample_characteristics_ch1\t"racegendersite: 5"\t"racegendersite: 1"\n'
            '!Sample_characteristics_ch1\t"bcell: 0.1"\t"bcell: 0.2"\n'
            '!Sample_characteristics_ch1\t"plaque: 4"\t"plaque: NA"\n'
            '!Sample_characteristics_ch1\t"cac: 73.2"\t"cac: 0"\n'
            "!series_matrix_table_begin\n!series_matrix_table_end\n"
        )
    p = read_people(path)
    assert p["id"].to_list() == ["100001", "100004"]
    assert p["age"].to_list() == [59.0, 63.0]
    assert p["plaque"].to_list() == [4.0, None]
    assert p["cac"].to_list() == [73.2, 0.0]
