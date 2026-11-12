"""RM28: the Cox fit recovers a known hazard ratio, handles ties, and the people table reads GEO's fields."""

import gzip

import numpy as np
import pytest

from entropy_age.experiments.survival_shield import cox, read_people


def test_cox_recovers_a_simulated_hazard_ratio_with_censoring():
    rng = np.random.default_rng(5)
    n = 4000
    x = rng.normal(size=(n, 2))
    true = np.array([0.7, -0.3])
    t_event = rng.exponential(1 / np.exp(x @ true))
    t_censor = rng.exponential(1.5, n)
    time, event = np.minimum(t_event, t_censor), t_event <= t_censor
    beta, se = cox(time, event, x)
    assert beta == pytest.approx(true, abs=4 * se.max())
    assert (se > 0).all() and se.max() < 0.05


def test_cox_handles_tied_times_and_null_effect():
    rng = np.random.default_rng(6)
    n = 2000
    x = rng.normal(size=(n, 1))
    time = np.ceil(rng.exponential(1.0, n) * 4)  # heavy ties
    event = rng.uniform(size=n) < 0.7
    beta, se = cox(time, event, x)
    assert abs(beta[0]) < 4 * se[0]


def test_people_read_status_follow_up_and_mmse(tmp_path):
    path = tmp_path / "m.txt.gz"
    with gzip.open(path, "wt") as h:
        h.write(
            '!Sample_geo_accession\t"GSM1"\t"GSM2"\n'
            '!Sample_characteristics_ch1\t"sex (1=woman; 0=man): 1"\t"sex (1=woman; 0=man): 0"\n'
            '!Sample_characteristics_ch1\t"methylation batch: 0"\t"methylation batch: 1"\n'
            '!Sample_characteristics_ch1\t"status (0=alive; 1=died): 1"\t"status (0=alive; 1=died): 0"\n'
            '!Sample_characteristics_ch1\t"follow-up time in days: 701"\t"follow-up time in days: 975"\n'
            '!Sample_characteristics_ch1\t"body mass index (bmi): 19.9"\t"body mass index (bmi): 24.0"\n'
            '!Sample_characteristics_ch1\t""\t"mini mental state examination (mmse): 18"\n'
            '!Sample_characteristics_ch1\t"cd8+cd28- cell proportion: 65.4"\t"cd8+cd28- cell proportion: 34.3"\n'
            "!series_matrix_table_begin\n!series_matrix_table_end\n"
        )
    p = read_people(path)
    assert p["died"].to_list() == [1.0, 0.0]
    assert p["days"].to_list() == [701.0, 975.0]
    assert p["mmse"].to_list() == [None, 18.0]
    assert p["woman"].to_list() == [1.0, 0.0]
