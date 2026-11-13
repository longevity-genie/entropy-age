"""RM17: transient reprogramming of human fibroblasts (Gill 2022)."""

import gzip

import numpy as np

from entropy_age.experiments.gill_reprogramming import donor_permutation_p, read_samples, within_donor_effect


def test_read_samples(tmp_path):
    path = tmp_path / "m.txt.gz"
    with gzip.open(path, "wt") as h:
        h.write('!Sample_title\t"O1_transiently_reprogrammed_15days_exp1"\t"O2_negative_control_15days_exp1"\t"O2_x"\n')
        h.write('!Sample_geo_accession\t"GSM1"\t"GSM2"\t"GSM3"\n')
        h.write(
            '!Sample_characteristics_ch1\t"cell type: Transiently reprogrammed fibroblast"\t'
            '"cell type: Negative control fibroblast"\t"cell type: Failed to transiently reprogram fibroblast"\n'
        )
        h.write(
            '!Sample_characteristics_ch1\t"length of reprogramming phase (days): 15"\t"length of reprogramming '
            'phase (days): 15"\t"length of reprogramming phase (days): 13"\n'
        )
        h.write(
            '!Sample_characteristics_ch1\t"donor age (years): 53"\t"donor age (years): 38"\t"donor age (years): 38"\n'
        )
    s = read_samples(path)
    assert s["arm"].to_list() == ["reprogrammed", "control", "failed"]
    assert s["donor"].to_list() == ["O1", "O2", "O2"] and s["donor_age"].to_list() == [53.0, 38.0, 38.0]


def test_within_donor_effect_ignores_donor_offsets():
    donor = np.array(["a", "a", "b", "b"])
    treated = np.array([True, False, True, False])
    y = np.array([11.0, 10.0, 101.0, 100.0])
    assert within_donor_effect(y, treated, donor) == 1.0


def test_permutation_p_small_for_consistent_effect():
    rng = np.random.default_rng(0)
    donor = np.repeat(["a", "b", "c"], 8)
    treated = np.tile([True] * 4 + [False] * 4, 3)
    y = treated * 2.0 + rng.normal(0, 0.1, donor.size) + np.repeat([0.0, 5.0, 10.0], 8)
    assert donor_permutation_p(y, treated, donor, 200, rng) < 0.02
