"""RM16 step 4: a second CR blood series."""

import gzip

import numpy as np

from entropy_age.experiments.cr_blood_replication import permutation_p, read_groups


def test_read_groups_splits_titles(tmp_path):
    path = tmp_path / "m.txt.gz"
    with gzip.open(path, "wt") as h:
        h.write('!Sample_title\t"CR1_7855__TLS+CR"\t"CR_T01__Time 0"\n')
    g = read_groups(path)
    assert g["sample"].to_list() == ["CR1_7855", "CR_T01"] and g["group"].to_list() == ["TLS+CR", "Time 0"]


def test_permutation_p_separates_shift_from_noise():
    rng = np.random.default_rng(0)
    treated = np.array([True] * 8 + [False] * 8)
    y = rng.normal(size=16)
    assert permutation_p(y + 3 * treated, treated, 500, rng) < 0.01
    same = np.concatenate([y[:8], y[:8][::-1]])  # equal means in both arms
    assert permutation_p(same, treated, 500, rng) > 0.5
