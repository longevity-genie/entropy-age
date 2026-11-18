"""RM16 step 4: a second CR blood series."""

import gzip

import numpy as np

from entropy_age.experiments.cr_blood_replication import assignments, overlap_chips, read_groups
from entropy_age.experiments.reprogramming_rrbs import permutation_p


def test_read_groups_splits_titles(tmp_path):
    path = tmp_path / "m.txt.gz"
    with gzip.open(path, "wt") as h:
        h.write('!Sample_title\t"CR1_7855__TLS+CR"\t"CR_T01__Time 0"\n')
        h.write(
            '!Sample_supplementary_file\t"ftp://x/GSM1_207211160013_R04C01_Grn.idat.gz"\t"ftp://x/GSM2_207211160010_R01C01_Grn.idat.gz"\n'
        )
    g = read_groups(path)
    assert g["sample"].to_list() == ["CR1_7855", "CR_T01"] and g["group"].to_list() == ["TLS+CR", "Time 0"]
    assert g["chip"].to_list() == ["207211160013", "207211160010"]


def _cr_layout():
    # GSE282661, CR against controls: chip 011 holds 7 CR and 3 controls, chip 013 1 CR and 5 controls
    chips = np.array(["011"] * 10 + ["013"] * 6)
    treated = np.array([True] * 7 + [False] * 3 + [True] + [False] * 5)
    return treated, chips


def test_within_chip_assignments_enumerate_all_720_and_keep_chip_counts():
    treated, chips = _cr_layout()
    draws, total, exact = assignments(treated, chips, 1000, np.random.default_rng(0))
    assert total == 720 and exact and len(draws) == 720
    assert np.array_equal(draws[0], treated)  # the observed assignment first
    assert len({d.tobytes() for d in draws}) == 720
    for d in draws:
        assert d[chips == "011"].sum() == 7 and d[chips == "013"].sum() == 1


def test_sampled_assignments_when_too_many():
    treated, chips = _cr_layout()
    draws, total, exact = assignments(treated, chips, 100, np.random.default_rng(0))
    assert total == 720 and not exact and len(draws) == 100 and np.array_equal(draws[0], treated)


def test_no_within_chip_overlap_is_flagged():
    # controls after 12 weeks against Time 0: Time 0 alone on one chip, controls on two others
    chips = np.array(["010"] * 8 + ["011"] * 3 + ["013"] * 5)
    treated = np.array([False] * 8 + [True] * 8)
    assert overlap_chips(treated, chips) == 0
    t, c = _cr_layout()
    assert overlap_chips(t, c) == 2


def test_permutation_p_separates_shift_from_noise():
    rng = np.random.default_rng(0)
    treated = np.array([True] * 8 + [False] * 8)
    y = rng.normal(size=16)
    assert permutation_p(y + 3 * treated, treated, 500, rng) < 0.01
    same = np.concatenate([y[:8], y[:8][::-1]])  # equal means in both arms
    assert permutation_p(same, treated, 500, rng) > 0.5
