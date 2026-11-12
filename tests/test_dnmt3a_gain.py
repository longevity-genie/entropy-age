"""RM12: Dnmt3a gain-of-function samples parsed from titles, and age centring on the controls."""

import gzip

import numpy as np

from entropy_age.experiments.dnmt3a_gain import age_centred, read_mice


def test_titles_give_column_tissue_age_and_genotype(tmp_path):
    path = tmp_path / "m.txt.gz"
    with gzip.open(path, "wt") as h:
        h.write('!Sample_title\t"k546_1_bm_m_4d: Bone marrow 4-day old W326R/+"\t"k1_li_w_1y: Liver 1-year old +/+"\n')
    m = read_mice(path)
    assert m["column"].to_list() == ["k546_1_bm_m_4d", "k1_li_w_1y"]
    assert m["tissue"].to_list() == ["Bone marrow", "Liver"]
    assert m["age"].to_list() == ["4-day", "1-year"] and m["mutant"].to_list() == [True, False]


def test_age_centring_removes_age_means_of_controls_and_keeps_genotype_differences():
    age = np.array(["a", "a", "a", "b", "b", "b"])
    control = np.array([True, True, False, True, True, False])
    x = np.array([[0.2], [0.2], [0.3], [0.6], [0.6], [0.7]])  # age b is +0.4; mutants +0.1 at both ages
    out = age_centred(x, age, control)
    assert np.allclose(out[control], 0.4)
    assert np.allclose(out[~control], 0.5)
