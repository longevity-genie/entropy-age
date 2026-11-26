import numpy as np

from entropy_age.experiments.identity_loss import (
    classify_fragments,
    fe_ols,
    fidelity,
    overlap_pieces,
    young_old_split,
)


def _calls(strings):
    meth = np.array([s.count("C") for s in strings])
    unmeth = np.array([s.count("T") for s in strings])
    return meth, unmeth


def test_classify_fragments_strict_and_uxm_patterns():
    meth, unmeth = _calls(["TTT", "TTC", "TTTTC", "CCCC", "TT", "T.T.T", "CCCCT"])
    f = classify_fragments(meth, unmeth, min_cpgs=3, u_max=0.25)
    assert f["informative"].tolist() == [True, True, True, True, False, True, True]
    assert f["u_full"].tolist() == [True, False, False, False, False, True, False]
    assert f["m_full"].tolist() == [False, False, False, True, False, False, False]
    # 1 of 5 methylated (0.2) passes UXM's 25 % rule, 1 of 3 (0.33) does not
    assert f["u_uxm"].tolist() == [True, False, True, False, False, True, False]
    assert f["m_uxm"].tolist() == [False, False, False, True, False, False, True]


def test_overlap_pieces_cuts_reads_to_the_marker():
    m_start, m_end = np.array([10, 20]), np.array([15, 25])
    # read 0 inside marker 0; read 1 straddles marker 0's end; read 2 spans both markers; read 3 between markers
    starts, lens = np.array([10, 13, 12, 16]), np.array([3, 5, 12, 3])
    rows, marks, offset, length = overlap_pieces(starts, lens, m_start, m_end)
    got = sorted(zip(rows.tolist(), marks.tolist(), offset.tolist(), length.tolist(), strict=True))
    assert got == [(0, 0, 0, 3), (1, 0, 0, 2), (2, 0, 0, 3), (2, 1, 8, 4)]


def test_fidelity_skips_thin_markers_and_other_types():
    u = np.array([[8, 1, 0], [5, 5, 5]])
    frags = np.array([[10, 2, 10], [10, 10, 10]])
    use = np.array([[True, True, False], [False, True, True]])
    assert np.allclose(fidelity(u, frags, use, min_frags=5), [0.8, 0.5])


def test_fe_ols_recovers_within_group_slope():
    rng = np.random.default_rng(1)
    group = np.repeat([0, 1, 2], 10)
    age = rng.uniform(20, 80, 30)
    y = 5.0 * group - 0.01 * age + rng.normal(0, 1e-4, 30)
    coef, se, df = fe_ols(y, age[:, None], group)
    assert abs(coef[0] + 0.01) < 1e-4 and df == 26


def test_young_old_split_trains_below_the_class_median():
    label = np.array(["a", "a", "a", "b", "b", "b", "b"])
    age = np.array([20, 40, 60, 30, 30, 50, 70])
    assert young_old_split(label, age).tolist() == [True, False, False, True, True, False, False]
