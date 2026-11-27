import numpy as np

from entropy_age.experiments.identity_loss import (
    classify_fragments,
    fe_ols,
    fidelity,
    overlap_pieces,
    select_u_markers,
    u_shares,
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


def test_overlap_pieces_follows_a_long_read_over_many_markers():
    m_start = np.array([0, 5, 10, 15, 20, 25])
    m_end = m_start + 3
    rows, marks, offset, length = overlap_pieces(np.array([1]), np.array([26]), m_start, m_end)
    got = sorted(zip(marks.tolist(), offset.tolist(), length.tolist(), strict=True))
    assert got == [(0, 0, 2), (1, 4, 3), (2, 9, 3), (3, 14, 3), (4, 19, 3), (5, 24, 2)]
    assert rows.tolist() == [0] * 6


def test_u_shares_masks_thin_regions():
    s = u_shares(np.array([[4, 1]]), np.array([[8, 2]]), min_frags=5)
    assert s[0, 0] == 0.5 and np.isnan(s[0, 1])


def test_select_u_markers_ranks_by_margin_and_ignores_held_out_rows():
    # regions: 0 strong marker of A, 1 weaker marker of A, 2 unmethylated everywhere, 3 A marker only in the held-out
    # A sample, 4 thin in one A training sample
    share = np.array(
        [
            [0.9, 0.6, 0.9, 0.0, np.nan],  # A, train
            [0.8, 0.7, 0.9, 0.0, 0.9],  # A, train
            [0.0, 0.0, 0.9, 0.9, 0.9],  # A, held out
            [0.0, 0.1, 0.8, 0.0, 0.0],  # B group 1
            [0.1, 0.1, 0.9, 0.0, 0.0],  # C (no type), group 2
        ],
        dtype=np.float32,
    )
    type_of = np.array(["A", "A", "A", "B", None], dtype=object)
    group_of = np.array(["a", "a", "a", "b", "c"], dtype=object)
    train = np.array([True, True, False, True, True])
    sel = select_u_markers(share, train, type_of, group_of, "A", n_top=5, min_delta=0.4, tg_quant=0.25)
    assert sel.tolist() == [0, 1]
    assert select_u_markers(share, train, type_of, group_of, "A", n_top=1, min_delta=0.4).tolist() == [0]
    assert select_u_markers(share, train & (type_of != "A"), type_of, group_of, "A").size == 0


def test_loso_summary_names_the_largest_change():
    import polars as pl

    from entropy_age.experiments.identity_loss import _loso_summary

    loso = pl.DataFrame(
        {
            "outcome": ["f"] * 3,
            "term": ["age"] * 3,
            "title": ["x", "y", "z"],
            "purity": [0.99, 0.7, None],
            "estimate_full": [-0.01] * 3,
            "estimate_without": [-0.011, 0.002, -0.009],
            "change": [-0.001, 0.012, 0.001],
        }
    )
    row = _loso_summary(loso).row(0, named=True)
    assert row["min_without_one"] == -0.011 and row["max_without_one"] == 0.002
    assert row["most_influential"].startswith("y (purity 0.7")


def test_marker_counts_are_the_same_across_read_chunks(tmp_path, monkeypatch):
    import gzip

    from entropy_age.experiments import identity_loss as il

    # CpG indices; regions [10, 15) and [20, 25), end exclusive
    pat = tmp_path / "x.pat.gz"
    lines = [
        "chr1\t8\tTTTTT\t2",  # 8-12: 3 CpGs in region 0, all U
        "chr1\t12\tTTCTTTTTTTTT\t1",  # 12-23: 3 in region 0 (T T C), 4 in region 1 (all U)
        "chr1\t15\tCCCCC\t3",  # 15-19: between regions
        "chr1\t21\tTTTTTTT\t1",  # 21-27: 4 in region 1
        "chr1\t24\tCCC\t5",  # 24-26: only 1 CpG in region 1, not informative
    ]
    with gzip.open(pat, "wt") as fh:
        fh.write("\n".join(lines) + "\n")
    m_start, m_end = np.array([10, 20]), np.array([15, 25])
    whole = il.marker_counts(pat, m_start, m_end, light=True)
    monkeypatch.setattr(il, "READ_CHUNK", 2)
    chunked = il.marker_counts(pat, m_start, m_end, light=True)
    assert whole.equals(chunked)
    assert whole["informative"].to_list() == [3, 2]
    assert whole["u_full"].to_list() == [2, 2]
    full = il.marker_counts(pat, m_start, m_end)
    assert full["calls"].to_list() == [2 * 3 + 3, 4 + 4 + 5 * 1]
    assert full["meth_calls"].to_list() == [1, 5]


class _Ctx:
    def rng(self, stream):
        return np.random.default_rng(stream)


def _frame(n=30, seed=0):
    import polars as pl

    rng = np.random.default_rng(seed)
    group = np.repeat(["a", "b", "c"], n // 3)
    age = rng.uniform(20, 80, n)
    hk, matched = rng.normal(size=n), rng.normal(size=n)
    y = 0.5 * hk + 0.1 * matched + rng.normal(0, 0.05, n)
    return pl.DataFrame(
        {
            "sample": [f"s{i}" for i in range(n)],
            "title": [f"t{i}" for i in range(n)],
            "group": group,
            "patient": [f"p{i}" for i in range(n)],
            "age": age,
            "female": rng.random(n) < 0.5,
            "purity": rng.uniform(0.8, 1.0, n),
            "entropy_housekeeping": hk,
            "entropy_matched": matched,
            "entropy_all": hk + matched,
            "y": y,
        }
    )


def test_side_rows_difference_comes_from_the_same_draws():
    from entropy_age.experiments import identity_loss as il

    il_n = il.N_BOOT
    il.N_BOOT = 200
    try:
        rows = il._side_rows(_frame(), "y", _Ctx())
    finally:
        il.N_BOOT = il_n
    link = {r["windows"]: r for r in rows if r["quantity"] == "fidelity on entropy"}
    hk, m, d = link["housekeeping"], link["matched"], link["housekeeping minus matched"]
    assert abs(d["estimate"] - (hk["estimate"] - m["estimate"])) < 1e-12
    assert d["lo_donors"] > 0  # housekeeping carries the link by construction


def test_pooled_rows_skip_nan_rows_and_empty_subsets():
    import polars as pl

    from entropy_age.experiments import identity_loss as il

    f = _frame().with_columns(pl.when(pl.col("sample") == "s0").then(None).otherwise(pl.col("y")).alias("y"))
    il_n = il.N_BOOT
    il.N_BOOT = 20
    try:
        rows = il._pooled_rows(f, ["y"], _Ctx(), controls=("purity",))
        assert rows and all(r["samples"] == 29 for r in rows) and np.isfinite(rows[0]["estimate"])
        assert il._pooled_rows(f.head(0), ["y"], _Ctx()) == []
    finally:
        il.N_BOOT = il_n
