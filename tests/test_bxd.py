"""RM3: the BXD liver series matrix, probe names, strain covariates and the OLS term."""

import gzip

import numpy as np
import polars as pl
import pytest

from entropy_age.experiments.bxd_liver import ols, probes_without_suffix, read_bxd, strain_dummies


def test_the_series_matrix_gives_samples_and_betas(tmp_path):
    path = tmp_path / "GSE199979_series_matrix.txt.gz"
    lines = [
        '!Sample_geo_accession\t"GSM1"\t"GSM2"',
        '!Sample_characteristics_ch1\t"age (days): 365"\t"age (days): 730"',
        '!Sample_characteristics_ch1\t"strain: BXD29"\t"strain: BXD29"',
        '!Sample_characteristics_ch1\t"Sex: F"\t"Sex: F"',
        '!Sample_characteristics_ch1\t"diet: CD"\t"diet: HFD"',
        "!series_matrix_table_begin",
        '"ID_REF"\t"GSM1"\t"GSM2"',
        '"cg01"\t0.5\t0.25',
        '"cg02"\t0.1\t',
        "!series_matrix_table_end",
    ]
    with gzip.open(path, "wt") as h:
        h.write("\n".join(lines) + "\n")
    samples, betas = read_bxd(path)
    assert samples["age_months"].to_list() == pytest.approx([365 / 30.44, 730 / 30.44])
    assert samples["diet"].to_list() == ["CD", "HFD"]
    assert betas["cpg"].to_list() == ["cg01", "cg02"]
    assert betas["GSM2"].to_list() == [0.25, None]


def test_probe_names_lose_the_array_suffix_and_keep_the_first():
    probes = pl.DataFrame({"cpg": ["cg01_BC21", "cg01_TC21", "cg02_BC11"], "entrez": [1, 2, 3]})
    out = probes_without_suffix(probes)
    assert out.rows() == [("cg01", 1), ("cg02", 3)]


def test_strains_with_one_mouse_are_pooled():
    d = strain_dummies(np.array(["A", "A", "B", "C", "D", "D"]))
    # levels: A, D, other (B and C pooled); one column per level but the first
    assert d.shape == (6, 2)
    assert d[2].tolist() == d[3].tolist()


def test_ols_returns_the_named_coefficient():
    rng = np.random.default_rng(0)
    a = rng.normal(0, 1, 200)
    b = rng.normal(0, 1, 200)
    y = 1 + 2 * a - 0.5 * b + rng.normal(0, 0.1, 200)
    coef, se, p = ols(y, np.column_stack([np.ones(200), a, b]), 2)
    assert coef == pytest.approx(-0.5, abs=0.03) and p < 1e-10
