import polars as pl
import pytest

from entropy_age.data import cohorts, reprocessed
from entropy_age.paths import ProjectPaths


def test_fibre_titles_become_the_authors_sample_names():
    assert reprocessed._fibre_id("s203_pre_WM_rep1") == "WM.s203-preWM"
    assert reprocessed._fibre_id("f29_12w_TII_rep1") == "T2.f29-12wT2"
    assert reprocessed._fibre_id("f28_pre_WM_rep2") == "f28_pre_WM_rep2"  # a technical replicate matches nothing
    ids = {"WM.f28-preWM-1", "WM.f28-preWM-2", "T1.f28-preT1"}
    assert reprocessed._fibre_id("f28_pre_WM_rep2", ids) == "WM.f28-preWM-2"  # replicates kept apart in the table
    assert reprocessed._fibre_id("f28_pre_WM_rep1", ids) == "WM.f28-preWM-1"
    assert reprocessed._fibre_id("f28_pre_TI_rep1", ids) == "T1.f28-preT1"


def test_switch_keeps_cohorts_without_raw_idats_on_deposited_betas(tmp_path, monkeypatch):
    monkeypatch.setenv("ENTROPY_AGE_BETAS", "sesame")
    paths = ProjectPaths(tmp_path)
    assert cohorts.layout("hannum", paths).betas.name == "hannum_betas.parquet"
    with pytest.raises(FileNotFoundError, match="no R SeSAMe rebuild"):  # never a quiet fall back to deposited betas
        _ = cohorts.layout("endicott", paths).betas
    monkeypatch.setenv("ENTROPY_AGE_BETAS", "deposited")
    assert cohorts.layout("endicott", paths).betas.name == "endicott_betas.parquet"


def test_cohort_table_renames_by_gsm_or_position(tmp_path):
    paths = ProjectPaths(tmp_path)
    d = reprocessed.series_dir(paths, "GSE1")
    d.mkdir(parents=True)
    pl.DataFrame({"probe": ["cg1", "cg2"], "GSM1": [0.1, 0.2], "GSM2": [0.3, 0.4], "GSM3": [0.5, 0.6]}).write_parquet(
        d / "GSE1.betas.chunk001.parquet"
    )
    pl.DataFrame(
        {
            "sample": ["GSM1", "GSM2", "GSM3"],
            "chip": [1, 1, 1],
            "position": ["R01C01", "R02C01", "R03C01"],
            "frac_dt": [0.99, 0.99, 0.5],
            "platform": ["HM450"] * 3,
        }
    ).write_csv(d / "GSE1.qc.csv")
    samples = tmp_path / "s.parquet"
    pl.DataFrame({"sample": ["GSM1", "1_R02C01", "GSM3"]}).write_parquet(samples)
    out = pl.read_parquet(reprocessed.cohort_betas(paths, "GSE1", "x", samples))
    assert out.columns == ["cpg", "GSM1", "1_R02C01"]  # GSM3 failed detection


def test_rebuilt_runs_write_beside_the_deposited_outputs(tmp_path, monkeypatch):
    from entropy_age.experiments.common import source_suffix

    monkeypatch.setenv("ENTROPY_AGE_BETAS", "sesame")
    assert source_suffix() == "_rsesame"
    monkeypatch.delenv("ENTROPY_AGE_BETAS")
    assert source_suffix() == ""


def test_processing_label_follows_the_switch(tmp_path, monkeypatch):
    from entropy_age.data.cohorts import layout
    from entropy_age.paths import ProjectPaths

    paths = ProjectPaths(tmp_path)
    monkeypatch.delenv("ENTROPY_AGE_BETAS", raising=False)
    assert layout("jaffe", paths).processing == "deposited"
    monkeypatch.setenv("ENTROPY_AGE_BETAS", "sesame")
    assert layout("jaffe", paths).processing == "rebuilt"
    assert layout("hannum", paths).processing == "deposited"  # no raw IDATs


def test_subset_keeps_rows_and_drops_constant_covariates():
    import numpy as np

    from entropy_age.experiments.tissue_contrast import Prepared, _subset

    p = Prepared("t", np.arange(12.0).reshape(4, 3), np.array([20.0, 30, 40, 90]),
                 np.array([[1.0], [0.0], [1.0], [0.0]]), ["a", "b", "c"], np.array(list("wxyz")), None)
    keep = p.ages <= 40
    q = _subset(p, keep, np.array([[0.2], [0.2], [0.2], [0.5]]))
    assert q.x.shape == (3, 3) and q.covariates.shape == (3, 1)  # the extra column is constant on the kept rows
