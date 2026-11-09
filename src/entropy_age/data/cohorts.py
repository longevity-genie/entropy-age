"""450k methylation cohorts for the cell-type test (prediction 2 of the selection hypothesis).

* ``hannum`` - GSE40279, whole blood, 656 people aged 19-101 (Hannum et al. 2013). The series
  matrix carries the beta values.
* ``jaffe`` - GSE74193, dorsolateral prefrontal cortex across the lifespan (Jaffe et al. 2016);
  betas in a supplementary CSV whose ``SampleN_Beta`` columns map to GEO samples through the
  sample titles ``SampleN_<group>``. Ships per-sample neuron-fraction estimates.
* ``lunnon`` - GSE59685, whole blood and four brain regions of the same ~117 donors (Lunnon et
  al. 2014), ages 40-105, most with Alzheimer's disease; betas in a supplementary CSV.

Every cohort is normalised to ``data/interim/cohorts/<name>_samples.parquet`` (``sample``,
``donor``, ``tissue``, ``age`` and the cohort's covariate fields) and ``<name>_betas.parquet``
(a matrix frame: ``cpg`` plus one Float32 column per GEO sample).
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

import polars as pl

from entropy_age.data.download import download
from entropy_age.data.geo import read_series_matrix, read_supplementary_betas, series_matrix_url
from entropy_age.paths import ProjectPaths

_SUPPL = "https://ftp.ncbi.nlm.nih.gov/geo/series/{prefix}nnn/{gse}/suppl/{name}"


@dataclass(frozen=True, slots=True)
class CohortSpec:
    name: str
    gse: str
    supplementary: str | None  # beta file; None when the series matrix carries the values
    description: str


COHORTS: dict[str, CohortSpec] = {
    "hannum": CohortSpec("hannum", "GSE40279", None, "whole blood, 656 people, 19-101 y"),
    "jaffe": CohortSpec("jaffe", "GSE74193", "GSE74193_GEO_procData.csv.gz", "prefrontal cortex across the lifespan"),
    "lunnon": CohortSpec("lunnon", "GSE59685", "GSE59685_betas.csv.gz", "blood and 4 brain regions, same donors"),
}


@dataclass(frozen=True, slots=True)
class CohortLayout:
    paths: ProjectPaths
    spec: CohortSpec

    @property
    def input(self) -> Path:
        return self.paths.input / "geo" / self.spec.gse

    @property
    def matrix(self) -> Path:
        return self.input / f"{self.spec.gse}_series_matrix.txt.gz"

    @property
    def supplementary(self) -> Path | None:
        return self.input / self.spec.supplementary if self.spec.supplementary else None

    @property
    def samples(self) -> Path:
        return self.paths.interim / "cohorts" / f"{self.spec.name}_samples.parquet"

    @property
    def betas(self) -> Path:
        return self.paths.interim / "cohorts" / f"{self.spec.name}_betas.parquet"


def layout(name: str, paths: ProjectPaths) -> CohortLayout:
    if name not in COHORTS:
        raise KeyError(f"unknown cohort {name!r}; known: {', '.join(COHORTS)}")
    return CohortLayout(paths, COHORTS[name])


def fetch_cohort(name: str, paths: ProjectPaths, *, overwrite: bool = False) -> list[Path]:
    lay = layout(name, paths)
    gse = lay.spec.gse
    files = [download(series_matrix_url(gse), lay.matrix, overwrite=overwrite, timeout=3600)]
    if lay.supplementary is not None:
        url = _SUPPL.format(prefix=gse[:-3], gse=gse, name=lay.spec.supplementary)
        files.append(download(url, lay.supplementary, overwrite=overwrite, timeout=7200))
    return files


def _float(column: str) -> pl.Expr:
    return pl.col(column).cast(pl.Float64, strict=False)


def _hannum(matrix: Path) -> tuple[pl.DataFrame, pl.DataFrame]:
    parsed = read_series_matrix(matrix)
    samples = parsed.samples.select(
        "sample",
        pl.col("sample").alias("donor"),
        pl.lit("whole blood").alias("tissue"),
        _float("age").alias("age"),
        pl.col("plate"),
        pl.col("gender").alias("sex"),
        pl.col("ethnicity"),
        pl.col("source"),
    )
    return samples, parsed.values.rename({"feature": "cpg"})


def _jaffe(matrix: Path, supplementary: Path) -> tuple[pl.DataFrame, pl.DataFrame]:
    meta = read_series_matrix(matrix).samples
    samples = meta.select(
        "sample",
        pl.col("brnum").alias("donor"),
        pl.lit("prefrontal cortex").alias("tissue"),
        _float("age").alias("age"),
        pl.col("plate"),
        pl.col("sex"),
        pl.col("race"),
        pl.col("group"),
        (pl.col("dropsample").str.to_uppercase() == "TRUE").alias("drop_sample"),
        (pl.col("bestqc").str.to_uppercase() == "TRUE").alias("best_qc"),
        _float("comp_neun_pos").alias("neuron_fraction"),
        pl.col("title").str.split("_").list.first().alias("column"),
    )
    betas = read_supplementary_betas(supplementary, suffix="_Beta")
    rename = dict(zip(samples["column"].to_list(), samples["sample"].to_list(), strict=True))
    missing = [c for c in betas.columns[1:] if c not in rename]
    if missing:
        raise ValueError(f"{len(missing)} beta columns have no sample title, e.g. {missing[:3]}")
    return samples.drop("column"), betas.rename(rename)


def _lunnon(matrix: Path, supplementary: Path) -> tuple[pl.DataFrame, pl.DataFrame]:
    meta = read_series_matrix(matrix).samples
    is_blood = pl.col("source_tissue") == "whole blood"
    samples = meta.select(
        "sample",
        pl.col("subjectid").alias("donor"),
        pl.col("source_tissue").alias("tissue"),
        pl.when(is_blood).then(_float("age.blood")).otherwise(_float("age.brain")).alias("age"),
        pl.col("sex"),
        pl.col("ad.disease.status").alias("ad_status"),
        pl.col("braak.stage").alias("braak_stage"),
        pl.col("barcode").str.split("_").list.first().alias("slide"),
        pl.col("barcode"),
        pl.col("title"),
    )
    betas = read_supplementary_betas(supplementary)
    columns = betas.columns[1:]
    for key in ("sample", "barcode", "title"):
        lookup = dict(zip(samples[key].to_list(), samples["sample"].to_list(), strict=True))
        if all(c in lookup for c in columns):
            return samples.drop("barcode", "title"), betas.rename(lookup)
    raise ValueError(
        f"cannot map the {len(columns)} beta columns of {supplementary.name} (e.g. {columns[:3]}) "
        "to GEO samples by accession, barcode or title"
    )


def prepare_cohort(name: str, paths: ProjectPaths) -> CohortLayout:
    lay = layout(name, paths)
    match name:
        case "hannum":
            samples, betas = _hannum(lay.matrix)
        case "jaffe":
            samples, betas = _jaffe(lay.matrix, lay.supplementary)
        case "lunnon":
            samples, betas = _lunnon(lay.matrix, lay.supplementary)
    lay.samples.parent.mkdir(parents=True, exist_ok=True)
    samples.write_parquet(lay.samples)
    betas.with_columns(pl.exclude("cpg").cast(pl.Float32)).write_parquet(lay.betas)
    return lay
