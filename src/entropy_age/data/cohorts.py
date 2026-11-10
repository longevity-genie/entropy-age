"""450k methylation cohorts for the cell-type test (prediction 2 of the selection hypothesis).

* ``hannum`` - GSE40279, whole blood, 656 people aged 19-101 (Hannum et al. 2013). The series
  matrix carries the beta values.
* ``jaffe`` - GSE74193, dorsolateral prefrontal cortex across the lifespan (Jaffe et al. 2016);
  betas in a supplementary CSV whose ``SampleN_Beta`` columns map to GEO samples through the
  sample titles ``SampleN_<group>``. Ships per-sample neuron-fraction estimates.
* ``lunnon`` - GSE59685, whole blood and four brain regions of the same ~117 donors (Lunnon et
  al. 2014), ages 40-105, most with Alzheimer's disease; betas in a supplementary CSV.
* ``huang`` - GSE279437, male mice, 11 tissues of the same animals on the mammalian array
  (HorvathMammalMethylChip320), young and aged, aged ones treated with exenatide (a GLP-1
  receptor agonist), rapamycin or vehicle (Huang et al.). The betas and a sample sheet
  (barcode, tissue, treatment, animal) are supplementary files; GEO gives no numeric ages.
* ``gtex`` - GSE213478, nine GTEx tissues (blood, colon, testis, lung, kidney, breast, prostate,
  ovary, skeletal muscle) of 424 post-mortem donors on the EPIC array (Oliva et al. 2023). Many
  donors gave several tissues. GEO gives ages as decade bins only (``20-29``); ``age`` is the bin
  midpoint and ``age_bin`` keeps the label. The supplementary betas (noob + BMIQ) are named by
  sample title.
* ``fibres`` - GSE304045, vastus lateralis of about 12 adults on EPIC v2: pooled type I fibres,
  pooled type II fibres and whole muscle of each person, before and after 12 weeks of interval
  training. The processed matrix names its columns ``<type>.<person>-<time><type>[-<rep>]``, which
  become the sample IDs here (GEO titles use another spelling). EPIC v2 probe IDs carry a suffix
  (``cg00002271_BC21``); replicate probes of one CpG are averaged. GEO gives no ages.

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
    series_matrix: bool = True  # False when a supplementary sample sheet replaces it
    metadata: str | None = None  # supplementary sample sheet


COHORTS: dict[str, CohortSpec] = {
    "hannum": CohortSpec("hannum", "GSE40279", None, "whole blood, 656 people, 19-101 y"),
    "jaffe": CohortSpec("jaffe", "GSE74193", "GSE74193_GEO_procData.csv.gz", "prefrontal cortex across the lifespan"),
    "lunnon": CohortSpec("lunnon", "GSE59685", "GSE59685_betas.csv.gz", "blood and 4 brain regions, same donors"),
    "huang": CohortSpec(
        "huang",
        "GSE279437",
        "GSE279437_Methylation_betas.csv.gz",
        "mouse, 11 tissues of the same animals, young/aged, exenatide/rapamycin",
        series_matrix=False,
        metadata="GSE279437_Methylation_metadata.csv.gz",
    ),
    "gtex": CohortSpec(
        "gtex",
        "GSE213478",
        "GSE213478_methylation_DNAm_noob_final_BMIQ_all_tissues_987.txt.gz",
        "human GTEx, 9 tissues of 424 donors, EPIC, decade age bins",
    ),
    "fibres": CohortSpec(
        "fibres",
        "GSE304045",
        "GSE304045_dasen_normalised_filtered_matrix_processed.csv.gz",
        "human muscle: type I, type II fibres and whole muscle of the same people, EPIC v2",
    ),
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
    def metadata(self) -> Path | None:
        return self.input / self.spec.metadata if self.spec.metadata else None

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
    files = []
    if lay.spec.series_matrix:
        files.append(download(series_matrix_url(gse), lay.matrix, overwrite=overwrite, timeout=3600))
    for target in (lay.metadata, lay.supplementary):
        if target is not None:
            url = _SUPPL.format(prefix=gse[:-3], gse=gse, name=target.name)
            files.append(download(url, target, overwrite=overwrite, timeout=7200))
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
    # 3 comment lines and a blank line, then a barcode header, then rows of new and old GSM accessions
    betas = read_supplementary_betas(supplementary, preamble=4, after_header=2)
    columns = betas.columns[1:]
    for key in ("sample", "barcode", "title"):
        lookup = dict(zip(samples[key].to_list(), samples["sample"].to_list(), strict=True))
        if all(c in lookup for c in columns):
            return samples.drop("barcode", "title"), betas.rename(lookup)
    raise ValueError(
        f"cannot map the {len(columns)} beta columns of {supplementary.name} (e.g. {columns[:3]}) "
        "to GEO samples by accession, barcode or title"
    )


#: GSE279437 tissue codes -> names
HUANG_TISSUES = {
    "WBC": "white blood cells",
    "Colon": "colon",
    "Spleen": "spleen",
    "Liver": "liver",
    "Kidney": "kidney",
    "Adipose": "adipose",
    "SkeletalMuscle": "skeletal muscle",
    "CardiacMuscle": "heart",
    "FCTX": "frontal cortex",
    "Hippo": "hippocampus",
    "HTH": "hypothalamus",
}


def _huang(metadata: Path, supplementary: Path) -> tuple[pl.DataFrame, pl.DataFrame]:
    sheet = pl.read_csv(metadata)
    unknown = sorted(set(sheet["tissue"].to_list()) - set(HUANG_TISSUES))
    if unknown:
        raise ValueError(f"unknown tissue codes in {metadata.name}: {unknown}")
    group = pl.col("treatment")
    samples = sheet.select(
        pl.col("sname").alias("sample"),
        pl.col("animal").alias("donor"),
        pl.col("tissue").replace_strict(HUANG_TISSUES).alias("tissue"),
        pl.lit(None, dtype=pl.Float64).alias("age"),
        group.alias("group"),
        group.str.split("_").list.first().str.to_lowercase().alias("age_group"),
        group.str.split("_").list.get(1).alias("treatment"),
        group.str.ends_with("_KD").alias("knockdown"),
        pl.col("sname").str.split("_").list.first().alias("chip"),
    )
    if not set(samples["age_group"].to_list()) <= {"young", "aged"}:
        raise ValueError(f"unexpected groups in {metadata.name}: {sorted(set(samples['group'].to_list()))}")
    betas = read_supplementary_betas(supplementary)
    missing = sorted(set(samples["sample"].to_list()) - set(betas.columns[1:]))
    if missing:
        raise ValueError(f"{len(missing)} samples of the sheet have no beta column, e.g. {missing[:3]}")
    return samples, betas.select("cpg", *samples["sample"].to_list())


#: GSE213478 tissue labels -> names
GTEX_TISSUES = {
    "Whole Blood": "whole blood",
    "Colon - Transverse": "colon",
    "Testis": "testis",
    "Lung": "lung",
    "Kidney - Cortex": "kidney",
    "Breast - Mammary Tissue": "breast",
    "Prostate": "prostate",
    "Ovary": "ovary",
    "Muscle - Skeletal": "skeletal muscle",
}


def _gtex(matrix: Path, supplementary: Path) -> tuple[pl.DataFrame, pl.DataFrame]:
    meta = read_series_matrix(matrix).samples
    unknown = sorted(set(meta["tissue"].to_list()) - set(GTEX_TISSUES))
    if unknown:
        raise ValueError(f"unknown GTEx tissues in {matrix.name}: {unknown}")
    bounds = pl.col("age").str.extract_all(r"\d+")
    samples = meta.select(
        "sample",
        pl.col("collaborator_participant_id").alias("donor"),
        pl.col("tissue").replace_strict(GTEX_TISSUES).alias("tissue"),
        ((bounds.list.first().cast(pl.Float64) + bounds.list.last().cast(pl.Float64) + 1) / 2).alias("age"),
        pl.col("age").alias("age_bin"),
        pl.col("sex").replace_strict({"1": "male", "2": "female"}).alias("sex"),
        pl.col("autolysis_score").alias("autolysis"),
        pl.col("smoker_status").alias("smoker"),
        pl.col("sample_plate").alias("plate"),
        pl.col("title"),
    )
    if samples["age"].null_count():
        raise ValueError(f"age bins without two numbers: {samples.filter(pl.col('age').is_null())['age_bin'].unique()}")
    betas = read_supplementary_betas(supplementary)
    rename = dict(zip(samples["title"].to_list(), samples["sample"].to_list(), strict=True))
    missing = [c for c in betas.columns[1:] if c not in rename]
    if missing:
        raise ValueError(f"{len(missing)} beta columns are not GEO sample titles, e.g. {missing[:3]}")
    return samples.drop("title"), betas.rename(rename)


_FIBRE_COLUMN = r"^(?<type>T1|T2|WM)\.(?<person>[sf]\d+)-(?<time>pre|4w|12w)(?:T1|T2|WM)(?:-(?<rep>\d+))?$"
_FIBRE_POOLED = r"^(?:T2X\.T2X|T[12]\.Pooled)-"  # pools of fibres from every person


def _fibres(matrix: Path, supplementary: Path) -> tuple[pl.DataFrame, pl.DataFrame]:
    sex = (
        read_series_matrix(matrix)
        .samples.select(pl.col("sample_id").alias("donor"), "sex")
        .filter(pl.col("sex").is_in(["male", "female"]))
        .unique()
    )
    if sex["donor"].is_duplicated().any():
        raise ValueError(f"{matrix.name} gives two sexes for one person")
    betas = read_supplementary_betas(supplementary)
    columns = pl.DataFrame({"sample": betas.columns[1:]}).filter(~pl.col("sample").str.contains(_FIBRE_POOLED))
    parts = columns.with_columns(parsed=pl.col("sample").str.extract_groups(_FIBRE_COLUMN)).unnest("parsed")
    bad = parts.filter(pl.col("type").is_null())["sample"].to_list()
    if bad:
        raise ValueError(f"{len(bad)} columns of {supplementary.name} do not parse, e.g. {bad[:3]}")
    samples = parts.select(
        "sample",
        pl.col("person").alias("donor"),
        pl.lit("skeletal muscle").alias("tissue"),
        pl.lit(None, dtype=pl.Float64).alias("age"),
        pl.col("type").alias("fibre"),
        pl.col("time").replace_strict({"pre": "pre", "4w": "4 weeks", "12w": "12 weeks"}).alias("timepoint"),
        pl.col("rep").cast(pl.Int64).fill_null(1).alias("replicate"),
    ).join(sex, on="donor", how="left")
    probe = pl.col("cpg").str.replace(r"_[A-Z]{2}\d+$", "")
    collapsed = (
        betas.select(probe, *samples["sample"].to_list())
        .group_by("cpg", maintain_order=True)
        .agg(pl.exclude("cpg").mean())
        .sort("cpg")
    )
    return samples, collapsed


def prepare_cohort(name: str, paths: ProjectPaths) -> CohortLayout:
    lay = layout(name, paths)
    match name:
        case "hannum":
            samples, betas = _hannum(lay.matrix)
        case "jaffe":
            samples, betas = _jaffe(lay.matrix, lay.supplementary)
        case "lunnon":
            samples, betas = _lunnon(lay.matrix, lay.supplementary)
        case "huang":
            samples, betas = _huang(lay.metadata, lay.supplementary)
        case "gtex":
            samples, betas = _gtex(lay.matrix, lay.supplementary)
        case "fibres":
            samples, betas = _fibres(lay.matrix, lay.supplementary)
    lay.samples.parent.mkdir(parents=True, exist_ok=True)
    samples.write_parquet(lay.samples)
    betas.with_columns(pl.exclude("cpg").cast(pl.Float32)).write_parquet(lay.betas)
    return lay
