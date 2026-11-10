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
* ``sturm`` - GSE179847, the Cellular Lifespan Study (Sturm et al. 2022): primary human dermal fibroblasts
  of 6 healthy donors and 3 SURF1 patients passaged for up to 265 days on EPIC, with mitochondrial,
  metabolic, glucocorticoid, contact-inhibition and 3 % oxygen arms. No sample is measured twice in one arm;
  the flasks split from one culture at a treatment's start are, and serve as replicates (``culture``).
* ``endicott`` - GSE197512, primary human cells (fetal, neonatal and adult fibroblasts, endothelial
  and smooth muscle cells, keratinocytes) passaged to replicative senescence on EPIC (Endicott et al.
  2022), with arms under mitomycin C growth arrest, serum withdrawal, 3 % oxygen and TERT
  immortalisation. Every sample carries its population doublings and days in culture; ``donor`` is
  the cell line, ``age`` is null (donor age stays in ``donor_age``). The series matrix carries the betas.

Every cohort is normalised to ``data/interim/cohorts/<name>_samples.parquet`` (``sample``,
``donor``, ``tissue``, ``age`` and the cohort's covariate fields) and ``<name>_betas.parquet``
(a matrix frame: ``cpg`` plus one Float32 column per GEO sample).
"""

from __future__ import annotations

import gzip
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
    "endicott": CohortSpec(
        "endicott",
        "GSE197512",
        None,
        "human primary cells passaged to senescence, with growth-arrested arms, EPIC",
    ),
    "sturm": CohortSpec(
        "sturm",
        "GSE179847",
        "GSE179847_Cell_lifespan_DNAm_processed_matrix.csv.gz",
        "human dermal fibroblasts passaged for up to 265 days, with treatment and contact-inhibition arms, EPIC",
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


def _endicott(matrix: Path) -> tuple[pl.DataFrame, pl.DataFrame]:
    parsed = read_series_matrix(matrix)
    meta = parsed.samples
    samples = meta.select(
        "sample",
        pl.col("cell_line").alias("donor"),
        pl.col("source_name_ch1").alias("tissue"),
        pl.lit(None, dtype=pl.Float64).alias("age"),
        "title",
        "subexperiment",
        "donor_age",
        pl.col("donor_sex").alias("sex"),
        pl.col("culture_oxygen").alias("oxygen"),
        pl.col("culture_agent").alias("agent"),
        pl.col("culture_pctfbs").alias("serum"),
        pl.col("expression_vector").alias("vector"),
        "subculture",
        _float("population_doublings").alias("doublings"),
        _float("days_in_culture").alias("days"),
        pl.col("batch_date").alias("batch"),
    )
    bad = samples.filter(pl.col("doublings").is_null() | pl.col("days").is_null())
    if bad.height:
        raise ValueError(
            f"{bad.height} samples without population doublings or days, e.g. {bad['sample'][:3].to_list()}"
        )
    return samples, parsed.values.rename({"feature": "cpg"})


def _characteristics(matrix: Path) -> pl.DataFrame:
    """Every ``key: value`` characteristic per sample, whatever row of the series matrix it sits in.

    GSE179847 gives some samples an extra characteristic, so one key sits in different rows for
    different samples and the row-wise reading of ``read_series_matrix`` names columns wrongly.
    """
    accessions: list[str] = []
    rows: list[list[str]] = []
    with gzip.open(matrix, "rt") as handle:
        for line in handle:
            cells = [c.strip().strip('"') for c in line.rstrip("\n").split("\t")[1:]]
            if line.startswith("!Sample_geo_accession"):
                accessions = cells
            elif line.startswith("!Sample_characteristics_ch1"):
                rows.append(cells)
    records = []
    for i, sample in enumerate(accessions):
        record = {"sample": sample}
        for row in rows:
            key, sep, value = row[i].partition(": ")
            if sep:
                if key in record:
                    raise ValueError(f"{matrix.name}: {sample} gives {key!r} twice")
                record[key] = value
        records.append(record)
    return pl.DataFrame(records, infer_schema_length=None)


_STURM_BASELINE_PARTS = ("1", "2")  # the lifespan parts; part 3 is the 3 % oxygen and contact study, part 4 metabolic


def _sturm(matrix: Path, supplementary: Path) -> tuple[pl.DataFrame, pl.DataFrame]:
    meta = _characteristics(matrix)
    treatment = pl.col("treatments")
    normal = pl.col("clinical_condition") == "Normal"
    baseline = normal & (treatment == "Control") & (pl.col("percent_oxygen") == "21")
    baseline = baseline & pl.col("study_part").is_in(_STURM_BASELINE_PARTS)
    before = _float("days_after_treatment_udays") <= 0
    samples = meta.select(
        "sample",
        pl.col("cell_line").alias("donor"),
        pl.lit("dermal fibroblast").alias("tissue"),
        pl.lit(None, dtype=pl.Float64).alias("age"),
        pl.when(baseline)
        .then(pl.lit("Baseline profiling"))
        .when(normal)
        .then(pl.lit("Study part ") + pl.col("study_part"))
        .otherwise(pl.lit("SURF1 mutation"))
        .alias("subexperiment"),
        _float("donor_age").alias("donor_age"),
        pl.col("Sex").alias("sex"),
        pl.when(pl.col("percent_oxygen") == "21")
        .then(pl.lit("Ambient"))
        .otherwise(pl.col("percent_oxygen") + "%")
        .alias("oxygen"),
        pl.when(treatment == "Control").then(pl.lit("control")).otherwise(treatment).alias("agent"),
        pl.lit("Full").alias("serum"),
        pl.lit("NA").alias("vector"),
        pl.col("passage").alias("subculture"),
        _float("population_doublings_utotal_divisions").alias("doublings"),
        _float("days_grown_udays").alias("days"),
        _float("days_after_treatment_udays").alias("days_treated"),
        pl.when(before)
        .then(pl.concat_str("cell_line", "passage", "days_grown_udays", separator="|"))
        .otherwise(pl.col("cell_line_group"))
        .alias("culture"),
        pl.col("sample_plate").alias("batch"),
        pl.concat_str("slide", "array", separator="_").alias("barcode"),
    )
    bad = samples.filter(pl.col("doublings").is_null() | pl.col("days").is_null())
    if bad.height:
        raise ValueError(
            f"{bad.height} samples without population doublings or days, e.g. {bad['sample'][:3].to_list()}"
        )
    betas = read_supplementary_betas(supplementary, suffix=" beta")
    rename = dict(zip(samples["barcode"].to_list(), samples["sample"].to_list(), strict=True))
    missing = [c for c in betas.columns[1:] if c not in rename]
    absent = sorted(set(rename) - set(betas.columns[1:]))
    if missing or absent:
        raise ValueError(
            f"{supplementary.name}: {len(missing)} columns name no sample (e.g. {missing[:3]}), "
            f"{len(absent)} samples have no column (e.g. {absent[:3]})"
        )
    return samples.drop("barcode"), betas.rename(rename)


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
        case "endicott":
            samples, betas = _endicott(lay.matrix)
        case "sturm":
            samples, betas = _sturm(lay.matrix, lay.supplementary)
    lay.samples.parent.mkdir(parents=True, exist_ok=True)
    samples.write_parquet(lay.samples)
    betas.with_columns(pl.exclude("cpg").cast(pl.Float32)).write_parquet(lay.betas)
    return lay
