"""Data shipped with the upstream repository (Meyer-DH/StochasticAgingClock, MIT licence).

Files are fetched from a pinned commit (or copied from a local clone) into
``data/input/upstream`` and normalised to parquet in ``data/interim/upstream``.
"""

from __future__ import annotations

import csv
import shutil
from dataclasses import dataclass
from pathlib import Path

import polars as pl

from entropy_age.data.download import download
from entropy_age.paths import ProjectPaths

UPSTREAM_REPO = "Meyer-DH/StochasticAgingClock"
UPSTREAM_COMMIT = "caf31fd37fc4ce67016229ae1dc1e985ef9e4745"
RAW_URL = f"https://raw.githubusercontent.com/{UPSTREAM_REPO}/{UPSTREAM_COMMIT}/Data/{{name}}"

UPSTREAM_FILES: tuple[str, ...] = (
    "horvath_clock_coef.csv",
    "phenoage_clock_coef.csv",
    "Lin_Weights.csv",
    "VidalBralo_weights.csv",
    "WeidnerAge.csv",
    "BitAgeClock.csv",
    "GSE41037_meta.csv",
    "GSE41037_estCellTypes.csv",
    "SupplementaryTable1.csv",
    "mydata_GitHub_SampleAnnotation.csv",
    "mydata_GitHub_meth_betas.csv",
    "taxonomy_anage_meta.csv",
)


@dataclass(frozen=True, slots=True)
class UpstreamLayout:
    paths: ProjectPaths

    @property
    def input(self) -> Path:
        return self.paths.input / "upstream"

    @property
    def interim(self) -> Path:
        return self.paths.interim / "upstream"

    def clock(self, name: str) -> Path:
        return self.interim / "clocks" / f"{name}.parquet"

    @property
    def dolphin_betas(self) -> Path:
        return self.interim / "dolphin_betas.parquet"

    @property
    def dolphin_samples(self) -> Path:
        return self.interim / "dolphin_samples.parquet"

    @property
    def gse41037_samples(self) -> Path:
        return self.interim / "gse41037_samples.parquet"

    @property
    def gse41037_cell_types(self) -> Path:
        return self.interim / "gse41037_cell_types.parquet"

    @property
    def celegans_samples(self) -> Path:
        return self.interim / "celegans_rnaseq_samples.parquet"

    @property
    def taxonomy(self) -> Path:
        return self.interim / "anage_taxonomy.parquet"


def fetch_upstream(paths: ProjectPaths, *, from_dir: Path | None = None, overwrite: bool = False) -> list[Path]:
    """Download (or copy from a local ``Data/`` folder) every upstream data file."""
    layout = UpstreamLayout(paths)
    layout.input.mkdir(parents=True, exist_ok=True)
    fetched: list[Path] = []
    for name in UPSTREAM_FILES:
        target = layout.input / name
        if from_dir is not None:
            source = Path(from_dir) / name
            if overwrite or not target.exists():
                shutil.copyfile(source, target)
            fetched.append(target)
        else:
            fetched.append(download(RAW_URL.format(name=name), target, overwrite=overwrite))
    return fetched


# ----------------------------------------------------------------------------- clocks
#: name -> (file, separator, feature column, coefficient column, explicit intercept)
_CLOCK_FILES: dict[str, tuple[str, str, str, str, float | None]] = {
    "horvath": ("horvath_clock_coef.csv", ",", "CpGmarker", "CoefficientTraining", None),
    "phenoage": ("phenoage_clock_coef.csv", ",", "CpGmarker", "CoefficientTraining", None),
    "lin": ("Lin_Weights.csv", ",", "CpGmarker", "CoefficientTraining", None),
    "vidal_bralo": ("VidalBralo_weights.csv", ",", "CpGmarker", "CoefficientTraining", None),
    "weidner": ("WeidnerAge.csv", "\t", "CpGmarker", "CoefficientTraining", None),
    # BiT age (Meyer & Schumacher 2021): intercept hard-coded upstream, output in hours
    "bitage": ("BitAgeClock.csv", "\t", "Genes", "Coef", 103.55),
}


def _normalise_clock(
    raw_dir: Path, file: str, sep: str, feature: str, coef: str, intercept: float | None
) -> pl.DataFrame:
    frame = pl.read_csv(raw_dir / file, separator=sep, columns=[feature, coef], infer_schema_length=0)
    frame = frame.select(
        pl.col(feature).str.strip_chars().alias("feature"),
        pl.col(coef).cast(pl.Float64).alias("coef"),
    ).with_columns(
        pl.when(pl.col("feature").is_in(["Intercept", "(Intercept)"]))
        .then(pl.lit("(Intercept)"))
        .otherwise(pl.col("feature"))
        .alias("feature")
    )
    if intercept is not None:
        frame = pl.concat([pl.DataFrame({"feature": ["(Intercept)"], "coef": [intercept]}), frame])
    if frame.filter(pl.col("feature") == "(Intercept)").height == 0:
        frame = pl.concat([pl.DataFrame({"feature": ["(Intercept)"], "coef": [0.0]}), frame])
    return frame


# ----------------------------------------------------------------------------- tables
def _gse41037_samples(raw_dir: Path) -> pl.DataFrame:
    """The transposed GEO header kept upstream: source, sex, age and disease status per GSM."""
    with (raw_dir / "GSE41037_meta.csv").open(newline="") as handle:
        rows = list(csv.reader(handle, delimiter="\t"))
    header, source, sex, age, disease = rows[0], rows[1], rows[2], rows[3], rows[4]
    return pl.DataFrame(
        {
            "sample": header[1:],
            "source": source[1:],
            "sex": sex[1:],
            "age": [float(a) for a in age[1:]],
            "disease_status": disease[1:],
        }
    ).with_columns(control=pl.col("source").str.contains("control"))


def prepare_upstream(paths: ProjectPaths) -> UpstreamLayout:
    """Normalise every raw upstream file into parquet."""
    layout = UpstreamLayout(paths)
    raw = layout.input
    (layout.interim / "clocks").mkdir(parents=True, exist_ok=True)

    for name, spec in _CLOCK_FILES.items():
        _normalise_clock(raw, *spec).write_parquet(layout.clock(name))

    betas = pl.read_csv(raw / "mydata_GitHub_meth_betas.csv")
    betas = betas.rename({betas.columns[0]: "cpg"}).with_columns(pl.exclude("cpg").cast(pl.Float32))
    betas.write_parquet(layout.dolphin_betas)

    annotation = pl.read_csv(raw / "mydata_GitHub_SampleAnnotation.csv").select(
        pl.col("Basename").alias("sample"),
        pl.col("SpeciesLatinName").alias("species"),
        pl.col("SpeciesCommonName").alias("common_name"),
        pl.col("Age").cast(pl.Float64).alias("age"),
        pl.col("Female").cast(pl.Int8).alias("female"),
        pl.col("Tissue").alias("tissue"),
    )
    annotation.write_parquet(layout.dolphin_samples)

    _gse41037_samples(raw).write_parquet(layout.gse41037_samples)
    cell_types = pl.read_csv(raw / "GSE41037_estCellTypes.csv")
    cell_types.rename({cell_types.columns[0]: "sample"}).write_parquet(layout.gse41037_cell_types)

    pl.read_csv(raw / "SupplementaryTable1.csv", infer_schema_length=10_000).write_parquet(layout.celegans_samples)
    pl.read_csv(raw / "taxonomy_anage_meta.csv", infer_schema_length=10_000).write_parquet(layout.taxonomy)
    return layout


def clock_names() -> list[str]:
    return list(_CLOCK_FILES)
