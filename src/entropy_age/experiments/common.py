"""Shared settings and data access for the figure scripts."""

from __future__ import annotations

import sys
from dataclasses import dataclass, field
from enum import StrEnum
from pathlib import Path

import numpy as np
import polars as pl

from entropy_age.clock import ClockConfig
from entropy_age.data.matrix import FeatureMatrix
from entropy_age.data.upstream import UpstreamLayout
from entropy_age.paths import ProjectPaths

#: Samples named in the paper's methods (GSE41037, whole blood).
HUMAN_YOUNG = "GSM1007467"  # 16-year-old control: ground state
HUMAN_OLD = "GSM1007832"  # 88-year-old control: equilibrium for empirical maintenance rates
#: Pan-mammalian example data (bottlenose dolphin blood).
DOLPHIN_YOUNG = "204529320059_R04C02"  # 0.57 years
DOLPHIN_OLD = "203877620075_R03C02"  # 57.87 years
DOLPHIN_MAX_AGE = 67  # simulated steps = maximum age used upstream
#: CpGs missing in more than this fraction of samples are dropped, the rest median-imputed.
MAX_MISSING = 0.05


class ScaleName(StrEnum):
    QUICK = "quick"
    PAPER = "paper"


@dataclass(frozen=True, slots=True)
class Scale:
    """Sizes of a run. ``paper`` mirrors the publication, ``quick`` keeps a laptop cool."""

    name: ScaleName
    replicates: int
    n_features: int
    per_age: int
    n_cells: int
    feature_sizes: tuple[int, ...]
    clock: ClockConfig = field(default_factory=ClockConfig.fast)

    @classmethod
    def of(cls, name: ScaleName | str) -> Scale:
        match ScaleName(name):
            case ScaleName.QUICK:
                return cls(ScaleName.QUICK, 2, 1000, 3, 1000, (2, 8, 32, 128, 512, 2048), ClockConfig.fast())
            case ScaleName.PAPER:
                return cls(
                    ScaleName.PAPER,
                    3,
                    2000,
                    3,
                    1000,
                    tuple(2**i for i in range(14)),
                    ClockConfig.paper(),
                )
        raise ValueError(name)


@dataclass(frozen=True, slots=True)
class Context:
    paths: ProjectPaths
    scale: Scale
    seed: int = 42

    def rng(self, stream: int) -> np.random.Generator:
        """Independent, reproducible generator per figure panel."""
        return np.random.default_rng([self.seed, stream])

    def out_dir(self, figure: str) -> Path:
        path = self.paths.output / figure
        path.mkdir(parents=True, exist_ok=True)
        return path

    def figure_path(self, name: str) -> Path:
        return self.paths.figures / f"{name}.png"


# ----------------------------------------------------------------------------- datasets
def human_blood(paths: ProjectPaths, *, controls_only: bool = True) -> tuple[FeatureMatrix, pl.DataFrame] | None:
    """GSE41037 beta values (27k array) and sample table, if fetched; samples sorted by age."""
    betas_path = paths.interim / "geo" / "gse41037_betas.parquet"
    samples_path = UpstreamLayout(paths).gse41037_samples
    if not (betas_path.exists() and samples_path.exists()):
        return None
    samples = pl.read_parquet(samples_path)
    if controls_only:
        samples = samples.filter(pl.col("control"))
    betas = pl.read_parquet(betas_path)
    present = set(betas.columns)
    samples = samples.filter(pl.col("sample").is_in(list(present))).sort("age", "sample")
    matrix = FeatureMatrix.from_frame(
        betas, feature_col="cpg", samples=samples["sample"].to_list(), max_missing=MAX_MISSING
    )
    return matrix, samples


def human_batches(paths: ProjectPaths, sample_ids: list[str]) -> np.ndarray:
    """GSE41037 processing dataset (1, 2, 3) per sample, from the GEO series matrix.

    Batch and age are confounded in GSE41037 (dataset 1 is mostly under 50, dataset 3
    mostly over 50), so age slopes across the cohort need these as covariates.
    """
    table = pl.read_parquet(paths.interim / "geo" / "gse41037_samples.parquet").select("sample", "dataset")
    joined = pl.DataFrame({"sample": sample_ids}).join(table, on="sample", how="left", maintain_order="left")
    if joined["dataset"].null_count():
        raise ValueError("some samples have no processing dataset in the GEO sample sheet")
    return joined["dataset"].to_numpy()


def dummies(labels: np.ndarray) -> np.ndarray:
    """One-hot columns for every level but the first (for use next to an intercept)."""
    levels = np.unique(labels)
    if levels.size < 2:
        return np.zeros((len(labels), 0))
    return np.column_stack([(labels == level).astype(float) for level in levels[1:]])


def dolphin_blood(paths: ProjectPaths) -> tuple[FeatureMatrix, pl.DataFrame] | None:
    layout = UpstreamLayout(paths)
    if not (layout.dolphin_betas.exists() and layout.dolphin_samples.exists()):
        return None
    samples = pl.read_parquet(layout.dolphin_samples).sort("age", "sample")
    matrix = FeatureMatrix.from_frame(
        pl.read_parquet(layout.dolphin_betas),
        feature_col="cpg",
        samples=samples["sample"].to_list(),
        max_missing=MAX_MISSING,
    )
    return matrix, samples


@dataclass(frozen=True, slots=True)
class GroundState:
    """A biological (or random) starting methylome and, optionally, an old sample for equilibria."""

    label: str
    young: np.ndarray
    old: np.ndarray | None
    features: list[str]

    def subset(self, index: np.ndarray) -> GroundState:
        return GroundState(
            self.label,
            self.young[index],
            None if self.old is None else self.old[index],
            [self.features[i] for i in index],
        )


def ground_state(paths: ProjectPaths, rng: np.random.Generator, n_random: int = 27_000) -> GroundState:
    """Human 16-year-old blood if GSE41037 is available, then dolphin, then a uniform random ground."""
    human = human_blood(paths)
    if human is not None:
        matrix, _ = human
        return GroundState(
            "human blood, 16 y (GSE41037)", matrix.sample(HUMAN_YOUNG), matrix.sample(HUMAN_OLD), matrix.features
        )
    dolphin = dolphin_blood(paths)
    if dolphin is not None:
        print(
            "ground state: GSE41037 is not prepared (`entropy-age fetch geo GSE41037`); "
            "SUBSTITUTING the 0.57-year-old dolphin blood sample",
            file=sys.stderr,
        )
        matrix, _ = dolphin
        return GroundState(
            "dolphin blood, 0.57 y", matrix.sample(DOLPHIN_YOUNG), matrix.sample(DOLPHIN_OLD), matrix.features
        )
    print(
        "ground state: no methylation data prepared (`entropy-age fetch upstream`, `entropy-age fetch geo GSE41037`); "
        "SUBSTITUTING a uniform random ground state",
        file=sys.stderr,
    )
    return GroundState("uniform random", rng.uniform(0, 1, n_random), None, [f"f{i}" for i in range(n_random)])


#: Tables up to this many rows are also written as CSV, which is tracked in git.
CSV_MAX_ROWS = 1_000


def write_table(frame: pl.DataFrame, path: Path) -> Path:
    """Parquet always (gitignored); a CSV copy for summaries small enough to keep in git."""
    path.parent.mkdir(parents=True, exist_ok=True)
    frame.write_parquet(path.with_suffix(".parquet"))
    if frame.height <= CSV_MAX_ROWS:
        frame.write_csv(path.with_suffix(".csv"))
    return path.with_suffix(".parquet")


def scaled_intervals(point: pl.DataFrame, boot: pl.DataFrame, keys: list[str]) -> pl.DataFrame:
    """``point`` with ``gap_sd`` and its 95 % bootstrap interval, all in SDs of the full sample's metric.

    Each row of ``point`` and ``boot`` carries the unscaled ``matched_diff`` and the metric's spread ``sd``
    over features. A resample of people adds sampling noise to every feature's metric, which inflates its
    spread, so dividing each resample by its own spread shrinks every resampled gap toward zero and puts the
    point estimate outside its own interval. The interval is therefore taken on ``matched_diff`` and divided
    by the full sample's ``sd``.
    """
    q = boot.group_by(keys).agg(
        pl.col("matched_diff").quantile(0.025).alias("diff_low"),
        pl.col("matched_diff").quantile(0.975).alias("diff_high"),
    )
    return (
        point.join(q, on=keys)
        .with_columns(
            gap_sd=pl.col("matched_diff") / pl.col("sd"),
            ci_low=pl.col("diff_low") / pl.col("sd"),
            ci_high=pl.col("diff_high") / pl.col("sd"),
        )
        .drop("diff_low", "diff_high")
    )
