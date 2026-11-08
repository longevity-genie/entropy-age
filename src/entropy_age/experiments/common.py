"""Shared settings and data access for the figure scripts."""

from __future__ import annotations

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
        matrix, _ = dolphin
        return GroundState(
            "dolphin blood, 0.57 y", matrix.sample(DOLPHIN_YOUNG), matrix.sample(DOLPHIN_OLD), matrix.features
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
