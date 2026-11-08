"""Feature-by-sample matrices stored as wide polars frames.

Convention used across the package: a *matrix frame* has one string column
naming the feature (``"cpg"`` or ``"gene"``) and one numeric column per
sample. Models want the transpose as numpy, so helpers convert both ways.
"""

from __future__ import annotations

from collections.abc import Sequence
from dataclasses import dataclass
from pathlib import Path

import numpy as np
import polars as pl


@dataclass(frozen=True, slots=True)
class FeatureMatrix:
    """``values[i, j]`` is feature ``features[j]`` measured in sample ``samples[i]``."""

    values: np.ndarray
    samples: list[str]
    features: list[str]

    def __post_init__(self) -> None:
        if self.values.shape != (len(self.samples), len(self.features)):
            raise ValueError(f"values {self.values.shape} != ({len(self.samples)}, {len(self.features)})")

    @classmethod
    def from_frame(
        cls,
        frame: pl.DataFrame,
        *,
        feature_col: str | None = None,
        samples: Sequence[str] | None = None,
        features: Sequence[str] | None = None,
        max_missing: float = 0.0,
    ) -> FeatureMatrix:
        """Features missing in more than ``max_missing`` of samples are dropped; the rest are
        median-imputed per feature (array data has sporadic failed probes, and dropping every
        probe with one gap would remove ~13 % of Horvath's CpGs from GSE41037)."""
        feature_col = feature_col or frame.columns[0]
        sample_cols = list(samples) if samples is not None else [c for c in frame.columns if c != feature_col]
        sub = frame.select(feature_col, *sample_cols)
        if features is not None:
            wanted = pl.DataFrame({feature_col: list(features)})
            sub = wanted.join(sub, on=feature_col, how="left", maintain_order="left")
        values = sub.select(sample_cols).to_numpy().astype(float)  # features x samples, nulls -> nan
        missing = np.isnan(values).mean(axis=1)
        keep = missing <= max_missing
        values = values[keep]
        if np.isnan(values).any():
            medians = np.nanmedian(values, axis=1, keepdims=True)
            values = np.where(np.isnan(values), medians, values)
        names = [f for f, k in zip(sub[feature_col].to_list(), keep, strict=True) if k]
        return cls(values=values.T, samples=sample_cols, features=names)

    def to_frame(self, feature_col: str = "cpg") -> pl.DataFrame:
        data = {feature_col: self.features}
        data.update({s: self.values[i] for i, s in enumerate(self.samples)})
        return pl.DataFrame(data)

    def select_features(self, features: Sequence[str]) -> FeatureMatrix:
        index = {f: j for j, f in enumerate(self.features)}
        missing = [f for f in features if f not in index]
        if missing:
            raise KeyError(f"{len(missing)} features missing, e.g. {missing[:5]}")
        cols = np.fromiter((index[f] for f in features), dtype=np.int64, count=len(features))
        return FeatureMatrix(self.values[:, cols], self.samples, list(features))

    def select_samples(self, samples: Sequence[str]) -> FeatureMatrix:
        index = {s: i for i, s in enumerate(self.samples)}
        rows = np.fromiter((index[s] for s in samples), dtype=np.int64, count=len(samples))
        return FeatureMatrix(self.values[rows], list(samples), self.features)

    def sample(self, name: str) -> np.ndarray:
        return self.values[self.samples.index(name)]


def read_matrix(path: Path, **kwargs) -> FeatureMatrix:
    return FeatureMatrix.from_frame(pl.read_parquet(path), **kwargs)
