"""Published linear clocks (coefficients from the upstream repository)."""

from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path

import numpy as np
import polars as pl

from entropy_age.clock import AgeTransform, HorvathAge, Identity
from entropy_age.data.matrix import FeatureMatrix

#: Output transform of each clock; BiT age predicts hours and is reported in days.
_TRANSFORMS: dict[str, AgeTransform] = {"horvath": HorvathAge()}
_OUTPUT_SCALE: dict[str, float] = {"bitage": 1 / 24}

CLOCK_LABELS: dict[str, str] = {
    "horvath": "Horvath 2013 (pan-tissue, 353 CpGs)",
    "phenoage": "Levine 2018 PhenoAge (513 CpGs)",
    "lin": "Lin 2016 (99 CpGs)",
    "vidal_bralo": "Vidal-Bralo 2016 (8 CpGs)",
    "weidner": "Weidner 2014 (3 CpGs)",
    "bitage": "BiT age, C. elegans (binarised transcriptome)",
}


@dataclass(frozen=True, slots=True)
class LinearClock:
    name: str
    intercept: float
    features: list[str]
    weights: np.ndarray
    transform: AgeTransform = field(default_factory=Identity)
    output_scale: float = 1.0

    @classmethod
    def read(cls, path: Path, name: str | None = None) -> LinearClock:
        name = name or path.stem
        frame = pl.read_parquet(path)
        intercept = float(frame.filter(pl.col("feature") == "(Intercept)")["coef"].sum())
        body = frame.filter(pl.col("feature") != "(Intercept)").sort("feature")
        return cls(
            name=name,
            intercept=intercept,
            features=body["feature"].to_list(),
            weights=body["coef"].to_numpy().astype(float),
            transform=_TRANSFORMS.get(name, Identity()),
            output_scale=_OUTPUT_SCALE.get(name, 1.0),
        )

    def predict(self, x: np.ndarray) -> np.ndarray:
        """``x`` is ``(n_samples, len(features))`` in the order of :attr:`features`."""
        raw = np.asarray(x, dtype=float) @ self.weights + self.intercept
        return self.transform.inverse(raw) * self.output_scale

    def predict_matrix(self, matrix: FeatureMatrix) -> np.ndarray:
        return self.predict(matrix.select_features(self.features).values)

    def coverage(self, available: set[str]) -> float:
        return sum(f in available for f in self.features) / len(self.features)
