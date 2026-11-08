"""Training and evaluating linear aging clocks on simulated or measured data."""

from __future__ import annotations

import warnings
from collections.abc import Sequence
from dataclasses import dataclass, field
from typing import Protocol

import numpy as np
import polars as pl
from scipy import stats
from sklearn.exceptions import ConvergenceWarning
from sklearn.linear_model import ElasticNet, ElasticNetCV
from sklearn.metrics import r2_score

PAPER_L1_RATIOS: tuple[float, ...] = (0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9)


# --------------------------------------------------------------------------- age transforms
class AgeTransform(Protocol):
    def forward(self, age: np.ndarray) -> np.ndarray: ...
    def inverse(self, value: np.ndarray) -> np.ndarray: ...


@dataclass(frozen=True, slots=True)
class Identity:
    def forward(self, age: np.ndarray) -> np.ndarray:
        return np.asarray(age, dtype=float)

    def inverse(self, value: np.ndarray) -> np.ndarray:
        return np.asarray(value, dtype=float)


@dataclass(frozen=True, slots=True)
class LinearRescale:
    """Map ``[in_min, in_max]`` onto ``[out_min, out_max]``.

    The human stochastic data-based clock (Fig. 5) rescales simulated steps
    1..73 to ``[-120, 280]`` (min-max, ``* 400 - 120``); the worm clock doubles them.
    """

    in_min: float
    in_max: float
    out_min: float = -120.0
    out_max: float = 280.0

    def forward(self, age: np.ndarray) -> np.ndarray:
        unit = (np.asarray(age, dtype=float) - self.in_min) / (self.in_max - self.in_min)
        return unit * (self.out_max - self.out_min) + self.out_min

    def inverse(self, value: np.ndarray) -> np.ndarray:
        """Back to simulated-step units (rarely needed: predictions are read on the output scale)."""
        unit = (np.asarray(value, dtype=float) - self.out_min) / (self.out_max - self.out_min)
        return unit * (self.in_max - self.in_min) + self.in_min


@dataclass(frozen=True, slots=True)
class LogLogRelativeAge:
    """Pan-mammalian transform: train on ``-log(-log(age / max_age))``, report relative age.

    ``inverse`` returns ``exp(-exp(-y))``, the relative age in ``(0, 1)`` used by
    Lu et al. (2023) and by the upstream pan-mammalian stochastic clocks.
    """

    max_age: float

    def forward(self, age: np.ndarray) -> np.ndarray:
        rel = np.clip(np.asarray(age, dtype=float) / self.max_age, 1e-9, 1 - 1e-9)
        return -np.log(-np.log(rel))

    def inverse(self, value: np.ndarray) -> np.ndarray:
        return np.exp(-np.exp(-np.asarray(value, dtype=float)))


@dataclass(frozen=True, slots=True)
class HorvathAge:
    """Horvath (2013) log-linear age transform with adult age 20."""

    adult_age: float = 20.0

    def forward(self, age: np.ndarray) -> np.ndarray:
        age = np.asarray(age, dtype=float)
        a = self.adult_age
        return np.where(age <= a, np.log(age + 1) - np.log(a + 1), (age - a) / (a + 1))

    def inverse(self, value: np.ndarray) -> np.ndarray:
        value = np.asarray(value, dtype=float)
        a = self.adult_age
        return np.where(value < 0, (1 + a) * np.exp(value) - 1, (1 + a) * value + a)


# --------------------------------------------------------------------------- metrics
@dataclass(frozen=True, slots=True)
class ClockMetrics:
    n: int
    pearson_r: float
    pearson_p: float
    spearman_r: float
    spearman_p: float
    r2: float
    median_abs_error: float

    def as_dict(self) -> dict[str, float | int]:
        return {f: getattr(self, f) for f in self.__slots__}


def evaluate(y_true: Sequence[float] | np.ndarray, y_pred: Sequence[float] | np.ndarray) -> ClockMetrics:
    """Pearson, Spearman, R^2 and median absolute error (the statistics reported upstream)."""
    t = np.asarray(y_true, dtype=float)
    p = np.asarray(y_pred, dtype=float)
    keep = np.isfinite(t) & np.isfinite(p)
    t, p = t[keep], p[keep]
    if t.size < 3 or np.ptp(p) == 0 or np.ptp(t) == 0:
        nan = float("nan")
        return ClockMetrics(int(t.size), nan, nan, nan, nan, nan, nan)
    pr = stats.pearsonr(t, p)
    sr = stats.spearmanr(t, p)
    return ClockMetrics(
        n=int(t.size),
        pearson_r=float(pr.statistic),
        pearson_p=float(pr.pvalue),
        spearman_r=float(sr.statistic),
        spearman_p=float(sr.pvalue),
        r2=float(r2_score(t, p)),
        median_abs_error=float(np.median(np.abs(t - p))),
    )


# --------------------------------------------------------------------------- elastic net clocks
@dataclass(frozen=True, slots=True)
class ClockConfig:
    """Elastic net settings. ``alphas=None`` searches the regularisation path by CV."""

    l1_ratio: tuple[float, ...] = PAPER_L1_RATIOS
    alphas: tuple[float, ...] | None = None
    n_alphas: int = 100
    cv: int = 5
    max_iter: int = 1000
    tol: float = 1e-4
    n_jobs: int | None = 1
    random_state: int = 0

    @classmethod
    def paper(cls) -> ClockConfig:
        """Upstream ``ElasticNetCV(l1_ratio=[0.1..0.9])`` used for the simulation figures."""
        return cls()

    @classmethod
    def fixed_alpha(cls, alpha: float = 1.0, l1_ratio: Sequence[float] = PAPER_L1_RATIOS) -> ClockConfig:
        """Upstream real-data clocks: ``alphas=[1]``, CV only over ``l1_ratio``."""
        return cls(l1_ratio=tuple(l1_ratio), alphas=(alpha,))

    @classmethod
    def fast(cls) -> ClockConfig:
        """Cheaper search for sweeps and tests: three l1 ratios, 30 alphas, 3 folds."""
        return cls(l1_ratio=(0.1, 0.5, 0.9), n_alphas=30, cv=3)


@dataclass(slots=True)
class StochasticClock:
    """A fitted linear clock: ``transform.inverse(X @ coef + intercept)``."""

    coef: np.ndarray
    intercept: float
    transform: AgeTransform = field(default_factory=Identity)
    features: list[str] | None = None
    alpha: float | None = None
    l1_ratio: float | None = None

    def predict(self, x: np.ndarray) -> np.ndarray:
        return self.transform.inverse(np.asarray(x, dtype=float) @ self.coef + self.intercept)

    def predict_raw(self, x: np.ndarray) -> np.ndarray:
        return np.asarray(x, dtype=float) @ self.coef + self.intercept

    @property
    def n_nonzero(self) -> int:
        return int(np.count_nonzero(self.coef))

    def coefficients(self) -> pl.DataFrame:
        names = self.features if self.features is not None else [f"f{i}" for i in range(self.coef.size)]
        return pl.DataFrame({"feature": ["(Intercept)", *names], "coef": [self.intercept, *self.coef.tolist()]})


def train_clock(
    x: np.ndarray,
    ages: Sequence[float] | np.ndarray,
    *,
    config: ClockConfig | None = None,
    transform: AgeTransform | None = None,
    features: Sequence[str] | None = None,
) -> StochasticClock:
    """Fit an elastic net on ``transform.forward(ages)``; convergence warnings are silenced as upstream."""
    config = config or ClockConfig()
    transform = transform or Identity()
    y = transform.forward(np.asarray(ages, dtype=float))
    x = np.asarray(x, dtype=float)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", category=ConvergenceWarning)
        if config.alphas is not None and len(config.alphas) == 1 and len(config.l1_ratio) == 1:
            model = ElasticNet(
                alpha=config.alphas[0],
                l1_ratio=config.l1_ratio[0],
                max_iter=config.max_iter,
                tol=config.tol,
                random_state=config.random_state,
            ).fit(x, y)
            alpha, l1 = config.alphas[0], config.l1_ratio[0]
        else:
            model = ElasticNetCV(
                l1_ratio=list(config.l1_ratio),
                alphas=list(config.alphas) if config.alphas is not None else config.n_alphas,
                cv=config.cv,
                max_iter=config.max_iter,
                tol=config.tol,
                n_jobs=config.n_jobs,
                random_state=config.random_state,
            ).fit(x, y)
            alpha, l1 = float(model.alpha_), float(model.l1_ratio_)
    return StochasticClock(
        coef=np.asarray(model.coef_, dtype=float),
        intercept=float(model.intercept_),
        transform=transform,
        features=list(features) if features is not None else None,
        alpha=alpha,
        l1_ratio=l1,
    )
