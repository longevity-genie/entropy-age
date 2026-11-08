"""Extension: measure how fast each feature (CpG, gene) loses its ground state.

The paper treats stochasticity as one knob per dataset. To test whether
selection protects some genes or cell types from entropy, stochasticity has to
be a property of each feature in each cell type. This module estimates it
from cross-sectional data (one sample per individual, ages known):

* :func:`drift_metrics` - model-free per-feature statistics: drift slope,
  rank correlation with age, slope of the binary (Shannon) entropy, and the
  age slope of between-individual dispersion (residual spread, scaled by the
  binomial ``sqrt(m (1 - m))`` so bounded features are comparable).
* :func:`fit_relaxation` - fits the two-state Markov expectation
  ``m(t) = z + (m0 - z) exp(-lambda t)`` per feature by profiling ``lambda``
  on a grid (closed-form least squares for ``z, m0`` at each ``lambda``).
  ``lambda`` is the entropy rate and maps back to maintenance efficiencies.
* :func:`compare_groups` / :func:`paired_comparison` - group (gene set) and
  cell-type contrasts with permutation tests stratified by baseline level,
  because every metric depends on where a feature starts.
* :func:`regress_out` - the paper's cell-composition correction (per-feature
  linear model on cell fractions; residuals plus the feature mean).
"""

from __future__ import annotations

from collections.abc import Sequence

import numpy as np
import polars as pl
from scipy import stats

from entropy_age.simulate.methylation import MaintenanceRates

_EPS = 1e-9


def binary_entropy(p: np.ndarray) -> np.ndarray:
    """Shannon entropy (bits) of a fraction: 0 at 0 or 1, maximal (1 bit) at 0.5."""
    p = np.clip(np.asarray(p, dtype=float), _EPS, 1 - _EPS)
    return -(p * np.log2(p) + (1 - p) * np.log2(1 - p))


def sample_entropy(x: np.ndarray) -> np.ndarray:
    """Mean binary entropy across features for every sample (Hannum-style methylome entropy)."""
    return binary_entropy(x).mean(axis=1)


def _slopes(t: np.ndarray, y: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """OLS slope and intercept of every column of ``y`` on ``t``."""
    tc = t - t.mean()
    slope = tc @ (y - y.mean(axis=0)) / (tc @ tc)
    return slope, y.mean(axis=0) - slope * t.mean()


def drift_metrics(
    x: np.ndarray,
    ages: Sequence[float] | np.ndarray,
    features: Sequence[str] | None = None,
    *,
    young_quantile: float = 0.2,
) -> pl.DataFrame:
    """Per-feature drift statistics from an ``(n_samples, n_features)`` matrix of values in [0, 1]."""
    x = np.asarray(x, dtype=float)
    t = np.asarray(ages, dtype=float)
    if x.shape[0] != t.size:
        raise ValueError("one age per sample is required")
    names = list(features) if features is not None else [f"f{j}" for j in range(x.shape[1])]

    slope, intercept = _slopes(t, x)
    ranked_t = stats.rankdata(t)
    ranked_x = stats.rankdata(x, axis=0)
    rc = ranked_t - ranked_t.mean()
    rx = ranked_x - ranked_x.mean(axis=0)
    spearman = (rc @ rx) / np.sqrt((rc @ rc) * (rx**2).sum(axis=0) + _EPS)

    entropy_slope, _ = _slopes(t, binary_entropy(x))

    residual = x - (intercept + np.outer(t, slope))
    level = np.clip(x.mean(axis=0), 0.01, 0.99)
    scaled_abs = np.abs(residual) / np.sqrt(level * (1 - level))
    dispersion_slope, _ = _slopes(t, scaled_abs)

    young = t <= np.quantile(t, young_quantile)
    old = t >= np.quantile(t, 1 - young_quantile)
    var_young = residual[young].var(axis=0, ddof=1)
    var_old = residual[old].var(axis=0, ddof=1)

    return pl.DataFrame(
        {
            "feature": names,
            "baseline": x[young].mean(axis=0),
            "late": x[old].mean(axis=0),
            "slope": slope,
            "spearman": spearman,
            "entropy_slope": entropy_slope,
            "dispersion_slope": dispersion_slope,
            "log2_dispersion_ratio": np.log2((var_old + _EPS) / (var_young + _EPS)),
        }
    )


def fit_relaxation(
    x: np.ndarray,
    ages: Sequence[float] | np.ndarray,
    features: Sequence[str] | None = None,
    *,
    t0: float | None = None,
    rates: np.ndarray | None = None,
) -> pl.DataFrame:
    """Fit ``m(t) = z + (m0 - z) exp(-lambda (t - t0))`` to every feature.

    ``lambda`` is profiled over ``rates`` (per age unit, default 1e-4..1 on a log
    grid); ``z`` and ``m0`` are solved in closed form for each candidate.
    Estimates at the lower edge of the grid mean the data only support a
    straight line (``z`` is then an extrapolation); ``z`` outside [0, 1] is flagged.

    ``window_change`` is the fitted change between the youngest and oldest age,
    the robust summary of how much a feature drifted. Read ``rate`` together with
    ``equilibrium``: selection protects a feature by holding it at a nearby
    equilibrium (mutation-selection balance), which it reaches quickly, so a
    protected feature can show a *higher* rate and a *smaller* change.
    """
    x = np.asarray(x, dtype=float)
    t_abs = np.asarray(ages, dtype=float)
    t = t_abs - (t_abs.min() if t0 is None else t0)
    grid = np.geomspace(1e-4, 1.0, 64) if rates is None else np.asarray(rates, dtype=float)
    names = list(features) if features is not None else [f"f{j}" for j in range(x.shape[1])]

    n_features = x.shape[1]
    best_sse = np.full(n_features, np.inf)
    best = np.zeros((3, n_features))  # rate, z, m0
    for rate in grid:
        decay = np.exp(-rate * t)
        design = np.column_stack([1.0 - decay, decay])
        coef, *_ = np.linalg.lstsq(design, x, rcond=None)
        sse = ((x - design @ coef) ** 2).sum(axis=0)
        better = sse < best_sse
        best_sse = np.where(better, sse, best_sse)
        best[0, better] = rate
        best[1:, better] = coef[:, better]

    total = ((x - x.mean(axis=0)) ** 2).sum(axis=0)
    rate, z, m0 = best
    in_range = (z >= 0) & (z <= 1)
    maintenance = MaintenanceRates.from_rate_and_equilibrium(rate, np.clip(z, 0, 1))
    window_change = (z - m0) * (np.exp(-rate * t.min()) - np.exp(-rate * t.max()))
    return pl.DataFrame(
        {
            "feature": names,
            "rate": rate,
            "half_life": np.log(2) / rate,
            "equilibrium": z,
            "start": m0,
            "window_change": window_change,
            "r2": 1 - best_sse / np.where(total > 0, total, np.nan),
            "rate_at_lower_bound": rate <= grid[0],
            "rate_at_upper_bound": rate >= grid[-1],
            "equilibrium_in_range": in_range,
            "em": np.where(in_range, maintenance.em, np.nan),
            "ed": np.where(in_range, maintenance.ed, np.nan),
        }
    )


def baseline_strata(values: np.ndarray, n_bins: int = 10) -> np.ndarray:
    """Quantile bin index (0..n_bins-1) of every value, used to match features by start level."""
    edges = np.unique(np.quantile(values, np.linspace(0, 1, n_bins + 1)))
    return np.clip(np.searchsorted(edges, values, side="right") - 1, 0, max(len(edges) - 2, 0))


def _matched_difference(values: np.ndarray, labels: np.ndarray, strata: np.ndarray) -> float:
    """In-group minus rest, averaged over strata weighted by in-group size (a CMH-style contrast)."""
    total, weight = 0.0, 0
    for s in np.unique(strata):
        here = strata == s
        inside, outside = values[here & labels], values[here & ~labels]
        total += inside.size * (inside.mean() - outside.mean())
        weight += inside.size
    return total / weight


def _stratified_permutation(
    values: np.ndarray, labels: np.ndarray, strata: np.ndarray, n_permutations: int, rng: np.random.Generator
) -> tuple[float, float]:
    """Matched difference and its two-sided p with labels shuffled within strata.

    Callers pass only strata that contain both groups; strata holding one group
    carry no within-level information and would make the test degenerate.
    """
    observed = _matched_difference(values, labels, strata)
    groups = [np.flatnonzero(strata == s) for s in np.unique(strata)]
    permuted = labels.copy()
    exceed = 0
    for _ in range(n_permutations):
        for idx in groups:
            permuted[idx] = rng.permutation(labels[idx])
        exceed += abs(_matched_difference(values, permuted, strata)) >= abs(observed) - 1e-12
    return observed, (exceed + 1) / (n_permutations + 1)


def cliffs_delta(a: np.ndarray, b: np.ndarray) -> float:
    """P(a > b) - P(a < b), computed from the Mann-Whitney U statistic."""
    if a.size == 0 or b.size == 0:
        return float("nan")
    u = stats.mannwhitneyu(a, b, alternative="two-sided").statistic
    return float(2 * u / (a.size * b.size) - 1)


def compare_groups(
    metrics: pl.DataFrame,
    membership: pl.DataFrame,
    *,
    metric: str,
    stratify_by: str | None = "baseline",
    n_bins: int = 10,
    n_permutations: int = 2000,
    seed: int = 0,
) -> pl.DataFrame:
    """Each group (gene set, region class) against all other features.

    ``membership`` has columns ``feature`` and ``group``; a feature may be in
    several groups. With ``stratify_by``, features are binned into ``n_bins``
    quantiles of that column, the contrast is the in-group minus rest difference
    within bins (weighted by group size) and labels are permuted within bins.
    Bins holding only one group are dropped and ``n_matched`` says how many
    group members remained.

    Choose the stratifier with care: it should be measured independently of the
    drift being compared (for example a separate young or fetal reference). A
    baseline taken from the cohort itself already contains early drift, so
    fast-drifting features sit in different bins from slow ones and matching
    removes part of the very effect being tested.
    """
    rng = np.random.default_rng(seed)
    base = metrics.drop_nulls(subset=[metric]).filter(pl.col(metric).is_finite())
    values = base[metric].to_numpy()
    feature_ids = base["feature"].to_numpy()
    strata = baseline_strata(base[stratify_by].to_numpy(), n_bins) if stratify_by else np.zeros(values.size, dtype=int)
    rows = []
    for group, members in membership.group_by("group", maintain_order=True):
        label = group[0] if isinstance(group, tuple) else group
        in_group = np.isin(feature_ids, members["feature"].to_numpy())
        if in_group.sum() < 2 or (~in_group).sum() < 2:
            continue
        mixed = [s for s in np.unique(strata) if in_group[strata == s].any() and (~in_group[strata == s]).any()]
        matched = np.isin(strata, mixed)
        if in_group[matched].sum() >= 2 and (~in_group[matched]).sum() >= 2:
            difference, p_value = _stratified_permutation(
                values[matched], in_group[matched], strata[matched], n_permutations, rng
            )
        else:
            difference, p_value = float("nan"), float("nan")
        rows.append(
            {
                "group": label,
                "n": int(in_group.sum()),
                "n_matched": int(in_group[matched].sum()),
                "median": float(np.median(values[in_group])),
                "median_rest": float(np.median(values[~in_group])),
                "mean_diff": float(values[in_group].mean() - values[~in_group].mean()),
                "matched_diff": float(difference),
                "cliffs_delta": cliffs_delta(values[in_group], values[~in_group]),
                "p_permutation": float(p_value),
            }
        )
    return pl.DataFrame(rows)


def paired_comparison(a: pl.DataFrame, b: pl.DataFrame, *, metric: str, names: tuple[str, str] = ("a", "b")) -> dict:
    """Same features measured in two cell types: Wilcoxon signed-rank on the per-feature difference."""
    joined = a.select("feature", pl.col(metric).alias("x")).join(
        b.select("feature", pl.col(metric).alias("y")), on="feature"
    )
    joined = joined.filter(pl.col("x").is_finite() & pl.col("y").is_finite())
    diff = (joined["x"] - joined["y"]).to_numpy()
    test = stats.wilcoxon(diff) if diff.size >= 10 and np.any(diff != 0) else None
    return {
        "n_features": int(diff.size),
        f"median_{names[0]}": float(joined["x"].median()) if diff.size else float("nan"),
        f"median_{names[1]}": float(joined["y"].median()) if diff.size else float("nan"),
        "median_difference": float(np.median(diff)) if diff.size else float("nan"),
        f"fraction_{names[0]}_greater": float(np.mean(diff > 0)) if diff.size else float("nan"),
        "p_wilcoxon": float(test.pvalue) if test is not None else float("nan"),
    }


def regress_out(x: np.ndarray, covariates: np.ndarray) -> np.ndarray:
    """Remove the part of every feature explained by ``covariates`` (e.g. cell fractions).

    Port of the paper's R correction: ``lm(x ~ B + NK + CD4T + CD8T + Mono +
    Neutro + Eosino)`` per CpG, residuals plus the CpG mean.
    """
    x = np.asarray(x, dtype=float)
    design = np.column_stack([np.ones(x.shape[0]), np.asarray(covariates, dtype=float)])
    coef, *_ = np.linalg.lstsq(design, x, rcond=None)
    return x - design @ coef + x.mean(axis=0)
