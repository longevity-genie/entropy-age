"""Extension: measure how fast each feature (CpG, gene) loses its ground state.

The paper treats stochasticity as one knob per dataset. To test whether
selection protects some genes or cell types from entropy, stochasticity has to
be a property of each feature in each cell type. This module estimates it
from cross-sectional data (one sample per individual, ages known):

* :func:`drift_metrics` - model-free per-feature statistics: drift slope,
  rank correlation with age, slope of the binary (Shannon) entropy, and the
  age slope of between-individual dispersion (residual spread, scaled by the
  binomial ``sqrt(m (1 - m))`` so bounded features are comparable).
* :func:`expression_drift` - the same idea for unbounded values such as log CPM: the age slope
  of the mean and of the absolute residual (between-individual spread).
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


def relative_entropy(x: np.ndarray, reference: np.ndarray) -> np.ndarray:
    """``KL(Bernoulli(x) || Bernoulli(m))`` in bits: how far ``x`` has moved from the reference level ``m``.

    Unlike binary entropy, which measures distance from 0.5 and so counts a move toward 0 or 1 as a gain in order,
    it counts a move up or down alike and needs no assumption about where a site drifts to; a move away from an
    extreme level weighs more, as binomial sampling says it should. The reference is clipped to [0.01, 0.99].
    """
    p = np.clip(np.asarray(x, dtype=float), _EPS, 1 - _EPS)
    m = np.clip(np.asarray(reference, dtype=float), 0.01, 0.99)  # as the per-person deviations clip the level
    return p * np.log2(p / m) + (1 - p) * np.log2((1 - p) / (1 - m))


def sample_entropy(x: np.ndarray) -> np.ndarray:
    """Mean binary entropy across features for every sample (Hannum-style methylome entropy)."""
    return binary_entropy(x).mean(axis=1)


def _age_fit(t: np.ndarray, y: np.ndarray, covariates: np.ndarray | None) -> tuple[np.ndarray, np.ndarray]:
    """Age coefficient of every column of ``y`` (adjusted for ``covariates``) and the fitted values."""
    columns = [np.ones_like(t), t] if covariates is None else [np.ones_like(t), t, covariates]
    design = np.column_stack(columns)
    coef, *_ = np.linalg.lstsq(design, y, rcond=None)
    return coef[1], design @ coef


def _partial_spearman(t: np.ndarray, x: np.ndarray, covariates: np.ndarray | None) -> np.ndarray:
    """Rank correlation of every column with ``t``, partialling out ``covariates`` when given."""
    ranked_t = stats.rankdata(t)
    ranked_x = stats.rankdata(x, axis=0)
    if covariates is not None:
        base = np.column_stack([np.ones_like(t), covariates])
        ranked_t = ranked_t - base @ np.linalg.lstsq(base, ranked_t, rcond=None)[0]
        ranked_x = ranked_x - base @ np.linalg.lstsq(base, ranked_x, rcond=None)[0]
    rc = ranked_t - ranked_t.mean()
    rx = ranked_x - ranked_x.mean(axis=0)
    return (rc @ rx) / np.sqrt((rc @ rc) * (rx**2).sum(axis=0) + _EPS)


def dispersion_levels(
    x: np.ndarray, features: Sequence[str] | None = None, *, covariates: np.ndarray | None = None
) -> pl.DataFrame:
    """Between-individual spread of every feature, for a cohort of similar age.

    No age slope is needed: in people of about the same age, how far individuals have diverged
    at a site is the drift it has accumulated (plus genetics and measurement noise, which a
    comparison of the same CpGs in the same donors' tissues holds fixed). ``dispersion`` is the
    residual standard deviation after regressing out ``covariates``, divided by the binomial
    ``sqrt(m (1 - m))`` so that bounded features starting at different levels are comparable.
    """
    x = np.asarray(x, dtype=float)
    names = list(features) if features is not None else [f"f{j}" for j in range(x.shape[1])]
    ones = np.ones((x.shape[0], 1))
    design = ones if covariates is None else np.column_stack([ones, np.asarray(covariates, dtype=float)])
    coef, *_ = np.linalg.lstsq(design, x, rcond=None)
    residual = x - design @ coef
    dof = max(x.shape[0] - design.shape[1], 1)
    level = x.mean(axis=0)
    scale = np.sqrt(np.clip(level, 0.01, 0.99) * (1 - np.clip(level, 0.01, 0.99)))
    return pl.DataFrame(
        {
            "feature": names,
            "baseline": level,
            "residual_sd": np.sqrt((residual**2).sum(axis=0) / dof),
            "dispersion": np.sqrt((residual**2).sum(axis=0) / dof) / scale,
        }
    )


def drift_metrics(
    x: np.ndarray,
    ages: Sequence[float] | np.ndarray,
    features: Sequence[str] | None = None,
    *,
    covariates: np.ndarray | None = None,
    young_quantile: float = 0.2,
    spearman: bool = True,
) -> pl.DataFrame:
    """Per-feature drift statistics from an ``(n_samples, n_features)`` matrix of values in [0, 1].

    ``covariates`` (``(n_samples, k)``, e.g. processing-batch dummies or cell fractions)
    turn every age slope into a partial slope, so only age variation *within* the
    covariates counts. Use it whenever batch and age are confounded: across a
    cohort whose old donors were run on other plates than its young ones, an
    unadjusted slope partly measures the batch. ``spearman`` is then a partial
    rank correlation. ``log2_dispersion_ratio`` compares the oldest with the
    youngest samples and is not adjusted.
    """
    x = np.asarray(x, dtype=float)
    t = np.asarray(ages, dtype=float)
    if x.shape[0] != t.size:
        raise ValueError("one age per sample is required")
    if covariates is not None:
        covariates = np.asarray(covariates, dtype=float).reshape(t.size, -1)
    names = list(features) if features is not None else [f"f{j}" for j in range(x.shape[1])]

    slope, fitted = _age_fit(t, x, covariates)
    rank_corr = _partial_spearman(t, x, covariates) if spearman else np.full(x.shape[1], np.nan)

    entropy_slope, _ = _age_fit(t, binary_entropy(x), covariates)
    young = t <= np.quantile(t, young_quantile)
    kl_slope, _ = _age_fit(t, relative_entropy(x, x[young].mean(axis=0)), covariates)

    residual = x - fitted
    level = np.clip(x.mean(axis=0), 0.01, 0.99)
    scaled_abs = np.abs(residual) / np.sqrt(level * (1 - level))
    dispersion_slope, _ = _age_fit(t, scaled_abs, covariates)

    old = t >= np.quantile(t, 1 - young_quantile)
    var_young = residual[young].var(axis=0, ddof=1)
    var_old = residual[old].var(axis=0, ddof=1)

    return pl.DataFrame(
        {
            "feature": names,
            "baseline": x[young].mean(axis=0),
            "late": x[old].mean(axis=0),
            "slope": slope,
            "abs_slope": np.abs(slope),
            "spearman": rank_corr,
            "entropy_slope": entropy_slope,
            "kl_slope": kl_slope,
            "dispersion_slope": dispersion_slope,
            "log2_dispersion_ratio": np.log2((var_old + _EPS) / (var_young + _EPS)),
        }
    )


def expression_drift(
    x: np.ndarray,
    ages: Sequence[float] | np.ndarray,
    features: Sequence[str] | None = None,
    *,
    covariates: np.ndarray | None = None,
) -> pl.DataFrame:
    """Per-gene drift of an ``(n_samples, n_genes)`` matrix on a log scale (log2 CPM).

    ``dispersion_slope`` is the age slope of the absolute residual of ``x ~ age + covariates``,
    adjusted for the same covariates: how fast individuals spread apart at that gene. On a log scale
    no binomial rescaling applies, but RNA spread still depends strongly on the mean, so compare
    genes only within strata of ``baseline`` (the gene's mean over all samples, which age barely
    moves). ``abs_slope`` is the absolute age slope of the mean and ``residual_sd`` the residual
    spread over all ages.
    """
    x = np.asarray(x, dtype=float)
    t = np.asarray(ages, dtype=float)
    if x.shape[0] != t.size:
        raise ValueError("one age per sample is required")
    if covariates is not None:
        covariates = np.asarray(covariates, dtype=float).reshape(t.size, -1)
    names = list(features) if features is not None else [f"f{j}" for j in range(x.shape[1])]
    slope, fitted = _age_fit(t, x, covariates)
    residual = x - fitted
    dispersion_slope, _ = _age_fit(t, np.abs(residual), covariates)
    dof = max(t.size - 2 - (0 if covariates is None else covariates.shape[1]), 1)
    return pl.DataFrame(
        {
            "feature": names,
            "baseline": x.mean(axis=0),
            "slope": slope,
            "abs_slope": np.abs(slope),
            "residual_sd": np.sqrt((residual**2).sum(axis=0) / dof),
            "dispersion_slope": dispersion_slope,
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


def _matched_difference(values: np.ndarray, labels: np.ndarray, codes: np.ndarray, n_strata: int) -> float:
    """In-group minus rest, averaged over strata weighted by in-group size (a CMH-style contrast).

    ``codes`` are stratum indices ``0..n_strata-1``; every stratum must hold both groups.
    """
    weight = labels.astype(float)
    n_in = np.bincount(codes, weights=weight, minlength=n_strata)
    n_all = np.bincount(codes, minlength=n_strata).astype(float)
    sum_in = np.bincount(codes, weights=values * weight, minlength=n_strata)
    sum_all = np.bincount(codes, weights=values, minlength=n_strata)
    difference = sum_in / n_in - (sum_all - sum_in) / (n_all - n_in)
    return float((n_in * difference).sum() / n_in.sum())


def _stratified_permutation(
    values: np.ndarray, labels: np.ndarray, strata: np.ndarray, n_permutations: int, rng: np.random.Generator
) -> tuple[float, float]:
    """Matched difference and its two-sided p with labels shuffled within strata.

    Callers pass only strata that contain both groups; strata holding one group
    carry no within-level information and would make the test degenerate.
    """
    _, codes = np.unique(strata, return_inverse=True)
    n_strata = int(codes.max()) + 1
    observed = _matched_difference(values, labels, codes, n_strata)
    if n_permutations == 0:
        return observed, float("nan")  # the matched difference alone (bootstrap resamples)
    exceed = 0
    for _ in range(n_permutations):
        permuted = _permute_within_strata(labels, codes, rng)
        exceed += abs(_matched_difference(values, permuted, codes, n_strata)) >= abs(observed) - 1e-12
    return observed, (exceed + 1) / (n_permutations + 1)


def _permute_within_strata(labels: np.ndarray, codes: np.ndarray, rng: np.random.Generator) -> np.ndarray:
    """Shuffle ``labels`` only among positions that share a stratum code."""
    grouped = np.argsort(codes, kind="stable")  # positions, grouped by stratum
    shuffled = np.lexsort((rng.random(codes.size), codes))  # also grouped by stratum, random within
    permuted = np.empty_like(labels)
    permuted[grouped] = labels[shuffled]
    return permuted


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

    A numeric float column is binned into quantiles; an integer, string or
    boolean column is used as given (e.g. CpG-island status x baseline decile).

    Choose the stratifier with care: it should be measured independently of the
    drift being compared (for example a separate young or fetal reference). A
    baseline taken from the cohort itself already contains early drift, so
    fast-drifting features sit in different bins from slow ones and matching
    removes part of the very effect being tested.
    """
    rng = np.random.default_rng(seed)
    base = metrics.drop_nulls(subset=[metric]).filter(pl.col(metric).is_finite())
    values = base[metric].to_numpy()
    if stratify_by is None:
        strata = np.zeros(values.size, dtype=int)
    elif base.schema[stratify_by].is_float():
        strata = baseline_strata(base[stratify_by].to_numpy(), n_bins)
    else:
        strata = np.unique(base[stratify_by].cast(pl.String).to_numpy(), return_inverse=True)[1]
    rows = []
    for group, members in membership.group_by("group", maintain_order=True):
        label = group[0] if isinstance(group, tuple) else group
        # a hash lookup: np.isin on string arrays took seconds per group at 10^5 features
        in_group = base["feature"].is_in(members["feature"].unique().implode()).to_numpy()
        if in_group.sum() < 2 or (~in_group).sum() < 2:
            continue
        n_in = np.bincount(strata, weights=in_group.astype(float))
        n_all = np.bincount(strata)
        mixed = np.flatnonzero((n_in > 0) & (n_in < n_all))
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


def composition_components(
    x: np.ndarray, k: int, *, covariates: np.ndarray | None = None
) -> tuple[np.ndarray, np.ndarray]:
    """Leading ``k`` principal components of ``x`` after regressing out ``covariates``.

    A reference-free proxy for cell composition, used when no sorted-cell reference is at hand:
    when samples differ in their mix of cell types, thousands of cell-type-specific CpGs move
    together and that shared movement dominates the leading components, while independent
    per-site drift does not. Pass background features (not the ones under test) as ``x``, and
    pass ``ages`` among ``covariates`` to get components that carry no age trend.

    Returns the ``(n_samples, k)`` component scores (orthogonal to ``covariates`` and to the
    intercept) and the fraction of the standardised residual variance each component explains.
    """
    x = np.asarray(x, dtype=float)
    n = x.shape[0]
    if not 0 < k < n:
        raise ValueError(f"k must be between 1 and {n - 1} for {n} samples, got {k}")
    ones = np.ones((n, 1))
    design = ones if covariates is None else np.column_stack([ones, np.asarray(covariates, dtype=float)])
    coef, *_ = np.linalg.lstsq(design, x, rcond=None)
    residual = x - design @ coef
    sd = residual.std(axis=0)
    z = residual[:, sd > _EPS] / sd[sd > _EPS]
    values, vectors = np.linalg.eigh(z @ z.T)
    order = np.argsort(values)[::-1][:k]
    scores = vectors[:, order] * np.sqrt(np.clip(values[order], 0.0, None))
    return scores, values[order] / values.sum()
