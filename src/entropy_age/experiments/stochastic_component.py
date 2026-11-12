"""How much of a clock is stochastic, and how much is knowing each CpG's age effect? (roadmap RM23)

Tong et al. 2024 (Nature Aging 4:886) simulate DNA methylation change from a young ground state and retrain a
clock's CpGs on the simulation; the share of the clock's accuracy the stochastic clock keeps
(RR2 = R²(stochastic) / R²(original), age against prediction) is its "stochastic component". In their model each
CpG moves, at each step, with probability ``1 - exp(-gamma |effect|)``, always in the direction of its observed
old-minus-young effect, by a half-normal step on the probit scale. The rate and direction are therefore set by
the CpG's real age effect; only the timing is random.

This module repeats their procedure in Hannum blood (ground state: the youngest donors; effect: oldest minus
youngest) for Horvath, PhenoAge and Zhang 2019, and strips information step by step (``MODES``):

* ``tong``: rate and direction from each CpG's effect (their model);
* ``direction``: one common rate, each CpG's own direction;
* ``noise``: one common rate, a random direction at every step (Meyer and Schumacher's noise);
* ``shield``: random direction, a common rate scaled down at the CpGs of essential and housekeeping promoters by
  ``SHIELD`` (the measured gap): the shield without any CpG's own age effect.

(gamma, sigma) are fitted per mode by grid search to the oldest donors' mean levels. Each mode trains an elastic
net (alpha 0.5 as in the paper) on a simulated cohort and is scored on the real donors. Intrinsic acceleration
removes 11 of 12 blood-cell fractions estimated by constrained projection (Houseman) on EpiDISH's cent12CT450k
reference (Tong et al. used EpiDISH's RPC).

``entropy-age figures stochastic_component`` (after ``fetch clocks``, ``fetch epidish`` and the Hannum cohort).
"""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import polars as pl
from scipy import optimize, stats
from sklearn.linear_model import ElasticNet

from entropy_age import plots
from entropy_age.data.clocks import ClockLayout, read_clock
from entropy_age.experiments.common import Context, write_table
from entropy_age.experiments.shield_clocks import linear_score
from entropy_age.experiments.tissue_contrast import _features, _inputs_missing, _matrix

PAIRS = {"Horvath 2013": "StocH", "PhenoAge": "StocP", "Zhang 2019": "StocZ"}
MODES = ("tong", "direction", "noise", "shield")
STEPS_PER_YEAR = 35
YOUNG_MAX, OLD_MIN = 35.0, 75.0
SHIELD = 0.7
GAMMAS = np.geomspace(0.05, 50, 9)
SIGMAS = np.geomspace(1e-4, 5e-2, 7)
SIM_AGES = np.repeat(np.arange(36, 101, 2), 6)
EPIDISH = "epidish/cent12CT450k.m.rda"
EPIDISH_URL = "https://github.com/sjczheng/EpiDISH/raw/devel/data/cent12CT450k.m.rda"


def step_probability(mode: str, effect: np.ndarray, gamma: float, shielded: np.ndarray) -> np.ndarray:
    """Per CpG, the probability of a change in one step."""
    mean_effect = float(np.mean(np.abs(effect)))
    if mode == "tong":
        rate = np.abs(effect)
    elif mode in ("direction", "noise"):
        rate = np.full(effect.size, mean_effect)
    elif mode == "shield":
        rate = np.where(shielded, SHIELD, 1.0) * mean_effect
    else:
        raise ValueError(f"unknown mode {mode!r}")
    return 1.0 - np.exp(-gamma * rate)


def simulate(
    ground: np.ndarray,
    ages: np.ndarray,
    start_age: float,
    effect: np.ndarray,
    shielded: np.ndarray,
    mode: str,
    gamma: float,
    sigma: float,
    rng,
) -> np.ndarray:
    """Samples x CpGs: each sample starts from a random ground profile and takes ``STEPS_PER_YEAR`` steps a year."""
    x = stats.norm.ppf(np.clip(ground[rng.integers(0, ground.shape[0], ages.size)], 1e-4, 1 - 1e-4))
    steps = np.round((ages - start_age) * STEPS_PER_YEAR).astype(int)
    prob = step_probability(mode, effect, gamma, shielded)
    direction = np.sign(effect)
    for t in range(int(steps.max(initial=0))):
        active = steps > t
        if not active.any():
            break
        n = int(active.sum())
        move = rng.random((n, effect.size)) < prob
        size = np.abs(rng.normal(0.0, sigma, (n, effect.size)))
        sign = rng.choice((-1.0, 1.0), (n, effect.size)) if mode in ("noise", "shield") else direction
        x[active] += move * sign * size
    return stats.norm.cdf(x)


def calibrate(ground, old_mean, old_age, start_age, effect, shielded, mode, rng) -> tuple[float, float, float]:
    """(gamma, sigma, MAE) minimising the MAE of simulated old-age means against the observed ones."""
    best = (np.nan, np.nan, np.inf)
    ages = np.full(10, old_age)
    for gamma in GAMMAS:
        for sigma in SIGMAS:
            sim = simulate(ground, ages, start_age, effect, shielded, mode, gamma, sigma, rng).mean(axis=0)
            mae = float(np.mean(np.abs(sim - old_mean)))
            if mae < best[2]:
                best = (float(gamma), float(sigma), mae)
    return best


def train_stochastic_clock(x: np.ndarray, ages: np.ndarray):
    mu, sd = x.mean(axis=0), x.std(axis=0)
    sd = np.where(sd > 0, sd, 1.0)
    model = ElasticNet(alpha=1e-3, l1_ratio=0.5, max_iter=20000).fit((x - mu) / sd, ages)
    return lambda z: model.predict((z - mu) / sd)


def deconvolve(x: np.ndarray, reference: np.ndarray) -> np.ndarray:
    """Samples x cell types: non-negative least squares per sample, rows summing to one (Houseman projection)."""
    out = np.vstack([optimize.nnls(reference, row)[0] for row in x])
    totals = out.sum(axis=1, keepdims=True)
    return out / np.where(totals > 0, totals, 1.0)


def r2(a: np.ndarray, b: np.ndarray) -> float:
    return float(np.corrcoef(a, b)[0, 1] ** 2)


def residual(y: np.ndarray, design: np.ndarray) -> np.ndarray:
    coef, *_ = np.linalg.lstsq(design, y, rcond=None)
    return y - design @ coef


def run(ctx: Context) -> list[Path]:
    import rdata

    if _inputs_missing(ctx, "stochastic_component"):
        return []
    ref_path = ctx.paths.input / EPIDISH
    clocks = ClockLayout(ctx.paths)
    if not ref_path.exists() or not all(clocks.file(c).exists() for c in PAIRS):
        print("stochastic_component: run `entropy-age fetch clocks` and `fetch epidish` first.", file=sys.stderr)
        return []
    plots.apply_style()
    out = ctx.out_dir("stochastic_component")
    rng = ctx.rng(170_000)
    ref = rdata.read_rda(str(ref_path))
    ref = next(iter(ref.values()))
    ref_cpgs, cell_types = [str(c) for c in ref.coords["dim_0"].values], [str(c) for c in ref.coords["dim_1"].values]
    _, _, membership = _features(ctx.paths, ctx.rng(120_000))
    shielded_cpgs = set(membership.filter(pl.col("group").is_in(["housekeeping", "core_essential"]))["feature"])
    from entropy_age.data.cohorts import layout as cohort_layout

    samples = pl.read_parquet(cohort_layout("hannum", ctx.paths).samples).drop_nulls("age").sort("age", "sample")
    ages = samples["age"].to_numpy()
    weights = {c: read_clock(clocks.file(c)) for c in [*PAIRS, *PAIRS.values()]}
    m_ref = _matrix("hannum", ctx.paths, samples["sample"].to_list(), ref_cpgs)
    keep_ref = [i for i, c in enumerate(ref_cpgs) if c in set(m_ref.features)]
    x_ref = np.where(np.isnan(m_ref.values), np.nanmean(m_ref.values, axis=0), m_ref.values)
    ref_matrix = ref.values[keep_ref][:, :]
    order = [m_ref.features.index(ref_cpgs[i]) for i in keep_ref]
    fractions = deconvolve(x_ref[:, order], ref_matrix)
    young, old = ages <= YOUNG_MAX, ages >= OLD_MIN
    start_age, old_age = float(np.median(ages[young])), float(np.median(ages[old]))
    test = ~young
    rows = []
    for original, stoc in PAIRS.items():
        cpgs = weights[original]["cpg"].to_list()
        m = _matrix("hannum", ctx.paths, samples["sample"].to_list(), cpgs)
        x = np.where(np.isnan(m.values), np.nanmean(m.values, axis=0), m.values)
        feats = m.features
        score_orig, _ = linear_score(x, feats, weights[original])
        score_pub, _ = linear_score(x, feats, weights[stoc])
        base = {"clock": original, "cpgs": len(feats), "test_n": int(test.sum())}
        design_age = np.column_stack([np.ones(test.sum()), ages[test]])
        design_iaa = np.column_stack([design_age, fractions[test][:, :-1]])
        eaa_orig, iaa_orig = residual(score_orig[test], design_age), residual(score_orig[test], design_iaa)
        r2_orig = r2(score_orig[test], ages[test])
        rows.append(
            {
                **base,
                "model": "published Stoc",
                "gamma": np.nan,
                "sigma": np.nan,
                "fit_mae": np.nan,
                "r2": r2(score_pub[test], ages[test]),
                "r2_original": r2_orig,
                "eaa_r": float(np.corrcoef(residual(score_pub[test], design_age), eaa_orig)[0, 1]),
                "iaa_r": float(np.corrcoef(residual(score_pub[test], design_iaa), iaa_orig)[0, 1]),
            }
        )
        effect = x[old].mean(axis=0) - x[young].mean(axis=0)
        shielded = np.array([f in shielded_cpgs for f in feats])
        for mode in MODES:
            gamma, sigma, mae = calibrate(
                x[young], x[old].mean(axis=0), old_age, start_age, effect, shielded, mode, rng
            )
            sim = simulate(x[young], SIM_AGES.astype(float), start_age, effect, shielded, mode, gamma, sigma, rng)
            predict = train_stochastic_clock(sim, SIM_AGES.astype(float))
            pred = predict(x[test])
            rows.append(
                {
                    **base,
                    "model": mode,
                    "gamma": gamma,
                    "sigma": sigma,
                    "fit_mae": mae,
                    "r2": r2(pred, ages[test]),
                    "r2_original": r2_orig,
                    "eaa_r": float(np.corrcoef(residual(pred, design_age), eaa_orig)[0, 1]),
                    "iaa_r": float(np.corrcoef(residual(pred, design_iaa), iaa_orig)[0, 1]),
                }
            )
            print(f"stochastic_component: {original}, {mode} done", file=sys.stderr, flush=True)
    table = pl.DataFrame(rows).with_columns(rr2=pl.col("r2") / pl.col("r2_original"))
    write_table(table, out / "fractions")
    write_table(pl.DataFrame(fractions, schema=cell_types).with_columns(age=pl.Series(ages)), out / "cell_fractions")
    return [_figure(ctx, table)]


def _figure(ctx: Context, table: pl.DataFrame) -> Path:
    fig, axes = plots.panels(2, ncols=2, size=4.4)
    models = ["published Stoc", *MODES]
    for i, (col, label) in enumerate(
        (("rr2", "RR2: stochastic R² / original R²"), ("iaa_r", "r with the original's intrinsic acceleration"))
    ):
        ax = axes[i]
        for c, clock in enumerate(PAIRS):
            sub = table.filter(pl.col("clock") == clock)
            v = [sub.filter(pl.col("model") == m)[col][0] for m in models]
            ax.plot(range(len(models)), v, "o-", color=plots.PALETTE[c], label=clock)
        ax.set_xticks(range(len(models)), models, rotation=30, ha="right", fontsize=7)
        ax.set(ylabel=label, title=label.split(":")[0])
        if i == 0:
            ax.legend(fontsize=7, frameon=False)
        plots.label_panel(ax, "ab"[i])
    fig.suptitle(
        "Stochastic component of clocks in Hannum blood, as information is stripped from the model",
        x=0.02,
        y=1.04,
        ha="left",
        fontsize=11,
        fontweight="bold",
    )
    return plots.save(fig, ctx.figure_path("stochastic_component"))
