"""Figure 1 / Extended Data Fig. 1: clocks from purely random, bounded data."""

from __future__ import annotations

from pathlib import Path

import numpy as np
import polars as pl

from entropy_age import plots
from entropy_age.clock import ClockMetrics, StochasticClock, evaluate, train_clock
from entropy_age.experiments.common import Context, write_table
from entropy_age.simulate.bulk import Bound, age_grid, random_ground_state, simulate_bulk

SIGMAS = (0.01, 0.025, 0.05, 0.075, 0.1, 0.15, 0.2, 0.3, 0.5)
TRANSFER_SIGMAS = (0.1, 0.2, 0.4, 0.8)


def _train_validate(
    ctx: Context, ground: np.ndarray, rng: np.random.Generator, *, sigma: float, bound: Bound
) -> tuple[StochasticClock, ClockMetrics, np.ndarray, np.ndarray]:
    ages = age_grid(100, ctx.scale.per_age)
    train = simulate_bulk(ground, ages, sigma=sigma, bound=bound, rng=rng)
    valid = simulate_bulk(ground, ages, sigma=sigma, bound=bound, rng=rng)
    clock = train_clock(train, ages, config=ctx.scale.clock)
    pred = clock.predict(valid)
    return clock, evaluate(ages, pred), ages, pred


def run(ctx: Context) -> list[Path]:
    plots.apply_style()
    n = ctx.scale.n_features
    out = ctx.out_dir("fig1_bulk")
    fig, axes = plots.panels(8, ncols=3, size=3.1)

    # a-c: the boundary is what makes a clock
    examples = [
        (Bound.NONE, 0.05, "No bound, σ = 0.05"),
        (Bound.CLIP, 0.05, "Clipped to [0, 1], σ = 0.05"),
        (Bound.LOGIT, 0.2, "Logit scale, σ = 0.2"),
    ]
    example_rows = []
    ground = random_ground_state(n, ctx.rng(1))
    logit_clock: StochasticClock | None = None
    for i, (bound, sigma, title) in enumerate(examples):
        clock, metrics, ages, pred = _train_validate(ctx, ground, ctx.rng(10 + i), sigma=sigma, bound=bound)
        plots.prediction_scatter(axes[i], ages, pred, metrics=metrics, title=title)
        plots.label_panel(axes[i], "abc"[i])
        example_rows.append({"bound": str(bound), "sigma": sigma, **metrics.as_dict()})
        if bound is Bound.LOGIT:
            logit_clock = clock
    write_table(pl.DataFrame(example_rows), out / "examples")

    # d: accuracy across noise levels (logit model)
    rows = []
    for sigma in SIGMAS:
        for rep in range(ctx.scale.replicates):
            rng = ctx.rng(100 + rep * 31 + int(sigma * 1000))
            g = random_ground_state(n, rng)
            _, metrics, *_ = _train_validate(ctx, g, rng, sigma=sigma, bound=Bound.LOGIT)
            rows.append({"sigma": sigma, "replicate": rep, **metrics.as_dict()})
    sweep = pl.DataFrame(rows)
    write_table(sweep, out / "noise_sweep")
    plots.sweep_box(
        axes[3],
        sweep,
        x="sigma",
        y="r2",
        xlabel="Noise SD per step (σ)",
        ylabel="R² (validation)",
        title="Accuracy vs noise level",
        ylim=(-0.05, 1.05),
    )
    plots.label_panel(axes[3], "d")

    # e-f: coefficients are reproducible and follow regression to the mean
    assert logit_clock is not None
    second, *_ = _train_validate(ctx, ground, ctx.rng(20), sigma=0.2, bound=Bound.LOGIT)
    r = float(np.corrcoef(logit_clock.coef, second.coef)[0, 1])
    axes[4].scatter(logit_clock.coef, second.coef, s=4, color=plots.POINT, alpha=0.5, linewidths=0)
    axes[4].set(
        xlabel="Coefficient, run 1", ylabel="Coefficient, run 2", title=f"Reproducible coefficients (r = {r:.2f})"
    )
    plots.label_panel(axes[4], "e")
    axes[5].scatter(ground, logit_clock.coef, s=4, color=plots.POINT, alpha=0.5, linewidths=0)
    axes[5].axhline(0, color=plots.MUTED, lw=0.8)
    axes[5].set(xlabel="Ground state value", ylabel="Elastic net coefficient", title="Regression to the mean")
    plots.label_panel(axes[5], "f")
    write_table(
        pl.DataFrame({"ground": ground, "coef_run1": logit_clock.coef, "coef_run2": second.coef}), out / "coefficients"
    )

    # g: number of features
    rows = []
    for size in ctx.scale.feature_sizes:
        for rep in range(ctx.scale.replicates):
            rng = ctx.rng(200 + size * 7 + rep)
            g = random_ground_state(size, rng)
            _, metrics, *_ = _train_validate(ctx, g, rng, sigma=0.2, bound=Bound.LOGIT)
            rows.append({"n_features": size, "replicate": rep, **metrics.as_dict()})
    sizes = pl.DataFrame(rows)
    write_table(sizes, out / "feature_sweep")
    plots.sweep_box(
        axes[6],
        sizes,
        x="n_features",
        y="r2",
        xlabel="Number of features",
        ylabel="R² (validation)",
        title="Accuracy vs feature count",
        ylim=(-0.05, 1.05),
    )
    plots.label_panel(axes[6], "g")

    # h: more noise per step = faster predicted aging
    rng = ctx.rng(300)
    ages = age_grid(100, ctx.scale.per_age)
    series = {
        f"σ = {s}" + (" (trained)" if s == 0.2 else ""): logit_clock.predict(
            simulate_bulk(ground, ages, sigma=s, bound=Bound.LOGIT, rng=rng)
        )
        for s in TRANSFER_SIGMAS
    }
    plots.series_scatter(axes[7], ages, series, title="Clock trained at σ = 0.2")
    plots.label_panel(axes[7], "h")
    write_table(pl.DataFrame({"age": ages, **{k: v for k, v in series.items()}}), out / "noise_transfer")

    fig.suptitle(
        "Stochastic variation on bounded data is enough to build an aging clock",
        x=0.02,
        ha="left",
        fontsize=12,
        fontweight="bold",
    )
    fig.tight_layout()
    return [plots.save(fig, ctx.figure_path("fig1_bulk"))]
