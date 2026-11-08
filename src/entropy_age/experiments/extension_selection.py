"""Extension: gene- and cell-type-specific stochasticity under selection.

Virtual experiment with known ground truth. One methylome is aged in four cell
types that share the same epimutation rates and differ only in turnover,
selection and replication errors. A quarter of the sites are "constrained"
(epimutations there cost fitness); the rest are neutral passengers.

Questions the figure answers:
  a. Are constrained sites protected from entropy, and in which cell types?
  b. Does the per-site estimator (:func:`fit_relaxation`) recover that from
     cross-sectional bulk data, as it would have to on real data?
  c. What do neutral sites show? (clonal drift / hitchhiking -> dispersion)
  d. Does a clock read protected tissues as younger?
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import polars as pl

from entropy_age import plots
from entropy_age.clock import ClockConfig, evaluate, train_clock
from entropy_age.experiments.common import Context, ground_state, write_table
from entropy_age.simulate.methylation import MaintenanceRates
from entropy_age.simulate.selection import CellPopulation, SelectionModel, simulate_with_selection
from entropy_age.stochasticity import compare_groups, drift_metrics, fit_relaxation

N_SITES = 400
N_CONSTRAINED = 100
MAX_AGE = 80
CELL_TYPES: tuple[tuple[str, float, float, float], ...] = (
    # name, turnover per step, selection strength, replication error multiplier
    ("post-mitotic", 0.0, 0.5, 0.0),
    ("renewing, neutral", 0.2, 0.0, 0.0),
    ("renewing + selection", 0.2, 0.5, 0.0),
    ("renewing + selection + replication errors", 0.2, 0.5, 1.0),
)


def _design(ctx: Context) -> tuple[np.ndarray, np.ndarray]:
    """Ground levels and constraint labels with matched baseline distributions.

    Sites are drawn from a biological ground state when available, restricted to
    near-bimodal levels (as promoters and gene bodies are); constrained and neutral
    sites are sampled from the same pool, so they start from the same levels.
    """
    rng = ctx.rng(90_000)
    full = ground_state(ctx.paths, rng)
    pool = np.flatnonzero((full.young < 0.1) | (full.young > 0.9))
    if pool.size < N_SITES:
        pool = np.arange(len(full.features))
    ground = full.young[rng.choice(pool, N_SITES, replace=False)]
    constrained = np.zeros(N_SITES, dtype=bool)
    constrained[rng.choice(N_SITES, N_CONSTRAINED, replace=False)] = True
    return ground, constrained


def run(ctx: Context) -> list[Path]:
    plots.apply_style()
    out = ctx.out_dir("extension_selection")
    ground, constrained = _design(ctx)
    weights = constrained.astype(float)
    rates = MaintenanceRates.universal(N_SITES, 0.99)
    ages = np.repeat(np.arange(0, MAX_AGE + 1, 2), 3)
    features = [f"{'constrained' if c else 'neutral'}_{i}" for i, c in enumerate(constrained)]
    membership = pl.DataFrame({"feature": features, "group": np.where(constrained, "constrained", "neutral")})

    data: dict[str, np.ndarray] = {}
    metric_frames, contrast_rows = [], []
    for i, (name, turnover, strength, division_error) in enumerate(CELL_TYPES):
        population = CellPopulation(
            name,
            rates,
            SelectionModel(weights, strength, turnover, division_error),
            n_cells=200,
            meta={"turnover": turnover, "strength": strength, "division_error": division_error},
        )
        x = simulate_with_selection(ground, ages, population, rng=ctx.rng(91_000 + i))
        data[name] = x
        relax = fit_relaxation(x, ages, features, t0=0.0)
        drift = drift_metrics(x, ages, features)
        metrics = relax.join(drift, on="feature").with_columns(
            pl.col("window_change").abs(), cell_type=pl.lit(name), constrained=pl.Series(constrained)
        )
        metric_frames.append(metrics)
        for metric in ("window_change", "rate", "entropy_slope", "dispersion_slope"):
            table = compare_groups(
                metrics, membership.filter(pl.col("group") == "constrained"), metric=metric, n_permutations=1000, seed=i
            )
            contrast_rows.extend({"cell_type": name, "metric": metric, **row} for row in table.to_dicts())
    metrics = pl.concat(metric_frames)
    contrasts = pl.DataFrame(contrast_rows)
    write_table(metrics, out / "site_metrics")
    write_table(contrasts, out / "constrained_vs_neutral")

    codes = dict(zip([c[0] for c in CELL_TYPES], ("PM", "RN", "RS", "RS+E"), strict=True))
    fig, axes = plots.panels(4, ncols=2, size=3.6)
    # a: loss of the ground state over time
    for (name, x), color in zip(data.items(), plots.PALETTE, strict=False):
        lost = np.abs(x - ground)
        plots.mean_band(axes[0], ages, lost[:, constrained].mean(axis=1), label=f"{codes[name]}: {name}", color=color)
        axes[0].plot(
            np.unique(ages),
            [lost[ages == a][:, ~constrained].mean() for a in np.unique(ages)],
            color=color,
            lw=1,
            ls="--",
        )
    axes[0].set(xlabel="Simulated age (steps)", ylabel="Mean |m − ground|", title="Loss of the ground state")
    axes[0].text(
        0.98,
        0.04,
        "solid: constrained\ndashed: neutral",
        transform=axes[0].transAxes,
        ha="right",
        fontsize=7.5,
        color=plots.MUTED,
    )
    plots.label_panel(axes[0], "a")

    # b: what the per-site estimator recovers from cross-sectional bulk data
    box = metrics.with_columns(
        level=pl.format(
            "{} | {}", pl.col("cell_type"), pl.when(pl.col("constrained")).then(pl.lit("C")).otherwise(pl.lit("N"))
        ),
        change=pl.col("window_change").abs(),
    )
    order = [f"{n} | {g}" for n in codes for g in ("N", "C")]
    labels = [f"{codes[n]} {g}" for n in codes for g in ("N", "C")]
    plots.sweep_box(
        axes[1],
        box,
        x="level",
        y="change",
        order=order,
        labels=labels,
        xlabel="Cell type, N = neutral / C = constrained sites",
        ylabel="Fitted change over life, |m(80) − m(0)|",
        title="Per-site estimate from bulk data",
    )
    plots.label_panel(axes[1], "b")

    # c: neutral sites diverge between individuals when clones expand
    for x, color in zip(data.values(), plots.PALETTE, strict=False):
        levels = np.unique(ages)
        sd = [x[ages == a][:, ~constrained].std(axis=0, ddof=1).mean() for a in levels]
        axes[2].plot(levels, sd, color=color, lw=1.5)
    axes[2].set(
        xlabel="Simulated age (steps)",
        ylabel="Between-individual SD (neutral sites)",
        title="Clonal drift and hitchhiking",
    )
    plots.label_panel(axes[2], "c")

    # d: one clock, trained on the post-mitotic type, reads the others
    reference = data["post-mitotic"]
    clock = train_clock(reference, ages, config=ClockConfig.fast())
    series, clock_rows = {}, []
    for name, x in data.items():
        pred = clock.predict(x)
        series[name] = pred
        clock_rows.append(
            {"cell_type": name, **evaluate(ages, pred).as_dict(), "slope": float(np.polyfit(ages, pred, 1)[0])}
        )
    plots.series_scatter(
        axes[3], ages, {codes[k]: v for k, v in series.items()}, title="Clock trained on post-mitotic cells", size=3
    )
    plots.label_panel(axes[3], "d")
    write_table(pl.DataFrame(clock_rows), out / "clock_transfer")

    handles, legend_labels = axes[0].get_legend_handles_labels()
    fig.legend(handles, legend_labels, loc="lower center", ncol=2, fontsize=8, bbox_to_anchor=(0.5, -0.06))
    fig.suptitle(
        "Selection among cells protects constrained genes from entropy only where cells turn over",
        x=0.02,
        ha="left",
        fontsize=12,
        fontweight="bold",
    )
    fig.tight_layout()
    return [plots.save(fig, ctx.figure_path("extension_selection"))]
