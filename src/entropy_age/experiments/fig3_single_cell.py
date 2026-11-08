"""Figure 3 / Extended Data Fig. 4-5: single-cell methylation with maintenance efficiencies."""

from __future__ import annotations

from pathlib import Path

import numpy as np
import polars as pl

from entropy_age import plots
from entropy_age.clock import evaluate, train_clock
from entropy_age.experiments.common import Context, GroundState, ground_state, write_table
from entropy_age.simulate.bulk import age_grid
from entropy_age.simulate.methylation import MaintenanceRates, simulate_methylation

EMS = (0.9, 0.95, 0.99, 0.999, 0.9999, 0.99995)
TRANSFER_EMS = (0.95, 0.99, 0.999, 0.9999)


def _pick(ground: GroundState, size: int, rng: np.random.Generator) -> GroundState:
    index = np.sort(rng.choice(len(ground.features), size=min(size, len(ground.features)), replace=False))
    return ground.subset(index)


def run(ctx: Context) -> list[Path]:
    plots.apply_style()
    out = ctx.out_dir("fig3_single_cell")
    full = ground_state(ctx.paths, ctx.rng(0))
    ages = age_grid(100, ctx.scale.per_age, start=1)
    n_cells = ctx.scale.n_cells

    def simulate(g: np.ndarray, rates: MaintenanceRates, rng: np.random.Generator) -> np.ndarray:
        return simulate_methylation(g, ages, rates, n_cells=n_cells, rng=rng)

    def train_validate(g: np.ndarray, rates: MaintenanceRates, rng: np.random.Generator):
        clock = train_clock(simulate(g, rates, rng), ages, config=ctx.scale.clock)
        pred = clock.predict(simulate(g, rates, rng))
        return clock, evaluate(ages, pred), pred

    fig, axes = plots.panels(7, ncols=4, size=3.0)

    # a: universal maintenance efficiency
    rows = []
    for em in EMS:
        for rep in range(ctx.scale.replicates):
            rng = ctx.rng(1000 + rep * 97 + int(em * 1e5))
            g = _pick(full, 500, rng)
            _, metrics, _ = train_validate(g.young, MaintenanceRates.universal(500, em), rng)
            rows.append({"em": em, "replicate": rep, **metrics.as_dict()})
    sweep = pl.DataFrame(rows)
    write_table(sweep, out / "em_sweep")
    plots.sweep_box(
        axes[0],
        sweep,
        x="em",
        y="r2",
        labels=[f"{100 * e:g}" for e in EMS],
        xlabel="Maintenance efficiency Em (%)",
        ylabel="R² (validation)",
        title="Universal Em, Ed = 1 − Em",
        ylim=(-0.05, 1.05),
    )
    plots.label_panel(axes[0], "a")

    # b: example at Em = 99.9 %
    rng = ctx.rng(2000)
    g = _pick(full, ctx.scale.n_features, rng)
    universal = MaintenanceRates.universal(len(g.features), 0.999)
    clock, metrics, pred = train_validate(g.young, universal, rng)
    plots.prediction_scatter(axes[1], ages, pred, metrics=metrics, title="Em = 99.9 %")
    plots.label_panel(axes[1], "b")

    # c: number of sites
    rows = []
    for size in ctx.scale.feature_sizes:
        for rep in range(ctx.scale.replicates):
            rng = ctx.rng(3000 + size * 13 + rep)
            sub = _pick(full, size, rng)
            _, m, _ = train_validate(sub.young, MaintenanceRates.universal(len(sub.features), 0.999), rng)
            rows.append({"n_features": size, "replicate": rep, **m.as_dict()})
    sizes = pl.DataFrame(rows)
    write_table(sizes, out / "feature_sweep")
    plots.sweep_box(
        axes[2],
        sizes,
        x="n_features",
        y="r2",
        xlabel="Number of CpG sites",
        ylabel="R² (validation)",
        title="CpG count (Em = 99.9 %)",
        ylim=(-0.05, 1.05),
    )
    plots.label_panel(axes[2], "c")

    # d: lower maintenance = faster predicted aging
    rng = ctx.rng(4000)
    series = {
        f"Em = {100 * em:g} %" + (" (trained)" if em == 0.999 else ""): clock.predict(
            simulate(g.young, MaintenanceRates.universal(len(g.features), em), rng)
        )
        for em in TRANSFER_EMS
    }
    plots.series_scatter(axes[3], ages, series, title="Clock trained at Em = 99.9 %")
    plots.label_panel(axes[3], "d")
    write_table(pl.DataFrame({"age": ages, **series}), out / "em_transfer")

    # e: empirical site-specific rates from the old sample
    rng = ctx.rng(5000)
    if g.old is not None:
        empirical = MaintenanceRates.from_equilibrium(g.old, rng, em_limit=0.95, ed_limit=0.23)
        label = "Empirical Em/Ed (old sample)"
    else:
        empirical = MaintenanceRates.random(len(g.features), rng, em_limit=0.95, ed_limit=0.23)
        label = "Random Em/Ed"
    _, m_emp, pred_emp = train_validate(g.young, empirical, rng)
    plots.prediction_scatter(axes[4], ages, pred_emp, metrics=m_emp, title=label)
    plots.label_panel(axes[4], "e")

    # f: no regression to equilibrium when starting at it
    rows, series = [], {}
    for start in (0.5, 0.51):
        rng = ctx.rng(6000 + int(start * 100))
        flat = np.full(500, start)
        rates = MaintenanceRates.universal(500, 0.999)
        trained = train_clock(
            simulate_methylation(flat, ages, rates, n_cells=n_cells, jitter=0, rng=rng), ages, config=ctx.scale.clock
        )
        pred = trained.predict(simulate_methylation(flat, ages, rates, n_cells=n_cells, jitter=0, rng=rng))
        m = evaluate(ages, pred)
        series[f"start {start} (r = {m.pearson_r:.2f})"] = pred
        rows.append({"start": start, **m.as_dict()})
    write_table(pl.DataFrame(rows), out / "equilibrium_start")
    plots.series_scatter(axes[5], ages, series, title="Start at the equilibrium (0.5)")
    plots.label_panel(axes[5], "f")

    # g: the binomial simulation follows the closed-form expectation
    rng = ctx.rng(7000)
    pick = rng.choice(len(g.features), 5, replace=False)
    steps = np.arange(0, 400, 10)
    sim_rates = MaintenanceRates.universal(5, 0.99)
    sim = simulate_methylation(g.young[pick], np.repeat(steps, 20), sim_rates, n_cells=n_cells, jitter=0, rng=rng)
    expected = sim_rates.expected(np.floor(n_cells * g.young[pick]) / n_cells, steps)
    for j, color in zip(range(5), plots.PALETTE, strict=False):
        axes[6].scatter(np.repeat(steps, 20), sim[:, j], s=2, color=color, alpha=0.25, linewidths=0)
        axes[6].plot(steps, expected[:, j], color=color, lw=1.4)
    axes[6].axhline(0.5, ls=":", color=plots.MUTED, lw=1)
    axes[6].set(xlabel="Steps", ylabel="Bulk methylation", title="Simulated vs closed form")
    plots.label_panel(axes[6], "g")

    fig.suptitle(
        f"Single-cell methylation simulations (ground state: {full.label})",
        x=0.02,
        ha="left",
        fontsize=12,
        fontweight="bold",
    )
    fig.tight_layout()
    return [plots.save(fig, ctx.figure_path("fig3_single_cell"))]
