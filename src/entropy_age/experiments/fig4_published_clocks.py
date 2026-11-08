"""Figure 4 / Extended Data Fig. 6-7: published clocks read simulated stochastic variation.

Every published clock is applied to methylomes simulated from a 16-year-old
blood sample (GSE41037), with simulated age starting at 16 so the axes are in
years. Needs ``entropy-age fetch geo GSE41037``.
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import polars as pl

from entropy_age import plots
from entropy_age.clock import evaluate
from entropy_age.data.upstream import UpstreamLayout
from entropy_age.experiments.common import HUMAN_OLD, HUMAN_YOUNG, Context, human_blood, write_table
from entropy_age.published import LinearClock
from entropy_age.simulate.methylation import MaintenanceRates, simulate_methylation

LIMITS_EM = (0.95, 0.97, 0.99)
LIMITS_ED = (0.01, 0.05, 0.23)
AGE_STEPS = 74  # simulated steps 1..73, shifted to ages 16..88 (the ground and equilibrium donors)
START_AGE = 15
LABEL = {"horvath": "Horvath", "phenoage": "PhenoAge"}


def run(ctx: Context) -> list[Path]:
    blood = human_blood(ctx.paths)
    if blood is None:
        print("fig4: GSE41037 not prepared; run `entropy-age fetch geo GSE41037` first - skipped")
        return []
    plots.apply_style()
    matrix, samples = blood
    out = ctx.out_dir("fig4_published_clocks")
    layout = UpstreamLayout(ctx.paths)
    clocks = {
        name: LinearClock.read(layout.clock(name), name)
        for name in ("horvath", "phenoage", "lin", "vidal_bralo", "weidner")
    }
    steps = np.arange(1, AGE_STEPS)
    ages = steps + START_AGE
    n_cells = ctx.scale.n_cells

    def simulate(clock: LinearClock, rates_for, young: str, rng: np.random.Generator) -> np.ndarray:
        sub = matrix.select_features(clock.features)
        ground = sub.sample(young)
        return clock.predict(simulate_methylation(ground, steps, rates_for(sub), n_cells=n_cells, rng=rng))

    def empirical(em_limit: float, ed_limit: float, rng: np.random.Generator):
        return lambda sub: MaintenanceRates.from_equilibrium(
            sub.sample(HUMAN_OLD), rng, em_limit=em_limit, ed_limit=ed_limit
        )

    def universal(em: float):
        return lambda sub: MaintenanceRates.universal(len(sub.features), em)

    fig, axes = plots.panels(8, ncols=4, size=3.0)
    rows = []
    for c_i, name in enumerate(("horvath", "phenoage")):
        clock = clocks[name]
        sweep_rows = []
        for rep in range(ctx.scale.replicates):
            for em in LIMITS_EM:
                for ed in LIMITS_ED:
                    rng = ctx.rng(10_000 + c_i * 1000 + rep * 100 + int(em * 100) + int(ed * 100))
                    pred = simulate(clock, empirical(em, ed, rng), HUMAN_YOUNG, rng)
                    m = evaluate(ages, pred)
                    sweep_rows.append(
                        {
                            "clock": name,
                            "em_limit": em,
                            "ed_limit": ed,
                            "limits": f"{em}/{ed}",
                            "replicate": rep,
                            **m.as_dict(),
                        }
                    )
        sweep = pl.DataFrame(sweep_rows)
        rows.append(sweep)
        base = c_i * 4
        plots.sweep_box(
            axes[base],
            sweep,
            x="limits",
            y="r2",
            xlabel="Em lower limit / Ed upper limit",
            ylabel="R² (simulated vs predicted)",
            title=f"{LABEL[name]}: rate limits",
            ylim=(-1.0, 1.0),
        )
        plots.label_panel(axes[base], "ae"[c_i])

        rng = ctx.rng(20_000 + c_i)
        pred = simulate(clock, empirical(0.97, 0.05, rng), HUMAN_YOUNG, rng)
        plots.prediction_scatter(
            axes[base + 1],
            ages,
            pred,
            metrics=evaluate(ages, pred),
            identity=False,
            title=f"{LABEL[name]}: empirical, 97 % / 5 %",
            ylabel=f"Predicted age ({name})",
        )
        plots.label_panel(axes[base + 1], "bf"[c_i])

        pred = simulate(clock, universal(0.99), HUMAN_YOUNG, rng)
        plots.prediction_scatter(
            axes[base + 2],
            ages,
            pred,
            metrics=evaluate(ages, pred),
            identity=False,
            title=f"{LABEL[name]}: universal Em = 99 %",
            ylabel=f"Predicted age ({name})",
        )
        plots.label_panel(axes[base + 2], "cg"[c_i])
    write_table(pl.concat(rows), out / "limit_sweep")

    # d: other first-generation clocks under a universal 99 % maintenance
    rng = ctx.rng(30_000)
    series, other_rows = {}, []
    for name in ("lin", "vidal_bralo", "weidner"):
        pred = simulate(clocks[name], universal(0.99), HUMAN_YOUNG, rng)
        m = evaluate(ages, pred)
        series[f"{name} (r = {m.pearson_r:.2f})"] = pred
        other_rows.append({"clock": name, **m.as_dict()})
    plots.series_scatter(axes[3], ages, series, title="Other clocks, Em = 99 %")
    plots.label_panel(axes[3], "d")

    # h: the ground state's own age shifts the curve (Horvath, empirical rates)
    ground_rows, series = [], {}
    for target_age in (16, 37, 81):
        donor = samples.filter(pl.col("age") == target_age)["sample"].to_list()
        donor_id = HUMAN_YOUNG if target_age == 16 else donor[0]
        rng = ctx.rng(40_000 + target_age)
        pred = simulate(clocks["horvath"], empirical(0.97, 0.05, rng), donor_id, rng)
        m = evaluate(steps, pred)
        series[f"ground {target_age} y"] = pred
        ground_rows.append({"ground_age": target_age, "sample": donor_id, **m.as_dict()})
    plots.series_scatter(
        axes[7],
        steps,
        series,
        xlabel="Simulated steps",
        ylabel="Predicted age (Horvath)",
        title="Ground state from older donors",
    )
    plots.label_panel(axes[7], "h")
    write_table(pl.DataFrame(other_rows), out / "other_clocks")
    write_table(pl.DataFrame(ground_rows), out / "ground_state_age")

    fig.suptitle(
        "Published clocks predict the amount of simulated stochastic variation",
        x=0.02,
        ha="left",
        fontsize=12,
        fontweight="bold",
    )
    fig.tight_layout()
    return [plots.save(fig, ctx.figure_path("fig4_published_clocks"))]
