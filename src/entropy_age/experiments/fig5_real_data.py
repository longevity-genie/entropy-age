"""Figure 5-6: a clock trained only on simulated data predicts the age of real samples.

* Human blood (GSE41037 controls): ground state = the 16-year-old, empirical
  maintenance rates from the 88-year-old, Horvath's 353 CpGs; with and without
  regression-based cell-composition correction; random CpG sets and a permuted-age null.
* Bottlenose dolphin blood (pan-mammalian example data shipped upstream):
  ground state = the 0.57-year-old, all 37,554 CpGs of the mammalian array, the
  log-log relative-age transform (analogues of the paper's clocks 3 and 4).
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import polars as pl

from entropy_age import plots
from entropy_age.clock import ClockConfig, LinearRescale, LogLogRelativeAge, evaluate, train_clock
from entropy_age.data.matrix import FeatureMatrix
from entropy_age.data.upstream import UpstreamLayout
from entropy_age.experiments.common import (
    DOLPHIN_MAX_AGE,
    DOLPHIN_OLD,
    DOLPHIN_YOUNG,
    HUMAN_OLD,
    HUMAN_YOUNG,
    Context,
    dolphin_blood,
    human_blood,
    write_table,
)
from entropy_age.published import LinearClock
from entropy_age.simulate.methylation import MaintenanceRates, simulate_methylation
from entropy_age.stochasticity import regress_out, sample_entropy

HUMAN_STEPS = np.arange(1, 74)
CELL_TYPES = ("B", "NK", "CD4T", "CD8T", "Mono", "Neutro", "Eosino")


def stochastic_clock_predictions(
    matrix: FeatureMatrix,
    *,
    young: str,
    old: str,
    steps: np.ndarray,
    transform,
    config: ClockConfig,
    rng: np.random.Generator,
    n_cells: int,
    em_limit: float = 0.97,
    ed_limit: float = 0.05,
    universal_em: float | None = None,
) -> tuple[np.ndarray, int]:
    """Simulate from ``young``, train on simulated steps, predict every sample of ``matrix``."""
    ground = matrix.sample(young)
    if universal_em is None:
        rates = MaintenanceRates.from_equilibrium(matrix.sample(old), rng, em_limit=em_limit, ed_limit=ed_limit)
    else:
        rates = MaintenanceRates.universal(len(matrix.features), universal_em)
    simulated = simulate_methylation(ground, steps, rates, n_cells=n_cells, rng=rng)
    clock = train_clock(simulated, steps, config=config, transform=transform, features=matrix.features)
    return clock.predict(matrix.values), clock.n_nonzero


def _human(ctx: Context, axes, out: Path) -> list[dict]:
    blood = human_blood(ctx.paths)
    if blood is None:
        for ax in axes:
            ax.set_visible(False)
        print("fig5: GSE41037 not prepared; human panels skipped")
        return []
    matrix, samples = blood
    horvath = LinearClock.read(UpstreamLayout(ctx.paths).clock("horvath"), "horvath")
    keep = ~samples["sample"].is_in([HUMAN_YOUNG, HUMAN_OLD]).to_numpy()
    ages = samples["age"].to_numpy()
    transform = LinearRescale(HUMAN_STEPS.min(), HUMAN_STEPS.max())
    config = ClockConfig.fixed_alpha(1.0)
    rows: list[dict] = []
    prediction_rows: list[pl.DataFrame] = []

    def predict(m: FeatureMatrix, stream: int) -> tuple[np.ndarray, int]:
        return stochastic_clock_predictions(
            m,
            young=HUMAN_YOUNG,
            old=HUMAN_OLD,
            steps=HUMAN_STEPS,
            transform=transform,
            config=config,
            rng=ctx.rng(stream),
            n_cells=ctx.scale.n_cells,
        )

    clock_sites = matrix.select_features(horvath.features)
    pred, n_nonzero = predict(clock_sites, 50_000)
    m = evaluate(ages[keep], pred[keep])
    rows.append(
        {
            "dataset": "GSE41037",
            "variant": "horvath_cpgs",
            "n_sites": len(clock_sites.features),
            "n_nonzero": n_nonzero,
            **m.as_dict(),
        }
    )
    prediction_rows.append(
        pl.DataFrame({"sample": samples["sample"], "age": ages, "predicted": pred, "variant": "horvath_cpgs"})
    )
    plots.prediction_scatter(
        axes[0],
        ages[keep],
        pred[keep],
        metrics=m,
        identity=False,
        xlabel="Chronological age (years)",
        ylabel="Predicted (stochastic clock)",
        title="Human blood, Horvath CpGs",
    )
    plots.label_panel(axes[0], "a")

    cell_path = UpstreamLayout(ctx.paths).gse41037_cell_types
    cells = pl.read_parquet(cell_path)
    order = pl.DataFrame({"sample": matrix.samples}).join(cells, on="sample", how="left", maintain_order="left")
    fractions = order.select(CELL_TYPES).to_numpy()
    corrected = FeatureMatrix(regress_out(clock_sites.values, fractions), clock_sites.samples, clock_sites.features)
    pred_c, n_nonzero = predict(corrected, 50_001)
    m = evaluate(ages[keep], pred_c[keep])
    rows.append(
        {
            "dataset": "GSE41037",
            "variant": "horvath_cpgs_cell_corrected",
            "n_sites": len(corrected.features),
            "n_nonzero": n_nonzero,
            **m.as_dict(),
        }
    )
    prediction_rows.append(
        pl.DataFrame(
            {"sample": samples["sample"], "age": ages, "predicted": pred_c, "variant": "horvath_cpgs_cell_corrected"}
        )
    )
    plots.prediction_scatter(
        axes[1],
        ages[keep],
        pred_c[keep],
        metrics=m,
        identity=False,
        xlabel="Chronological age (years)",
        ylabel="Predicted (stochastic clock)",
        title="… after cell-composition correction",
    )
    plots.label_panel(axes[1], "b")

    size_rows = []
    permuted = ctx.rng(50_100).permutation(ages[keep])
    for size in (100, 500, 2000, 5000):
        for rep in range(ctx.scale.replicates):
            rng = ctx.rng(51_000 + size + rep)
            index = np.sort(rng.choice(len(matrix.features), size=size, replace=False))
            sub = matrix.select_features([matrix.features[i] for i in index])
            p, _ = predict(sub, 52_000 + size + rep)
            size_rows.append(
                {"n_sites": size, "replicate": rep, "ages": "true", **evaluate(ages[keep], p[keep]).as_dict()}
            )
            size_rows.append(
                {"n_sites": size, "replicate": rep, "ages": "permuted", **evaluate(permuted, p[keep]).as_dict()}
            )
    sizes = pl.DataFrame(size_rows)
    write_table(sizes, out / "human_random_cpgs")
    plot_frame = sizes.with_columns(level=pl.format("{} {}", pl.col("n_sites"), pl.col("ages")))
    order_levels = [f"{s} {a}" for s in (100, 500, 2000, 5000) for a in ("true", "permuted")]
    plots.sweep_box(
        axes[2],
        plot_frame,
        x="level",
        y="pearson_r",
        order=order_levels,
        xlabel="Random CpGs (true / permuted ages)",
        ylabel="Pearson r",
        title="Random CpG sets",
        ylim=(-0.5, 1.0),
    )
    plots.label_panel(axes[2], "c")

    entropy = sample_entropy(matrix.values)
    axes[3].scatter(ages, entropy, s=5, color=plots.POINT, alpha=0.5, linewidths=0)
    r = np.corrcoef(ages, entropy)[0, 1]
    axes[3].set(
        xlabel="Chronological age (years)",
        ylabel="Mean CpG entropy (bits)",
        title=f"Methylome entropy rises with age (r = {r:.2f})",
    )
    plots.label_panel(axes[3], "d")
    write_table(pl.concat(prediction_rows), out / "human_predictions")
    write_table(pl.DataFrame({"sample": samples["sample"], "age": ages, "entropy": entropy}), out / "human_entropy")
    return rows


def _dolphin(ctx: Context, axes, out: Path) -> list[dict]:
    dolphin = dolphin_blood(ctx.paths)
    if dolphin is None:
        for ax in axes:
            ax.set_visible(False)
        print("fig5: upstream data not prepared; dolphin panels skipped")
        return []
    matrix, samples = dolphin
    ages = samples["age"].to_numpy()
    keep = ~samples["sample"].is_in([DOLPHIN_YOUNG, DOLPHIN_OLD]).to_numpy()
    steps = np.arange(1, DOLPHIN_MAX_AGE)
    transform = LogLogRelativeAge(DOLPHIN_MAX_AGE)
    config = ClockConfig.fixed_alpha(1.0, l1_ratio=(0.01, 0.001))
    rows: list[dict] = []
    frames = []
    for i, (label, universal) in enumerate(
        (("clock 3: empirical Em/Ed", None), ("clock 4: universal Em = 99 %", 0.99))
    ):
        pred, n_nonzero = stochastic_clock_predictions(
            matrix,
            young=DOLPHIN_YOUNG,
            old=DOLPHIN_OLD,
            steps=steps,
            transform=transform,
            config=config,
            rng=ctx.rng(60_000 + i),
            n_cells=ctx.scale.n_cells,
            universal_em=universal,
        )
        rel = ages / DOLPHIN_MAX_AGE
        m = evaluate(rel[keep], pred[keep])
        rows.append(
            {
                "dataset": "dolphin",
                "variant": label,
                "n_sites": len(matrix.features),
                "n_nonzero": n_nonzero,
                **m.as_dict(),
            }
        )
        frames.append(
            pl.DataFrame({"sample": samples["sample"], "age": ages, "predicted_relative_age": pred, "variant": label})
        )
        plots.prediction_scatter(
            axes[i],
            rel[keep],
            pred[keep],
            metrics=m,
            identity=False,
            xlabel="Relative age (age / 67 y)",
            ylabel="Predicted relative age",
            title=f"Dolphin blood, {label}",
        )
        plots.label_panel(axes[i], "ef"[i])
    write_table(pl.concat(frames), out / "dolphin_predictions")

    entropy = sample_entropy(matrix.values)
    axes[2].scatter(ages, entropy, s=8, color=plots.POINT, alpha=0.7, linewidths=0)
    r = np.corrcoef(ages, entropy)[0, 1]
    axes[2].set(
        xlabel="Age (years)", ylabel="Mean CpG entropy (bits)", title=f"Dolphin methylome entropy (r = {r:.2f})"
    )
    plots.label_panel(axes[2], "g")
    return rows


def run(ctx: Context) -> list[Path]:
    plots.apply_style()
    out = ctx.out_dir("fig5_real_data")
    fig, axes = plots.panels(7, ncols=4, size=3.1)
    rows = _human(ctx, axes[:4], out) + _dolphin(ctx, axes[4:7], out)
    if not rows:
        return []
    write_table(pl.DataFrame(rows), out / "summary")
    fig.suptitle(
        "Clocks trained only on simulated stochastic variation predict real ages",
        x=0.02,
        ha="left",
        fontsize=12,
        fontweight="bold",
    )
    fig.tight_layout()
    return [plots.save(fig, ctx.figure_path("fig5_real_data"))]
