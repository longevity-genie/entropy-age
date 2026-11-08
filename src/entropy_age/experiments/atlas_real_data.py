"""Extension on real data: a per-CpG stochasticity atlas.

For every CpG of a cross-sectional cohort this estimates drift, entropy slope,
dispersion slope and a fitted relaxation rate, then asks whether CpGs used by
published clocks are the fast-drifting ("poorly maintained") ones, compared
only with CpGs that start at a similar level. The paper argues clock CpGs
change fastest because they matter least; this is the first test of that
claim, and the same machinery takes any gene set or cell type.
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import polars as pl

from entropy_age import plots
from entropy_age.data.upstream import UpstreamLayout
from entropy_age.experiments.common import Context, dolphin_blood, human_blood, write_table
from entropy_age.published import LinearClock
from entropy_age.stochasticity import baseline_strata, compare_groups, drift_metrics, fit_relaxation

METRICS = ("rate", "entropy_slope", "dispersion_slope")


def _atlas(matrix, ages: np.ndarray, *, t0: float) -> pl.DataFrame:
    relax = fit_relaxation(matrix.values, ages, matrix.features, t0=t0)
    drift = drift_metrics(matrix.values, ages, matrix.features)
    return relax.join(drift, on="feature")


def _within_stratum_percentile(values: np.ndarray, baseline: np.ndarray, n_bins: int = 10) -> np.ndarray:
    """Percentile of each value among features of the same baseline decile (0.5 = typical for its level)."""
    strata = baseline_strata(baseline, n_bins)
    pct = np.empty(values.size)
    for s in np.unique(strata):
        idx = np.flatnonzero(strata == s)
        pct[idx] = (np.argsort(np.argsort(values[idx])) + 0.5) / idx.size
    return pct


def run(ctx: Context) -> list[Path]:
    plots.apply_style()
    out = ctx.out_dir("atlas_real_data")
    layout = UpstreamLayout(ctx.paths)
    clock_sets = {name: set(LinearClock.read(layout.clock(name), name).features) for name in ("horvath", "phenoage")}

    # (label, matrix, ages, compare clock CpGs?) - human clock probes do not hybridise on the
    # dolphin (mammalian array) samples: they sit at ~0.5 with no variance, so no clock contrast there
    cohorts = []
    human = human_blood(ctx.paths)
    if human is not None:
        cohorts.append(("Human blood", human[0], human[1]["age"].to_numpy(), True))
    dolphin = dolphin_blood(ctx.paths)
    if dolphin is not None:
        cohorts.append(("Dolphin blood", dolphin[0], dolphin[1]["age"].to_numpy(), False))
    if not cohorts:
        print("atlas: no cohort prepared - skipped")
        return []

    fig, axes = plots.panels(3 * len(cohorts), ncols=3, size=3.3)
    letters = iter("abcdef")
    contrast_rows = []
    for c, (label, matrix, ages, compare_clocks) in enumerate(cohorts):
        atlas = _atlas(matrix, ages, t0=0.0)
        write_table(atlas, out / f"{label.split()[0].lower()}_cpg_metrics")

        ax = axes[3 * c]
        hb = ax.hexbin(atlas["baseline"], atlas["entropy_slope"], gridsize=45, bins="log", cmap="Greys", mincnt=1)
        ax.axhline(0, color=plots.ACCENT, lw=0.8)
        ax.set(
            xlabel="Baseline methylation (youngest 20 %)",
            ylabel="Entropy slope (bits / year)",
            title=f"{label}: drift vs start level",
        )
        fig.colorbar(hb, ax=ax, shrink=0.7, label="CpGs (log)")
        plots.label_panel(ax, next(letters))

        ax = axes[3 * c + 1]
        fitted = atlas.filter(
            ~pl.col("rate_at_lower_bound")
            & ~pl.col("rate_at_upper_bound")
            & pl.col("equilibrium_in_range")
            & (pl.col("r2") > 0.05)
        )
        ax.hist(np.log10(fitted["rate"].to_numpy()), bins=40, color=plots.POINT, alpha=0.8)
        ax.set(
            xlabel="log10 relaxation rate λ (per year)",
            ylabel="CpGs",
            title=f"{label}: λ identifiable for {fitted.height:,} / {atlas.height:,}",
        )
        plots.label_panel(ax, next(letters))

        ax = axes[3 * c + 2]
        if not compare_clocks:
            ax.set_visible(False)
            continue
        membership = pl.concat(
            [
                pl.DataFrame({"feature": sorted(s & set(matrix.features)), "group": name})
                for name, s in clock_sets.items()
            ]
        )
        for metric in METRICS:
            table = compare_groups(atlas, membership, metric=metric, n_permutations=2000, seed=c)
            contrast_rows.extend({"cohort": label, "metric": metric, **row} for row in table.to_dicts())
        pct = _within_stratum_percentile(atlas["entropy_slope"].to_numpy(), atlas["baseline"].to_numpy())
        features = atlas["feature"].to_numpy()
        groups = {
            "other CpGs": pct[~np.isin(features, list(clock_sets["horvath"] | clock_sets["phenoage"]))],
            "Horvath": pct[np.isin(features, list(clock_sets["horvath"]))],
            "PhenoAge": pct[np.isin(features, list(clock_sets["phenoage"]))],
        }
        frame = pl.concat(
            [
                pl.DataFrame({"group": k, "pct": v if k != "other CpGs" else ctx.rng(7).choice(v, 600, replace=False)})
                for k, v in groups.items()
            ]
        )
        plots.sweep_box(
            ax,
            frame,
            x="group",
            y="pct",
            order=list(groups),
            xlabel="",
            ylabel="Percentile within start level",
            title="Clock vs level-matched CpGs",
            ylim=(0, 1),
        )
        ax.axhline(0.5, ls=":", color=plots.MUTED, lw=1)
        plots.label_panel(ax, next(letters))

    if contrast_rows:
        write_table(pl.DataFrame(contrast_rows), out / "clock_cpgs_vs_rest")
    fig.suptitle(
        "Per-CpG stochasticity atlas from cross-sectional cohorts", x=0.02, ha="left", fontsize=12, fontweight="bold"
    )
    fig.tight_layout()
    return [plots.save(fig, ctx.figure_path("atlas_real_data"))]
