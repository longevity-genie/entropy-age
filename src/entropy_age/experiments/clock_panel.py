"""The clock panel in one figure: how the published clocks relate to each other, to drift and to the shield, and how
each fares beside the shield in the breast-cancer models (figures review 2026-11-24, section 2; Anton: "I want to see
multiple clocks").

Drawn from ``clock_increment``'s tables, no new fitting: ``per_woman`` (each clock's acceleration, the residual on age
among the EPIC-Italy women, standardised; B, F and P per woman) and ``nested_models``.

a. Correlation between the clocks' accelerations, clocks ordered by hierarchical clustering (average linkage on
   1 - r), a strip giving each clock's generation.
b. Each clock's acceleration against B (drift at housekeeping look-alike promoters), F (the shield share) and P (the
   absolute gap), in the same order.
c. Odds ratio per SD of each clock alone (M0) and of F once that clock and B are in the model (M2), same chip first.
d. Per-woman disagreement between the adult blood clocks against B and against F.

``entropy-age figures clock_panel`` after ``figures clock_increment``.
"""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import polars as pl
from scipy.cluster import hierarchy
from scipy.spatial.distance import squareform

from entropy_age import plots
from entropy_age.experiments.clock_increment import GENERATION, PRIMARY, B, F, P
from entropy_age.experiments.common import Context

GEN_COLOURS = {
    "1st: chronological": plots.PALETTE[0],
    "2nd: phenotypic and mortality": plots.PALETTE[1],
    "3rd: pace of ageing": plots.PALETTE[2],
    "mitotic and stochastic": plots.PALETTE[3],
    "other": "#cbd2d9",
}
DESIGN = "same chip first"


def generation(clock: str) -> str:
    return GENERATION.get(clock, "other")


def order_clocks(acc: np.ndarray) -> np.ndarray:
    r = np.corrcoef(acc, rowvar=False)
    dist = np.clip(1 - r, 0, 2)
    np.fill_diagonal(dist, 0)
    link = hierarchy.linkage(squareform(dist, checks=False), "average")
    return hierarchy.leaves_list(link)


def run(ctx: Context) -> list[Path]:
    src = ctx.paths.output / "clock_increment"
    if not (src / "per_woman.parquet").exists():
        print("clock_panel: needs figures clock_increment (per_woman)", file=sys.stderr)
        return []
    women = pl.read_parquet(src / "per_woman.parquet")
    models = pl.read_parquet(src / "nested_models.parquet")
    clocks = [c.removeprefix("acc ") for c in women.columns if c.startswith("acc ")]
    acc = women.select([f"acc {c}" for c in clocks]).to_numpy()
    order = order_clocks(acc)
    clocks = [clocks[i] for i in order]
    acc = acc[:, order]
    r = np.corrcoef(acc, rowvar=False)
    n = len(clocks)

    plots.apply_style()
    fig = plots.plt.figure(figsize=(15, 0.22 * n + 3.2))
    grid = fig.add_gridspec(1, 5, width_ratios=[0.25, n * 0.22, 2.2, 3.2, 3.4], wspace=0.08)
    ax_gen = fig.add_subplot(grid[0])
    ax_r = fig.add_subplot(grid[1], sharey=ax_gen)
    ax_d = fig.add_subplot(grid[2], sharey=ax_gen)
    ax_or = fig.add_subplot(grid[3], sharey=ax_gen)
    ax_dis = fig.add_subplot(grid[4])

    y = np.arange(n)
    ax_gen.barh(y, 1, color=[GEN_COLOURS[generation(c)] for c in clocks], height=1)
    ax_gen.set(yticks=y, yticklabels=clocks, xticks=[], ylim=(n - 0.5, -0.5))
    ax_gen.tick_params(axis="y", labelsize=7)
    ax_gen.set_title("a", loc="left", fontweight="bold")

    im = ax_r.imshow(r, cmap="RdBu_r", vmin=-1, vmax=1, aspect="auto")
    ax_r.set_xticks(y, clocks, rotation=90, fontsize=6)
    ax_r.tick_params(axis="y", left=False, labelleft=False)
    ax_r.set_title("Clocks against each other (r of accelerations)", fontsize=9)
    fig.colorbar(im, ax=ax_r, fraction=0.03, pad=0.01, shrink=0.5)

    measures = [(B, "B drift"), (F, "F shield"), (P, "P gap")]
    strip = np.array(
        [[np.corrcoef(women[f"acc {c}"].to_numpy(), women[m].to_numpy())[0, 1] for m, _ in measures] for c in clocks]
    )
    ax_d.imshow(strip, cmap="RdBu_r", vmin=-0.6, vmax=0.6, aspect="auto")
    for i in range(n):
        for j in range(len(measures)):
            ax_d.text(
                j,
                i,
                f"{strip[i, j]:+.2f}".replace("+0.", "+.").replace("-0.", "−."),
                ha="center",
                va="center",
                fontsize=5.5,
                color="white" if abs(strip[i, j]) > 0.35 else plots.INK,
            )
    ax_d.set_xticks(range(len(measures)), [name for _, name in measures], rotation=90, fontsize=7)
    ax_d.tick_params(axis="y", left=False, labelleft=False)
    ax_d.set_title("b  with drift (r)", fontsize=9, loc="left")

    prim = models.filter(pl.col("design") == DESIGN)
    rows = {row["clock"]: row for row in prim.iter_rows(named=True)}
    for i, c in enumerate(clocks):
        row = rows.get(c)
        if row is None or row.get("clock_or") is None:
            continue
        colour = GEN_COLOURS[generation(c)]
        ax_or.plot([row["clock_low"], row["clock_high"]], [i - 0.17] * 2, color=colour, lw=1)
        ax_or.plot(row["clock_or"], i - 0.17, "o", color=colour, ms=3.5)
        ax_or.plot([row["F_low"], row["F_high"]], [i + 0.17] * 2, color=plots.INK, lw=1)
        ax_or.plot(
            row["F_or_given_B"], i + 0.17, "s", color=plots.INK, ms=3, mfc="white" if c != PRIMARY else plots.INK
        )
    ax_or.axvline(1, color=plots.MUTED, lw=0.8, ls="--")
    ax_or.set_xscale("log")
    ax_or.tick_params(axis="y", left=False, labelleft=False)
    ax_or.set(xlabel="odds ratio per SD, breast cancer > 2 y after the draw")
    ax_or.set_title("c  ● clock alone   ■ shield F given clock and B", fontsize=9, loc="left")

    controls = women
    ax_dis.scatter(
        controls[B].to_numpy(),
        controls["disagreement"].to_numpy(),
        s=5,
        alpha=0.45,
        color=plots.PALETTE[1],
        label=f"B drift, r {np.corrcoef(controls[B].to_numpy(), controls['disagreement'].to_numpy())[0, 1]:+.2f}",
    )
    ax_dis.scatter(
        controls[F].to_numpy(),
        controls["disagreement"].to_numpy(),
        s=5,
        alpha=0.45,
        color=plots.PALETTE[0],
        label=f"F shield, r {np.corrcoef(controls[F].to_numpy(), controls['disagreement'].to_numpy())[0, 1]:+.2f}",
    )
    ax_dis.set(xlabel="drift measure (SD)", ylabel="disagreement between clocks (SD)")
    ax_dis.set_title("d  where clocks disagree", fontsize=9, loc="left")
    ax_dis.legend(fontsize=7, frameon=False)

    handles = [plots.plt.Rectangle((0, 0), 1, 1, color=v) for v in GEN_COLOURS.values()]
    fig.legend(
        handles, list(GEN_COLOURS), loc="lower left", ncol=5, fontsize=7, frameon=False, bbox_to_anchor=(0.02, 0)
    )
    fig.subplots_adjust(bottom=0.2, top=0.93, left=0.09, right=0.98)
    return [plots.save(fig, ctx.figure_path("clock_panel"))]
