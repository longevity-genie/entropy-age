"""Two main figures drawn from tables other experiments wrote: the shield follows expression measured in the same
samples, and the shield goes with cell division. Nothing is computed here beyond reading those tables; every value is
the housekeeping gap on relative entropy (``kl_slope``) with the interval its experiment gave.

``entropy-age figures headline`` (needs the outputs of mesa_expression, sturm_expression, epidermis_expression,
jena_rrbs, nuclei_paired, neuron_series, bj_telomerase, astrocyte_passage, astrocyte_htert and arrest_timecourse; a
missing table drops its rows and is named on stderr).
"""

from __future__ import annotations

import sys
from pathlib import Path

import matplotlib.pyplot as plt
import polars as pl

from entropy_age import plots
from entropy_age.experiments.common import Context

WIDTH = 7.0  # inches: the printed width of a two-column figure, so the font sizes below are the printed ones
BEFORE, AFTER = "island x level", "+ own expression decile"
EXPRESSION = (  # label, folder, filter column, filter value
    ("Monocytes, 1,202 people (age)", "mesa_expression", "cell", "CD14+ monocytes"),
    ("Fibroblasts in culture (doublings)", "sturm_expression_rsesame", None, None),
    ("Epidermis, 61 women, EPIC (age)", "epidermis_expression", "arm", "EPIC"),
    ("Epidermis, 56 women, 450k (age)", "epidermis_expression", "arm", "450k"),
)


def _read(ctx: Context, folder: str, name: str) -> pl.DataFrame | None:
    path = ctx.paths.output / folder / f"{name}.csv"
    if not path.exists():
        print(f"headline: {path} missing; its rows are left out", file=sys.stderr)
        return None
    return pl.read_csv(path)


def _housekeeping(frame: pl.DataFrame, column: str | None, value: str | None) -> pl.DataFrame:
    frame = frame.filter((pl.col("metric") == "kl_slope") & (pl.col("group") == "housekeeping"))
    return frame if column is None else frame.filter(pl.col(column) == value)


def expression_rows(ctx: Context) -> tuple[list[dict], dict[str, pl.DataFrame]]:
    """Per series: the gap before and after the expression decile joins the matching, and the gradient table."""
    rows, gradients = [], {}
    for label, folder, column, value in EXPRESSION:
        gaps = _read(ctx, folder, "matched_gaps")
        if gaps is None:
            continue
        gaps = _housekeeping(gaps, column, value)
        pick = {m: gaps.filter(pl.col("matching") == m).row(0, named=True) for m in (BEFORE, AFTER)}
        rows.append({"label": label, "before": pick[BEFORE], "after": pick[AFTER]})
        gradient = _read(ctx, folder, "gradient")
        if gradient is not None:
            gradient = gradient.filter(pl.col("metric") == "kl_slope")
            gradients[label] = gradient if column is None else gradient.filter(pl.col(column) == value)
    return rows, gradients


def division_rows(ctx: Context) -> list[dict]:
    """Block, label, gap and interval of every division comparison, in reading order."""
    rows: list[dict] = []

    def add(block: str, label: str, gap: float, low: float, high: float, divides: bool) -> None:
        # intervals resample animals, donors or cultures; untreated fibroblasts and the confluent cultures keep the
        # gene-resampling interval (no culture bootstrap was run for them), as the caption says
        rows.append({"block": block, "label": label, "gap": gap, "low": low, "high": high, "divides": divides})

    jena = _read(ctx, "jena_rrbs", "gaps")
    if jena is not None:
        for tissue, divides in (("colon", True), ("hippocampus", False)):
            r = _housekeeping(jena, "tissue", tissue).row(0, named=True)
            add("Mouse tissues, one cohort", tissue.capitalize(), r["gap_sd"], r["ci_low"], r["ci_high"], divides)
    paired = _read(ctx, "nuclei_paired", "paired")
    if paired is not None:
        r = _housekeeping(paired, "glia", "other glia (NeuN-/SOX10-)").row(0, named=True)
        add("Human cortex, sorted nuclei", "Other glia (NeuN-/SOX10-)", r["gap_glia"], r["gap_glia_low"], r["gap_glia_high"], True)
        add(
            "Human cortex, sorted nuclei",
            "Neurons, same donors",
            r["gap_neurons"],
            r["gap_neurons_low"],
            r["gap_neurons_high"],
            False,
        )
    neurons = _read(ctx, "series_neurons", "gaps")  # the shared runner: donor-resampling interval
    if neurons is not None:
        r = neurons.filter(
            (pl.col("metric") == "kl_slope")
            & (pl.col("group") == "housekeeping")
            & (pl.col("band") == "all")
            & (pl.col("covariates") == "held + chip")
            & (pl.col("cohort") == "neurons, every donor")
        ).row(0, named=True)
        add(
            "Human cortex, sorted nuclei",
            "Neurons, second series",
            r["gap_sd"],
            r["donor_ci_low"],
            r["donor_ci_high"],
            False,
        )
    wgbs = _read(ctx, "wgbs_nuclei", "paired")
    if wgbs is not None:
        r = _housekeeping(wgbs, "donors_used", "every donor (diagnosis held)").row(0, named=True)
        block = "Human cortex, sorted nuclei, by sequencing"
        add(block, "Oligodendrocyte lineage", r["gap_oligodendrocytes"], r["gap_oligodendrocytes_low"],
            r["gap_oligodendrocytes_high"], True)  # fmt: skip
        add(block, "Neurons, same donors", r["gap_neurons"], r["gap_neurons_low"], r["gap_neurons_high"], False)
    culture = "Human cells in culture, per doubling"
    for label, folder, column, value in (
        ("Fibroblasts", "sturm_expression_rsesame", "matching", BEFORE),
        ("Fibroblasts with telomerase", "bj_telomerase", "cells", "hTERT"),
        ("Astrocytes", "astrocyte_passage", None, None),
        ("Astrocytes with telomerase", "astrocyte_htert", None, None),
    ):
        table = _read(ctx, folder, "matched_gaps")
        if table is not None:
            r = _housekeeping(table, column, value).row(0, named=True)
            low, high = r["ci_low"], r["ci_high"]
            boot = _read(ctx, folder, "culture_bootstrap")  # cultures resampled within line, where it was run
            if boot is not None:
                b = _housekeeping(boot, column if column == "cells" else None, value).row(0, named=True)
                low, high = b["ci_low"], b["ci_high"]
            add(culture, label, r["gap_sd"], low, high, True)
    arrest = _read(ctx, "arrest_timecourse", "gaps")
    if arrest is not None:
        for donor in ("Donor #1", "Donor #4"):  # the two donors whose arrays are not ordered by time along the chip
            r = _housekeeping(arrest, "donor", donor).row(0, named=True)
            add(
                "Human fibroblasts kept from dividing, per month",
                f"Confluent, {donor.lower()}",
                r["gap_sd"],
                r["ci_low"],
                r["ci_high"],
                False,
            )
    return rows


def _style() -> None:
    plots.apply_style()
    plt.rcParams.update(
        {"font.size": 7.5, "axes.labelsize": 7.5, "xtick.labelsize": 7, "ytick.labelsize": 7, "legend.fontsize": 7}
    )
    plt.rcParams["savefig.dpi"] = 400


def expression_figure(ctx: Context) -> Path | None:
    rows, gradients = expression_rows(ctx)
    if not rows:
        return None
    _style()
    fig, (left, right) = plt.subplots(1, 2, figsize=(WIDTH, 2.9), gridspec_kw={"width_ratios": [1.15, 1]})
    for i, row in enumerate(rows):
        y = len(rows) - 1 - i
        for offset, key, color, name in (
            (0.14, "before", plots.POINT, "island and level"),
            (-0.14, "after", plots.ACCENT, "plus expression decile"),
        ):
            r = row[key]
            left.errorbar(
                r["gap_sd"],
                y + offset,
                xerr=[[r["gap_sd"] - r["ci_low"]], [r["ci_high"] - r["gap_sd"]]],
                fmt="o",
                ms=3.5,
                lw=1,
                color=color,
                label=name if i == 0 else None,
            )
        left.text(
            0.004,
            y,
            f"{100 * row['after']['share_removed']:.0f} % removed",
            va="center",
            fontsize=6.5,
            color=plots.MUTED,
        )
    left.axvline(0, color=plots.MUTED, lw=0.6)
    left.set_yticks(range(len(rows)), [r["label"] for r in reversed(rows)])
    left.set_ylim(-0.6, len(rows) - 0.4)
    left.set_xlim(right=0.05)
    left.set_xlabel("Housekeeping minus matched CpGs (SD)")
    left.set_xticks([-0.10, -0.05, 0.0])
    left.legend(title="CpGs matched on", title_fontsize=7, loc="lower center", bbox_to_anchor=(0.45, 1.0), ncol=2)
    for (label, table), color in zip(gradients.items(), plots.PALETTE, strict=False):
        table = table.sort("bin")
        right.errorbar(
            table["bin"], table["relative_slope_sd"], yerr=1.96 * table["se"], fmt="o-", ms=2.5, lw=0.8, color=color,
            label=label.split(",")[0] + (", 450k" if "450k" in label else ", EPIC" if "EPIC" in label else ""),
        )  # fmt: skip
    right.axhline(0, color=plots.MUTED, lw=0.6)
    right.set_xticks(range(1, 11))
    right.set_xlabel("Expression decile in the same samples (1 = lowest)")
    right.set_ylabel("Drift relative to matched CpGs (SD)")
    right.legend(loc="lower left")
    for ax, letter in ((left, "a"), (right, "b")):
        ax.text(-0.02, 1.16, letter, transform=ax.transAxes, fontsize=10, fontweight="bold", color=plots.INK)
    fig.tight_layout()
    return plots.save(fig, ctx.paths.figures / "expression_shield.png")


def division_figure(ctx: Context) -> Path | None:
    rows = division_rows(ctx)
    if not rows:
        return None
    _style()
    fig, ax = plt.subplots(figsize=(WIDTH * 0.72, 0.21 * len(rows) + 1.5))
    y, ticks, labels, block = 0.0, [], [], None
    for row in rows:
        if row["block"] != block:
            block = row["block"]
            y -= 0.9
            ax.text(-0.335, y, block, fontsize=7.5, fontweight="bold", color=plots.INK, va="center", ha="left")
            y -= 0.9
        color = plots.PALETTE[0] if row["divides"] else plots.ACCENT
        ax.errorbar(
            row["gap"],
            y,
            xerr=[[row["gap"] - row["low"]], [row["high"] - row["gap"]]],
            fmt="o",
            ms=3.5,
            lw=1,
            color=color,
        )
        ticks.append(y)
        labels.append(row["label"])
        y -= 1.0
    ax.axvline(0, color=plots.MUTED, lw=0.6)
    ax.set_yticks(ticks, labels)
    ax.set_xlim(-0.34, 0.06)
    ax.set_ylim(y, 0.2)
    ax.set_xlabel("Housekeeping minus matched CpGs (SD); below zero = shielded")
    ax.plot([], [], "o", ms=3.5, color=plots.PALETTE[0], label="cells divide")
    ax.plot([], [], "o", ms=3.5, color=plots.ACCENT, label="cells do not divide")
    ax.legend(loc="upper right")
    fig.tight_layout()
    return plots.save(fig, ctx.paths.figures / "division_shield.png")


def run(ctx: Context) -> list[Path]:
    return [p for p in (expression_figure(ctx), division_figure(ctx)) if p is not None]
