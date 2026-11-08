"""Extension, prediction 1 on real data: do essential-gene CpGs drift less in blood?

GSE41037 controls (394 people, 16-88 y, whole blood, a renewing tissue), all
27k-array CpGs. Each gene set is compared with every other CpG *at the same
start*: strata are CpG-island status (sequence-defined, so untouched by drift)
crossed with the decile of young-adult methylation. P values are stratified
permutations, BH-adjusted over every test in the figure.

GSE41037 was processed in three datasets whose ages barely overlap (dataset 1
mostly under 50, dataset 3 mostly over 50), so an unadjusted age slope partly
measures the batch. Every metric is therefore also estimated as a partial
slope within batches, and with blood-cell fractions added as covariates; the
contrast is repeated inside each batch separately.

What this can and cannot show: lower drift for essential genes in blood is what
selection predicts, but active, open promoters are also shielded by chromatin
regardless of selection. Separating the two needs the cell-type contrast
(prediction 2: the protection should vanish in post-mitotic cells).
"""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import polars as pl
from scipy.stats import false_discovery_control

from entropy_age import plots
from entropy_age.data.annotation import AnnotationLayout, probe_membership
from entropy_age.data.upstream import UpstreamLayout
from entropy_age.experiments.common import Context, dummies, human_batches, human_blood, write_table
from entropy_age.experiments.fig5_real_data import CELL_TYPES
from entropy_age.stochasticity import baseline_strata, compare_groups, drift_metrics

SETS = ("core_essential", "housekeeping", "non_essential")
SET_LABELS = {"core_essential": "core essential", "housekeeping": "housekeeping", "non_essential": "non-essential"}
METRICS = {
    "abs_slope": "|age slope|",
    "entropy_slope": "entropy slope",
    "dispersion_slope": "dispersion slope",
}
VARIANTS = ("unadjusted", "batch-adjusted", "batch + cells")
VARIANT_STYLE = {
    "unadjusted": ("o", plots.PALETTE[0]),
    "batch-adjusted": ("s", plots.PALETTE[1]),
    "batch + cells": ("D", plots.PALETTE[2]),
}
N_PERMUTATIONS = 1000


def _with_strata(metrics: pl.DataFrame, probes: pl.DataFrame) -> pl.DataFrame:
    metrics = metrics.join(
        probes.select(pl.col("cpg").alias("feature"), "symbol", "cpg_island"), on="feature", how="left"
    )
    decile = baseline_strata(metrics["baseline"].to_numpy(), 10)
    return metrics.with_columns(stratum=pl.format("{}|{}", pl.col("cpg_island").fill_null(False), pl.Series(decile)))


def _gene_metrics(cpg: pl.DataFrame) -> pl.DataFrame:
    genes = (
        cpg.filter(pl.col("symbol").is_not_null() & (pl.col("symbol") != ""))
        .group_by("symbol")
        .agg(*[pl.col(m).median() for m in METRICS], pl.col("baseline").mean(), pl.col("cpg_island").any())
        .rename({"symbol": "feature"})
    )
    decile = baseline_strata(genes["baseline"].to_numpy(), 10)
    return genes.with_columns(stratum=pl.format("{}|{}", pl.col("cpg_island"), pl.Series(decile)))


def _contrasts(metrics: pl.DataFrame, membership: pl.DataFrame, metrics_used, seed: int, **labels) -> list[dict]:
    rows = []
    for metric in metrics_used:
        table = compare_groups(
            metrics, membership, metric=metric, stratify_by="stratum", n_permutations=N_PERMUTATIONS, seed=seed
        )
        spread = float(metrics[metric].std())
        rows.extend(
            {**labels, "metric": metric, **row, "matched_diff_sd": row["matched_diff"] / spread}
            for row in table.to_dicts()
        )
    return rows


def _forest(ax, contrasts: pl.DataFrame, title: str) -> None:
    y, ticks, tick_labels = 0.0, [], []
    for group in SETS:
        for metric in METRICS:
            for k, variant in enumerate(VARIANTS):
                row = contrasts.filter(
                    (pl.col("group") == group) & (pl.col("metric") == metric) & (pl.col("variant") == variant)
                )
                if row.height == 0:
                    continue
                r = row.row(0, named=True)
                marker, color = VARIANT_STYLE[variant]
                ax.scatter(
                    r["matched_diff_sd"],
                    y + (k - 1) * 0.22,
                    marker=marker,
                    color=color,
                    s=20,
                    alpha=1.0 if r["p_bh"] < 0.05 else 0.3,
                    label=variant if (group, metric) == (SETS[0], "abs_slope") else None,
                )
            ticks.append(y)
            tick_labels.append(f"{SET_LABELS[group]} · {METRICS[metric]}")
            y += 1
        y += 0.6
    ax.axvline(0, color=plots.MUTED, lw=0.8)
    ax.set_yticks(ticks, tick_labels, fontsize=7)
    ax.invert_yaxis()
    ax.set(xlabel="Set minus matched CpGs (SD of metric)", title=title)


def run(ctx: Context) -> list[Path]:
    blood = human_blood(ctx.paths)
    annotation = AnnotationLayout(ctx.paths)
    if blood is None or not annotation.gene_sets.exists():
        print(
            "genesets: needs `entropy-age fetch geo GSE41037` and `entropy-age fetch annotation`; nothing was run",
            file=sys.stderr,
        )
        return []
    plots.apply_style()
    out = ctx.out_dir("gene_sets")
    matrix, samples = blood
    ages = samples["age"].to_numpy()
    batch = human_batches(ctx.paths, matrix.samples)
    cells = (
        pl.DataFrame({"sample": matrix.samples})
        .join(
            pl.read_parquet(UpstreamLayout(ctx.paths).gse41037_cell_types),
            on="sample",
            how="left",
            maintain_order="left",
        )
        .select(CELL_TYPES)
        .to_numpy()
    )
    covariates = {"unadjusted": None, "batch-adjusted": dummies(batch), "batch + cells": np.c_[dummies(batch), cells]}

    probes = pl.read_parquet(annotation.probes)
    membership = probe_membership(probes, pl.read_parquet(annotation.gene_sets)).filter(
        pl.col("feature").is_in(matrix.features) & pl.col("group").is_in(SETS)
    )
    gene_membership = membership.select(pl.col("gene").alias("feature"), "group").unique(maintain_order=True)
    write_table(
        membership.group_by("group").agg(pl.len().alias("cpgs"), pl.col("gene").n_unique().alias("genes")),
        out / "set_sizes",
    )

    rows, per_variant = [], {}
    for v, variant in enumerate(VARIANTS):
        cpg = _with_strata(drift_metrics(matrix.values, ages, matrix.features, covariates=covariates[variant]), probes)
        per_variant[variant] = cpg
        rows += _contrasts(cpg, membership, METRICS, 10 + v, level="CpG", variant=variant, batch="all")
        rows += _contrasts(
            _gene_metrics(cpg), gene_membership, METRICS, 20 + v, level="gene", variant=variant, batch="all"
        )
    write_table(per_variant["batch-adjusted"], out / "cpg_metrics_batch_adjusted")

    # replication inside each processing batch (no batch term needed there)
    for b, level in enumerate(np.unique(batch)):
        inside = batch == level
        cpg = _with_strata(
            drift_metrics(matrix.values[inside], ages[inside], matrix.features, covariates=cells[inside]), probes
        )
        rows += _contrasts(
            cpg,
            membership,
            ("entropy_slope", "dispersion_slope"),
            30 + b,
            level="CpG",
            variant="within batch",
            batch=str(level),
        )

    contrasts = pl.DataFrame(rows)
    contrasts = contrasts.with_columns(
        p_bh=pl.Series(false_discovery_control(contrasts["p_permutation"].fill_nan(1.0).to_numpy(), method="bh"))
    )
    write_table(contrasts, out / "contrasts")

    adjusted = per_variant["batch-adjusted"].join(membership.select("feature", "group"), on="feature", how="left")
    adjusted = adjusted.with_columns(pl.col("group").fill_null("other CpGs"))
    order = [*SETS, "other CpGs"]
    labels = [SET_LABELS.get(g, g) for g in order]

    def sample_rest(frame: pl.DataFrame) -> pl.DataFrame:
        return frame.filter(pl.col("group") != "other CpGs").vstack(
            frame.filter(pl.col("group") == "other CpGs").sample(n=1500, seed=ctx.seed)
        )

    fig, axes = plots.panels(6, ncols=3, size=3.5)

    composition = adjusted.group_by("group").agg(pl.col("cpg_island").mean().alias("cgi"))
    cgi = dict(zip(composition["group"].to_list(), composition["cgi"].to_list(), strict=True))
    plots.sweep_box(
        axes[0],
        sample_rest(adjusted),
        x="group",
        y="baseline",
        order=order,
        labels=[f"{lab}\n{100 * cgi[g]:.0f} % island" for g, lab in zip(order, labels, strict=True)],
        xlabel="",
        ylabel="Methylation, youngest 20 %",
        title="Gene sets start in different places",
    )
    plots.label_panel(axes[0], "a")

    for i, metric in enumerate(("entropy_slope", "dispersion_slope")):
        pct = adjusted.with_columns(
            pct=(pl.col(metric).rank("average").over("stratum") - 0.5) / pl.len().over("stratum")
        )
        ax = axes[1 + i]
        plots.sweep_box(
            ax,
            sample_rest(pct),
            x="group",
            y="pct",
            order=order,
            labels=labels,
            xlabel="",
            ylabel="Percentile within island × start stratum",
            ylim=(0, 1),
            title=f"{METRICS[metric].capitalize()} (batch-adj.)",
        )
        ax.axhline(0.5, ls=":", color=plots.MUTED, lw=1)
        plots.label_panel(ax, "bc"[i])

    for k, level in enumerate(("CpG", "gene")):
        _forest(
            axes[3 + k],
            contrasts.filter((pl.col("level") == level) & (pl.col("batch") == "all")),
            f"{level} level (faded: BH p ≥ 0.05)",
        )
        plots.label_panel(axes[3 + k], "de"[k])
    handles, legend_labels = axes[3].get_legend_handles_labels()
    fig.legend(handles, legend_labels, loc="lower center", ncol=3, fontsize=8, bbox_to_anchor=(0.5, -0.03))

    ax = axes[5]
    within = contrasts.filter(pl.col("variant") == "within batch")
    for k, metric in enumerate(("entropy_slope", "dispersion_slope")):
        for j, group in enumerate(SETS):
            sub = within.filter((pl.col("metric") == metric) & (pl.col("group") == group)).sort("batch")
            xs = j + (k - 0.5) * 0.3 + np.linspace(-0.08, 0.08, sub.height)
            ax.scatter(
                xs,
                sub["matched_diff_sd"],
                color=plots.PALETTE[k + 3],
                s=22,
                alpha=[1.0 if p < 0.05 else 0.3 for p in sub["p_bh"]],
                label=METRICS[metric] if j == 0 else None,
            )
    ax.axhline(0, color=plots.MUTED, lw=0.8)
    ax.set_xticks(range(len(SETS)), [SET_LABELS[g] for g in SETS], rotation=20, ha="right")
    ax.set(ylabel="Set minus matched CpGs (SD)", title="Inside each of the 3 batches")
    ax.legend(fontsize=7, loc="upper left")
    plots.label_panel(ax, "f")

    fig.suptitle(
        "Do essential-gene CpGs drift less in blood? (GSE41037, 394 people, 27k CpGs)",
        x=0.02,
        ha="left",
        fontsize=12,
        fontweight="bold",
    )
    fig.tight_layout()
    return [plots.save(fig, ctx.figure_path("gene_sets"))]
