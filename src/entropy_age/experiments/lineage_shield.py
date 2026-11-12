"""The shield by blood lineage, in silico (roadmap RM31, second part)

Sorted cells (R44) showed the shield inside monocytes and CD4 T cells at about half of whole blood's size. Here the
mixed blood of Hannum (656 people) is split by lineage without sorting, as CellDMC does for mean effects (Zheng et
al. 2018): with each person's 12 EpiDISH fractions from run R38 grouped into lineages, every CpG is fitted as

    beta = sum_k f_k (a_k + c_k age)

so ``c_k`` is the CpG's age slope within lineage k. Per lineage, essential and housekeeping CpGs are compared with
look-alikes matched on island status and level (``tissue_contrast``'s strata), on |c_k|: a shield in lineage k is a
smaller |c_k| for the set. Lineages: myeloid (neutrophils, eosinophils, basophils, monocytes), T cells (CD4 and CD8,
naive and memory, Tregs) and B and NK cells.

``entropy-age figures lineage_shield`` (after R38; seconds to minutes).
"""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import polars as pl

from entropy_age import plots
from entropy_age.data.cohorts import layout as cohort_layout
from entropy_age.experiments.common import Context, write_table
from entropy_age.experiments.tissue_contrast import _features, _matched_sd, _matrix, _with_strata

LINEAGES = {
    "myeloid": ("Neu", "Eos", "Baso", "Mono"),
    "T cells": ("CD4Tnv", "CD4Tmem", "CD8Tnv", "CD8Tmem", "Treg"),
    "B and NK": ("Bnv", "Bmem", "NK"),
}
PERMUTATIONS = 1000


def lineage_slopes(x: np.ndarray, age: np.ndarray, fractions: np.ndarray) -> np.ndarray:
    """Per lineage (columns of ``fractions``) and CpG: the age slope within that lineage (lineages x CpGs)."""
    design = np.column_stack([fractions, fractions * (age - age.mean())[:, None]])
    coef, *_ = np.linalg.lstsq(design, x, rcond=None)
    return coef[fractions.shape[1] :]


def run(ctx: Context) -> list[Path]:
    fpath = ctx.paths.output / "stochastic_component" / "cell_fractions.parquet"
    if not fpath.exists():
        print("lineage_shield: run `figures stochastic_component` (R38) first. Nothing run.", file=sys.stderr)
        return []
    plots.apply_style()
    out = ctx.out_dir("lineage_shield")
    samples = pl.read_parquet(cohort_layout("hannum", ctx.paths).samples).drop_nulls("age").sort("age", "sample")
    fr = pl.read_parquet(fpath)
    if fr.height != samples.height or not np.allclose(fr["age"].to_numpy(), samples["age"].to_numpy()):
        raise ValueError("cell fractions do not line up with the Hannum samples")
    lineage = np.column_stack([fr.select(list(cols)).sum_horizontal().to_numpy() for cols in LINEAGES.values()])
    features, probes, membership = _features(ctx.paths, ctx.rng(140_000))
    m = _matrix("hannum", ctx.paths, samples["sample"].to_list(), features)
    age = samples["age"].to_numpy()
    slopes = lineage_slopes(m.values, age, lineage)
    young = age <= np.quantile(age, 0.2)
    rows = []
    for k, name in enumerate(LINEAGES):
        metrics = pl.DataFrame(
            {"feature": m.features, "baseline": m.values[young].mean(axis=0), "lineage_abs_slope": np.abs(slopes[k])}
        )
        table = _matched_sd(_with_strata(metrics, probes), membership, ["lineage_abs_slope"], PERMUTATIONS, ctx.seed)
        rows.append(table.with_columns(lineage=pl.lit(name), mean_fraction=pl.lit(float(lineage[:, k].mean()))))
        print(f"lineage_shield: {name} done", file=sys.stderr, flush=True)
    result = pl.concat(rows)
    write_table(result, out / "gaps")
    return [_figure(ctx, result)]


def _figure(ctx: Context, table: pl.DataFrame) -> Path:
    fig, axes = plots.panels(1, ncols=1, size=4.6)
    ax = axes[0]
    groups = [g for g in ("core_essential", "housekeeping", "non_essential") if g in table["group"].unique()]
    for i, group in enumerate(groups):
        t = table.filter(pl.col("group") == group)
        ax.barh(
            np.arange(t.height) + 0.27 * i,
            t["matched_diff_sd"],
            height=0.25,
            color=plots.PALETTE[i],
            label=group.replace("_", " "),
        )
    names = table.filter(pl.col("group") == groups[0])["lineage"].to_list()
    ax.set_yticks(np.arange(len(names)) + 0.27, names)
    ax.axvline(0, color=plots.MUTED, lw=0.8)
    ax.invert_yaxis()
    ax.legend(fontsize=7)
    ax.set(xlabel="set minus look-alikes, |age slope within lineage| (SD)", title="The shield by lineage, Hannum")
    return plots.save(fig, ctx.figure_path("lineage_shield"))
