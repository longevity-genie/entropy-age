"""Which tissues age fastest, and is the order shared across species and layers? (roadmap RM25)

From docs/BEYOND_SHIELD.md (idea 5), asked by antonkulaga (2026-11-12) for questions beyond the shield. Reads
result tables earlier runs wrote; computes nothing new from raw data. A tissue's pace is the median per-feature age
slope of the features outside the gene sets (``median_rest`` of the housekeeping contrast), in two kinds:

- ``abs_slope``: how far features move in mean with age (directed change);
- ``dispersion_slope``: how much more variable they become between individuals (noise).

Tissue orders are compared by Spearman correlation over shared tissues: mouse against rat RNA (``rna_rodents``,
one atlas), human methylation against human RNA (GTEx, same tissues), human against mouse methylation. Units differ
(per year in people, per month in rodents), so only orders are compared.

``entropy-age figures tissue_pace`` (seconds; needs the tables of gtex_tissues, mouse_tissues, rna_gtex,
rna_rodents).
"""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import polars as pl
from scipy.stats import spearmanr

from entropy_age import plots
from entropy_age.experiments.common import Context, write_table

KINDS = ("abs_slope", "dispersion_slope")
#: GTEx RNA tissue names -> the methylation tables' names
GTEX_RNA = {
    "Whole Blood": "whole blood",
    "Colon - Transverse": "colon",
    "Testis": "testis",
    "Kidney - Cortex": "kidney",
    "Prostate": "prostate",
    "Ovary": "ovary",
    "Lung": "lung",
    "Breast - Mammary Tissue": "breast",
    "Muscle - Skeletal": "skeletal muscle",
}
#: mouse methylation tissue names -> human methylation names
MOUSE_METH = {
    "white blood cells": "whole blood",
    "colon": "colon",
    "kidney": "kidney",
    "skeletal muscle": "skeletal muscle",
}
SOURCES = {
    "human methylation": ("gtex_tissues", {"adjustment": "covariates"}),
    "mouse methylation": ("mouse_tissues", {"measure": "aging"}),
    "human RNA": ("rna_gtex", {"adjustment": "covariates"}),
    "mouse RNA": ("rna_rodents", {"adjustment": "covariates", "species": "mouse"}),
    "rat RNA": ("rna_rodents", {"adjustment": "covariates", "species": "rat"}),
}
PAIRS = (("mouse RNA", "rat RNA"), ("human methylation", "human RNA"), ("human methylation", "mouse methylation"))


def paces(table: pl.DataFrame, where: dict[str, str]) -> pl.DataFrame:
    """``tissue``, ``kind``, ``pace``: the median slope of features outside the sets, per tissue and kind."""
    t = table.filter(pl.col("metric").is_in(KINDS) & (pl.col("group") == "housekeeping"))
    for column, value in where.items():
        t = t.filter(pl.col(column) == value)
    return t.select("tissue", pl.col("metric").alias("kind"), pl.col("median_rest").alias("pace"))


def agreement(a: pl.DataFrame, b: pl.DataFrame) -> list[dict]:
    """Per kind: Spearman correlation of the two tissue orders over the shared tissues."""
    rows = []
    for kind in KINDS:
        j = a.filter(pl.col("kind") == kind).join(b.filter(pl.col("kind") == kind), on="tissue", suffix="_b")
        if j.height < 4:
            continue
        rho, p = spearmanr(j["pace"], j["pace_b"])
        rows.append({"kind": kind, "tissues": j.height, "rho": float(rho), "p": float(p)})
    return rows


def _rename(frame: pl.DataFrame, mapping: dict[str, str]) -> pl.DataFrame:
    return frame.filter(pl.col("tissue").is_in(list(mapping))).with_columns(pl.col("tissue").replace_strict(mapping))


def run(ctx: Context) -> list[Path]:
    tables = {}
    for name, (folder, where) in SOURCES.items():
        path = ctx.paths.output / folder / "point_estimates.parquet"
        if not path.exists():
            print(f"tissue_pace: {path} missing; run `figures {folder}` first. Nothing run.", file=sys.stderr)
            return []
        tables[name] = paces(pl.read_parquet(path), where)
    plots.apply_style()
    out = ctx.out_dir("tissue_pace")
    aligned = dict(tables)
    aligned["human RNA"] = _rename(tables["human RNA"], GTEX_RNA)
    aligned["mouse methylation"] = _rename(tables["mouse methylation"], MOUSE_METH)
    rows = [{"a": a, "b": b, **r} for a, b in PAIRS for r in agreement(aligned[a], aligned[b])]
    table = pl.DataFrame(rows)
    write_table(table, out / "agreement")
    write_table(pl.concat([t.with_columns(source=pl.lit(n)) for n, t in tables.items()]), out / "paces")
    return [_figure(ctx, aligned, table)]


def _figure(ctx: Context, aligned: dict[str, pl.DataFrame], table: pl.DataFrame) -> Path:
    fig, axes = plots.panels(4, ncols=4, size=3.4)
    for i, kind in enumerate(KINDS):
        for k, (a, b) in enumerate(PAIRS[:2]):
            ax = axes[2 * i + k]
            j = (
                aligned[a]
                .filter(pl.col("kind") == kind)
                .join(aligned[b].filter(pl.col("kind") == kind), on="tissue", suffix="_b")
            )
            ax.scatter(j["pace"], j["pace_b"], s=14, color=plots.PALETTE[i])
            ax.locator_params(nbins=4)
            for t, x, y in j.select("tissue", "pace", "pace_b").iter_rows():
                ax.annotate(t, (x, y), fontsize=5.5, alpha=0.8)
            r = table.filter((pl.col("a") == a) & (pl.col("b") == b) & (pl.col("kind") == kind))
            rho = r["rho"][0] if r.height else np.nan
            ax.set(xlabel=a, ylabel=b, title=f"{kind.replace('_', ' ')}: rho {rho:.2f}")
            plots.label_panel(ax, "abcd"[2 * i + k])
    fig.suptitle(
        "Tissue aging pace: directed change (a, b) and noise (c, d), median slope outside the gene sets",
        x=0.02,
        y=1.04,
        ha="left",
        fontsize=11,
        fontweight="bold",
    )
    return plots.save(fig, ctx.figure_path("tissue_pace"))
