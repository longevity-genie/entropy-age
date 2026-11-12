"""A second methylation-machinery model: does losing Tet1 remove the shield in cortex? (roadmap RM12)

Tet3 heterozygous mice lost the shield in cortex and striatum against their littermates (runs R15, R16). Foong
et al. 2024 (GSE271375) measured adult cortex of Tet1 knockout (``KO``), Tet1 catalytic-dead (``Hxd``) and wild
type mice, 8 each, on the Infinium mouse array (MM285), with bisulfite (``BS``: 5mC + 5hmC) and bACE (5hmC) per
animal. The bisulfite arrays are used, as everywhere else in this project. Each mutant arm is scored against the
wild types as in ``consortium_interventions`` (controls leave themselves out): drift ``B``, shield share ``F``,
signed shift ``S``; raw and with 3 composition components.

``entropy-age fetch tet1`` then ``entropy-age figures tet1``.
"""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import polars as pl

from entropy_age import plots
from entropy_age.data.mouse import MouseAnnotationLayout, mouse_membership
from entropy_age.experiments.bxd_liver import ols
from entropy_age.experiments.common import Context, write_table
from entropy_age.experiments.consortium_interventions import arm_scores
from entropy_age.experiments.interventions import _sets
from entropy_age.experiments.tissue_contrast import SET_LABELS, SETS

GSE = "GSE271375"
BETAS = "Tet1_array_betas.csv.gz"
ARMS = ("KO", "Hxd")
CONTROL = "WT"
COMPONENTS = (0, 3)


def layout(paths) -> Path:
    return paths.input / "geo" / GSE / BETAS


def wide_bisulfite(long: pl.DataFrame) -> pl.DataFrame:
    """The long ``ProbeID, Betas, Sample`` table as probes x bisulfite samples (``WT1_BS`` -> ``WT1``)."""
    bs = long.filter(pl.col("Sample").str.ends_with("_BS")).with_columns(
        animal=pl.col("Sample").str.strip_suffix("_BS")
    )
    return bs.pivot(on="animal", index="ProbeID", values="Betas").rename({"ProbeID": "cpg"})


def genotype(animal: str) -> str:
    """``WT3`` -> ``WT``, ``Hxd7`` -> ``Hxd``."""
    return animal.rstrip("0123456789")


def run(ctx: Context) -> list[Path]:
    path = layout(ctx.paths)
    lay = MouseAnnotationLayout(ctx.paths)
    if not (path.exists() and lay.probes.exists()):
        print("tet1: run `entropy-age fetch tet1` (and mouse-annotation) first. Nothing run.", file=sys.stderr)
        return []
    plots.apply_style()
    out = ctx.out_dir("tet1_cortex")
    probes = pl.read_parquet(lay.probes).filter(pl.col("promoter") & ~pl.col("chr").is_in(["X", "Y", "MT"]))
    wide = wide_bisulfite(pl.read_csv(path, null_values=["NA"], schema_overrides={"Betas": pl.Float64})).filter(
        pl.col("cpg").is_in(probes["cpg"].implode())
    )
    animals = [c for c in wide.columns if c != "cpg"]
    features = wide["cpg"].to_list()
    x_all = wide.select(animals).to_numpy().T.astype(float)
    ok = ~np.isnan(x_all).any(axis=0)
    x_all, features = x_all[:, ok], [f for f, k in zip(features, ok, strict=True) if k]
    members = mouse_membership(probes.filter(pl.col("cpg").is_in(features)), pl.read_parquet(lay.gene_sets)).filter(
        pl.col("group").is_in(SETS)
    )
    island_of = dict(probes.select("cpg", "cpg_island").iter_rows())
    island = np.array([bool(island_of.get(f, False)) for f in features])
    sets = _sets(features, members)
    geno = np.array([genotype(a) for a in animals])
    print(f"tet1: {len(animals)} animals, {len(features)} promoter CpGs", file=sys.stderr, flush=True)
    rows = []
    for arm in ARMS:
        keep = (geno == CONTROL) | (geno == arm)
        treated = geno[keep] == arm
        design = np.column_stack([np.ones(treated.size), treated.astype(float)])
        for k in COMPONENTS:
            scores = arm_scores(x_all[keep], treated, island, sets, components=k)
            for group in sets:
                for measure in ("B", "F", "S"):
                    y = scores[f"{measure}_{group}"]
                    coef, se, p = ols(y, design, 1)
                    rows.append(
                        {
                            "arm": arm,
                            "components": k,
                            "group": group,
                            "measure": measure,
                            "effect": coef,
                            "se": se,
                            "p": p,
                            "control_mean": float(y[~treated].mean()),
                            "treated_n": int(treated.sum()),
                            "control_n": int((~treated).sum()),
                        }
                    )
    table = pl.DataFrame(rows)
    write_table(table, out / "effects")
    return [_figure(ctx, table)]


def _figure(ctx: Context, table: pl.DataFrame) -> Path:
    fig, axes = plots.panels(3, ncols=3, size=3.8)
    groups = [g for g in SETS if g in table["group"].unique()]
    for i, measure in enumerate(("F", "B", "S")):
        ax = axes[i]
        labels, y = [], 0
        for arm in ARMS:
            for group in groups:
                for j, k in enumerate(COMPONENTS):
                    r = table.filter(
                        (pl.col("arm") == arm)
                        & (pl.col("group") == group)
                        & (pl.col("measure") == measure)
                        & (pl.col("components") == k)
                    )
                    ax.errorbar(
                        r["effect"][0],
                        y + 0.25 * j,
                        xerr=1.96 * r["se"][0],
                        fmt="os"[j],
                        ms=4,
                        color=plots.PALETTE[groups.index(group)],
                        alpha=1.0 if k else 0.5,
                    )
                labels.append((y + 0.12, f"{arm} · {SET_LABELS[group]}"))
                y += 1
            y += 0.5
        ax.axvline(0, color=plots.MUTED, lw=0.8)
        ax.set_yticks([p for p, _ in labels], [t for _, t in labels] if i == 0 else [""] * len(labels), fontsize=7)
        ax.invert_yaxis()
        ax.set(
            xlabel="Mutant minus wild type (95 % CI)",
            title={"F": "Shield share F (> 0: weaker)", "B": "Drift B", "S": "Signed shift S"}[measure],
        )
        plots.label_panel(ax, "abc"[i])
    fig.suptitle(
        "Tet1 knockout and catalytic-dead cortex (GSE271375), 8 vs 8: circles raw, squares 3 components",
        x=0.02,
        y=1.04,
        ha="left",
        fontsize=11,
        fontweight="bold",
    )
    return plots.save(fig, ctx.figure_path("tet1_cortex"))
