"""The opposite direction: does a Dnmt3a gain-of-function mutation move the shield? (roadmap RM12)

Tet3 heterozygous brain lost the shield; Tet1 loss strengthened it. Sarni et al. 2026 (GSE324160) measured bone
marrow, intestine and liver of mice carrying Dnmt3a W326R/+ (a gain of function that causes a progeroid
syndrome) and wild-type littermates at 4 days, 23 days, 3 months and 1 year, 3 per genotype per age, on the
mouse methylation array. Age dominates these data, so each CpG is first centred on the wild types of its age
group (``age_centred``); then, per tissue, the 12 mutants are scored against the 12 wild types as in
``consortium_interventions`` (controls leave themselves out) and the effect carries age group as a covariate.

``entropy-age fetch dnmt3a`` then ``entropy-age figures dnmt3a``.
"""

from __future__ import annotations

import gzip
import sys
from pathlib import Path

import numpy as np
import polars as pl

from entropy_age import plots
from entropy_age.data.mouse import MouseAnnotationLayout, mouse_membership
from entropy_age.experiments.bxd_liver import ols
from entropy_age.experiments.common import Context, dummies, write_table
from entropy_age.experiments.consortium_interventions import arm_scores
from entropy_age.experiments.interventions import _sets
from entropy_age.experiments.tissue_contrast import SET_LABELS, SETS

GSE = "GSE324160"
BETAS = f"{GSE}_sarni_et_al_mouse_infinium_beta.tsv.gz"
AGES = ("4-day", "23-day", "3-month", "1-year")


def layout(paths) -> tuple[Path, Path]:
    base = paths.input / "geo" / GSE
    return base / f"{GSE}_series_matrix.txt.gz", base / BETAS


def read_mice(matrix: Path) -> pl.DataFrame:
    """``column`` (the beta table's name), ``tissue``, ``age`` (group) and ``mutant`` from the sample titles."""
    with gzip.open(matrix, "rt") as handle:
        titles = next(
            [c.strip().strip('"') for c in line.rstrip("\n").split("\t")[1:]]
            for line in handle
            if line.startswith("!Sample_title")
        )
    rows = []
    for t in titles:
        column, desc = (part.strip() for part in t.split(":", 1))
        tissue = desc.split(" ")[0] if not desc.startswith("Bone marrow") else "Bone marrow"
        age = next(a for a in AGES if a in desc)
        rows.append({"column": column, "tissue": tissue, "age": age, "mutant": "W326R" in desc})
    return pl.DataFrame(rows)


def age_centred(x: np.ndarray, age: np.ndarray, control: np.ndarray) -> np.ndarray:
    """``x`` with each age group's control mean replaced by the overall control mean, CpG by CpG."""
    out = x.copy()
    overall = x[control].mean(axis=0)
    for a in np.unique(age):
        rows = age == a
        out[rows] = x[rows] - x[rows & control].mean(axis=0) + overall
    return np.clip(out, 0.0, 1.0)


def run(ctx: Context) -> list[Path]:
    matrix, betas = layout(ctx.paths)
    lay = MouseAnnotationLayout(ctx.paths)
    if not (matrix.exists() and betas.exists() and lay.probes.exists()):
        print("dnmt3a: run `entropy-age fetch dnmt3a` (and mouse-annotation) first. Nothing run.", file=sys.stderr)
        return []
    plots.apply_style()
    out = ctx.out_dir("dnmt3a_gain")
    mice = read_mice(matrix)
    probes = pl.read_parquet(lay.probes).filter(pl.col("promoter") & ~pl.col("chr").is_in(["X", "Y", "MT"]))
    frame = (
        pl.scan_csv(betas, separator="\t", null_values=["NA", ""])
        .select(["ID_REF", *mice["column"].to_list()])
        .filter(pl.col("ID_REF").is_in(probes["cpg"].implode()))
        .collect()
    )
    features = frame["ID_REF"].to_list()
    members = mouse_membership(probes.filter(pl.col("cpg").is_in(features)), pl.read_parquet(lay.gene_sets)).filter(
        pl.col("group").is_in(SETS)
    )
    island_of = dict(probes.select("cpg", "cpg_island").iter_rows())
    rows = []
    for (tissue,), sel in mice.group_by("tissue", maintain_order=True):
        x = frame.select(sel["column"].to_list()).to_numpy().T.astype(float)
        ok = ~np.isnan(x).any(axis=0)
        feats = [f for f, k in zip(features, ok, strict=True) if k]
        sets = _sets(feats, members)
        island = np.array([bool(island_of.get(f, False)) for f in feats])
        mutant, age = sel["mutant"].to_numpy(), sel["age"].to_numpy()
        xc = age_centred(x[:, ok], age, ~mutant)
        scores = arm_scores(xc, mutant, island, sets)
        design = np.column_stack([np.ones(mutant.size), mutant.astype(float), dummies(age)])
        for group in sets:
            for measure in ("B", "F", "S"):
                y = scores[f"{measure}_{group}"]
                coef, se, p = ols(y, design, 1)
                rows.append(
                    {
                        "tissue": tissue,
                        "group": group,
                        "measure": measure,
                        "effect": coef,
                        "se": se,
                        "p": p,
                        "control_mean": float(y[~mutant].mean()),
                        "mutants": int(mutant.sum()),
                        "controls": int((~mutant).sum()),
                    }
                )
        print(f"dnmt3a: {tissue} done", file=sys.stderr, flush=True)
    table = pl.DataFrame(rows)
    write_table(table, out / "effects")
    return [_figure(ctx, table)]


def _figure(ctx: Context, table: pl.DataFrame) -> Path:
    fig, axes = plots.panels(3, ncols=3, size=3.6)
    groups = [g for g in SETS if g in table["group"].unique()]
    tissues = table["tissue"].unique(maintain_order=True).to_list()
    for i, measure in enumerate(("F", "B", "S")):
        ax = axes[i]
        labels, y = [], 0
        for tissue in tissues:
            for g, group in enumerate(groups):
                r = table.filter(
                    (pl.col("tissue") == tissue) & (pl.col("group") == group) & (pl.col("measure") == measure)
                )
                ax.errorbar(r["effect"][0], y, xerr=1.96 * r["se"][0], fmt="o", ms=4, color=plots.PALETTE[g])
                labels.append((y, f"{tissue} · {SET_LABELS[group]}"))
                y += 1
            y += 0.5
        ax.axvline(0, color=plots.MUTED, lw=0.8)
        ax.set_yticks([p for p, _ in labels], [t for _, t in labels] if i == 0 else [""] * len(labels), fontsize=7)
        ax.invert_yaxis()
        ax.set(
            xlabel="W326R/+ minus wild type (95 % CI)",
            title={"F": "Shield share F (> 0: weaker)", "B": "Drift B", "S": "Signed shift S"}[measure],
        )
        plots.label_panel(ax, "abc"[i])
    fig.suptitle(
        "Dnmt3a gain of function (W326R/+), 12 vs 12 per tissue over four ages (GSE324160)",
        x=0.02,
        y=1.04,
        ha="left",
        fontsize=11,
        fontweight="bold",
    )
    return plots.save(fig, ctx.figure_path("dnmt3a_gain"))
