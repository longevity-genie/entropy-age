"""Does partial reprogramming restore the shield? (roadmap RM17)

Browder et al. 2022 (GSE190665) induced the four Yamaguchi factors (OSKM, "4F") cyclically in vivo with
doxycycline for 1, 7 or 10 months and measured six tissues (kidney, liver, lung, muscle, skin, spleen) on the
mammalian array, with untreated old 4F mice, old wild types and young 4F mice as references. Per tissue, the
treated mice (all three durations pooled; ``TREATED``) are scored against the untreated old mice (``OLD``) as in
``consortium_interventions``: controls leave themselves out, drift ``B``, the shield as a share of drift ``F``,
the signed shift ``S``. The young mice against the same old controls are the positive control (``YOUNG``).

Chromatin maintenance predicts a reset that restores the shield (``F`` more negative, ``B`` lower in treated);
selection predicts nothing specific, since no cells are removed.

``entropy-age fetch reprogramming`` then ``entropy-age figures reprogramming``.
"""

from __future__ import annotations

import gzip
import sys
from pathlib import Path

import numpy as np
import polars as pl

from entropy_age import plots
from entropy_age.data.cohorts import _characteristics
from entropy_age.data.mouse import MouseAnnotationLayout, mouse_membership
from entropy_age.experiments.bxd_liver import ols, probes_without_suffix
from entropy_age.experiments.common import Context, dummies, write_table
from entropy_age.experiments.consortium_interventions import arm_scores
from entropy_age.experiments.interventions import _sets
from entropy_age.experiments.tissue_contrast import SET_LABELS, SETS

GSE = "GSE190665"
BETAS = f"{GSE}_datBetaNormalized.csv.gz"
TREATED = {"2.4F+1mDox", "3.4F+7mDox", "4.4F+10mDox"}
OLD = {"5.4F-old", "6.B6-old"}
YOUNG = {"7.4F-young"}
MIN_ARM = 4


def layout(paths) -> tuple[Path, Path]:
    base = paths.input / "geo" / GSE
    return base / f"{GSE}_series_matrix.txt.gz", base / BETAS


def read_mice(matrix: Path) -> pl.DataFrame:
    """``sample``, ``barcode`` (the beta column, ``X`` + array position), ``tissue``, ``group``, ``sex``, ``age``."""
    chars = _characteristics(matrix)
    titles = accessions = sources = None
    with gzip.open(matrix, "rt") as handle:
        for line in handle:
            cells = [c.strip().strip('"') for c in line.rstrip("\n").split("\t")[1:]]
            if line.startswith("!Sample_title"):
                titles = cells
            elif line.startswith("!Sample_geo_accession"):
                accessions = cells
            elif line.startswith("!Sample_source_name_ch1"):
                sources = cells
    frame = pl.DataFrame({"sample": accessions, "title": titles, "tissue": sources}).with_columns(
        barcode=pl.lit("X") + pl.col("title").str.extract(r"(\d+_R\d+C\d+)$")
    )
    return frame.join(
        chars.select("sample", "group", pl.col("Sex").alias("sex"), pl.col("age").cast(pl.Float64, strict=False)),
        on="sample",
    )


def arm_of(group: str) -> str | None:
    if group in TREATED:
        return "treated"
    if group in OLD:
        return "old"
    if group in YOUNG:
        return "young"
    return None


def run(ctx: Context) -> list[Path]:
    matrix, betas = layout(ctx.paths)
    lay = MouseAnnotationLayout(ctx.paths)
    if not (matrix.exists() and betas.exists() and lay.probes.exists()):
        print("reprogramming: run `entropy-age fetch reprogramming` (and mouse-annotation) first.", file=sys.stderr)
        return []
    plots.apply_style()
    out = ctx.out_dir("reprogramming")
    mice = read_mice(matrix).with_columns(arm=pl.col("group").map_elements(arm_of, return_dtype=pl.String))
    mice = mice.filter(pl.col("arm").is_not_null())
    probes = probes_without_suffix(
        pl.read_parquet(lay.probes).filter(pl.col("promoter") & ~pl.col("chr").is_in(["X", "Y", "MT"]))
    )
    header = pl.read_csv(betas, n_rows=0).columns
    mice = mice.filter(pl.col("barcode").is_in(header))
    frame = (
        pl.scan_csv(betas, null_values=["NA", ""])
        .select([header[0], *mice["barcode"].to_list()])
        .filter(pl.col(header[0]).is_in(probes["cpg"].implode()))
        .collect()
        .rename({header[0]: "cpg"})
    )
    features = frame["cpg"].to_list()
    members = mouse_membership(probes.filter(pl.col("cpg").is_in(features)), pl.read_parquet(lay.gene_sets)).filter(
        pl.col("group").is_in(SETS)
    )
    island_of = dict(probes.select("cpg", "cpg_island").iter_rows())
    island = np.array([bool(island_of.get(f, False)) for f in features])
    sets = _sets(features, members)
    rows = []
    for (tissue,), sel in mice.group_by("tissue", maintain_order=True):
        x_all = frame.select(sel["barcode"].to_list()).to_numpy().T.astype(float)
        ok = ~np.isnan(x_all).any(axis=0)
        x_all = x_all[:, ok]
        sub_sets = {
            g: (s[ok], [np.searchsorted(np.flatnonzero(ok), gi[ok[gi]]) for gi in genes if ok[gi].any()])
            for g, (s, genes) in sets.items()
        }
        is_old = (sel["arm"] == "old").to_numpy()
        for arm in ("treated", "young"):
            is_arm = (sel["arm"] == arm).to_numpy()
            if is_arm.sum() < MIN_ARM or is_old.sum() < MIN_ARM:
                continue
            keep = is_old | is_arm
            treated = is_arm[keep]
            scores = arm_scores(x_all[keep], treated, island[ok], sub_sets)
            design = np.column_stack(
                [np.ones(treated.size), treated.astype(float), dummies(sel["sex"].to_numpy()[keep])]
            )
            for group in sub_sets:
                for measure in ("B", "F", "S"):
                    y = scores[f"{measure}_{group}"]
                    coef, se, p = ols(y, design, 1)
                    rows.append(
                        {
                            "tissue": tissue,
                            "arm": arm,
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
        print(f"reprogramming: {tissue} done", file=sys.stderr, flush=True)
    table = pl.DataFrame(rows)
    write_table(table, out / "effects")
    return [_figure(ctx, table)]


def _figure(ctx: Context, table: pl.DataFrame) -> Path:
    fig, axes = plots.panels(2, ncols=2, size=4.4)
    for i, measure in enumerate(("F", "B")):
        ax = axes[i]
        sub = table.filter((pl.col("measure") == measure) & (pl.col("group") == "housekeeping")).sort("tissue", "arm")
        labels = [
            f"{r['tissue']} · {r['arm']} ({r['treated_n']} vs {r['control_n']})" for r in sub.iter_rows(named=True)
        ]
        colors = [plots.PALETTE[0] if a == "treated" else plots.MUTED for a in sub["arm"]]
        ys = np.arange(sub.height)
        ax.errorbar(sub["effect"], ys, xerr=1.96 * sub["se"].to_numpy(), fmt="none", ecolor=colors)
        ax.scatter(sub["effect"], ys, c=colors, s=18)
        ax.axvline(0, color=plots.MUTED, lw=0.8)
        ax.set_yticks(ys, labels if i == 0 else [""] * len(labels), fontsize=7)
        ax.invert_yaxis()
        ax.set(
            xlabel="Against untreated old mice (95 % CI)",
            title={"F": "Shield share F (< 0: stronger shield)", "B": "Drift B (< 0: younger)"}[measure],
        )
        plots.label_panel(ax, "ab"[i])
    fig.suptitle(
        f"Cyclic partial reprogramming in vivo (GSE190665), {SET_LABELS['housekeeping']}: treated (blue), young (grey)",
        x=0.02,
        y=1.04,
        ha="left",
        fontsize=11,
        fontweight="bold",
    )
    return plots.save(fig, ctx.figure_path("reprogramming"))
