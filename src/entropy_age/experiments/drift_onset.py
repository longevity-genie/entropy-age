"""When does drift begin? A change point in methylation noise against age at sexual maturity (roadmap RM29)

Idea 6 of docs/BEYOND_SHIELD.md; open in docs/LITERATURE_COVERAGE.md (the drift rate scales with lifespan,
Bertucci-Richter 2023, and clocks scale with maturity, Lu 2023, but no one has looked for a change point). If
noise accumulates from birth at a steady rate, the noise of an animal grows linearly with age; if maintenance
holds until reproduction and lapses after (the disposable-soma reading), it stays flat until about sexual maturity
and rises after.

Data: the Mammalian Methylation Consortium (GSE223748), which gives each species' age at sexual maturity and
maximum lifespan. For every species and tissue with animals both before maturity and well after it:

1. **Noise per animal, development removed.** Each promoter CpG's directed change with age is fitted with a cubic
   spline in log age, 10 knots evenly spaced (``trend_residuals``); residuals are rescaled by their leverage, so
   that where the spline bends to few young animals their noise is not understated. What is left is noise. An
   animal's noise ``D`` is the mean over CpGs of |residual| / sqrt(m (1 - m)), separately for look-alike
   (background) CpGs and for essential-gene CpGs (does the shield begin at the same time?). A test checks that
   steep development with steady noise yields no bend.
2. **The change point.** ``D ~ a + b1 age + b2 max(0, age - tau)``, with ``tau`` on a grid of the animals' ages
   (``hinge_fit``); the slope before ``tau`` is ``b1``, after it ``b1 + b2``; ``delta_bic`` against a straight
   line says whether a bend is supported at all.
3. **Across species.** ``log tau`` regressed on log age at maturity, and separately on log lifespan, with tissue
   held: a slope near 1 on maturity, with the onset near maturity, says noise starts when reproduction does.

``entropy-age figures drift_onset`` (as a Slurm job; the consortium's beta table, about 4.5 GB).
"""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import polars as pl

from entropy_age import plots
from entropy_age.data.cohorts import _characteristics
from entropy_age.data.mouse import MouseAnnotationLayout, mouse_membership
from entropy_age.experiments.bxd_liver import ols, probes_without_suffix
from entropy_age.experiments.common import Context, dummies, write_table
from entropy_age.experiments.species_shield import layout, read_samples
from entropy_age.experiments.tissue_contrast import SETS

MIN_ANIMALS = 20
MIN_SIDE = 4  # animals before maturity, and after twice maturity
GRID = np.linspace(0.1, 0.9, 33)  # quantiles of the group's ages tried as the change point
SET = "housekeeping"


def trend_residuals(x: np.ndarray, age: np.ndarray, maturity: float, knots: int = 10) -> np.ndarray:
    """``x`` minus each CpG's directed trend: a cubic spline in log(age + maturity / 20) with ``knots`` interior
    knots evenly spaced over log age, flexible enough to absorb development's steep early change."""
    t = np.log(age + maturity / 20.0)
    inner = np.linspace(t.min(), t.max(), knots + 2)[1:-1]
    design = np.column_stack([np.ones(t.size), t, t**2, t**3, *(np.maximum(0.0, t - k) ** 3 for k in inner)])
    coef, *_ = np.linalg.lstsq(design, x, rcond=None)
    # a residual is smaller where the spline bends to few animals: rescale by its leverage
    q, _ = np.linalg.qr(design)
    leverage = np.clip((q**2).sum(axis=1), 0.0, 0.95)
    return (x - design @ coef) / np.sqrt(1.0 - leverage)[:, None]


def noise_per_animal(x: np.ndarray, residual: np.ndarray, columns: np.ndarray) -> np.ndarray:
    """Mean over ``columns`` of |residual| scaled by the CpG's binomial spread at its mean level."""
    m = np.clip(x[:, columns].mean(axis=0), 0.01, 0.99)
    return (np.abs(residual[:, columns]) / np.sqrt(m * (1 - m))).mean(axis=1)


def hinge_fit(age: np.ndarray, d: np.ndarray, grid=GRID) -> dict:
    """The best continuous two-slope fit of ``d`` on ``age`` over change points at ``grid`` quantiles of age."""
    n = age.size
    line = np.column_stack([np.ones(n), age])
    sse_line = float(((d - line @ np.linalg.lstsq(line, d, rcond=None)[0]) ** 2).sum())
    best = None
    for tau in np.unique(np.quantile(age, grid)):
        design = np.column_stack([np.ones(n), age, np.maximum(0.0, age - tau)])
        coef, *_ = np.linalg.lstsq(design, d, rcond=None)
        sse = float(((d - design @ coef) ** 2).sum())
        if best is None or sse < best["sse"]:
            best = {
                "tau": float(tau),
                "slope_before": float(coef[1]),
                "slope_after": float(coef[1] + coef[2]),
                "sse": sse,
            }
    # two extra parameters (the bend and its place) against the line
    best["delta_bic"] = float(n * np.log(sse_line / best["sse"]) - 2 * np.log(n))
    # an onset: a supported bend from (nearly) flat to rising, not a fall that levels off
    best["onset"] = bool(
        best["delta_bic"] > 2 and best["slope_after"] > 0 and best["slope_after"] > 2 * abs(best["slope_before"])
    )
    return best


def read_traits(matrix: Path) -> pl.DataFrame:
    """``sample`` and ``maturity`` (age at sexual maturity, years) from the series matrix."""
    chars = _characteristics(matrix)
    col = "ageatsexualmaturityyears"
    if col not in chars.columns:
        raise ValueError(f"{matrix.name}: no {col}")
    return chars.select("sample", pl.col(col).cast(pl.Float64, strict=False).alias("maturity"))


def across_species(onsets: pl.DataFrame) -> list[dict]:
    """log tau on log maturity, and on log lifespan, with tissue dummies; groups with an onset (``hinge_fit``)."""
    sel = onsets.filter(pl.col("onset") & (pl.col("cpgs") == "background"))
    rows = []
    if sel.height < 8:
        return rows
    y = np.log(sel["tau"].to_numpy())
    tissue = dummies(sel["tissue"].to_numpy()) if sel["tissue"].n_unique() > 1 else np.zeros((sel.height, 0))
    for trait in ("maturity", "max_lifespan"):
        design = np.column_stack([np.ones(sel.height), np.log(sel[trait].to_numpy()), tissue])
        coef, se, p = ols(y, design, 1)
        rows.append(
            {
                "predictor": f"log {trait}",
                "slope": coef,
                "se": se,
                "p": p,
                "groups": sel.height,
                "median_tau_over_maturity": float((sel["tau"] / sel["maturity"]).median()),
            }
        )
    return rows


def run(ctx: Context) -> list[Path]:
    matrix, betas_path = layout(ctx.paths)
    if not betas_path.exists():
        print("drift_onset: run `entropy-age fetch mammals` first. Nothing run.", file=sys.stderr)
        return []
    plots.apply_style()
    out = ctx.out_dir("drift_onset")
    samples = (
        read_samples(matrix)
        .join(read_traits(matrix), on="sample")
        .filter(
            pl.col("age").is_not_null()
            & pl.col("maturity").is_not_null()
            & pl.col("barcode").is_not_null()
            & (pl.col("intervention").is_null() | pl.col("intervention").is_in(["WT", "Control", ""]))
        )
    )
    header = pl.read_csv(betas_path, n_rows=0).columns
    samples = samples.filter(pl.col("barcode").is_in(header))
    groups = (
        samples.group_by("species", "tissue")
        .agg(
            pl.len().alias("n"),
            (pl.col("age") < pl.col("maturity")).sum().alias("before"),
            (pl.col("age") > 2 * pl.col("maturity")).sum().alias("after"),
        )
        .filter((pl.col("n") >= MIN_ANIMALS) & (pl.col("before") >= MIN_SIDE) & (pl.col("after") >= MIN_SIDE))
    )
    keep = samples.join(groups.select("species", "tissue"), on=["species", "tissue"])
    print(f"drift_onset: {groups.height} species-tissue groups, {keep.height} animals", file=sys.stderr, flush=True)
    lay = MouseAnnotationLayout(ctx.paths)
    probes = probes_without_suffix(
        pl.read_parquet(lay.probes).filter(pl.col("promoter") & ~pl.col("chr").is_in(["X", "Y", "MT"]))
    )
    id_col = header[0]
    frame = (
        pl.scan_csv(betas_path, null_values=["NA", ""])
        .select([id_col, *keep["barcode"].unique().to_list()])
        .filter(pl.col(id_col).is_in(probes["cpg"].implode()))
        .collect()
        .rename({id_col: "cpg"})
    )
    features = frame["cpg"].to_list()
    members = mouse_membership(probes.filter(pl.col("cpg").is_in(features)), pl.read_parquet(lay.gene_sets))
    in_any = np.isin(features, members.filter(pl.col("group").is_in(SETS))["feature"].unique().to_list())
    in_set = np.isin(features, members.filter(pl.col("group") == SET)["feature"].unique().to_list())
    rows, curves = [], []
    for (species, tissue), sel in keep.group_by("species", "tissue", maintain_order=True):
        x = frame.select(sel["barcode"].to_list()).to_numpy().T.astype(float)
        ok = ~np.isnan(x).any(axis=0)
        age, maturity = sel["age"].to_numpy(), float(sel["maturity"][0])
        residual = trend_residuals(x[:, ok], age, maturity)
        for label, cols in (("background", ~in_any[ok]), (SET, in_set[ok])):
            if cols.sum() < 50:
                continue
            d = noise_per_animal(x[:, ok], residual, np.flatnonzero(cols))
            fit = hinge_fit(age, d)
            rows.append(
                {
                    "species": species,
                    "tissue": tissue,
                    "cpgs": label,
                    "animals": sel.height,
                    "maturity": maturity,
                    "max_lifespan": float(sel["max_lifespan"][0]),
                    **fit,
                }
            )
            curves.append(
                pl.DataFrame(
                    {
                        "species": species,
                        "tissue": tissue,
                        "cpgs": label,
                        "rel_age": age / maturity,
                        "noise": (d - d.mean()) / d.std(),
                    }
                )
            )
    onsets = pl.DataFrame(rows)
    write_table(onsets, out / "onsets")
    write_table(pl.concat(curves), out / "curves")
    write_table(pl.DataFrame(across_species(onsets)), out / "across_species")
    return [_figure(ctx, onsets, pl.concat(curves))]


def _figure(ctx: Context, onsets: pl.DataFrame, curves: pl.DataFrame) -> Path:
    fig, axes = plots.panels(3, ncols=3, size=4.0)
    ax = axes[0]
    bg = curves.filter(pl.col("cpgs") == "background")
    edges = np.geomspace(0.05, 20, 25)
    mids = np.sqrt(edges[1:] * edges[:-1])
    binned = [
        bg.filter(pl.col("rel_age").is_between(lo, hi))["noise"] for lo, hi in zip(edges[:-1], edges[1:], strict=True)
    ]
    ax.scatter(bg["rel_age"], bg["noise"], s=2, alpha=0.15, color=plots.MUTED)
    ax.plot(mids, [b.median() if b.len() >= 10 else np.nan for b in binned], color=plots.PALETTE[0], lw=2)
    ax.axvline(1.0, color=plots.PALETTE[1], lw=1, ls="--")
    ax.set(
        xscale="log",
        xlabel="age / age at sexual maturity",
        ylabel="noise per animal (z within group)",
        title="Noise against relative age, all groups",
    )
    plots.label_panel(ax, "a")
    sup = onsets.filter((pl.col("cpgs") == "background") & pl.col("onset"))
    ax = axes[1]
    ax.scatter(sup["maturity"], sup["tau"], s=12, color=plots.PALETTE[0])
    lim = [
        min(sup["maturity"].min() or 0.1, sup["tau"].min() or 0.1),
        max(sup["maturity"].max() or 1, sup["tau"].max() or 1),
    ]
    ax.plot(lim, lim, color=plots.MUTED, lw=0.8, ls="--")
    ax.set(
        xscale="log",
        yscale="log",
        xlabel="age at sexual maturity (years)",
        ylabel="change point tau (years)",
        title=f"Onset against maturity ({sup.height} groups with an onset)",
    )
    plots.label_panel(ax, "b")
    ax = axes[2]
    for i, label in enumerate(("background", SET)):
        o = onsets.filter(pl.col("cpgs") == label)
        ax.scatter(
            o["slope_before"] * o["maturity"],
            o["slope_after"] * o["maturity"],
            s=10,
            color=plots.PALETTE[i],
            label=label,
            alpha=0.7,
        )
    ax.axhline(0, color=plots.MUTED, lw=0.8)
    ax.axvline(0, color=plots.MUTED, lw=0.8)
    ax.legend(fontsize=7)
    ax.set(xlabel="slope before tau (per maturity span)", ylabel="slope after tau", title="Before and after the bend")
    plots.label_panel(ax, "c")
    fig.suptitle(
        "When does methylation noise begin? Mammalian consortium, development removed",
        x=0.02,
        y=1.04,
        ha="left",
        fontsize=11,
        fontweight="bold",
    )
    return plots.save(fig, ctx.figure_path("drift_onset"))
