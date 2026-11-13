"""Is the species shield a floor? Essential promoters near 0 % cannot drift down (roadmap RM13 (4))

On the conserved mammalian array most essential-promoter CpGs sit near 0 % methylation, where drift is floored
(``figures species``, R12). Matching on the young level decile puts look-alikes in the same decile, but inside the
lowest decile set CpGs may still sit lower. This refits every species-tissue group's gap with only CpGs off the
floor and the ceiling (young-quarter mean between ``LOW`` and ``HIGH``), set and look-alikes alike, and refits the
lifespan regression (``gap ~ log lifespan + log weight + tissue``) on those gaps. A shield that is only a floor
vanishes off it; a lifespan result that rests on the floor changes.

``entropy-age figures species_floor`` as a Slurm job (the same input as ``figures species``).
"""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import polars as pl

from entropy_age import plots
from entropy_age.data.mouse import MouseAnnotationLayout, mouse_membership
from entropy_age.experiments.bxd_liver import probes_without_suffix
from entropy_age.experiments.common import Context, write_table
from entropy_age.experiments.interventions import _sets
from entropy_age.experiments.species_shield import (
    MIN_AGE_SPAN,
    MIN_ANIMALS,
    N_RESAMPLES,
    _lifespan_fit,
    group_gap,
    layout,
    read_samples,
)
from entropy_age.experiments.tissue_contrast import SETS

LOW, HIGH = 0.10, 0.90
VARIANTS = ("all CpGs", "off floor and ceiling")


def off_floor(x: np.ndarray, ages: np.ndarray, low: float = LOW, high: float = HIGH) -> np.ndarray:
    """Columns whose mean over the youngest quarter lies strictly between ``low`` and ``high``."""
    young = ages <= np.quantile(ages, 0.25)
    level = x[young].mean(axis=0)
    return (level > low) & (level < high)


def run(ctx: Context) -> list[Path]:
    matrix, betas_path = layout(ctx.paths)
    if not betas_path.exists():
        print("species_floor: run `entropy-age fetch mammals` first. Nothing run.", file=sys.stderr)
        return []
    plots.apply_style()
    out = ctx.out_dir("species_floor")
    samples = read_samples(matrix).filter(
        pl.col("age").is_not_null()
        & pl.col("max_lifespan").is_not_null()
        & pl.col("barcode").is_not_null()
        & (pl.col("intervention").is_null() | pl.col("intervention").is_in(["WT", "Control", ""]))
    )
    lay = MouseAnnotationLayout(ctx.paths)
    probes = probes_without_suffix(
        pl.read_parquet(lay.probes).filter(pl.col("promoter") & ~pl.col("chr").is_in(["X", "Y", "MT"]))
    )
    header = pl.read_csv(betas_path, n_rows=0).columns
    groups = (
        samples.filter(pl.col("barcode").is_in(header))
        .group_by("species", "tissue")
        .agg(
            pl.len().alias("n"),
            pl.col("age").min().alias("lo"),
            pl.col("age").max().alias("hi"),
            pl.col("max_lifespan").first(),
        )
        .filter((pl.col("n") >= MIN_ANIMALS) & ((pl.col("hi") - pl.col("lo")) >= MIN_AGE_SPAN * pl.col("max_lifespan")))
    )
    keep = samples.join(groups.select("species", "tissue"), on=["species", "tissue"]).filter(
        pl.col("barcode").is_in(header)
    )
    frame = (
        pl.scan_csv(betas_path, null_values=["NA", ""])
        .select([header[0], *keep["barcode"].unique().to_list()])
        .filter(pl.col(header[0]).is_in(probes["cpg"].implode()))
        .collect()
        .rename({header[0]: "cpg"})
    )
    features = frame["cpg"].to_list()
    members = mouse_membership(probes.filter(pl.col("cpg").is_in(features)), pl.read_parquet(lay.gene_sets)).filter(
        pl.col("group").is_in(SETS)
    )
    island = dict(probes.select("cpg", "cpg_island").iter_rows())
    isl = np.array([bool(island.get(f, False)) for f in features])
    sets = _sets(features, members)
    rng = ctx.rng(190_000)
    print(f"species_floor: {groups.height} groups, {len(features)} CpGs", file=sys.stderr, flush=True)
    rows, shares = [], []
    for (species, tissue), sel in keep.group_by("species", "tissue", maintain_order=True):
        x = frame.select(sel["barcode"].to_list()).to_numpy().T.astype(float)
        ok = ~np.isnan(x).any(axis=0)
        if ok.sum() < 0.9 * x.shape[1]:
            continue
        ages = sel["age"].to_numpy()
        mid = off_floor(x[:, ok], ages)
        first = sel.row(0, named=True)
        for g, (s, _) in sets.items():
            shares.append(
                {"species": species, "tissue": tissue, "group": g, "share_off_floor": float(mid[s[ok]].mean())}
            )
        for variant, cols in ((VARIANTS[0], np.flatnonzero(ok)), (VARIANTS[1], np.flatnonzero(ok)[mid])):
            sub_sets = {
                g: (s[cols], [np.searchsorted(cols, gi[np.isin(gi, cols)]) for gi in genes if np.isin(gi, cols).any()])
                for g, (s, genes) in sets.items()
            }
            sub_sets = {g: v for g, v in sub_sets.items() if v[0].sum() >= 10}
            gaps = group_gap(
                x[:, cols],
                ages,
                [features[c] for c in cols],
                isl[cols],
                sub_sets,
                N_RESAMPLES[str(ctx.scale.name)],
                rng,
            )
            rows += [
                {
                    "variant": variant,
                    "species": species,
                    "tissue": tissue,
                    "animals": sel.height,
                    "max_lifespan": first["max_lifespan"],
                    "adult_weight": first["adult_weight"],
                    **r,
                }
                for r in gaps
            ]
    gaps = pl.DataFrame(rows)
    write_table(gaps, out / "species_gaps")
    write_table(pl.DataFrame(shares), out / "share_off_floor")
    fits = []
    for variant in VARIANTS:
        for metric in ("entropy_slope", "dispersion_slope", "kl_slope"):
            for group in SETS:
                sel = gaps.filter(
                    (pl.col("variant") == variant)
                    & (pl.col("metric") == metric)
                    & (pl.col("group") == group)
                    & (pl.col("adult_weight") > 0)
                ).drop_nulls(["adult_weight", "gap_sd"])
                if sel.height < 10:
                    continue
                fits += [
                    {"variant": variant, **r}
                    for r in _lifespan_fit(sel, sel["gap_sd"].to_numpy(), metric, group, "gap_sd")
                ]
    fits = pl.DataFrame(fits)
    write_table(fits, out / "lifespan_fits")
    return [_figure(ctx, gaps)]


def _figure(ctx: Context, gaps: pl.DataFrame) -> Path:
    fig, axes = plots.panels(2, ncols=2, size=4.4)
    for ax, group, letter in zip(axes, ("housekeeping", "core_essential"), "ab", strict=True):
        wide = (
            gaps.filter((pl.col("metric") == "entropy_slope") & (pl.col("group") == group))
            .select("variant", "species", "tissue", "gap_sd")
            .pivot(on="variant", index=["species", "tissue"], values="gap_sd")
            .drop_nulls()
        )
        ax.scatter(wide[VARIANTS[0]], wide[VARIANTS[1]], s=10, alpha=0.6, color=plots.PALETTE[0])
        lim = [
            min(wide[VARIANTS[0]].min(), wide[VARIANTS[1]].min()),
            max(wide[VARIANTS[0]].max(), wide[VARIANTS[1]].max()),
        ]
        ax.plot(lim, lim, color=plots.MUTED, lw=0.8)
        ax.axhline(0, color=plots.MUTED, lw=0.6, ls=":")
        ax.axvline(0, color=plots.MUTED, lw=0.6, ls=":")
        ax.set(
            xlabel="entropy gap, all CpGs (SD)",
            ylabel=f"entropy gap, CpGs at {LOW:.0%}-{HIGH:.0%} (SD)",
            title=f"{group.replace('_', ' ')} ({wide.height} groups)",
        )
        plots.label_panel(ax, letter)
    fig.subplots_adjust(wspace=0.4)
    return plots.save(fig, ctx.figure_path("species_floor"))
