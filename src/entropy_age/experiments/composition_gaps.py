"""Blood's shield with reference-based cell fractions held (roadmap RM34 item 2).

The omics re-analysis practice asks for cell composition first: blood methylation shifts with the cell mix, and a
composition change with age reads as an epigenetic one. Hannum's headline gaps (``figures tissues``) hold plate,
source, sex and ethnicity; the composition checks so far used principal components (``tissues_cells``), sorted cells
(R44) and lineages (R49). Here 12 blood-cell fractions are estimated per donor by constrained projection (Houseman)
on EpiDISH's cent12CT450k reference, as in R38, and 11 of them (one left out, the fractions sum to one) are added to
the covariates of every age slope. The gaps on all four drift metrics are reported with and without them.

``entropy-age figures composition_gaps`` as a Slurm job (after ``fetch epidish`` and the Hannum cohort).
"""

from __future__ import annotations

import dataclasses
import sys
from pathlib import Path

import numpy as np
import polars as pl

from entropy_age import plots
from entropy_age.experiments.common import Context, write_table
from entropy_age.experiments.stochastic_component import EPIDISH, deconvolve
from entropy_age.experiments.tissue_contrast import (
    LIFESPAN_METRICS,
    _features,
    _hannum,
    _lifespan_metrics,
    _matched_sd,
    _matrix,
    _with_strata,
)

LABELS = {
    "entropy_slope": "entropy",
    "kl_slope": "relative entropy",
    "abs_slope": "|mean change|",
    "dispersion_slope": "dispersion",
}


def fractions(ctx: Context, samples: list[str]) -> tuple[np.ndarray, list[str]]:
    """Donors x 12 cell fractions (Houseman projection on EpiDISH cent12CT450k) and the cell-type names."""
    import rdata

    ref = next(iter(rdata.read_rda(str(ctx.paths.input / EPIDISH)).values()))
    ref_cpgs, cells = [str(c) for c in ref.coords["dim_0"].values], [str(c) for c in ref.coords["dim_1"].values]
    m = _matrix("hannum", ctx.paths, samples, ref_cpgs)
    keep = [i for i, c in enumerate(ref_cpgs) if c in set(m.features)]
    x = np.where(np.isnan(m.values), np.nanmean(m.values, axis=0), m.values)
    order = [m.features.index(ref_cpgs[i]) for i in keep]
    return deconvolve(x[:, order], ref.values[keep]), cells


def run(ctx: Context) -> list[Path]:
    if not (ctx.paths.input / EPIDISH).exists():
        print("composition_gaps: run `entropy-age fetch epidish` first; nothing run", file=sys.stderr)
        return []
    plots.apply_style()
    out = ctx.out_dir("composition_gaps")
    features, probes, membership = _features(ctx.paths, ctx.rng(120_000))
    p = _hannum(ctx.paths, features)
    from entropy_age.data.cohorts import layout as cohort_layout

    samples = pl.read_parquet(cohort_layout("hannum", ctx.paths).samples).drop_nulls("age").sort("age", "sample")
    frac, cells = fractions(ctx, samples["sample"].to_list())
    write_table(pl.DataFrame(frac, schema=cells).with_columns(sample=samples["sample"]), out / "fractions")
    with_cells = dataclasses.replace(p, covariates=np.column_stack([p.covariates, frac[:, :-1]]))
    rows = []
    for label, q in (("plate, source, sex, ethnicity", p), ("+ 12 cell fractions", with_cells)):
        t = _matched_sd(_with_strata(_lifespan_metrics(q), probes), membership, LIFESPAN_METRICS, 1000, ctx.seed)
        rows.append(t.with_columns(covariates=pl.lit(label)))
    table = pl.concat(rows)
    write_table(table, out / "gaps")
    age_r = {c: float(np.corrcoef(frac[:, i], p.ages)[0, 1]) for i, c in enumerate(cells)}
    write_table(pl.DataFrame({"cell": list(age_r), "r_with_age": list(age_r.values())}), out / "fraction_age")
    print(f"composition_gaps: {p.x.shape[0]} donors, {len(p.features)} CpGs", file=sys.stderr)
    return [_figure(ctx, table)]


def _figure(ctx: Context, table: pl.DataFrame) -> Path:
    fig, axes = plots.panels(1, ncols=1, size=5.0)
    ax = axes[0]
    keys = table.filter(pl.col("covariates") == table["covariates"][0]).select("metric", "group")
    y = np.arange(keys.height)
    for i, cov in enumerate(table["covariates"].unique(maintain_order=True)):
        t = keys.join(table.filter(pl.col("covariates") == cov), on=["metric", "group"], how="left")
        ax.barh(y + 0.38 * i, t["matched_diff_sd"], height=0.36, color=plots.PALETTE[i], label=cov)
    ax.set_yticks(y + 0.19, [f"{LABELS.get(m, m)} · {g.replace('_', ' ')}" for m, g in keys.iter_rows()], fontsize=7)
    ax.axvline(0, color=plots.MUTED, lw=0.8)
    ax.invert_yaxis()
    ax.set(xlabel="age gap, set minus look-alikes (SD)", title="Hannum blood: cell fractions held")
    ax.legend(fontsize=7, loc="upper center", bbox_to_anchor=(0.5, -0.1), ncol=2, frameon=False)
    return plots.save(fig, ctx.figure_path("composition_gaps"))
