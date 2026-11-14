"""Placebo gene sets: does the matching itself make a shield? (roadmap RM34 item 1)

Every shield gap compares a gene set's promoter CpGs with look-alikes of the same island status and young-level
decile. If promoter CpGs of *any* set of genes came out below their look-alikes, the shield would be a property of
the matching (genes' CpGs are clustered and correlated; promoters differ from random CpGs in ways the strata miss),
not of essential genes. This draws placebo sets: random genes from those linked to a promoter CpG but in none of the
three sets, as many genes as the real set, all their promoter CpGs, and runs them through the same contrast on the
same drift metrics (binary entropy, relative entropy, ``|mean change|``, dispersion; ``figures tissues``' samples and
covariates, on every autosomal promoter CpG). The real gap is then read against the placebo distribution: its
rank, and the placebo mean, which should be 0.

``entropy-age figures placebo_sets`` as a Slurm job.
"""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import polars as pl

from entropy_age import plots
from entropy_age.data.annotation import AnnotationLayout, link_membership
from entropy_age.experiments.common import Context, write_table
from entropy_age.experiments.tissue_contrast import (
    LIFESPAN_METRICS,
    SETS,
    _hannum,
    _jaffe,
    _lifespan_metrics,
    _matched_sd,
    _with_strata,
)

N_PLACEBO = {"paper": 200, "quick": 20}
LABELS = {
    "entropy_slope": "entropy",
    "kl_slope": "relative entropy",
    "abs_slope": "|mean change|",
    "dispersion_slope": "dispersion",
}


def placebo_membership(
    links: pl.DataFrame, real: pl.DataFrame, features: set[str], n_genes: int, name: str, rng
) -> pl.DataFrame:
    """``feature``/``group`` rows for ``n_genes`` random genes outside every real set, all their promoter CpGs."""
    taken = set(real["gene"].to_list())
    pool = links.filter(pl.col("promoter") & pl.col("cpg").is_in(list(features)) & ~pl.col("symbol").is_in(list(taken)))
    genes = np.array(sorted(set(pool["symbol"].drop_nulls().to_list())))
    pick = set(rng.choice(genes, size=min(n_genes, genes.size), replace=False).tolist())
    rows = pool.filter(pl.col("symbol").is_in(list(pick))).select(pl.col("cpg").alias("feature")).unique()
    return rows.with_columns(group=pl.lit(name))


def summarise(real: pl.DataFrame, placebo: pl.DataFrame) -> pl.DataFrame:
    """Per tissue, set and metric: the real gap, the placebo mean and 95 % range, and the real gap's rank."""
    rows = []
    for (tissue, group, metric), g in real.group_by("tissue", "group", "metric", maintain_order=True):
        p = placebo.filter((pl.col("tissue") == tissue) & (pl.col("set") == group) & (pl.col("metric") == metric))
        v = p["matched_diff_sd"].to_numpy()
        r = float(g["matched_diff_sd"][0])
        rows.append(
            {
                "tissue": tissue,
                "group": group,
                "metric": metric,
                "real": r,
                "placebo_mean": float(v.mean()),
                "placebo_low": float(np.quantile(v, 0.025)),
                "placebo_high": float(np.quantile(v, 0.975)),
                "share_placebo_below_real": float((v <= r).mean()),
                "placebos": int(v.size),
            }
        )
    return pl.DataFrame(rows)


def run(ctx: Context) -> list[Path]:
    plots.apply_style()
    out = ctx.out_dir("placebo_sets")
    ann = AnnotationLayout(ctx.paths)
    probes = pl.read_parquet(ann.probes_450k).filter(pl.col("promoter") & ~pl.col("chr").is_in(["X", "Y"]))
    links = pl.read_parquet(ann.links_450k)
    membership = link_membership(links, pl.read_parquet(ann.gene_sets))
    membership = membership.filter(pl.col("group").is_in(SETS) & pl.col("feature").is_in(probes["cpg"].implode()))
    features = sorted(probes["cpg"].to_list())
    n_genes = {g: membership.filter(pl.col("group") == g)["gene"].n_unique() for g in SETS}
    rng = ctx.rng(680_000)
    reals, placebos = [], []
    for p in (_hannum(ctx.paths, features), _jaffe(ctx.paths, features)):
        metrics = _with_strata(_lifespan_metrics(p), probes)
        have = set(p.features)
        real = _matched_sd(metrics, membership.select("feature", "group"), LIFESPAN_METRICS, 1000, ctx.seed)
        reals.append(real.with_columns(tissue=pl.lit(p.label)))
        for group in SETS:
            for b in range(N_PLACEBO[str(ctx.scale.name)]):
                mem = placebo_membership(links, membership, have, n_genes[group], f"placebo_{group}_{b}", rng)
                t = _matched_sd(metrics, mem, LIFESPAN_METRICS, 0, ctx.seed)
                placebos.append(t.with_columns(tissue=pl.lit(p.label), set=pl.lit(group), draw=pl.lit(b)))
        print(f"placebo_sets: {p.label}: {len(p.features)} CpGs, {p.x.shape[0]} donors", file=sys.stderr, flush=True)
    real = pl.concat(reals)
    placebo = pl.concat(placebos, how="diagonal")
    summary = summarise(real, placebo)
    write_table(summary, out / "summary")
    placebo.write_parquet(out / "placebo_draws.parquet")
    write_table(pl.DataFrame([{"set": g, "genes": n} for g, n in n_genes.items()]), out / "set_genes")
    return [_figure(ctx, summary)]


def _figure(ctx: Context, summary: pl.DataFrame) -> Path:
    tissues = summary["tissue"].unique(maintain_order=True).to_list()
    fig, axes = plots.panels(len(tissues), ncols=len(tissues), size=4.4)
    for ax, tissue, letter in zip(axes, tissues, "ab", strict=False):
        t = summary.filter(pl.col("tissue") == tissue)
        y = np.arange(t.height)
        ax.hlines(y, t["placebo_low"], t["placebo_high"], color=plots.MUTED, lw=3, label="placebo sets, 95 %")
        ax.scatter(t["placebo_mean"], y, color=plots.MUTED, s=12, zorder=3)
        ax.scatter(t["real"], y, color=plots.PALETTE[0], s=22, zorder=4, label="real set")
        labels = [f"{g.replace('_', ' ')} · {LABELS[m]}" for g, m in zip(t["group"], t["metric"], strict=True)]
        ax.set_yticks(y, labels if letter == "a" else [""] * len(labels), fontsize=6)
        ax.axvline(0, color=plots.MUTED, lw=0.8)
        ax.invert_yaxis()
        ax.set(xlabel="set minus look-alikes (SD)", title=tissue)
        if letter == "a":
            ax.legend(fontsize=7, loc="lower left")
        plots.label_panel(ax, letter)
    fig.suptitle("Placebo gene sets through the same matching", x=0.02, y=1.03, ha="left", fontsize=11)
    return plots.save(fig, ctx.figure_path("placebo_sets"))
