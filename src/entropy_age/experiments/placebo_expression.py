"""Placebo gene sets matched on expression (roadmap RM38 item 3; a simulated referee's concern on round v0.1).

R69 (``placebo_sets``) drew placebo genes at random from all genes outside the three sets, and they gave a gap of
0. Essential and housekeeping genes are highly expressed, and a promoter in use may drift less for that reason
alone, so a random draw does not answer whether expression, not essentiality, makes the shield. Here each placebo
set is drawn decile by decile of the tissue's GTEx v8 median expression (log2 TPM + 1; Whole Blood for Hannum,
Brain - Frontal Cortex BA9 for Jaffe): as many genes outside the three sets in each decile as the real set has
there. Everything else is R69's: all promoter CpGs of the drawn genes, the shield's matching (island x young-level
decile) and metrics, the real gap read against 200 placebo draws.

``entropy-age figures placebo_expression`` as a Slurm job.
"""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import polars as pl

from entropy_age import plots
from entropy_age.data.annotation import AnnotationLayout, link_membership
from entropy_age.experiments.common import Context, write_table
from entropy_age.experiments.placebo_sets import LABELS, N_PLACEBO, summarise
from entropy_age.experiments.tissue_contrast import (
    LIFESPAN_METRICS,
    SETS,
    _hannum,
    _jaffe,
    _lifespan_metrics,
    _matched_sd,
    _with_strata,
)

GTEX = "gtex/GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_median_tpm.gct.gz"
TISSUE_TPM = {"blood (Hannum)": "Whole Blood", "cortex (Jaffe)": "Brain - Frontal Cortex (BA9)"}
DECILES = 10


def expression(path: Path, column: str) -> pl.DataFrame:
    """``symbol`` and ``log_tpm`` (log2 median TPM + 1) of one GTEx tissue; a symbol seen twice keeps its maximum."""
    frame = pl.read_csv(path, separator="\t", skip_rows=2, columns=["Description", column])
    return (
        frame.select(
            pl.col("Description").alias("symbol"), (pl.col(column).cast(pl.Float64) + 1).log(2).alias("log_tpm")
        )
        .group_by("symbol")
        .agg(pl.col("log_tpm").max())
    )


def with_deciles(genes: pl.DataFrame, n: int = DECILES) -> pl.DataFrame:
    """Adds ``decile`` (0..n-1) of ``log_tpm`` over the genes given; ties share a decile by rank."""
    rank = genes["log_tpm"].rank("ordinal").to_numpy() - 1
    return genes.with_columns(decile=pl.Series((rank * n // max(genes.height, 1)).astype(np.int64)))


def matched_draw(real: pl.DataFrame, pool: pl.DataFrame, rng) -> tuple[list[str], int]:
    """Genes from ``pool`` with the decile counts of ``real``; returns them and the shortfall where a decile ran out."""
    need = dict(real.group_by("decile").len().iter_rows())
    picked, short = [], 0
    for d, k in sorted(need.items()):
        cand = np.array(sorted(pool.filter(pl.col("decile") == d)["symbol"].to_list()))
        take = min(k, cand.size)
        short += k - take
        if take:
            picked += rng.choice(cand, size=take, replace=False).tolist()
    return picked, short


def run(ctx: Context) -> list[Path]:
    path = ctx.paths.input / GTEX
    if not path.exists():
        print(f"placebo_expression: {path} missing", file=sys.stderr)
        return []
    plots.apply_style()
    out = ctx.out_dir("placebo_expression")
    ann = AnnotationLayout(ctx.paths)
    probes = pl.read_parquet(ann.probes_450k).filter(pl.col("promoter") & ~pl.col("chr").is_in(["X", "Y"]))
    links = pl.read_parquet(ann.links_450k)
    membership = link_membership(links, pl.read_parquet(ann.gene_sets))
    membership = membership.filter(pl.col("group").is_in(SETS) & pl.col("feature").is_in(probes["cpg"].implode()))
    features = sorted(probes["cpg"].to_list())
    rng = ctx.rng(380_003)
    reals, placebos, balance = [], [], []
    for p in (_hannum(ctx.paths, features), _jaffe(ctx.paths, features)):
        have = set(p.features)
        prom = links.filter(pl.col("promoter") & pl.col("cpg").is_in(list(have))).select("cpg", "symbol").drop_nulls()
        expr = with_deciles(expression(path, TISSUE_TPM[p.label]).join(prom.select("symbol").unique(), on="symbol"))
        taken = set(membership["gene"].to_list())
        pool = expr.filter(~pl.col("symbol").is_in(list(taken)))
        metrics = _with_strata(_lifespan_metrics(p), probes)
        real = _matched_sd(metrics, membership.select("feature", "group"), LIFESPAN_METRICS, 1000, ctx.seed)
        reals.append(real.with_columns(tissue=pl.lit(p.label)))
        for group in SETS:
            real_genes = expr.filter(
                pl.col("symbol").is_in(membership.filter(pl.col("group") == group)["gene"].implode())
            )
            for b in range(N_PLACEBO[str(ctx.scale.name)]):
                genes, short = matched_draw(real_genes, pool, rng)
                mem = (
                    prom.filter(pl.col("symbol").is_in(genes))
                    .select(pl.col("cpg").alias("feature"))
                    .unique()
                    .with_columns(group=pl.lit(f"placebo_{group}_{b}"))
                )
                t = _matched_sd(metrics, mem, LIFESPAN_METRICS, 0, ctx.seed)
                placebos.append(t.with_columns(tissue=pl.lit(p.label), set=pl.lit(group), draw=pl.lit(b)))
                if b == 0:
                    drawn = pool.filter(pl.col("symbol").is_in(genes))
                    balance.append(
                        {
                            "tissue": p.label,
                            "set": group,
                            "real_genes": real_genes.height,
                            "placebo_genes": drawn.height,
                            "shortfall": short,
                            "real_log_tpm": float(real_genes["log_tpm"].mean()),
                            "placebo_log_tpm": float(drawn["log_tpm"].mean()),
                        }
                    )
        print(f"placebo_expression: {p.label}: {len(p.features)} CpGs", file=sys.stderr, flush=True)
    summary = summarise(pl.concat(reals), pl.concat(placebos, how="diagonal"))
    write_table(summary, out / "summary")
    write_table(pl.DataFrame(balance), out / "balance")
    pl.concat(placebos, how="diagonal").write_parquet(out / "placebo_draws.parquet")
    return [_figure(ctx, summary)]


def _figure(ctx: Context, summary: pl.DataFrame) -> Path:
    tissues = summary["tissue"].unique(maintain_order=True).to_list()
    fig, axes = plots.panels(len(tissues), ncols=len(tissues), size=4.4)
    for ax, tissue, letter in zip(axes, tissues, "ab", strict=False):
        t = summary.filter(pl.col("tissue") == tissue)
        y = np.arange(t.height)
        ax.hlines(
            y, t["placebo_low"], t["placebo_high"], color=plots.MUTED, lw=3, label="expression-matched placebo, 95 %"
        )
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
    fig.suptitle("Placebo gene sets matched on expression", x=0.02, y=1.03, ha="left", fontsize=11)
    return plots.save(fig, ctx.figure_path("placebo_expression"))
