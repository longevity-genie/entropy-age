"""Does methylation noise reach transcription? The same GTEx donors, both layers (roadmap RM30)

Idea 7 of docs/BEYOND_SHIELD.md; open in docs/LITERATURE_COVERAGE.md for spread (mean-level coupling is published:
Oliva 2023). GTEx measured promoter methylation (450k-type EPIC subset, 9 tissues) and RNA in the same donors.
Per tissue, over the donors with both:

1. **Signed coupling (the positive control).** Per gene, the correlation across donors of the promoter methylation
   residual and the expression residual (both on age and technical covariates). Promoter methylation silences, so
   the mean over genes should be negative; if it is not, donors are mismatched or the layers do not meet.
2. **Spread coupling, per donor.** Per gene, the correlation across donors of the *absolute* residuals: is a donor
   who is unusual in methylation at a gene's promoter also unusual in its expression? Mean over genes, against a
   null from permuting donors in one layer (the same permutation for all genes keeps gene-gene structure).
3. **Spread coupling, with age.** Across genes, Spearman correlation of the age slope of methylation spread with
   that of expression spread, each first residualised on the gene's mean level (spread depends on level).

Methylation per gene: the mean over its promoter CpGs (TSS200, TSS1500) of the residual divided by
``sqrt(m (1 - m))``. Expression: expressed protein-coding genes as log2(CPM + 1), as in ``rna_gtex``.

``entropy-age figures meth_expr_spread`` (as a Slurm job; needs the GTEx RNA settings).
"""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import polars as pl
from scipy import stats

from entropy_age import plots
from entropy_age.data.annotation import AnnotationLayout
from entropy_age.data.cohorts import layout as cohort_layout
from entropy_age.data.gtex_rna import GtexRnaMissing, GtexRnaSource, read_counts, read_genes, read_samples
from entropy_age.experiments.common import Context, dummies, write_table
from entropy_age.experiments.rna_gtex import expressed, log_cpm
from entropy_age.experiments.tissue_pace import GTEX_RNA

PERMUTATIONS = 200
MIN_DONORS = 40
PROMOTER = ("TSS200", "TSS1500")


def residuals(y: np.ndarray, design: np.ndarray) -> np.ndarray:
    coef, *_ = np.linalg.lstsq(design, y, rcond=None)
    return y - design @ coef


def columnwise_r(a: np.ndarray, b: np.ndarray) -> np.ndarray:
    """Pearson r of each column of ``a`` with the same column of ``b`` (rows: donors)."""
    a = a - a.mean(axis=0)
    b = b - b.mean(axis=0)
    den = np.sqrt((a**2).sum(axis=0) * (b**2).sum(axis=0))
    return np.where(den > 0, (a * b).sum(axis=0) / np.where(den > 0, den, 1.0), np.nan)


def spread_coupling(
    meth_res: np.ndarray, expr_res: np.ndarray, rng: np.random.Generator, permutations: int = PERMUTATIONS
) -> dict:
    """Signed and absolute per-gene coupling across donors, with a donor-permutation null for the absolute one."""
    signed = columnwise_r(meth_res, expr_res)
    absolute = columnwise_r(np.abs(meth_res), np.abs(expr_res))
    observed = float(np.nanmean(absolute))
    null = np.array(
        [
            np.nanmean(columnwise_r(np.abs(meth_res[rng.permutation(meth_res.shape[0])]), np.abs(expr_res)))
            for _ in range(permutations)
        ]
    )
    return {
        "signed_mean_r": float(np.nanmean(signed)),
        "signed_share_negative": float(np.nanmean(signed < 0)),
        "abs_mean_r": observed,
        "abs_null_mean": float(null.mean()),
        "abs_null_sd": float(null.std(ddof=1)),
        "abs_p": float((1 + (null >= observed).sum()) / (1 + permutations)),
    }


def slope_coupling(meth_slope, meth_level, expr_slope, expr_level) -> tuple[float, float]:
    """Spearman r over genes of the two spread slopes, each residualised on its rank of mean level."""

    def resid(v, level):
        rv, rl = stats.rankdata(v), stats.rankdata(level)
        return residuals(rv, np.column_stack([np.ones(rv.size), rl]))

    r, p = stats.spearmanr(resid(meth_slope, meth_level), resid(expr_slope, expr_level))
    return float(r), float(p)


def _age_slope(abs_res: np.ndarray, age: np.ndarray, cov: np.ndarray) -> np.ndarray:
    design = np.column_stack([np.ones(age.size), age, cov])
    coef, *_ = np.linalg.lstsq(design, abs_res, rcond=None)
    return coef[1]


def _promoter_genes(annotation: AnnotationLayout) -> pl.DataFrame:
    """``cpg`` -> ``ensembl`` for promoter CpGs of protein-coding genes."""
    links = pl.read_parquet(annotation.links_450k).filter(pl.col("region").is_in(PROMOTER))
    hgnc = pl.read_parquet(annotation.hgnc).filter(pl.col("protein_coding") & pl.col("ensembl").is_not_null())
    return links.join(hgnc.select("symbol", "ensembl"), on="symbol").select("cpg", "ensembl").unique()


def run(ctx: Context) -> list[Path]:
    try:
        source = GtexRnaSource.from_environment()
    except GtexRnaMissing as error:
        print(f"meth_expr_spread: {error}. Nothing run.", file=sys.stderr)
        return []
    plots.apply_style()
    out = ctx.out_dir("meth_expr_spread")
    rng = ctx.rng(30_000)
    cpg_gene = _promoter_genes(AnnotationLayout(ctx.paths))
    meth_samples = pl.read_parquet(cohort_layout("gtex", ctx.paths).samples)
    rna_samples = read_samples(source)
    genes_all = read_genes(source)
    rows, gene_rows = [], []
    for rna_tissue, tissue in GTEX_RNA.items():
        m = meth_samples.filter(pl.col("tissue") == tissue).unique("donor", keep="first")
        r = rna_samples.filter(pl.col("tissue") == rna_tissue).unique("donor", keep="first")
        both = m.join(r.select("donor", pl.col("sample").alias("rna_sample"), "rin", "ischemic_minutes"), on="donor")
        both = both.drop_nulls(["age", "rin", "ischemic_minutes"]).sort("donor")
        if both.height < MIN_DONORS:
            print(f"meth_expr_spread: {tissue}: {both.height} donors with both; left out", file=sys.stderr)
            continue
        age = both["age"].to_numpy().astype(float)
        sex = dummies(both["sex"].to_numpy())
        # expression
        values, library = log_cpm(read_counts(source, both["rna_sample"].to_list()).astype(float))
        keep = expressed(values)
        genes = [g for g, k in zip(genes_all, keep, strict=True) if k]
        expr = values[:, keep]
        tech = np.column_stack([sex, both["rin"].to_numpy(), both["ischemic_minutes"].to_numpy(), np.log2(library)])
        expr_res = residuals(expr, np.column_stack([np.ones(age.size), age, tech]))
        # methylation, per gene
        links = cpg_gene.filter(pl.col("ensembl").is_in(genes))
        frame = (
            pl.scan_parquet(cohort_layout("gtex", ctx.paths).betas)
            .select("cpg", *both["sample"].to_list())
            .filter(pl.col("cpg").is_in(links["cpg"].unique().implode()))
            .collect()
        )
        cpgs = frame["cpg"].to_list()
        beta = frame.drop("cpg").to_numpy().T.astype(float)
        ok = ~np.isnan(beta).any(axis=0)
        beta, cpgs = beta[:, ok], [c for c, k in zip(cpgs, ok, strict=True) if k]
        level = beta.mean(axis=0)
        scale = np.sqrt(np.clip(level * (1 - level), 1e-4, None))
        cpg_res = residuals(beta, np.column_stack([np.ones(age.size), age, sex])) / scale
        col = {c: i for i, c in enumerate(cpgs)}
        by_gene = links.filter(pl.col("cpg").is_in(cpgs)).group_by("ensembl").agg("cpg")
        gene_idx = {g: i for i, g in enumerate(genes)}
        shared = [g for g in by_gene["ensembl"].to_list() if g in gene_idx]
        cpg_lists = dict(by_gene.iter_rows())
        meth_res = np.column_stack([cpg_res[:, [col[c] for c in cpg_lists[g]]].mean(axis=1) for g in shared])
        meth_level = np.array([level[[col[c] for c in cpg_lists[g]]].mean() for g in shared])
        e_idx = [gene_idx[g] for g in shared]
        expr_g = expr_res[:, e_idx]
        coupling = spread_coupling(meth_res, expr_g, rng)
        ms = _age_slope(np.abs(meth_res), age, sex)
        es = _age_slope(np.abs(expr_g), age, tech)
        slope_r, slope_p = slope_coupling(ms, meth_level, es, expr[:, e_idx].mean(axis=0))
        rows.append(
            {
                "tissue": tissue,
                "donors": both.height,
                "genes": len(shared),
                **coupling,
                "slope_r": slope_r,
                "slope_p": slope_p,
            }
        )
        gene_rows.append(
            pl.DataFrame(
                {
                    "tissue": tissue,
                    "gene": shared,
                    "signed_r": columnwise_r(meth_res, expr_g),
                    "abs_r": columnwise_r(np.abs(meth_res), np.abs(expr_g)),
                    "meth_spread_slope": ms,
                    "expr_spread_slope": es,
                }
            )
        )
        print(f"meth_expr_spread: {tissue}: {both.height} donors, {len(shared)} genes", file=sys.stderr, flush=True)
    table = pl.DataFrame(rows)
    write_table(table, out / "coupling")
    write_table(pl.concat(gene_rows), out / "genes")
    return [_figure(ctx, table)]


def _figure(ctx: Context, table: pl.DataFrame) -> Path:
    fig, axes = plots.panels(3, ncols=3, size=3.8)
    t = table.sort("donors", descending=True)
    y = np.arange(t.height)
    labels = [f"{a} ({n})" for a, n in t.select("tissue", "donors").iter_rows()]
    axes[0].barh(y, t["signed_mean_r"], color=plots.PALETTE[0])
    axes[0].set(title="Signed coupling (control: < 0)", xlabel="mean r over genes")
    axes[1].barh(y, t["abs_mean_r"] - t["abs_null_mean"], xerr=1.96 * t["abs_null_sd"], color=plots.PALETTE[1])
    axes[1].set(title="Spread coupling per donor", xlabel="mean r minus permutation null")
    axes[2].barh(y, t["slope_r"], color=plots.PALETTE[2])
    axes[2].set(title="Spread rise with age, across genes", xlabel="Spearman r (level held)")
    for i, ax in enumerate(axes):
        ax.axvline(0, color=plots.MUTED, lw=0.8)
        ax.set_yticks(y, labels if i == 0 else [""] * len(labels), fontsize=8)
        ax.invert_yaxis()
        plots.label_panel(ax, "abc"[i])
    fig.suptitle(
        "Does methylation noise reach expression? GTEx, same donors (donors in brackets)",
        x=0.02,
        y=1.04,
        ha="left",
        fontsize=11,
        fontweight="bold",
    )
    return plots.save(fig, ctx.figure_path("meth_expr_spread"))
