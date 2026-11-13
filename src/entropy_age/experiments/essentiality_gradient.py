"""Does a gene's essentiality predict its promoter's protection beyond chromatin and expression? (roadmap RM13 (1))

Prior work knows that essential and housekeeping promoters drift less (Teschendorff 2013, West 2013); what would be
new is essentiality predicting protection after the things that travel with it are held fixed. Among broadly
expressed genes only (expressed at 1 TPM or more in at least ``BROAD`` of GTEx v8's 54 tissues), each promoter
CpG's age slope of entropy and of dispersion (``figures tissues``' measures; Hannum blood, Jaffe cortex) is
regressed on its gene's DepMap essentiality: the mean CRISPR gene effect (Chronos) over the cell lines of DepMap
24Q4 (more negative = more essential). Held fixed, in steps:

1. island status × starting-level decile (the matching of every cohort test);
2. + the CpG's chromatin state in the tissue's own Roadmap epigenome (E062 blood, E073 cortex; ``figures
   chromatin``);
3. + the gene's expression in that tissue (log2 median TPM, GTEx whole blood or frontal cortex BA9) and its
   breadth (tissues at 1 TPM or more).

Strata enter as fixed effects (absorbed by demeaning); the metric is in SD; the interval comes from resampling
genes. A coefficient above 0 means more essential genes drift less beyond what is held. Each CpG takes its one gene
(CpGs at several genes' promoters are left out).

``entropy-age figures essentiality_gradient`` as a Slurm job, with ``data/input/depmap/CRISPRGeneEffect_24Q4.csv``
(Figshare 51064667) and ``data/input/gtex/GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_median_tpm.gct.gz``.
"""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import polars as pl

from entropy_age import plots
from entropy_age.data.annotation import AnnotationLayout
from entropy_age.experiments.chromatin_match import NO_STATE, OWN_AND_OTHER
from entropy_age.experiments.chromatin_match import _states as chromatin_states
from entropy_age.experiments.common import Context, write_table
from entropy_age.experiments.tissue_contrast import _hannum, _jaffe, _lifespan_metrics
from entropy_age.stochasticity import baseline_strata

DEPMAP = "depmap/CRISPRGeneEffect_24Q4.csv"
GTEX = "gtex/GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_median_tpm.gct.gz"
TISSUE_TPM = {"blood (Hannum)": "Whole Blood", "cortex (Jaffe)": "Brain - Frontal Cortex (BA9)"}
BROAD = 50
METRICS = ("entropy_slope", "dispersion_slope")
STEPS = ("island x level", "+ chromatin state", "+ expression and breadth")
N_BOOT = {"paper": 200, "quick": 30}


def depmap_effects(path: Path) -> pl.DataFrame:
    """``entrez`` and ``effect``, the mean gene effect over cell lines, from DepMap's ``SYMBOL (ENTREZ)`` columns."""
    frame = pl.read_csv(path)
    genes = frame.columns[1:]
    values = frame.select(genes).to_numpy().astype(float)
    entrez = [int(g.rsplit("(", 1)[1].rstrip(")")) for g in genes]
    return pl.DataFrame({"entrez": entrez, "effect": np.nanmean(values, axis=0)})


def gtex_expression(path: Path, hgnc: pl.DataFrame, tissues: list[str]) -> pl.DataFrame:
    """``entrez``, ``breadth`` (tissues at 1 TPM or more) and log2(TPM + 1) of each of ``tissues``."""
    frame = pl.read_csv(path, separator="\t", skip_rows=2)
    columns = [c for c in frame.columns if c not in ("Name", "Description")]
    tpm = frame.select(columns).to_numpy().astype(float)
    out = pl.DataFrame(
        {
            "ensembl": frame["Name"].str.replace(r"\.\d+$", ""),
            "breadth": (tpm >= 1).sum(axis=1),
            **{t: np.log2(frame[t].to_numpy().astype(float) + 1) for t in tissues},
        }
    )
    ids = hgnc.select("ensembl", "entrez").drop_nulls()
    return out.join(ids, on="ensembl").drop("ensembl").unique(subset="entrez", keep="first")


def cpg_entrez(links: pl.DataFrame, hgnc: pl.DataFrame) -> pl.DataFrame:
    """``cpg`` and ``entrez`` for promoter CpGs linked to exactly one gene (symbols through HGNC, old ones too)."""
    current = hgnc.select("symbol", "entrez")
    old = hgnc.select(pl.col("previous").alias("symbol"), "entrez").explode("symbol").drop_nulls()
    lookup = pl.concat([current, old.filter(~pl.col("symbol").is_in(current["symbol"].implode()))]).unique(
        subset="symbol", keep="first"
    )
    pairs = links.filter(pl.col("promoter")).select("cpg", "symbol").unique().join(lookup, on="symbol")
    single = pairs.group_by("cpg").agg(pl.col("entrez").n_unique().alias("n"), pl.col("entrez").first())
    return single.filter(pl.col("n") == 1).select("cpg", "entrez")


def demean(values: np.ndarray, groups: np.ndarray) -> np.ndarray:
    """``values`` (rows) minus their group's mean: fixed effects absorbed."""
    _, inv = np.unique(groups, return_inverse=True)
    v = np.asarray(values, dtype=float)
    v2 = v.reshape(v.shape[0], -1)
    sums = np.zeros((inv.max() + 1, v2.shape[1]))
    np.add.at(sums, inv, v2)
    counts = np.bincount(inv)[:, None]
    return (v2 - (sums / counts)[inv]).reshape(v.shape)


def slope(y: np.ndarray, x: np.ndarray, covariates: np.ndarray, strata: np.ndarray) -> float:
    """Coefficient of ``x`` in ``y ~ x + covariates``, strata as fixed effects."""
    design = demean(np.column_stack([x, covariates]) if covariates.size else x[:, None], strata)
    coef, *_ = np.linalg.lstsq(design, demean(y, strata), rcond=None)
    return float(coef[0])


def run(ctx: Context) -> list[Path]:
    paths = ctx.paths
    if not (paths.input / DEPMAP).exists() or not (paths.input / GTEX).exists():
        print("essentiality_gradient: DepMap or GTEx file missing; nothing run", file=sys.stderr)
        return []
    plots.apply_style()
    out = ctx.out_dir("essentiality_gradient")
    lay = AnnotationLayout(paths)
    probes = pl.read_parquet(lay.probes_450k).filter(pl.col("promoter") & ~pl.col("chr").is_in(["X", "Y"]))
    hgnc = pl.read_parquet(lay.hgnc)
    genes = (
        cpg_entrez(pl.read_parquet(lay.links_450k), hgnc)
        .join(depmap_effects(paths.input / DEPMAP), on="entrez")
        .join(gtex_expression(paths.input / GTEX, hgnc, list(TISSUE_TPM.values())), on="entrez")
        .filter(pl.col("breadth") >= BROAD)
    )
    features = sorted(set(genes["cpg"].to_list()) & set(probes["cpg"].to_list()))
    print(f"essentiality_gradient: {len(features)} CpGs of {genes['entrez'].n_unique()} broad genes", file=sys.stderr)
    island = dict(probes.select("cpg", "cpg_island").iter_rows())
    rng = ctx.rng(260_000)
    rows, deciles = [], []
    for label, prep in (("blood (Hannum)", _hannum), ("cortex (Jaffe)", _jaffe)):
        p = prep(paths, features)
        metrics = _lifespan_metrics(p)
        feats = metrics["feature"].to_list()
        g = pl.DataFrame({"cpg": feats}).join(genes, on="cpg", how="left", maintain_order="left")
        state = chromatin_states(feats, probes, paths)[OWN_AND_OTHER[label][0]]
        isl = np.array([bool(island.get(f, False)) for f in feats])
        base = np.char.add(isl.astype(str), baseline_strata(metrics["baseline"].to_numpy(), 10).astype(str))
        strata = {STEPS[0]: base, STEPS[1]: np.char.add(base, state), STEPS[2]: np.char.add(base, state)}
        expr = np.column_stack([g[TISSUE_TPM[label]].to_numpy(), g["breadth"].to_numpy().astype(float)])
        cov = {STEPS[0]: np.empty((len(feats), 0)), STEPS[1]: np.empty((len(feats), 0)), STEPS[2]: expr}
        effect = g["effect"].to_numpy()
        entrez = g["entrez"].to_numpy()
        unique_genes = np.unique(entrez)
        by_gene = {e: np.flatnonzero(entrez == e) for e in unique_genes}
        for metric in METRICS:
            y = metrics[metric].to_numpy()
            y = y / y.std()
            for step in STEPS:
                est = slope(y, effect, cov[step], strata[step])
                draws = []
                for _ in range(N_BOOT[str(ctx.scale.name)]):
                    pick = np.concatenate([by_gene[e] for e in rng.choice(unique_genes, unique_genes.size)])
                    draws.append(slope(y[pick], effect[pick], cov[step][pick], strata[step][pick]))
                draws = np.array(draws)
                rows.append(
                    {
                        "tissue": label,
                        "metric": metric,
                        "held": step,
                        "per_unit_effect": est,
                        "low": float(np.quantile(draws, 0.025)),
                        "high": float(np.quantile(draws, 0.975)),
                        "p_boot": float(2 * min((draws <= 0).mean(), (draws >= 0).mean())),
                        "cpgs": len(feats),
                        "genes": int(unique_genes.size),
                        "no_state": int((state == NO_STATE).sum()),
                    }
                )
            resid = demean(y, strata[STEPS[2]])
            cut = np.quantile(effect, np.linspace(0, 1, 11)[1:-1])
            dec = np.searchsorted(cut, effect)
            deciles += [
                {
                    "tissue": label,
                    "metric": metric,
                    "decile": int(d),
                    "effect_mid": float(np.median(effect[dec == d])),
                    "residual": float(resid[dec == d].mean()),
                }
                for d in range(10)
            ]
        print(f"essentiality_gradient: {label} done", file=sys.stderr, flush=True)
    table = pl.DataFrame(rows)
    write_table(table, out / "gradient")
    write_table(pl.DataFrame(deciles), out / "deciles")
    return [_figure(ctx, table, pl.DataFrame(deciles))]


def _figure(ctx: Context, table: pl.DataFrame, deciles: pl.DataFrame) -> Path:
    fig, axes = plots.panels(2, ncols=2, size=4.4)
    ax = axes[0]
    keys = table.select("tissue", "metric").unique(maintain_order=True)
    y = np.arange(keys.height)
    for i, step in enumerate(STEPS):
        t = keys.join(table.filter(pl.col("held") == step), on=["tissue", "metric"], how="left")
        ax.errorbar(
            t["per_unit_effect"],
            y + 0.2 * (i - 1),
            xerr=[t["per_unit_effect"] - t["low"], t["high"] - t["per_unit_effect"]],
            fmt="o",
            ms=4,
            lw=1,
            color=plots.PALETTE[i],
            label=step,
        )
    ax.set_yticks(y, [f"{a} · {m.replace('_slope', '')}" for a, m in keys.iter_rows()])
    ax.axvline(0, color=plots.MUTED, lw=0.8)
    ax.invert_yaxis()
    ax.set(xlabel="drift (SD) per unit of DepMap gene effect", title="Broadly expressed genes")
    ax.legend(fontsize=7, loc="upper center", bbox_to_anchor=(0.5, -0.15), ncol=1, frameon=False)
    plots.label_panel(ax, "a")
    ax = axes[1]
    for i, ((tissue, metric), d) in enumerate(deciles.group_by("tissue", "metric", maintain_order=True)):
        ax.plot(d["decile"] + 1, d["residual"], "o-", ms=3, color=plots.PALETTE[i], label=f"{tissue} · {metric[:4]}")
    ax.axhline(0, color=plots.MUTED, lw=0.8)
    ax.set(
        xlabel="DepMap essentiality decile (1 = most essential)",
        ylabel="drift left after all held (SD)",
        title="Gradient by decile",
    )
    ax.legend(fontsize=7, loc="upper center", bbox_to_anchor=(0.5, -0.15), ncol=2, frameon=False)
    plots.label_panel(ax, "b")
    fig.subplots_adjust(wspace=0.75)
    return plots.save(fig, ctx.figure_path("essentiality_gradient"))
