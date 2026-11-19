"""Promoter activity measured in the same cells: is the shield what expression predicts? (roadmap RM40 item 4)

The corrected reading of this project is that the shield is a property of active promoters, not of essential genes
(R108). Until now activity came from reference maps of other people's cells. MESA measured CD14+ monocytes of the
same 1,202 people twice: methylation on the 450k array (GSE56046, deposited values, no IDATs) and expression on the
HT-12 array (GSE56045); the samples join one to one on the participant id in the title.

1. **Matching on the gene's own expression.** Per-CpG age slopes as in ``sorted_cells`` (race-sex-site and
   contamination scores held). Set CpGs against matched CpGs under four matchings: island x level; + the decile of
   the gene's mean expression in these monocytes (a CpG whose gene has no probe on the expression array is its own
   class); + the chromatin state of blood mononuclear cells (Roadmap E062); + both. The share of the gap each
   removes says how much of the shield expression in the same cells accounts for.
2. **A gradient outside the gene sets.** Among promoter CpGs of genes in no set, each CpG's slope minus the mean of
   its island x level stratum, by expression decile, and its rank correlation with expression across genes. If the
   shield is promoter activity, well-expressed genes outside every set are shielded too.
3. **Positive control.** Per gene, the correlation across people of promoter methylation and expression, both
   residualised on age and the covariates. Promoter methylation silences, so the mean should be negative; if it
   is not, the two series are not joined correctly.

Expression per gene: the probe with the highest mean among the gene's probes. ``entropy-age figures
mesa_expression`` as a Slurm job.
"""

from __future__ import annotations

import gzip
import sys
from pathlib import Path

import numpy as np
import polars as pl
from scipy import stats

from entropy_age.data.annotation import AnnotationLayout
from entropy_age.experiments.chromatin_match import N_RESAMPLES, matched_gaps
from entropy_age.experiments.chromatin_match import _states as chromatin_states
from entropy_age.experiments.common import Context, write_table
from entropy_age.experiments.per_person import _sets
from entropy_age.experiments.shield_clocks import residual
from entropy_age.experiments.sorted_cells import covariates, normalized_path, read_betas, read_people
from entropy_age.experiments.syndromes import matrix_path
from entropy_age.experiments.tissue_contrast import Prepared, _features, _lifespan_metrics
from entropy_age.stochasticity import baseline_strata

METHYLATION, EXPRESSION = "GSE56046", "GSE56045"
PLATFORM = "GPL10558.annot.gz"
METRICS = ("kl_slope", "entropy_slope", "dispersion_slope")
NO_PROBE = "none"
BINS = 10
PERMUTATIONS = {"paper": 200, "quick": 20}


def read_expression(matrix: Path) -> tuple[list[str], list[str], np.ndarray]:
    """Probe ids, participant ids (from the sample titles) and the probes x samples table of a series matrix."""
    titles, skip = [], 0
    with gzip.open(matrix, "rt") as handle:
        for i, line in enumerate(handle):
            if line.startswith("!Sample_title"):
                titles = [c.strip().strip('"').split("_")[0] for c in line.rstrip("\n").split("\t")[1:]]
            if line.startswith("!series_matrix_table_begin"):
                skip = i + 1
                break
    frame = pl.read_csv(matrix, separator="\t", skip_rows=skip, comment_prefix="!", null_values=["", "null", "NA"])
    values = frame.select(frame.columns[1:]).to_numpy().astype(float)
    return frame[frame.columns[0]].to_list(), titles, values


def read_platform(path: Path) -> dict[str, str]:
    """Expression probe -> gene symbol from GEO's platform annotation; probes of no or several genes left out."""
    out = {}
    with gzip.open(path, "rt", errors="replace") as handle:
        started = False
        for line in handle:
            if not started:
                started = line.startswith("ID\t")
                continue
            cells = line.rstrip("\n").split("\t")
            if len(cells) > 2 and cells[2] and "///" not in cells[2]:
                out[cells[0]] = cells[2]
    return out


def gene_level(probes: list[str], values: np.ndarray, symbol_of: dict[str, str]) -> tuple[list[str], np.ndarray]:
    """Genes and their genes x samples expression: per gene the probe with the highest mean."""
    means = np.nanmean(values, axis=1)
    best: dict[str, int] = {}
    for i, probe in enumerate(probes):
        gene = symbol_of.get(probe)
        if gene is not None and (gene not in best or means[i] > means[best[gene]]):
            best[gene] = i
    genes = sorted(best)
    return genes, values[[best[g] for g in genes]]


def cpg_expression(features: list[str], links: pl.DataFrame, gene_mean: dict[str, float]) -> tuple[np.ndarray, list]:
    """Per CpG: the mean expression of its promoter's gene (the best expressed, when several) and that gene; NaN and
    ``None`` when no gene of the CpG has a probe on the expression array."""
    by_cpg: dict[str, list[str]] = {}
    for cpg, symbol in links.filter(pl.col("promoter")).select("cpg", "symbol").iter_rows():
        by_cpg.setdefault(cpg, []).append(symbol)
    level, gene = np.full(len(features), np.nan), [None] * len(features)
    for i, f in enumerate(features):
        known = [g for g in by_cpg.get(f, []) if g in gene_mean]
        if known:
            gene[i] = max(known, key=gene_mean.__getitem__)
            level[i] = gene_mean[gene[i]]
    return level, gene


def expression_bins(level: np.ndarray, n: int = BINS) -> np.ndarray:
    """``"none"`` for a CpG without an expression value, else its bin ``"1"``..``"n"`` of expression among the rest."""
    out = np.full(level.size, NO_PROBE, dtype="<U4")
    known = np.isfinite(level)
    if known.any():
        out[known] = (baseline_strata(level[known], n) + 1).astype(str)
    return out


def stratum_residual(values: np.ndarray, strata: np.ndarray) -> np.ndarray:
    """Each value minus the mean of its stratum."""
    out = np.empty(values.size)
    for s in np.unique(strata):
        sel = strata == s
        out[sel] = values[sel] - values[sel].mean()
    return out


def gradient(values: np.ndarray, base: np.ndarray, bins: np.ndarray, genes: list, eligible: np.ndarray) -> list[dict]:
    """Among ``eligible`` CpGs: the slope relative to its island x level stratum (strata means from eligible CpGs),
    in SDs of ``values``, per expression bin (mean over genes, SE over genes), and the rank correlation with the bin
    across genes."""
    idx = np.flatnonzero(eligible & (bins != NO_PROBE))
    resid = stratum_residual(values[idx], base[idx]) / values.std()
    frame = (
        pl.DataFrame({"gene": [genes[i] for i in idx], "bin": bins[idx].astype(int), "resid": resid})
        .group_by("gene")
        .agg(pl.col("bin").first(), pl.col("resid").mean())
    )
    rho, p = stats.spearmanr(frame["bin"].to_numpy(), frame["resid"].to_numpy())
    rows = [
        {
            "bin": int(b),
            "genes": sel.height,
            "relative_slope_sd": float(sel["resid"].mean()),
            "se": float(sel["resid"].std() / np.sqrt(sel.height)) if sel.height > 1 else float("nan"),
            "spearman": float(rho),
            "p_spearman": float(p),
        }
        for (b,), sel in frame.sort("bin").group_by("bin", maintain_order=True)
    ]
    return rows


def coupling(meth: np.ndarray, expr: np.ndarray, design: np.ndarray, n_perm: int, rng) -> dict:
    """Mean over genes of the across-people correlation of methylation and expression residuals (genes x people
    each), with a null from permuting the people of one layer (one permutation for all genes)."""

    def standardise(a: np.ndarray) -> np.ndarray:
        r = np.vstack([residual(row, design) for row in a])
        r -= r.mean(axis=1, keepdims=True)
        sd = r.std(axis=1, keepdims=True)
        return r / np.where(sd > 0, sd, np.nan)

    zm, ze = standardise(meth), standardise(expr)
    r = np.nanmean(zm * ze, axis=1)
    null = np.array([np.nanmean(np.nanmean(zm * ze[:, rng.permutation(ze.shape[1])], axis=1)) for _ in range(n_perm)])
    observed = float(np.nanmean(r))
    return {
        "genes": int(np.isfinite(r).sum()),
        "mean_r": observed,
        "share_negative": float(np.nanmean(r < 0)),
        "null_mean": float(null.mean()),
        "null_sd": float(null.std()),
        "p_permutation": float((1 + (np.abs(null) >= abs(observed)).sum()) / (1 + n_perm)),
    }


def run(ctx: Context) -> list[Path]:
    paths = ctx.paths
    expr_dir = paths.input / "geo" / EXPRESSION
    needed = [normalized_path(paths, METHYLATION), matrix_path(paths, EXPRESSION), expr_dir / PLATFORM]
    if not all(p.exists() for p in needed):
        print(
            f"mesa_expression: missing {[str(p.name) for p in needed if not p.exists()]}. Nothing run.", file=sys.stderr
        )
        return []
    out = ctx.out_dir("mesa_expression")
    scale = str(ctx.scale.name)
    features, probes, membership = _features(paths, ctx.rng(140_000))
    people = read_people(matrix_path(paths, METHYLATION)).drop_nulls("age")
    probe_ids, titles, table = read_expression(matrix_path(paths, EXPRESSION))
    column = {t: j for j, t in enumerate(titles)}
    people = people.filter(pl.col("id").is_in(titles))
    ids = people["id"].to_list()
    genes, expr = gene_level(probe_ids, table[:, [column[i] for i in ids]], read_platform(expr_dir / PLATFORM))
    gene_mean = dict(zip(genes, np.nanmean(expr, axis=1).tolist(), strict=True))
    print(f"mesa_expression: {len(ids)} people with both layers, {len(genes)} genes", file=sys.stderr, flush=True)

    x = read_betas(normalized_path(paths, METHYLATION), features, ids)
    ok = ~np.isnan(x).any(axis=0)
    x, feats = x[:, ok], [f for f, k in zip(features, ok, strict=True) if k]
    cov = covariates(people)
    ages = people["age"].to_numpy()
    metrics = _lifespan_metrics(Prepared("CD14+ monocytes", x, ages, cov, feats, np.array(ids), cov))
    island_of = dict(probes.select("cpg", "cpg_island").iter_rows())
    island = np.array([bool(island_of.get(f, False)) for f in feats])
    level, gene = cpg_expression(feats, pl.read_parquet(AnnotationLayout(paths).links_450k), gene_mean)
    bins = expression_bins(level)
    state = chromatin_states(feats, probes, paths)["E062"]
    base = np.char.add(
        np.char.add(island.astype(str), "|"), baseline_strata(metrics["baseline"].to_numpy(), 10).astype(str)
    )
    with_expr = np.char.add(np.char.add(base, "|"), bins)
    variants = {
        "island x level": base,
        "+ own expression decile": with_expr,
        "+ chromatin state (E062)": np.char.add(np.char.add(base, "|"), state),
        "+ expression + state": np.char.add(np.char.add(with_expr, "|"), state),
    }
    sets = _sets(feats, membership.filter(pl.col("feature").is_in(feats)))
    in_any = np.any([s for s, _ in sets.values()], axis=0)
    gaps, grad = [], []
    for metric in METRICS:
        values = metrics[metric].to_numpy()
        gaps += [{"metric": metric, **r} for r in matched_gaps(values, sets, variants, N_RESAMPLES[scale], ctx.seed)]
        grad += [{"metric": metric, **r} for r in gradient(values, base, bins, gene, ~in_any)]
    write_table(pl.DataFrame(gaps), out / "matched_gaps")
    write_table(pl.DataFrame(grad), out / "gradient")

    # set CpGs by expression bin: where the sets sit
    where = [
        {"group": g, "bin": b, "cpgs": int(((bins == b) & s).sum()), "share": float(((bins == b) & s).sum() / s.sum())}
        for g, (s, _) in {**sets, "outside the sets": (~in_any, [])}.items()
        for b in [NO_PROBE, *map(str, range(1, BINS + 1))]
    ]
    write_table(pl.DataFrame(where), out / "set_expression")

    # positive control: per gene, promoter methylation against expression across people
    row_of = {g: i for i, g in enumerate(genes)}
    by_gene: dict[str, list[int]] = {}
    for j, g in enumerate(gene):
        if g is not None:
            by_gene.setdefault(g, []).append(j)
    used = sorted(by_gene)
    meth = np.vstack([x[:, by_gene[g]].mean(axis=1) for g in used])
    design = np.column_stack([np.ones(len(ids)), ages, cov[:, np.ptp(cov, axis=0) > 1e-12]])
    row = coupling(meth, expr[[row_of[g] for g in used]], design, PERMUTATIONS[scale], ctx.rng(141_000))
    write_table(pl.DataFrame([{**row, "people": len(ids), "cpgs": len(feats)}]), out / "coupling")
    return []
