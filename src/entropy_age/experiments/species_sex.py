"""Drift and the shield by sex across mammals (roadmap RM39 item 7).

R78 found no sex difference in human tissues on relative entropy, and slightly weaker spread in women in two.
The literature check (docs/LITERATURE_COVERAGE.md, idea 10) left open drift by sex in other species. The
Mammalian Methylation Consortium (GSE223748) records each animal's sex (``female``). For every species and
tissue of ``species_shield`` (no intervention, an age span of at least 15 % of the species' lifespan), with at
least ``MIN_PER_SEX`` animals of each sex, ``species_shield.group_gap`` is run on females and on males apart: the
look-alike CpGs' drift rate per year (``background_rate``) and each set's matched gap. Across groups: the
female-to-male ratio of the look-alike rate (log, its median and sign share, a Wilcoxon signed-rank p) and the
difference of the gap, in SD of the group's pooled metric.

``entropy-age figures species_sex`` as a Slurm job (the consortium's beta table, about 4.5 GB).
"""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import polars as pl
from scipy import stats

from entropy_age.data.cohorts import _characteristics
from entropy_age.data.mouse import MouseAnnotationLayout, mouse_membership
from entropy_age.experiments.bxd_liver import probes_without_suffix
from entropy_age.experiments.common import Context, write_table
from entropy_age.experiments.interventions import _sets
from entropy_age.experiments.species_shield import MIN_AGE_SPAN, group_gap, layout, read_samples
from entropy_age.experiments.tissue_contrast import SETS
from entropy_age.stochasticity import drift_metrics

MIN_PER_SEX = 12
N_RESAMPLES = {"paper": 50, "quick": 10}
METRICS = ("dispersion_slope", "kl_slope")


def sex_of(value) -> str | None:
    """``F``/``M`` from the consortium's ``female`` field (1/0, or words); None when unknown."""
    v = str(value).strip().lower() if value is not None else ""
    if v in ("1", "1.0", "true", "female", "f"):
        return "F"
    if v in ("0", "0.0", "false", "male", "m"):
        return "M"
    return None


def summarise(gaps: pl.DataFrame) -> pl.DataFrame:
    """Per metric and set, across groups: female/male look-alike rate ratio and gap difference."""
    wide = gaps.pivot(
        on="sex", index=["species", "tissue", "metric", "group"], values=["background_rate", "gap_pooled_sd"]
    ).drop_nulls()
    rows = []
    for (metric, group), g in wide.group_by("metric", "group", maintain_order=True):
        f, m = g["background_rate_F"].to_numpy(), g["background_rate_M"].to_numpy()
        ok = (f > 0) & (m > 0)
        ratio = np.log(f[ok] / m[ok])
        diff = (g["gap_pooled_sd_F"] - g["gap_pooled_sd_M"]).to_numpy()
        rows.append(
            {
                "metric": metric,
                "group": group,
                "groups": g.height,
                "species": g["species"].n_unique(),
                "rate_groups": int(ok.sum()),
                "median_log_rate_ratio_f_m": float(np.median(ratio)) if ratio.size else None,
                "share_rate_higher_in_f": float((ratio > 0).mean()) if ratio.size else None,
                "p_rate": float(stats.wilcoxon(ratio).pvalue) if ratio.size >= 6 else None,
                "median_gap_diff_f_m": float(np.median(diff)),
                "share_gap_more_negative_in_f": float((diff < 0).mean()),
                "p_gap": float(stats.wilcoxon(diff).pvalue) if diff.size >= 6 else None,
            }
        )
    return pl.DataFrame(rows)


def run(ctx: Context) -> list[Path]:
    matrix, betas_path = layout(ctx.paths)
    if not betas_path.exists():
        print(f"species_sex: {betas_path} missing", file=sys.stderr)
        return []
    out = ctx.out_dir("species_sex")
    sexes = _characteristics(matrix).select("sample", pl.col("female").map_elements(sex_of, return_dtype=pl.String))
    samples = (
        read_samples(matrix)
        .join(sexes.rename({"female": "sex"}), on="sample", how="left")
        .filter(
            pl.col("age").is_not_null()
            & pl.col("max_lifespan").is_not_null()
            & pl.col("barcode").is_not_null()
            & pl.col("sex").is_not_null()
            & (pl.col("intervention").is_null() | pl.col("intervention").is_in(["WT", "Control", ""]))
        )
    )
    header = pl.read_csv(betas_path, n_rows=0).columns
    id_col = header[0]
    samples = samples.filter(pl.col("barcode").is_in(header))
    per_sex = (
        samples.group_by("species", "tissue", "sex")
        .agg(
            pl.len().alias("n"),
            (pl.col("age").max() - pl.col("age").min()).alias("span"),
            pl.col("max_lifespan").first(),
        )
        .filter((pl.col("n") >= MIN_PER_SEX) & (pl.col("span") >= MIN_AGE_SPAN * pl.col("max_lifespan")))
    )
    both = per_sex.group_by("species", "tissue").len().filter(pl.col("len") == 2).select("species", "tissue")
    keep = samples.join(both, on=["species", "tissue"])
    print(f"species_sex: {both.height} species-tissue groups with both sexes, {keep.height} animals", file=sys.stderr)
    lay = MouseAnnotationLayout(ctx.paths)
    probes = probes_without_suffix(
        pl.read_parquet(lay.probes).filter(pl.col("promoter") & ~pl.col("chr").is_in(["X", "Y", "MT"]))
    )
    frame = (
        pl.scan_csv(betas_path, null_values=["NA", ""])
        .select([id_col, *keep["barcode"].unique().to_list()])
        .filter(pl.col(id_col).is_in(probes["cpg"].implode()))
        .collect()
        .rename({id_col: "cpg"})
    )
    features = frame["cpg"].to_list()
    members = mouse_membership(probes.filter(pl.col("cpg").is_in(features)), pl.read_parquet(lay.gene_sets)).filter(
        pl.col("group").is_in(SETS)
    )
    island = dict(probes.select("cpg", "cpg_island").iter_rows())
    isl = np.array([bool(island.get(f, False)) for f in features])
    sets = _sets(features, members)
    rng = ctx.rng(390_007)
    rows = []
    for (species, tissue), sel in keep.group_by("species", "tissue", maintain_order=True):
        x_all = frame.select(sel["barcode"].to_list()).to_numpy().T.astype(np.float32)
        ok = ~np.isnan(x_all).any(axis=0)
        if ok.sum() < 0.9 * x_all.shape[1]:
            continue
        cols = np.flatnonzero(ok)
        sub_sets = {
            g: (s[cols], [np.searchsorted(cols, gi[np.isin(gi, cols)]) for gi in genes if np.isin(gi, cols).any()])
            for g, (s, genes) in sets.items()
        }
        feats = [features[c] for c in cols]
        pooled = drift_metrics(x_all[:, cols], sel["age"].to_numpy(), feats, spearman=False)
        pooled_sd = {m: float(pooled[m].to_numpy().std()) for m in METRICS}
        for sex in ("F", "M"):
            rows_sex = (sel["sex"] == sex).to_numpy()
            gaps = group_gap(
                x_all[rows_sex][:, cols],
                sel["age"].to_numpy()[rows_sex],
                feats,
                isl[cols],
                sub_sets,
                N_RESAMPLES[str(ctx.scale.name)],
                rng,
            )
            rows += [
                {
                    "species": species,
                    "tissue": tissue,
                    "sex": sex,
                    "animals": int(rows_sex.sum()),
                    **g,
                    "gap_pooled_sd": g["gap"] / pooled_sd[g["metric"]],
                }
                for g in gaps
                if g["metric"] in METRICS
            ]
        print(f"species_sex: {species}, {tissue}", file=sys.stderr, flush=True)
    gaps = pl.DataFrame(rows)
    write_table(gaps, out / "gaps_by_sex")
    write_table(summarise(gaps), out / "summary")
    return []
