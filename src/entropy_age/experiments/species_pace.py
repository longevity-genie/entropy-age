"""Directed pace against noise pace, tissue by tissue, across mammals (roadmap RM39 item 3).

R40 (``tissue_pace``) found that tissues rank alike in mouse and rat RNA by how far expression moves with age, but
not by how much animals come to differ, and that human methylation and RNA give unrelated orders. The literature
check (docs/LITERATURE_COVERAGE.md, idea 5) left open the same split in methylation across species. In the
Mammalian Methylation Consortium (GSE223748; ``species_shield``'s groups), per species and tissue, over the
look-alike CpGs (promoter CpGs in none of the three sets): the **directed pace**, the median absolute age slope of
the mean, and the **noise pace**, the median age slope of dispersion; both per year and per lifespan (times the
species' maximum lifespan). Then:

* within each species with at least ``MIN_TISSUES`` tissues: the Spearman correlation of the tissues' directed
  and noise paces (do the tissues that change fastest also become noisiest?);
* across species, for each tissue pair present in at least ``MIN_SPECIES`` species: how often the pair is ordered
  alike, on each pace (is a tissue order shared across species?).

``entropy-age figures species_pace`` as a Slurm job (the consortium's beta table, about 4.5 GB).
"""

from __future__ import annotations

import itertools
import sys
from pathlib import Path

import numpy as np
import polars as pl
from scipy import stats

from entropy_age.data.mouse import MouseAnnotationLayout, mouse_membership
from entropy_age.experiments.bxd_liver import probes_without_suffix
from entropy_age.experiments.common import Context, write_table
from entropy_age.experiments.species_shield import MIN_AGE_SPAN, MIN_ANIMALS, layout, read_samples
from entropy_age.experiments.tissue_contrast import SETS
from entropy_age.stochasticity import drift_metrics

MIN_TISSUES = 3
MIN_SPECIES = 5
PACES = ("directed", "noise")


def within_species(paces: pl.DataFrame, unit: str) -> pl.DataFrame:
    rows = []
    for (species,), g in paces.group_by("species", maintain_order=True):
        if g.height < MIN_TISSUES:
            continue
        r = stats.spearmanr(g[f"directed_{unit}"], g[f"noise_{unit}"])
        rows.append({"species": species, "unit": unit, "tissues": g.height, "rho": float(r.statistic)})
    return pl.DataFrame(rows)


def pair_agreement(paces: pl.DataFrame, unit: str) -> pl.DataFrame:
    """Per tissue pair seen in ``MIN_SPECIES`` species: share of species where the first tissue is faster."""
    rows = []
    tissues = sorted(paces["tissue"].unique().to_list())
    for a, b in itertools.combinations(tissues, 2):
        ga = paces.filter(pl.col("tissue") == a).select(
            "species", *(pl.col(f"{p}_{unit}").alias(f"{p}_a") for p in PACES)
        )
        gb = paces.filter(pl.col("tissue") == b).select(
            "species", *(pl.col(f"{p}_{unit}").alias(f"{p}_b") for p in PACES)
        )
        j = ga.join(gb, on="species")
        if j.height < MIN_SPECIES:
            continue
        row = {"tissue_a": a, "tissue_b": b, "unit": unit, "species": j.height}
        for p in PACES:
            k = int((j[f"{p}_a"] > j[f"{p}_b"]).sum())
            row[f"{p}_share_a_faster"] = k / j.height
            row[f"{p}_p"] = float(stats.binomtest(k, j.height).pvalue)
        rows.append(row)
    return pl.DataFrame(rows)


def run(ctx: Context) -> list[Path]:
    matrix, betas_path = layout(ctx.paths)
    if not betas_path.exists():
        print(f"species_pace: {betas_path} missing", file=sys.stderr)
        return []
    out = ctx.out_dir("species_pace")
    samples = read_samples(matrix).filter(
        pl.col("age").is_not_null()
        & pl.col("max_lifespan").is_not_null()
        & pl.col("barcode").is_not_null()
        & (pl.col("intervention").is_null() | pl.col("intervention").is_in(["WT", "Control", ""]))
    )
    header = pl.read_csv(betas_path, n_rows=0).columns
    id_col = header[0]
    samples = samples.filter(pl.col("barcode").is_in(header))
    groups = (
        samples.group_by("species", "tissue")
        .agg(
            pl.len().alias("n"),
            (pl.col("age").max() - pl.col("age").min()).alias("span"),
            pl.col("max_lifespan").first(),
        )
        .filter((pl.col("n") >= MIN_ANIMALS) & (pl.col("span") >= MIN_AGE_SPAN * pl.col("max_lifespan")))
    )
    keep = samples.join(groups.select("species", "tissue"), on=["species", "tissue"])
    lay = MouseAnnotationLayout(ctx.paths)
    probes = probes_without_suffix(
        pl.read_parquet(lay.probes).filter(pl.col("promoter") & ~pl.col("chr").is_in(["X", "Y", "MT"]))
    )
    members = mouse_membership(probes, pl.read_parquet(lay.gene_sets)).filter(pl.col("group").is_in(SETS))
    lookalike = probes.filter(~pl.col("cpg").is_in(members["feature"].implode()))["cpg"]
    frame = (
        pl.scan_csv(betas_path, null_values=["NA", ""])
        .select([id_col, *keep["barcode"].unique().to_list()])
        .filter(pl.col(id_col).is_in(lookalike.implode()))
        .collect()
        .rename({id_col: "cpg"})
    )
    features = frame["cpg"].to_list()
    print(f"species_pace: {groups.height} groups, {len(features)} look-alike CpGs", file=sys.stderr, flush=True)
    rows = []
    for (species, tissue), sel in keep.group_by("species", "tissue", maintain_order=True):
        x = frame.select(sel["barcode"].to_list()).to_numpy().T.astype(np.float32)
        ok = ~np.isnan(x).any(axis=0)
        if ok.sum() < 0.9 * x.shape[1]:
            continue
        m = drift_metrics(
            x[:, ok], sel["age"].to_numpy(), [f for f, k in zip(features, ok, strict=True) if k], spearman=False
        )
        life = float(sel["max_lifespan"][0])
        directed = float(m["abs_slope"].median())
        noise = float(m["dispersion_slope"].median())
        rows.append(
            {
                "species": species,
                "tissue": tissue,
                "animals": sel.height,
                "max_lifespan": life,
                "directed_year": directed,
                "noise_year": noise,
                "directed_life": directed * life,
                "noise_life": noise * life,
            }
        )
    paces = pl.DataFrame(rows)
    write_table(paces, out / "paces")
    within = pl.concat([within_species(paces, u) for u in ("year", "life")])
    write_table(within, out / "within_species")
    write_table(pl.concat([pair_agreement(paces, u) for u in ("year", "life")], how="diagonal"), out / "tissue_pairs")
    for u in ("year", "life"):
        w = within.filter(pl.col("unit") == u)["rho"].to_numpy()
        if w.size:
            print(f"species_pace: per {u}: median rho {np.median(w):.2f} over {w.size} species", file=sys.stderr)
    return []
