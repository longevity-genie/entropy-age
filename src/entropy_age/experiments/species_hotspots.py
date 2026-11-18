"""CpGs that become disordered in every mammal? (roadmap RM35 item 3; BEYOND_SHIELD idea 11)

R72 found that in human blood the CpGs whose spread rises fastest with age sit at Polycomb targets, a programme,
while in cortex they are spread broadly. If age drift has a shared programme across mammals, the same conserved
CpGs should be hotspots in many species and tissues; if it is noise on each genome, hotspots should overlap no
more than their starting levels make them.

On the Mammalian Methylation Consortium array (GSE223748; the groups of ``species_shield``: a species and tissue
with at least 15 animals spanning 15 % of the species' lifespan), each autosomal mouse-mapped CpG gets its age
slope of relative entropy from the young state and of spread (``drift_metrics``). The top 5 % per group are its
hotspots. For each CpG, the number of groups where it is a hotspot is set against a null that reshuffles each
group's hotspot labels within its own young-level deciles (keeping how hotspots depend on level). Reported: the
excess of CpGs hot in many groups, the mean pairwise rank correlation of the slopes between groups (within one
order and across orders), and, for the CpGs hot in at least half the groups, the enrichment in CpG islands,
promoters and the essential sets against a level-matched background (``drift_hotspots.enrichment``).

``entropy-age figures species_hotspots`` as a Slurm job.
"""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import polars as pl
from scipy.stats import spearmanr

from entropy_age import plots
from entropy_age.data.mouse import MouseAnnotationLayout, mouse_membership
from entropy_age.experiments.bxd_liver import probes_without_suffix
from entropy_age.experiments.common import Context, write_table
from entropy_age.experiments.drift_hotspots import enrichment
from entropy_age.experiments.species_shield import MIN_AGE_SPAN, MIN_ANIMALS, layout, read_samples
from entropy_age.experiments.tissue_contrast import SETS
from entropy_age.stochasticity import baseline_strata, drift_metrics

TOP = 0.05
METRICS = ("kl_slope", "dispersion_slope")
N_NULL = 200
UNIVERSAL = 0.5  # share of groups


def hot_labels(v: np.ndarray) -> np.ndarray:
    """Top ``TOP`` of the finite values."""
    ok = np.isfinite(v)
    hot = np.zeros(v.size, dtype=bool)
    hot[ok] = v[ok] >= np.quantile(v[ok], 1 - TOP)
    return hot


def shuffled(
    hot: np.ndarray, strata: np.ndarray, rng: np.random.Generator, measured: np.ndarray | None = None
) -> np.ndarray:
    """``hot`` with its labels permuted within each stratum, among the ``measured`` CpGs only."""
    out = hot.copy()
    measured = np.ones(hot.size, dtype=bool) if measured is None else measured
    for s in np.unique(strata):
        k = np.flatnonzero((strata == s) & measured)
        out[k] = hot[k][rng.permutation(k.size)]
    return out


def run(ctx: Context) -> list[Path]:
    matrix, betas_path = layout(ctx.paths)
    if not betas_path.exists():
        print(f"species_hotspots: {betas_path} missing", file=sys.stderr)
        return []
    plots.apply_style()
    out = ctx.out_dir("species_hotspots")
    samples = read_samples(matrix).filter(
        pl.col("age").is_not_null()
        & pl.col("max_lifespan").is_not_null()
        & pl.col("barcode").is_not_null()
        & (pl.col("intervention").is_null() | pl.col("intervention").is_in(["WT", "Control", ""]))
    )
    lay = MouseAnnotationLayout(ctx.paths)
    probes = probes_without_suffix(pl.read_parquet(lay.probes).filter(~pl.col("chr").is_in(["X", "Y", "MT"])))
    header = pl.read_csv(betas_path, n_rows=0).columns
    id_col = header[0]
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
        .select([id_col, *keep["barcode"].unique().to_list()])
        .filter(pl.col(id_col).is_in(probes["cpg"].implode()))
        .collect()
        .rename({id_col: "cpg"})
    )
    features = frame["cpg"].to_list()
    print(f"species_hotspots: {groups.height} groups, {len(features)} CpGs", file=sys.stderr, flush=True)

    slopes = {m: [] for m in METRICS}
    levels, labels = [], []
    for (species, tissue), sel in keep.group_by("species", "tissue", maintain_order=True):
        x = frame.select(sel["barcode"].to_list()).to_numpy().T.astype(float)
        ok = ~np.isnan(x).any(axis=0)
        if ok.sum() < 0.9 * x.shape[1]:
            continue
        x = np.where(np.isnan(x), np.nanmean(x, axis=0), x)
        met = drift_metrics(x, sel["age"].to_numpy(), features, spearman=False)
        for m in METRICS:
            v = met[m].to_numpy().copy()
            v[~ok] = np.nan
            slopes[m].append(v)
        levels.append(met["baseline"].to_numpy())
        labels.append((species, tissue))
    n_groups = len(labels)
    level = np.vstack(levels)
    strata = [baseline_strata(lv, 10) for lv in level]
    order = _orders(ctx, [s for s, _ in labels])
    rng = ctx.rng(35_003)

    summary, per_cpg, corr_rows = [], {"cpg": features}, []
    for m in METRICS:
        v = np.vstack(slopes[m])
        hot = np.vstack([hot_labels(r) for r in v])
        count = hot.sum(axis=0)
        per_cpg[f"hot_groups_{m}"] = count
        null_max = []
        thresholds = (5, 10, int(UNIVERSAL * n_groups))
        null_counts = {t: [] for t in thresholds}
        for _ in range(N_NULL):
            c = np.vstack([shuffled(hot[g], strata[g], rng, np.isfinite(v[g])) for g in range(n_groups)]).sum(axis=0)
            null_max.append(c.max())
            for t in thresholds:
                null_counts[t].append(int((c >= t).sum()))
        for t in thresholds:
            nc = np.array(null_counts[t])
            obs = int((count >= t).sum())
            summary.append(
                {
                    "metric": m,
                    "groups": n_groups,
                    "hot_in_at_least": t,
                    "observed_cpgs": obs,
                    "null_mean": float(nc.mean()),
                    "null_q95": float(np.quantile(nc, 0.95)),
                    "excess_ratio": obs / nc.mean() if nc.mean() > 0 else (None if obs == 0 else float("inf")),
                    "p_null": float(((nc >= obs).sum() + 1) / (N_NULL + 1)),
                }
            )
        # pairwise rank correlation of slopes, within and across orders
        r = np.full((n_groups, n_groups), np.nan)
        for i in range(n_groups):
            for j in range(i + 1, n_groups):
                ok = np.isfinite(v[i]) & np.isfinite(v[j])
                r[i, j] = spearmanr(v[i, ok], v[j, ok]).statistic
        same = np.array([[order[i] == order[j] for j in range(n_groups)] for i in range(n_groups)])
        same_tissue = np.array([[labels[i][1] == labels[j][1] for j in range(n_groups)] for i in range(n_groups)])
        iu = np.triu_indices(n_groups, 1)
        for name, mask in (
            ("all pairs", np.ones_like(same)),
            ("same order", same),
            ("different order", ~same),
            ("same tissue, different order", same_tissue & ~same),
        ):
            vals = r[iu][mask[iu]]
            vals = vals[np.isfinite(vals)]
            corr_rows.append(
                {
                    "metric": m,
                    "pairs": name,
                    "n": vals.size,
                    "mean_r": float(vals.mean()) if vals.size else None,
                    "q25": float(np.quantile(vals, 0.25)) if vals.size else None,
                    "q75": float(np.quantile(vals, 0.75)) if vals.size else None,
                }
            )
    write_table(pl.DataFrame(summary), out / "universality")
    write_table(pl.DataFrame(corr_rows), out / "slope_correlation")
    per = pl.DataFrame(per_cpg).with_columns(mean_level=pl.Series(level.mean(axis=0)))
    per.write_parquet(out / "per_cpg.parquet")

    # what the universal hotspots are, against a level-matched background
    rows = []
    ann = pl.DataFrame({"cpg": features}).join(
        probes.select("cpg", "promoter", "cpg_island").unique("cpg"), on="cpg", how="left"
    )
    members = mouse_membership(probes.filter(pl.col("cpg").is_in(features)), pl.read_parquet(lay.gene_sets))
    masks = {
        "CpG island": ann["cpg_island"].fill_null(False).to_numpy(),
        "promoter": ann["promoter"].fill_null(False).to_numpy(),
    }
    for g in SETS:
        s = set(members.filter(pl.col("group") == g)["feature"].to_list())
        masks[f"set: {g}"] = np.array([c in s for c in features])
    mstrata = baseline_strata(level.mean(axis=0), 10)
    for m in METRICS:
        universal = per[f"hot_groups_{m}"].to_numpy() >= UNIVERSAL * n_groups
        for name, mask in masks.items():
            ratio, z, obs, exp = enrichment(universal, mask, mstrata)
            rows.append(
                {
                    "metric": m,
                    "class": name,
                    "universal_cpgs": int(universal.sum()),
                    "observed": obs,
                    "expected": exp,
                    "ratio": ratio,
                    "z": z,
                }
            )
    write_table(pl.DataFrame(rows), out / "universal_enrichment")
    return [_figure(ctx, per, n_groups)]


MAMMAL_ORDERS = {
    "Primates",
    "Rodentia",
    "Lagomorpha",
    "Scandentia",
    "Carnivora",
    "Artiodactyla",
    "Perissodactyla",
    "Chiroptera",
    "Eulipotyphla",
    "Pholidota",
    "Proboscidea",
    "Sirenia",
    "Hyracoidea",
    "Afrosoricida",
    "Macroscelidea",
    "Tubulidentata",
    "Cingulata",
    "Pilosa",
    "Diprotodontia",
    "Dasyuromorphia",
    "Didelphimorphia",
    "Peramelemorphia",
    "Monotremata",
    "Dermoptera",
}


def _orders(ctx: Context, species: list[str]) -> list[str]:
    """Taxonomic order per species from the NCBI lineages on disk (``fetch taxonomy``); else the species itself."""
    from entropy_age.experiments.species_phylo import lineages_path, organisms

    path = lineages_path(ctx.paths)
    if not path.exists():
        print("species_hotspots: no lineages on disk; species stand in for orders", file=sys.stderr)
        return species
    names = dict(organisms(layout(ctx.paths)[0]).iter_rows())
    lineage = {r["organism"]: r["lineage"].split(";") for r in pl.read_csv(path).iter_rows(named=True)}
    out = []
    for s in species:
        ranks = [t.strip() for t in lineage.get(names.get(s, ""), [])]
        hit = [t for t in ranks if t in MAMMAL_ORDERS]
        out.append(hit[0] if hit else s)
    return out


def _figure(ctx: Context, per: pl.DataFrame, n_groups: int) -> Path:
    fig, axes = plots.panels(2, ncols=2, size=4.4)
    for ax, m, letter in zip(axes, METRICS, "ab", strict=False):
        c = per[f"hot_groups_{m}"].to_numpy()
        ax.hist(c, bins=np.arange(0, c.max() + 2) - 0.5, color=plots.PALETTE[0])
        ax.axvline(UNIVERSAL * n_groups, color=plots.MUTED, lw=0.8, ls="--")
        ax.set_yscale("log")
        ax.set(
            xlabel=f"groups where the CpG is a top-5 % hotspot (of {n_groups})",
            ylabel="CpGs",
            title=m.replace("_", " "),
        )
        plots.label_panel(ax, letter)
    fig.suptitle("Do the same CpGs become disordered across mammals?", x=0.02, y=1.02, ha="left", fontsize=11)
    return plots.save(fig, ctx.figure_path("species_hotspots"))
