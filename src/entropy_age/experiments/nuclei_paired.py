"""Neurons against glia in the same donors: is the shield larger where cells can divide? (roadmap RM40, RM44)

R109 found the matched gap largest in glial nuclei and smallest in neurons of GSE279509 (Walker et al. 2025), each
cell type analysed on its own donors and with intervals that resampled genes. Here the contrast is tested directly.
For each glial fraction, the donors who have both that fraction and NeuN+ nuclei passing QC are kept, the CpGs
measured in every array of both fractions are used, and the gap (set minus matched CpGs, age slope with sex held,
island status x young-level decile strata from that cell type's own young level) is computed in each cell type. The
statistic is the gap in glia minus the gap in neurons. Its interval comes from a donor bootstrap that draws donors
once and takes both of each donor's arrays, so everything a donor shares between fractions (age, sex, post-mortem
handling, most chips) is resampled together; each resample is scaled by the full sample's SD of the measure.

Negative means glia are more shielded than neurons. ``entropy-age figures nuclei_paired`` as a Slurm job.
"""

from __future__ import annotations

import os
import sys
from pathlib import Path

import numpy as np
import polars as pl

from entropy_age.data import reprocessed
from entropy_age.data.annotation import AnnotationLayout
from entropy_age.experiments.brain_replication import WALKER, complete_matrix, walker_groups
from entropy_age.experiments.common import Context, dummies, write_table
from entropy_age.experiments.per_person import _sets
from entropy_age.experiments.tissue_contrast import _features
from entropy_age.protection import matched_scores
from entropy_age.stochasticity import baseline_strata, drift_metrics

NEURONS = f"{WALKER} neurons (NeuN+)"
GLIA = (
    f"{WALKER} other glia (NeuN-/SOX10-)",
    f"{WALKER} oligodendrocyte lineage (SOX10+)",
)
METRICS = ("kl_slope", "entropy_slope", "abs_slope")
N_BOOT = {"paper": 400, "quick": 40}
MIN_DONORS = 40


def gaps(x: np.ndarray, ages: np.ndarray, sex: np.ndarray, feats: list[str], island: np.ndarray, sets: dict) -> dict:
    """(metric, group) -> (unscaled matched difference, SD of the metric over CpGs) for one cell type."""
    cov = sex[:, np.ptp(sex, axis=0) > 1e-12] if sex.shape[1] else sex
    metrics = drift_metrics(x, ages, feats, covariates=cov if cov.shape[1] else None, spearman=False)
    base = metrics["baseline"].to_numpy()
    strata = np.char.add(np.char.add(island.astype(str), "|"), baseline_strata(base, 10).astype(str))
    out = {}
    for metric in METRICS:
        values = metrics[metric].to_numpy()
        for group, (in_set, _) in sets.items():
            out[(metric, group)] = (float(matched_scores(values[None, :], in_set, strata)[0][0]), float(values.std()))
    return out


def paired_difference(point_glia: dict, point_neuron: dict, boots: list[tuple[dict, dict]]) -> list[dict]:
    """Per metric and set: both gaps in SDs of the full sample, their difference and its donor-bootstrap interval."""
    rows = []
    for key, (g, sd_g) in point_glia.items():
        n, sd_n = point_neuron[key]
        draws_g = np.array([b[0][key][0] for b in boots]) / sd_g
        draws_n = np.array([b[1][key][0] for b in boots]) / sd_n
        diff = draws_g - draws_n
        tail = min((diff <= 0).mean(), (diff >= 0).mean()) if diff.size else np.nan

        def interval(v: np.ndarray) -> tuple[float, float]:
            return (float(np.quantile(v, 0.025)), float(np.quantile(v, 0.975))) if v.size else (np.nan, np.nan)

        rows.append(
            {
                "metric": key[0],
                "group": key[1],
                "gap_glia": g / sd_g,
                "gap_glia_low": interval(draws_g)[0],
                "gap_glia_high": interval(draws_g)[1],
                "gap_neurons": n / sd_n,
                "gap_neurons_low": interval(draws_n)[0],
                "gap_neurons_high": interval(draws_n)[1],
                "difference": g / sd_g - n / sd_n,
                "ci_low": interval(diff)[0],
                "ci_high": interval(diff)[1],
                "p_bootstrap": float(min(1.0, max(2 * tail, 1 / max(diff.size, 1)))) if diff.size else np.nan,
                "resamples": int(diff.size),
            }
        )
    return rows


def run(ctx: Context) -> list[Path]:
    if not reprocessed.available(ctx.paths, WALKER):
        print(f"nuclei_paired: no SeSAMe rebuild for {WALKER}; nothing run", file=sys.stderr)
        return []
    out = ctx.out_dir("nuclei_paired")
    n_boot = int(os.environ.get("ENTROPY_AGE_BOOT", N_BOOT[str(ctx.scale.name)]))
    features, _, membership = _features(ctx.paths, ctx.rng(120_000))
    probes = pl.read_parquet(AnnotationLayout(ctx.paths).probes_450k).filter(pl.col("cpg").is_in(features))
    island_of = dict(probes.select("cpg", "cpg_island").iter_rows())
    parts = walker_groups(ctx.paths.input / "geo" / WALKER / f"{WALKER}_series_matrix.txt.gz")
    qc = reprocessed.read_qc(ctx.paths, WALKER)
    betas = reprocessed.by_accession(ctx.paths, WALKER, {s: s for s in qc["sample"].to_list()})
    measured = set(betas.columns)
    neurons = parts[NEURONS].filter(pl.col("sample").is_in(measured)).drop_nulls("age")
    rng = ctx.rng(310_000)
    rows, sizes = [], []
    for label in GLIA:
        glia = parts[label].filter(pl.col("sample").is_in(measured)).drop_nulls("age")
        both = sorted(set(glia["donor"].to_list()) & set(neurons["donor"].to_list()))
        if len(both) < MIN_DONORS:
            print(f"nuclei_paired: {label}: {len(both)} donors with both fractions, skipped", file=sys.stderr)
            continue
        g = glia.filter(pl.col("donor").is_in(both)).sort("donor")
        n = neurons.filter(pl.col("donor").is_in(both)).sort("donor")
        x_all, feats = complete_matrix(betas, features, [*g["sample"].to_list(), *n["sample"].to_list()])
        xg, xn = x_all[: g.height], x_all[g.height :]
        ages = g["age"].to_numpy()
        sex = dummies(g["sex"].fill_null("NA").to_numpy())
        island = np.array([bool(island_of.get(f, False)) for f in feats])
        sets = _sets(feats, membership)
        point = (gaps(xg, ages, sex, feats, island, sets), gaps(xn, ages, sex, feats, island, sets))
        boots = []
        for _ in range(n_boot):
            idx = rng.integers(0, len(both), len(both))
            boots.append(
                (
                    gaps(xg[idx], ages[idx], sex[idx], feats, island, sets),
                    gaps(xn[idx], ages[idx], sex[idx], feats, island, sets),
                )
            )
        name = label.removeprefix(f"{WALKER} ")
        rows += [{"glia": name, **r} for r in paired_difference(point[0], point[1], boots)]
        sizes.append(
            {
                "glia": name,
                "donors": len(both),
                "age_min": float(ages.min()),
                "age_max": float(ages.max()),
                "cpgs": len(feats),
                "resamples": n_boot,
            }
        )
        print(f"nuclei_paired: {name}: {len(both)} donors, {len(feats)} CpGs", file=sys.stderr, flush=True)
    write_table(pl.DataFrame(sizes), out / "sizes")
    write_table(pl.DataFrame(rows), out / "paired")
    return []


# ----------------------------------------------------------------------------- a second series (low power)
GASPARONI = "GSE66351"


def run_gasparoni(ctx: Context) -> list[Path]:
    """The same contrast in GSE66351 (Gasparoni et al. 2018; 450k): NeuN- glia minus NeuN+ neurons of the donors who
    have both, sex and diagnosis (Alzheimer's or control) held. A check of low power: 31 donors, 29 of them over 60,
    so the age slope rests on two young donors."""
    from entropy_age.data.cohorts import _characteristics

    if not reprocessed.available(ctx.paths, GASPARONI):
        print(f"nuclei_paired: no SeSAMe rebuild for {GASPARONI}; nothing run", file=sys.stderr)
        return []
    out = ctx.out_dir("nuclei_paired_gse66351")
    n_boot = int(os.environ.get("ENTROPY_AGE_BOOT", N_BOOT[str(ctx.scale.name)]))
    features, _, membership = _features(ctx.paths, ctx.rng(120_000))
    probes = pl.read_parquet(AnnotationLayout(ctx.paths).probes_450k).filter(pl.col("cpg").is_in(features))
    island_of = dict(probes.select("cpg", "cpg_island").iter_rows())
    qc = reprocessed.read_qc(ctx.paths, GASPARONI).filter(pl.col("passed"))
    chars = _characteristics(ctx.paths.input / "geo" / GASPARONI / f"{GASPARONI}_series_matrix.txt.gz")
    s = qc.select("sample").join(chars, on="sample").with_columns(age=pl.col("age").cast(pl.Float64, strict=False))
    s = s.filter(pl.col("cell type").is_in(["Neuron", "Glia"])).drop_nulls("age")
    both = sorted(
        s.group_by("donor_id").agg(pl.col("cell type").n_unique().alias("k")).filter(pl.col("k") == 2)["donor_id"]
    )
    g = s.filter((pl.col("cell type") == "Glia") & pl.col("donor_id").is_in(both)).sort("donor_id")
    n = s.filter((pl.col("cell type") == "Neuron") & pl.col("donor_id").is_in(both)).sort("donor_id")
    g, n = (d.unique(subset="donor_id", keep="first", maintain_order=True) for d in (g, n))
    betas = reprocessed.by_accession(ctx.paths, GASPARONI, {v: v for v in s["sample"].to_list()})
    x_all, feats = complete_matrix(betas, features, [*g["sample"].to_list(), *n["sample"].to_list()])
    xg, xn = x_all[: g.height], x_all[g.height :]
    ages = g["age"].to_numpy()
    cov = np.column_stack(
        [dummies(g["Sex" if "Sex" in g.columns else "sex"].to_numpy()), dummies(g["diagnosis"].to_numpy())]
    )
    island = np.array([bool(island_of.get(f, False)) for f in feats])
    sets = _sets(feats, membership)
    point = (gaps(xg, ages, cov, feats, island, sets), gaps(xn, ages, cov, feats, island, sets))
    rng = ctx.rng(311_000)
    boots = []
    for _ in range(n_boot):
        idx = rng.integers(0, len(both), len(both))
        if np.ptp(ages[idx]) < 20:  # a resample without a young donor carries no age slope
            continue
        boots.append(
            (
                gaps(xg[idx], ages[idx], cov[idx], feats, island, sets),
                gaps(xn[idx], ages[idx], cov[idx], feats, island, sets),
            )
        )
    rows = [{"glia": "glia (NeuN-)", **r} for r in paired_difference(point[0], point[1], boots)]
    write_table(pl.DataFrame(rows), out / "paired")
    write_table(
        pl.DataFrame(
            {
                "donors": [len(both)],
                "under_60": [int((ages < 60).sum())],
                "age_min": [float(ages.min())],
                "age_max": [float(ages.max())],
                "cpgs": [len(feats)],
                "resamples_used": [len(boots)],
            }
        ),
        out / "sizes",
    )
    print(f"nuclei_paired: {GASPARONI}: {len(both)} donors, {len(feats)} CpGs, {len(boots)} resamples", file=sys.stderr)
    return []
