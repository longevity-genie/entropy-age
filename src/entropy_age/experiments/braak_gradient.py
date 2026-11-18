"""Does disorder follow Alzheimer's regional order within the same brains? (roadmap RM39 item 2)

R39 (``health_states``) found no effect of Alzheimer's or Braak stage on drift in any of Lunnon's regions taken one
at a time. Braak staging says tangles reach the entorhinal cortex first, then the temporal neocortex, the frontal
cortex last, and the cerebellum hardly at all. If pathology drives local disorder, the regions should separate
within a donor, and more so at a higher stage, even when donor-level differences (age, agonal state, the person's
overall drift) hide it between donors.

Lunnon GSE59685, donors with entorhinal cortex, superior temporal gyrus, frontal cortex and cerebellum. Per region,
every donor's scores (``health_states.loo_scores``: drift ``B``, shield share ``F``, signed shift ``S`` per gene
set, each against a reference fitted on the other donors) are standardised over donors, without and with 5
principal components removed (Lunnon ships no neuron fraction; the components carry the neuron-glia balance).
Within each donor the scores are regressed on the region's place in the Braak order (cerebellum 0, frontal 1,
temporal 2, entorhinal 3): that slope is the donor's regional gradient, free of everything about the donor. The
gradient is then regressed on Braak stage with age and sex. Also each region against cerebellum, per stage.

``entropy-age figures braak_gradient`` as a Slurm job.
"""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import polars as pl

from entropy_age.data.cohorts import layout as cohort_layout
from entropy_age.experiments.bxd_liver import ols
from entropy_age.experiments.common import Context, dummies, write_table
from entropy_age.experiments.health_states import BRAAK, MEASURES, loo_scores
from entropy_age.experiments.per_person import _islands, _sets
from entropy_age.experiments.tissue_contrast import _features, _matrix

ORDER = {"cerebellum": 0, "frontal cortex": 1, "superior temporal gyrus": 2, "entorhinal cortex": 3}
COMPONENTS = (0, 5)


def gradients(scores: pl.DataFrame, keys: list[str]) -> pl.DataFrame:
    """Per donor and score: the slope of the region-standardised score on the region's Braak order."""
    z = scores.with_columns(
        [((pl.col(k) - pl.col(k).mean()) / pl.col(k).std()).over("region", "components").alias(k) for k in keys]
    ).with_columns(rank=pl.col("region").replace_strict(ORDER, return_dtype=pl.Float64))
    rows = []
    for (donor, comp), g in z.group_by("donor", "components", maintain_order=True):
        if g.height < 3:
            continue
        r = g["rank"].to_numpy()
        rc = r - r.mean()
        row = {"donor": donor, "components": comp, "regions": g.height}
        for k in keys:
            v = g[k].to_numpy()
            row[k] = float(rc @ (v - v.mean()) / (rc @ rc))
        rows.append(row)
    return pl.DataFrame(rows)


def stage_effects(grad: pl.DataFrame, donors: pl.DataFrame, keys: list[str]) -> list[dict]:
    """Per score: the gradient's mean (does disorder rise along the order at all?) and its slope on Braak stage."""
    t = grad.join(donors, on="donor")
    out = []
    for (comp,), g in t.group_by("components", maintain_order=True):
        design = np.column_stack(
            [np.ones(g.height), g["braak"].to_numpy(), g["age"].to_numpy(), dummies(g["sex"].to_numpy())]
        )
        for k in keys:
            y = g[k].to_numpy()
            coef, se, p = ols(y, design, 1)
            mean = float(y.mean())
            se_mean = float(y.std(ddof=1) / np.sqrt(y.size))
            measure, _, group = k.partition("_")
            out.append(
                {
                    "components": comp,
                    "measure": measure,
                    "group": group,
                    "donors": g.height,
                    "mean_gradient": mean,
                    "mean_gradient_se": se_mean,
                    "braak_slope": coef,
                    "braak_se": se,
                    "braak_p": p,
                }
            )
    return out


def run(ctx: Context) -> list[Path]:
    out = ctx.out_dir("braak_gradient")
    features, probes, membership = _features(ctx.paths, ctx.rng(140_000))
    s = (
        pl.read_parquet(cohort_layout("lunnon", ctx.paths).samples)
        .drop_nulls("age")
        .filter(pl.col("tissue").is_in(list(ORDER)) & pl.col("braak_stage").is_in(list(BRAAK)))
    )
    parts = []
    for (region,), t in s.group_by("tissue", maintain_order=True):
        t = t.sort("donor")
        m = _matrix("lunnon", ctx.paths, t["sample"].to_list(), features)
        island = _islands(m.features, probes)
        sets = _sets(m.features, membership)
        for k in COMPONENTS:
            sc = loo_scores(m.values, island, sets, components=k)
            parts.append(
                pl.DataFrame(
                    {
                        "donor": t["donor"],
                        "region": region,
                        "components": k,
                        **{key: v for key, v in sc.items() if key.split("_")[0] in MEASURES},
                    }
                )
            )
        print(f"braak_gradient: {region}: {t.height} donors", file=sys.stderr, flush=True)
    scores = pl.concat(parts, how="diagonal")
    keys = [c for c in scores.columns if c.split("_")[0] in MEASURES]
    donors = s.group_by("donor").agg(
        pl.col("braak_stage").first().replace_strict(BRAAK).cast(pl.Float64).alias("braak"),
        pl.col("age").mean(),
        pl.col("sex").first(),
    )
    grad = gradients(scores, keys)
    write_table(scores.join(donors, on="donor"), out / "region_scores")
    write_table(grad, out / "gradients")
    write_table(pl.DataFrame(stage_effects(grad, donors, keys)), out / "effects")
    return []
