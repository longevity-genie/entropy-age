"""The joint test in cortex: is anything left of the shield beyond promoter activity? (roadmap RM44)

Three checks each left a gap in Jaffe's prefrontal cortex when run alone: CpGs kept off the 0 % floor (R80),
gene sets matched on expression (R82), and matching on chromatin state and sorted accessibility (R68b). Here they are
held together. Within a band of young methylation level (all CpGs, 10-90 %, 20-80 %), a set's promoter CpGs are
compared with other promoter CpGs in the same stratum of

* CpG-island status x young-level decile (the base matching),
* plus the expression decile of the CpG's gene (GTEx v8 frontal cortex BA9 median TPM; a CpG of several genes takes
  the highest; CpGs of genes GTEx does not list form their own class),
* plus the chromatin state (Roadmap E073, prefrontal cortex),
* plus sorted neuron and glia accessibility (GSE96614, BA9),
* plus all three.

If the gap is gone under the last matching, what the data support is "active promoters drift less"; if it stays,
essential and housekeeping promoters are protected beyond activity as these references measure it. The expression,
state and accessibility tracks come from other people than the methylation (Supplementary Table S4).

``entropy-age figures joint_cortex`` (with ``ENTROPY_AGE_BETAS=sesame`` for the rebuilt arrays), as a Slurm job.
"""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import polars as pl

from entropy_age.data.annotation import AnnotationLayout
from entropy_age.data.cohorts import layout as cohort_layout
from entropy_age.experiments.chromatin_match import N_RESAMPLES, _states, matched_gaps
from entropy_age.experiments.common import Context, write_table
from entropy_age.experiments.per_person import _sets
from entropy_age.experiments.placebo_expression import GTEX, TISSUE_TPM, expression, with_deciles
from entropy_age.experiments.tissue_contrast import _features, _inputs_missing, _jaffe, _lifespan_metrics
from entropy_age.stochasticity import baseline_strata

LABEL = "cortex (Jaffe)"
#: band name -> (lowest, highest) young methylation level kept
BANDS = {"all": (0.0, 1.0), "10-90 %": (0.1, 0.9), "20-80 %": (0.2, 0.8)}
METRICS = ("kl_slope", "dispersion_slope", "entropy_slope")
NO_EXPRESSION = "na"


def expression_class(features: list[str], links: pl.DataFrame, genes: pl.DataFrame) -> np.ndarray:
    """Per CpG: the expression decile (``"0"``-``"9"``) of its most expressed promoter gene, ``"na"`` when none of its
    genes is in ``genes`` (``symbol``, ``decile``)."""
    best = (
        links.filter(pl.col("promoter"))
        .join(genes.select("symbol", "decile"), on="symbol")
        .group_by("cpg")
        .agg(pl.col("decile").max())
    )
    joined = pl.DataFrame({"cpg": features}).join(best, on="cpg", how="left")
    return joined["decile"].cast(pl.String).fill_null(NO_EXPRESSION).to_numpy().astype(str)


def joint_variants(
    island: np.ndarray, baseline: np.ndarray, expr: np.ndarray, state: np.ndarray, cells: np.ndarray
) -> dict[str, np.ndarray]:
    """The stratum label of every CpG under each matching, from the base one to all tracks together."""

    def add(*parts: np.ndarray) -> np.ndarray:
        out = parts[0]
        for part in parts[1:]:
            out = np.char.add(np.char.add(out, "|"), part)
        return out

    base = add(island.astype(str), baseline_strata(baseline, 10).astype(str))
    return {
        "island x level": base,
        "+ expression": add(base, expr),
        "+ chromatin state": add(base, state),
        "+ sorted ATAC": add(base, cells),
        "+ expression + state + sorted ATAC": add(base, expr, state, cells),
    }


def usable_share(in_set: np.ndarray, strata: np.ndarray) -> float:
    """Share of the set's CpGs sitting in a stratum that also holds control CpGs (the others drop out of the gap)."""
    codes = np.unique(strata, return_inverse=True)[1]
    has_set = np.bincount(codes, weights=in_set.astype(float)) > 0
    has_other = np.bincount(codes, weights=(~in_set).astype(float)) > 0
    return float(in_set[(has_set & has_other)[codes]].sum() / max(1, in_set.sum()))


def run(ctx: Context) -> list[Path]:
    gtex = ctx.paths.input / GTEX
    if _inputs_missing(ctx, "joint_cortex") or not gtex.exists():
        print("joint_cortex: inputs missing (cohorts, annotation or GTEx medians); nothing run", file=sys.stderr)
        return []
    out = ctx.out_dir("joint_cortex")
    n_resamples = N_RESAMPLES[str(ctx.scale.name)]
    features, _, membership = _features(ctx.paths, ctx.rng(120_000))
    annotation = AnnotationLayout(ctx.paths)
    probes = pl.read_parquet(annotation.probes_450k).filter(pl.col("cpg").is_in(features))
    island_of = dict(probes.select("cpg", "cpg_island").iter_rows())
    p = _jaffe(ctx.paths, features)
    metrics = _lifespan_metrics(p)
    feats = metrics["feature"].to_list()
    baseline = metrics["baseline"].to_numpy()
    states = _states(feats, probes, ctx.paths)
    genes = with_deciles(expression(gtex, TISSUE_TPM[LABEL]))
    expr = expression_class(feats, pl.read_parquet(annotation.links_450k), genes)
    island = np.array([bool(island_of.get(f, False)) for f in feats])
    cells = np.char.add(np.char.add(states["NeuN+"], "|"), states["NeuN-"])
    rows, sizes, coverage = [], [], []
    for band, (lo, hi) in BANDS.items():
        keep = (baseline >= lo) & (baseline <= hi)
        kept = [f for f, k in zip(feats, keep, strict=True) if k]
        sets = _sets(kept, membership)
        variants = joint_variants(island[keep], baseline[keep], expr[keep], states["E073"][keep], cells[keep])
        sizes.append(
            {
                "band": band,
                "cpgs": int(keep.sum()),
                **{f"n_{g}": int(s[0].sum()) for g, s in sets.items()},
                "strata_all_tracks": int(np.unique(variants["+ expression + state + sorted ATAC"]).size),
            }
        )
        coverage += [
            {"band": band, "group": g, "matching": v, "set_cpgs": int(s[0].sum()), "usable": usable_share(s[0] > 0, st)}
            for g, s in sets.items()
            for v, st in variants.items()
        ]
        for metric in METRICS:
            values = metrics[metric].to_numpy()[keep]
            rows += [
                {"band": band, "metric": metric, **r}
                for r in matched_gaps(values, sets, variants, n_resamples, seed=ctx.seed)
            ]
        print(f"joint_cortex: band {band}: {int(keep.sum())} CpGs", file=sys.stderr, flush=True)
    processing = cohort_layout("jaffe", ctx.paths).processing
    write_table(pl.DataFrame(rows).with_columns(processing=pl.lit(processing)), out / "gaps")
    write_table(pl.DataFrame(coverage), out / "coverage")
    write_table(
        pl.DataFrame(sizes).with_columns(donors=pl.lit(p.x.shape[0]), processing=pl.lit(processing)), out / "sizes"
    )
    return []
