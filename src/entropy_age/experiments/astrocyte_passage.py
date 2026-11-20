"""A second cell type in culture: fetal astrocytes passaged to exhaustion (roadmap RM40, GSE202554).

GSE202554 (found by the external dataset search) passaged primary human fetal astrocytes of two donors for up to 15
cumulative population doublings, on the EPIC array, with passages well spread over its 8 chips; it is rebuilt from
IDATs. In fibroblasts the housekeeping gap per doubling is about −0.10 SD (Sturm, R124). Here the same slopes are
taken in another cell type: unsorted (bulk) cultures only, per promoter CpG the slopes of relative entropy, entropy
and dispersion on doublings with a donor intercept (and chip intercepts under ``ENTROPY_AGE_CHIP=1``), then set
CpGs against CpGs matched on island status and starting level.

Limits: two donors, 37 of the 45 arrays from one; a span of 15 doublings against 79 in Sturm; no arrested control.
``entropy-age figures astrocyte_passage`` (needs the SeSAMe rebuild of GSE202554).
"""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import polars as pl

from entropy_age.data import reprocessed
from entropy_age.data.cohorts import _characteristics
from entropy_age.experiments.brain_replication import complete_matrix
from entropy_age.experiments.chromatin_match import N_RESAMPLES, matched_gaps
from entropy_age.experiments.common import Context, chip_terms, dummies, write_table
from entropy_age.experiments.per_person import _sets
from entropy_age.experiments.sturm_expression import chip_overlap
from entropy_age.experiments.tissue_contrast import Prepared, _features, _lifespan_metrics
from entropy_age.stochasticity import baseline_strata

GSE = "GSE202554"
METRICS = ("kl_slope", "entropy_slope", "dispersion_slope")


def bulk_cultures(chars: pl.DataFrame) -> pl.DataFrame:
    """``sample``, ``donor`` and ``doublings`` of the unsorted cultures with doublings recorded."""
    chars = chars.rename({c: c.lower() for c in chars.columns})
    return (
        chars.filter(pl.col("facs.code") == "Bulk")
        .select(
            "sample", pl.col("donor").cast(pl.String), pl.col("cpd").cast(pl.Float64, strict=False).alias("doublings")
        )
        .drop_nulls("doublings")
        .sort("donor", "doublings")
    )


def run(ctx: Context) -> list[Path]:
    paths = ctx.paths
    matrix = paths.input / "geo" / GSE / f"{GSE}_series_matrix.txt.gz"
    if not (reprocessed.available(paths, GSE) and matrix.exists()):
        print(f"astrocyte_passage: no SeSAMe rebuild or series matrix for {GSE}; nothing run", file=sys.stderr)
        return []
    out = ctx.out_dir("astrocyte_passage")
    features, probes, membership = _features(paths, ctx.rng(140_000))
    qc = reprocessed.read_qc(paths, GSE).filter(pl.col("passed"))
    cultures = bulk_cultures(_characteristics(matrix)).join(qc.select("sample", "chip"), on="sample")
    betas = reprocessed.by_accession(paths, GSE, {s: s for s in cultures["sample"].to_list()})
    x, feats = complete_matrix(betas, features, cultures["sample"].to_list())
    doublings = cultures["doublings"].to_numpy()
    donors, chips = cultures["donor"].to_numpy(), cultures["chip"].to_numpy()
    write_table(chip_overlap(doublings, donors, chips), out / "chip_overlap")
    cov = dummies(donors)
    if chip_terms():
        cov = np.column_stack([cov, dummies(chips)])
    metrics = _lifespan_metrics(Prepared("astrocytes", x, doublings, cov, feats, donors, cov))
    island_of = dict(probes.select("cpg", "cpg_island").iter_rows())
    island = np.array([bool(island_of.get(f, False)) for f in feats])
    levels = baseline_strata(metrics["baseline"].to_numpy(), 10).astype(str)
    variants = {"island x level": np.char.add(np.char.add(island.astype(str), "|"), levels)}
    sets = _sets(feats, membership.filter(pl.col("feature").is_in(feats)))
    gaps = []
    for metric in METRICS:
        values = metrics[metric].to_numpy()
        rows = matched_gaps(values, sets, variants, N_RESAMPLES[str(ctx.scale.name)], ctx.seed)
        gaps += [{"metric": metric, **r} for r in rows]
    write_table(pl.DataFrame(gaps), out / "matched_gaps")
    sizes = cultures.group_by("donor", maintain_order=True).agg(
        pl.len().alias("cultures"), pl.col("doublings").min().alias("min"), pl.col("doublings").max().alias("max")
    )
    write_table(sizes.with_columns(cpgs=pl.lit(len(feats))), out / "sizes")
    print(f"astrocyte_passage: {cultures.height} cultures, {len(feats)} CpGs", file=sys.stderr, flush=True)
    return []
