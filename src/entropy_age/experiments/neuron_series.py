"""A third series of sorted cortical neurons: is there a shield in cells that do not divide? (roadmap RM40)

GSE112179 (Pai et al. 2019; found by the external dataset search) sorted NeuN+ nuclei from the frontal cortex of 83
donors (schizophrenia, bipolar disorder and 28 controls; 23 to 77 years; 100 arrays with technical replicates, of
which one per donor is used) and ran them on the EPIC array; the series
is rebuilt from IDATs with the project's SeSAMe recipe. Neurons do not divide, so a shield here cannot come from
selection among dividing cells. The first sorted series (GSE279509) had a neuronal gap near zero within chips.

Two groups: every donor with the diagnosis held as a covariate, and the controls alone. Per group, as in
``brain_replication``: the age slope of each promoter CpG's drift measures (sex, and diagnosis where it varies, held;
a second fit adds chip intercepts), then the matched gap of each gene set against CpGs of the same island status and
young-level decile, on all CpGs and on those off the floor. Only CpGs measured in every array of the group are used.
The overlap of age with chip is tabulated per group.

Limits: one cell fraction and no glia of the same donors, so the series gives the size of the neuronal gap, not a
contrast; two thirds of the donors carry a psychiatric diagnosis; the controls alone are few.
``entropy-age figures neuron_series`` (needs the SeSAMe rebuild of GSE112179).
"""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import polars as pl

from entropy_age.data import reprocessed
from entropy_age.data.annotation import AnnotationLayout
from entropy_age.data.cohorts import _characteristics
from entropy_age.experiments.brain_replication import MIN_SAMPLES, complete_matrix, group_gaps
from entropy_age.experiments.chromatin_match import N_RESAMPLES
from entropy_age.experiments.common import Context, dummies, write_table
from entropy_age.experiments.sturm_expression import chip_overlap
from entropy_age.experiments.tissue_contrast import _features

GSE = "GSE112179"
ALL, CONTROLS = "neurons (NeuN+), every donor (diagnosis held)", "neurons (NeuN+), controls"


def neuron_groups(chars: pl.DataFrame, usable: list[str] | None = None) -> dict[str, pl.DataFrame]:
    """Group -> ``sample``, ``donor``, ``age``, ``sex``, ``diagnosis``: one array per donor (the series holds
    technical replicates of some donors; the first array by accession among
    ``usable`` is kept)."""
    frame = (
        chars.select(
            "sample",
            donor=pl.col("tissuebank.id"),
            age=pl.col("age").cast(pl.Float64, strict=False),
            sex=pl.col("Sex"),
            diagnosis=pl.col("dist.dx"),
        )
        .drop_nulls("age")
        .filter(pl.col("sample").is_in(usable) if usable is not None else pl.lit(True))
        .sort("sample")
        .unique(subset="donor", keep="first", maintain_order=True)
    )
    return {ALL: frame, CONTROLS: frame.filter(pl.col("diagnosis") == "Control")}


def run(ctx: Context) -> list[Path]:
    matrix = ctx.paths.input / "geo" / GSE / f"{GSE}_series_matrix.txt.gz"
    if not (reprocessed.available(ctx.paths, GSE) and matrix.exists()):
        print(f"neuron_series: no SeSAMe rebuild or series matrix for {GSE}; nothing run", file=sys.stderr)
        return []
    out = ctx.out_dir("neuron_series")
    n_resamples = N_RESAMPLES[str(ctx.scale.name)]
    features, _, membership = _features(ctx.paths, ctx.rng(120_000))
    probes = pl.read_parquet(AnnotationLayout(ctx.paths).probes_450k).filter(pl.col("cpg").is_in(features))
    island_of = dict(probes.select("cpg", "cpg_island").iter_rows())
    qc = reprocessed.read_qc(ctx.paths, GSE).filter(pl.col("passed"))
    chip_of = dict(qc.select("sample", "chip").iter_rows())
    betas = reprocessed.by_accession(ctx.paths, GSE, {s: s for s in qc["sample"].to_list()})
    rows, sizes, overlap = [], [], []
    for label, s in neuron_groups(_characteristics(matrix), betas.columns).items():
        s = s.filter(pl.col("sample").is_in(betas.columns)).sort("age", "sample")
        if s.height < MIN_SAMPLES:
            print(f"neuron_series: {label}: {s.height} arrays, skipped", file=sys.stderr)
            continue
        x, feats = complete_matrix(betas, features, s["sample"].to_list())
        ages = s["age"].to_numpy()
        chips = np.array([chip_of.get(v) or "NA" for v in s["sample"].to_list()])
        diagnosis = s["diagnosis"].fill_null("NA").to_numpy()
        held = np.column_stack([dummies(s["sex"].fill_null("NA").to_numpy()), dummies(diagnosis)])
        overlap += [{"group": label, **r} for r in chip_overlap(ages, diagnosis, chips).to_dicts()]
        sizes.append(
            {
                "group": label,
                "arrays": s.height,
                "age_min": float(ages.min()),
                "age_max": float(ages.max()),
                "cpgs": len(feats),
                "chips": int(np.unique(chips).size),
                "chips_with_2_or_more": sum(1 for c in np.unique(chips) if (chips == c).sum() >= 2),
            }
        )
        for design, cov in {"sex": held, "sex + chip": np.column_stack([held, dummies(chips)])}.items():
            cov = cov[:, np.ptp(cov, axis=0) > 1e-12] if cov.shape[1] else cov
            if s.height - cov.shape[1] - 2 < 10:  # too few residual degrees of freedom for a slope
                continue
            rows += [
                {"cohort": label, "covariates": design, **r}
                for r in group_gaps(x, ages, cov, feats, island_of, membership, n_resamples, ctx.seed)
            ]
        print(f"neuron_series: {label}: {s.height} arrays, {len(feats)} CpGs", file=sys.stderr, flush=True)
    for name, table in (("sizes", sizes), ("gaps", rows), ("chip_overlap", overlap)):
        if table:
            write_table(pl.DataFrame(table), out / name)
    return []
