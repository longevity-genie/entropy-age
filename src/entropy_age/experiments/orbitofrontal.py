"""A third cortex series: the age gap in bulk orbitofrontal cortex (roadmap RM40, GSE254293).

GSE254293 (found by the external dataset search) holds bulk orbitofrontal cortex of 91 people aged 24 to 84 on the
EPIC array, 60 with a psychiatric diagnosis and 36 without, age well spread over its 12 chips. It is rebuilt from
IDATs and analysed as the other cortex series (``brain_replication.group_gaps``): the matched gap of the age slopes
in each diagnosis group and in everyone with diagnosis held, with sex or sex and chip as covariates. Technical
duplicates are dropped. The record codes diagnosis as 0 and 1 without naming them; 0 is read as no diagnosis, the
smaller group, as in the paper's design, and both groups are reported so nothing rests on that reading.

Limits: bulk tissue, cell mix not held; post-mortem interval up to 72 h and smoking are not held.
``entropy-age figures orbitofrontal`` (needs the SeSAMe rebuild of GSE254293).
"""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import polars as pl

from entropy_age.data import reprocessed
from entropy_age.data.annotation import AnnotationLayout
from entropy_age.data.cohorts import _characteristics
from entropy_age.experiments.brain_replication import complete_matrix, group_gaps
from entropy_age.experiments.chromatin_match import N_RESAMPLES
from entropy_age.experiments.common import Context, dummies, write_table
from entropy_age.experiments.tissue_contrast import _features

GSE = "GSE254293"
MIN_SAMPLES = 20


def read_people(matrix: Path) -> pl.DataFrame:
    """``sample``, ``age``, ``sex``, ``status`` ("0" or "1"); technical duplicates left out."""
    chars = _characteristics(matrix)
    chars = chars.rename({c: c.lower() for c in chars.columns})
    return (
        chars.filter(pl.col("technical_duplicate") == "no")
        .select(
            "sample",
            pl.col("age").cast(pl.Float64, strict=False),
            pl.col("sex").fill_null("NA"),
            pl.col("disease_status").alias("status"),
        )
        .drop_nulls("age")
    )


def groups(people: pl.DataFrame) -> dict[str, pl.DataFrame]:
    return {
        "status 0": people.filter(pl.col("status") == "0"),
        "status 1": people.filter(pl.col("status") == "1"),
        "everyone, status held": people,
    }


def run(ctx: Context) -> list[Path]:
    matrix = ctx.paths.input / "geo" / GSE / f"{GSE}_series_matrix.txt.gz"
    if not (reprocessed.available(ctx.paths, GSE) and matrix.exists()):
        print(f"orbitofrontal: no SeSAMe rebuild or series matrix for {GSE}; nothing run", file=sys.stderr)
        return []
    out = ctx.out_dir("orbitofrontal")
    n_resamples = N_RESAMPLES[str(ctx.scale.name)]
    features, _, membership = _features(ctx.paths, ctx.rng(120_000))
    probes = pl.read_parquet(AnnotationLayout(ctx.paths).probes_450k).filter(pl.col("cpg").is_in(features))
    island_of = dict(probes.select("cpg", "cpg_island").iter_rows())
    qc = reprocessed.read_qc(ctx.paths, GSE).filter(pl.col("passed"))
    chip_of = dict(qc.select("sample", "chip").iter_rows())
    betas = reprocessed.by_accession(ctx.paths, GSE, {s: s for s in qc["sample"].to_list()})
    people = read_people(matrix).filter(pl.col("sample").is_in(betas.columns)).sort("age", "sample")
    rows, sizes = [], []
    for label, s in groups(people).items():
        if s.height < MIN_SAMPLES:
            print(f"orbitofrontal: {label}: {s.height} arrays, skipped", file=sys.stderr)
            continue
        x, feats = complete_matrix(betas, features, s["sample"].to_list())
        ages = s["age"].to_numpy()
        chips = np.array([chip_of[v] for v in s["sample"].to_list()])
        base = dummies(s["sex"].to_numpy())
        if s["status"].n_unique() > 1:
            base = np.column_stack([base, dummies(s["status"].to_numpy())])
        for design, cov in {"sex": base, "sex + chip": np.column_stack([base, dummies(chips)])}.items():
            cov = cov[:, np.ptp(cov, axis=0) > 1e-12]
            if s.height - cov.shape[1] - 2 < 10:
                continue
            rows += [
                {"cohort": label, "covariates": design, **r}
                for r in group_gaps(x, ages, cov, feats, island_of, membership, n_resamples, ctx.seed)
            ]
        within = sum(((ages[chips == c] - ages[chips == c].mean()) ** 2).sum() for c in np.unique(chips))
        sizes.append(
            {
                "group": label,
                "people": s.height,
                "age_min": float(ages.min()),
                "age_max": float(ages.max()),
                "cpgs": len(feats),
                "chips": int(np.unique(chips).size),
                "r2_age_by_chip": float(1 - within / ((ages - ages.mean()) ** 2).sum()),
            }
        )
        print(f"orbitofrontal: {label}: {s.height} people, {len(feats)} CpGs", file=sys.stderr, flush=True)
    if rows:
        write_table(pl.DataFrame(sizes), out / "sizes")
        write_table(pl.DataFrame(rows), out / "gaps")
    return []
