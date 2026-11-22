"""Three brain regions of the same donors, a fourth human brain series (roadmap RM40, GSE349373).

GSE349373 (found by the external dataset search) holds prefrontal cortex, amygdala and hippocampus of about 300
donors aged 15 to 80 on the EPIC array: controls and donors with major depression or post-traumatic stress disorder,
about 100 of each per region. It is rebuilt from IDATs. The brain bank is the one behind Jaffe's cohort, so donors
who are among Jaffe's adult controls are left out; what is left is independent of the discovery cohort.

Per region, two groups: the controls, and every donor with the diagnosis held. As in ``brain_replication``: the age
slope of each promoter CpG's drift measures (sex, and diagnosis where it varies, held; a second fit adds chip
intercepts), then the matched gap of each gene set against CpGs of the same island status and young-level decile, on
all CpGs and on those off the floor. Only CpGs measured in every array of the group are used. The overlap of age
with chip is tabulated per group.

Limits: bulk tissue, so the mix of neurons and glia is not held; regions of one donor are analysed apart, not as
pairs. ``entropy-age figures brain_regions`` (needs the SeSAMe rebuild of GSE349373).
"""

from __future__ import annotations

import gzip
import sys
from pathlib import Path

import numpy as np
import polars as pl

from entropy_age.data import reprocessed
from entropy_age.data.annotation import AnnotationLayout
from entropy_age.data.cohorts import _characteristics
from entropy_age.data.cohorts import layout as cohort_layout
from entropy_age.experiments.brain_replication import MIN_SAMPLES, complete_matrix, group_gaps
from entropy_age.experiments.chromatin_match import N_RESAMPLES
from entropy_age.experiments.common import Context, dummies, write_table
from entropy_age.experiments.sturm_expression import chip_overlap
from entropy_age.experiments.tissue_contrast import _features

GSE = "GSE349373"
REGIONS = {"PFC": "prefrontal cortex", "AMG": "amygdala", "HPC": "hippocampus"}


def region_groups(chars: pl.DataFrame, titles: list[str], exclude_donors: set[str]) -> dict[str, pl.DataFrame]:
    """Group -> ``sample``, ``donor`` (from the title ``<donor>_<region>``), ``age``, ``sex``, ``diagnosis``; one
    array per donor and region, the donors in ``exclude_donors`` left out."""
    frame = (
        chars.select(
            "sample",
            donor=pl.Series([t.split("_")[0] for t in titles]),
            region=pl.col("brain region"),
            age=pl.col("age at death").cast(pl.Float64, strict=False),
            sex=pl.col("Sex"),
            diagnosis=pl.col("primary dx"),
        )
        .drop_nulls("age")
        .filter(~pl.col("donor").is_in(list(exclude_donors)))
        .sort("sample")
        .unique(subset=["donor", "region"], keep="first", maintain_order=True)
    )
    out = {}
    for code, name in REGIONS.items():
        part = frame.filter(pl.col("region") == code).drop("region")
        out[f"{name}, controls"] = part.filter(pl.col("diagnosis") == "Control")
        out[f"{name}, every donor (diagnosis held)"] = part
    return out


def _titles(matrix: Path) -> list[str]:
    with gzip.open(matrix, "rt") as handle:
        for line in handle:
            if line.startswith("!Sample_title"):
                return [c.strip().strip('"') for c in line.rstrip("\n").split("\t")[1:]]
    return []


def run(ctx: Context) -> list[Path]:
    matrix = ctx.paths.input / "geo" / GSE / f"{GSE}_series_matrix.txt.gz"
    if not (reprocessed.available(ctx.paths, GSE) and matrix.exists()):
        print(f"brain_regions: no SeSAMe rebuild or series matrix for {GSE}; nothing run", file=sys.stderr)
        return []
    out = ctx.out_dir("brain_regions")
    n_resamples = N_RESAMPLES[str(ctx.scale.name)]
    features, _, membership = _features(ctx.paths, ctx.rng(120_000))
    probes = pl.read_parquet(AnnotationLayout(ctx.paths).probes_450k).filter(pl.col("cpg").is_in(features))
    island_of = dict(probes.select("cpg", "cpg_island").iter_rows())
    jaffe = pl.read_parquet(cohort_layout("jaffe", ctx.paths).samples)
    jaffe_donors = set(jaffe.filter((pl.col("group") == "Control") & (pl.col("age") >= 20))["donor"].to_list())
    qc = reprocessed.read_qc(ctx.paths, GSE).filter(pl.col("passed"))
    chip_of = dict(qc.select("sample", "chip").iter_rows())
    passed = set(qc["sample"].to_list())
    chars = _characteristics(matrix)
    groups = region_groups(
        chars.filter(pl.col("sample").is_in(list(passed))),
        [t for t, s in zip(_titles(matrix), chars["sample"].to_list(), strict=True) if s in passed],
        jaffe_donors,
    )
    rows, sizes, overlap = [], [], []
    for label, s in groups.items():
        if s.height < MIN_SAMPLES:
            print(f"brain_regions: {label}: {s.height} arrays, skipped", file=sys.stderr)
            continue
        s = s.sort("age", "sample")
        betas = reprocessed.by_accession(ctx.paths, GSE, {v: v for v in s["sample"].to_list()})
        s = s.filter(pl.col("sample").is_in(betas.columns))
        x, feats = complete_matrix(betas, features, s["sample"].to_list())
        ages = s["age"].to_numpy()
        chips = np.array([str(chip_of.get(v) or "NA") for v in s["sample"].to_list()])
        diagnosis = s["diagnosis"].fill_null("NA").to_numpy()
        held = np.column_stack([dummies(s["sex"].fill_null("NA").to_numpy()), dummies(diagnosis)])
        overlap += [{"group": label, **r} for r in chip_overlap(ages, diagnosis, chips).to_dicts()]
        sizes.append(
            {
                "group": label,
                "arrays": s.height,
                "donors": s["donor"].n_unique(),
                "age_min": float(ages.min()),
                "age_max": float(ages.max()),
                "cpgs": len(feats),
                "chips": int(np.unique(chips).size),
            }
        )
        for design, cov in {"sex": held, "sex + chip": np.column_stack([held, dummies(chips)])}.items():
            cov = cov[:, np.ptp(cov, axis=0) > 1e-12] if cov.shape[1] else cov
            if s.height - cov.shape[1] - 2 < 10:
                continue
            rows += [
                {"cohort": label, "covariates": design, **r}
                for r in group_gaps(x, ages, cov, feats, island_of, membership, n_resamples, ctx.seed)
            ]
        print(f"brain_regions: {label}: {s.height} arrays, {len(feats)} CpGs", file=sys.stderr, flush=True)
    for name, table in (("sizes", sizes), ("gaps", rows), ("chip_overlap", overlap)):
        if table:
            write_table(pl.DataFrame(table), out / name)
    return []
