"""Replication of the cortex shield in two new brain series, one of sorted nuclei (roadmap RM40, RM44).

Both series were chosen for raw IDATs, one array type (EPIC) and a design in which age is not bound to chip more than
chance allows (``docs/DATASET_SEARCH_2026-11-18.md``), and are rebuilt with the project's SeSAMe recipe.

* **GSE252501** (White et al. 2024): dorsolateral prefrontal cortex and nucleus accumbens of the same donors. Used:
  donors without alcohol use disorder, each region apart. The brain bank is the one behind Jaffe's cohort, so donors
  also among Jaffe's adult controls are left out; what is left is independent of the discovery cohort.
* **GSE279509** (Walker et al. 2025): nuclei sorted from cortex into neurons (NeuN+), oligodendrocyte lineage
  (SOX10+), other glia (double negative) and unsorted nuclei, 18-108 years. Neurons do not divide, so a shield in
  NeuN+ nuclei cannot come from selection among dividing cells. GEO gives no diagnosis for this series.

Per group: the age slope of each promoter CpG's drift measures (sex held; a second fit adds chip intercepts, which
estimates age within chips only), then the matched gap of each gene set against CpGs of the same island status and
young-level decile, on all CpGs and on those off the floor (10-90 %). Only CpGs measured in every array of the group
are used: nothing is filled. Intervals resample the set's genes.

``ENTROPY_AGE_BETAS`` is not consulted: these series exist only as rebuilds. ``entropy-age figures brain_replication``.
"""

from __future__ import annotations

import gzip
import re
import sys
from pathlib import Path

import numpy as np
import polars as pl

from entropy_age.data import reprocessed
from entropy_age.data.annotation import AnnotationLayout
from entropy_age.data.cohorts import _characteristics
from entropy_age.data.cohorts import layout as cohort_layout
from entropy_age.experiments.chromatin_match import N_RESAMPLES, matched_gaps
from entropy_age.experiments.common import Context, dummies, write_table
from entropy_age.experiments.per_person import _sets
from entropy_age.experiments.tissue_contrast import _features
from entropy_age.stochasticity import baseline_strata, drift_metrics

WHITE, WALKER = "GSE252501", "GSE279509"
CELL_TYPES = {"NeuNPos": "neurons (NeuN+)", "Sox10Pos": "oligodendrocyte lineage (SOX10+)",
              "DoubleNeg": "other glia (NeuN-/SOX10-)", "Total": "unsorted nuclei"}  # fmt: skip
BANDS = {"all": (0.0, 1.0), "10-90 %": (0.1, 0.9)}
METRICS = ("kl_slope", "dispersion_slope", "entropy_slope")
MIN_SAMPLES = 30


def _titles(matrix: Path) -> list[str]:
    with gzip.open(matrix, "rt") as handle:
        for line in handle:
            if line.startswith("!Sample_title"):
                return [c.strip().strip('"') for c in line.rstrip("\n").split("\t")[1:]]
    return []


def white_groups(matrix: Path, exclude_donors: set[str]) -> dict[str, pl.DataFrame]:
    """Region -> ``sample``, ``donor``, ``age``, ``sex`` of donors without alcohol use disorder, Jaffe's donors out."""
    chars = _characteristics(matrix)
    titles = _titles(matrix)
    frame = chars.with_columns(
        donor=pl.Series([m.group(0) if (m := re.search(r"Br\d+", t)) else t for t in titles]),
        region=pl.Series([t.split(",")[1].strip() if "," in t else "" for t in titles]),
        age=pl.col("agedeath").cast(pl.Float64, strict=False),
        sex=pl.col("Sex"),
    ).filter((pl.col("dsm-5 aud status") == "Control") & ~pl.col("donor").is_in(list(exclude_donors)))
    return {
        f"{WHITE} {region}, controls": part.select("sample", "donor", "age", "sex")
        for (region,), part in frame.group_by("region", maintain_order=True)
    }


def walker_groups(matrix: Path) -> dict[str, pl.DataFrame]:
    """Cell type -> ``sample``, ``donor``, ``age``, ``sex``, one array per individual and cell type."""
    chars = _characteristics(matrix).with_columns(
        donor=pl.col("individual id").cast(pl.String),
        age=pl.col("age").cast(pl.Float64, strict=False),
        sex=pl.col("Sex"),
    )
    out = {}
    for code, label in CELL_TYPES.items():
        part = chars.filter(pl.col("cell type") == code).unique(subset="donor", keep="first", maintain_order=True)
        out[f"{WALKER} {label}"] = part.select("sample", "donor", "age", "sex")
    return out


def complete_matrix(frame: pl.DataFrame, features: list[str], samples: list[str]) -> tuple[np.ndarray, list[str]]:
    """Samples x CpGs for the ``features`` measured in every one of ``samples`` (no value is filled)."""
    x = reprocessed.matrix(frame, features, samples)
    keep = ~np.isnan(x).any(axis=0)
    return x[:, keep], [f for f, k in zip(features, keep, strict=True) if k]


def group_gaps(x, ages, covariates, feats, island_of, membership, n_resamples: int, seed: int) -> list[dict]:
    """Matched gaps of every set, metric and band for one group and one covariate set."""
    metrics = drift_metrics(x, ages, feats, covariates=covariates if covariates.shape[1] else None, spearman=False)
    baseline = metrics["baseline"].to_numpy()
    island = np.array([bool(island_of.get(f, False)) for f in feats])
    rows = []
    for band, (lo, hi) in BANDS.items():
        keep = (baseline >= lo) & (baseline <= hi)
        kept = [f for f, k in zip(feats, keep, strict=True) if k]
        sets = _sets(kept, membership)
        strata = np.char.add(
            np.char.add(island[keep].astype(str), "|"), baseline_strata(baseline[keep], 10).astype(str)
        )
        for metric in METRICS:
            values = metrics[metric].to_numpy()[keep]
            rows += [
                {"band": band, "metric": metric, "cpgs": int(keep.sum()), **r}
                for r in matched_gaps(values, sets, {"island x level": strata}, n_resamples, seed=seed)
            ]
    return rows


def run(ctx: Context) -> list[Path]:
    missing = [g for g in (WHITE, WALKER) if not reprocessed.available(ctx.paths, g)]
    if missing:
        print(f"brain_replication: no SeSAMe rebuild for {missing}; nothing run", file=sys.stderr)
        return []
    out = ctx.out_dir("brain_replication")
    n_resamples = N_RESAMPLES[str(ctx.scale.name)]
    features, _, membership = _features(ctx.paths, ctx.rng(120_000))
    probes = pl.read_parquet(AnnotationLayout(ctx.paths).probes_450k).filter(pl.col("cpg").is_in(features))
    island_of = dict(probes.select("cpg", "cpg_island").iter_rows())
    jaffe = pl.read_parquet(cohort_layout("jaffe", ctx.paths).samples)
    jaffe_donors = set(jaffe.filter((pl.col("group") == "Control") & (pl.col("age") >= 20))["donor"].to_list())
    geo = ctx.paths.input / "geo"
    groups = {
        WHITE: white_groups(geo / WHITE / f"{WHITE}_series_matrix.txt.gz", jaffe_donors),
        WALKER: walker_groups(geo / WALKER / f"{WALKER}_series_matrix.txt.gz"),
    }
    rows, sizes = [], []
    for gse, parts in groups.items():
        qc = reprocessed.read_qc(ctx.paths, gse)
        chip_of = dict(qc.select("sample", "chip").iter_rows())
        betas = reprocessed.by_accession(ctx.paths, gse, {s: s for s in qc["sample"].to_list()})
        for label, s in parts.items():
            s = s.filter(pl.col("sample").is_in(betas.columns)).drop_nulls("age").sort("age", "sample")
            if s.height < MIN_SAMPLES:
                print(f"brain_replication: {label}: {s.height} arrays, skipped", file=sys.stderr)
                continue
            x, feats = complete_matrix(betas, features, s["sample"].to_list())
            ages = s["age"].to_numpy()
            chips = np.array([chip_of.get(v) or "NA" for v in s["sample"].to_list()])
            sex = dummies(s["sex"].fill_null("NA").to_numpy())
            designs = {"sex": sex, "sex + chip": np.column_stack([sex, dummies(chips)])}
            shared = sum(1 for c in np.unique(chips) if (chips == c).sum() >= 2)
            sizes.append(
                {
                    "group": label,
                    "arrays": s.height,
                    "donors": s["donor"].n_unique(),
                    "age_min": float(ages.min()),
                    "age_max": float(ages.max()),
                    "cpgs": len(feats),
                    "chips": int(np.unique(chips).size),
                    "chips_with_2_or_more": shared,
                }
            )
            for design, cov in designs.items():
                cov = cov[:, np.ptp(cov, axis=0) > 1e-12] if cov.shape[1] else cov
                if s.height - cov.shape[1] - 2 < 10:  # too few residual degrees of freedom for a slope
                    continue
                rows += [
                    {"cohort": label, "covariates": design, **r}
                    for r in group_gaps(x, ages, cov, feats, island_of, membership, n_resamples, ctx.seed)
                ]
            print(f"brain_replication: {label}: {s.height} arrays, {len(feats)} CpGs", file=sys.stderr, flush=True)
    write_table(pl.DataFrame(sizes), out / "sizes")
    write_table(pl.DataFrame(rows), out / "gaps")
    return []
