"""One runner for every array series analysed as an age slope, driven by a description file per dataset.

Each file ``config/series/<name>.toml`` names a rebuilt series and how to read its sample table; this module holds
every analytic decision, so that no series is processed its own way:

* arrays that pass the 90 % detection cut; **one array per donor** and group (the first by accession);
* only CpGs measured in every array of the group (nothing is filled);
* the age slope of each promoter CpG's drift measures with sex and the columns under ``hold`` as covariates, and a
  second fit that adds chip intercepts;
* the matched gap of each gene set against CpGs of the same island status and young-level decile
  (``brain_replication.group_gaps``), on all CpGs and on those starting between 10 and 90 %;
* **two intervals for every gap**: resampling the set's genes with the arrays fixed (stability over genes), and
  resampling donors with the strata fixed (sampling of people);
* the overlap of age with chip per group.

Description file::

    accession = "GSE349373"
    matrix = "GSE349373_series_matrix.txt.gz"   # optional; for series split by platform
    age = "age at death"          # column of the sample characteristics (or "source_name")
    sex = "Sex"
    donor = "title:^([^_]+)_"      # a column, or "title:<regex with one group>" applied to the sample title
    hold = ["primary dx"]         # further covariates, held where they vary
    annotation = "human"          # "mm285" (mouse MM285 ids), or "mouse" (GPL32543 table; sparse island flags)
    strip_probe_suffix = false    # true for the 40k mammalian array, which writes cg... without _BC21
    age_values = {Young = 0, Old = 1}   # optional: the age column holds group labels, not numbers
    exclude_cpgs = "config/series/lists/x.csv"   # optional: CpGs left out (a ``cpg`` column)
    [derive]                      # optional: columns read from the title, "title:<regex with one group>"
    condition = "title:\\((Myo[NFT])\\)"
    exclude_jaffe_donors = true   # donors of the discovery cortex cohort (same brain bank) left out
    [groups."prefrontal cortex, controls"]
    "brain region" = "PFC"
    "primary dx" = "Control"

``entropy-age figures series_<name>``; outputs in ``data/output/series_<name>``.
"""

from __future__ import annotations

import gzip
import re
import sys
import tomllib
from pathlib import Path

import numpy as np
import polars as pl

from entropy_age.data import reprocessed
from entropy_age.data.annotation import AnnotationLayout
from entropy_age.data.cohorts import _characteristics
from entropy_age.data.cohorts import layout as cohort_layout
from entropy_age.experiments.brain_replication import BANDS, METRICS, MIN_SAMPLES, complete_matrix, group_gaps
from entropy_age.experiments.chromatin_match import N_RESAMPLES
from entropy_age.experiments.common import Context, dummies, write_table
from entropy_age.experiments.per_person import _sets
from entropy_age.experiments.sturm_expression import chip_overlap
from entropy_age.experiments.tissue_contrast import _features
from entropy_age.protection import matched_scores
from entropy_age.stochasticity import baseline_strata, drift_metrics

CONFIG = Path(__file__).resolve().parents[3] / "config" / "series"
N_DONOR_BOOT = {"paper": 200, "quick": 10}


def specs() -> dict[str, dict]:
    """Name -> parsed description of every file in ``config/series``."""
    return {p.stem: tomllib.loads(p.read_text()) for p in sorted(CONFIG.glob("*.toml"))}


def titles(matrix: Path, row: str = "!Sample_title") -> list[str]:
    """One sample row of the series matrix (titles by default; ``!Sample_source_name_ch1`` for source names)."""
    with gzip.open(matrix, "rt") as handle:
        for line in handle:
            if line.startswith(row):
                return [c.strip().strip('"') for c in line.rstrip("\n").split("\t")[1:]]
    return []


def from_title(sample_titles: list[str], column: str, rule: str) -> pl.Series:
    """``column`` read from every sample title by ``rule``, ``"title:<regex with one group>"``; every title must match."""
    if not rule.startswith("title:"):
        raise ValueError(f"series: {column} must be read as 'title:<regex>', not {rule!r}")
    pattern = re.compile(rule.removeprefix("title:"))
    found = [pattern.search(t) for t in sample_titles]
    missed = [t for m, t in zip(found, sample_titles, strict=True) if m is None]
    if missed:
        raise ValueError(f"series: {len(missed)} sample titles do not match the {column} pattern, e.g. {missed[0]!r}")
    return pl.Series(column, [m.group(1) for m in found])


def sample_table(chars: pl.DataFrame, sample_titles: list[str], spec: dict) -> pl.DataFrame:
    """``sample``, ``donor``, ``age``, ``sex``, the held columns and the columns the groups filter on."""
    if "donor" not in spec:
        raise ValueError("series: the description file must name the donor (a column or 'title:<regex>')")
    donor = spec["donor"]
    derived = [from_title(sample_titles, column, rule) for column, rule in spec.get("derive", {}).items()]
    if derived:
        chars = chars.with_columns(derived)
    if donor.startswith("title:"):
        donors = from_title(sample_titles, "donor", donor)
    else:
        donors = chars[donor].cast(pl.String).alias("donor")
    extra = sorted({*spec.get("hold", []), *(c for g in spec["groups"].values() for c in g)})
    sex = pl.col(spec["sex"]).alias("sex") if spec.get("sex") else pl.lit("NA").alias("sex")
    if "age_values" in spec:  # no numeric ages: a label per age group, e.g. {Young = 0, Old = 1}
        mapping = {str(k): float(v) for k, v in spec["age_values"].items()}
        age = pl.col(spec["age"]).cast(pl.String).replace_strict(mapping, default=None, return_dtype=pl.Float64)
    else:
        age = pl.col(spec["age"]).cast(pl.Float64, strict=False)
    table = chars.select("sample", donors, age.alias("age"), sex, *extra)
    if table["donor"].null_count():
        raise ValueError(f"series: {table['donor'].null_count()} arrays have no donor id")
    return table.drop_nulls("age")


def groups_of(table: pl.DataFrame, spec: dict, usable: list[str], exclude: set[str]) -> dict[str, pl.DataFrame]:
    """Group name -> its arrays: passing detection, not excluded, matching the group's filters, one per donor."""
    table = table.filter(pl.col("sample").is_in(usable) & ~pl.col("donor").is_in(list(exclude))).sort("sample")
    out = {}
    for name, filters in spec["groups"].items():
        part = table
        for column, value in filters.items():
            values = value if isinstance(value, list) else [value]
            part = part.filter(pl.col(column).is_in(values))
        out[name] = part.unique(subset="donor", keep="first", maintain_order=True)
    return out


def donor_intervals(x, ages, cov, feats, island, sets, n: int, rng) -> dict[tuple[str, str, str], tuple[float, float]]:
    """(band, metric, group) -> the 2.5th and 97.5th percentiles of the gap over ``n`` resamples of donors, with the
    strata, the band and the SD that scales the gap all taken from the full sample. A resample whose distinct donors
    leave fewer than 10 residual degrees of freedom is left out."""
    full = drift_metrics(x, ages, feats, covariates=cov if cov.shape[1] else None, spearman=False)
    baseline = full["baseline"].to_numpy()
    plan = {}
    for band, (lo, hi) in BANDS.items():
        keep = (baseline >= lo) & (baseline <= hi)
        strata = np.char.add(
            np.char.add(island[keep].astype(str), "|"), baseline_strata(baseline[keep], 10).astype(str)
        )
        plan[band] = (keep, strata, {m: float(np.std(full[m].to_numpy()[keep])) for m in METRICS})
    draws: dict[tuple[str, str, str], list[float]] = {}
    for _ in range(n):
        rows = rng.integers(0, ages.size, ages.size)
        c = cov[rows]
        c = c[:, np.ptp(c, axis=0) > 1e-12] if c.shape[1] else c
        if np.unique(rows).size - c.shape[1] - 2 < 10 or np.ptp(ages[rows]) == 0:
            continue  # too few distinct donors for the design: repeated donors add no degrees of freedom
        boot = drift_metrics(x[rows], ages[rows], feats, covariates=c if c.shape[1] else None, spearman=False)
        for band, (keep, strata, full_sd) in plan.items():
            for metric in METRICS:
                values = boot[metric].to_numpy()[keep]
                for group, (in_set, _) in sets.items():
                    # divided by the full sample's SD: a resample's own SD is inflated by sampling noise and would
                    # pull every resampled gap toward zero (see common.scaled_intervals)
                    gap = float(matched_scores(values[None, :], in_set[keep], strata)[0][0]) / full_sd[metric]
                    draws.setdefault((band, metric, group), []).append(gap)
    # fewer than half the resamples usable: the survivors are the donor draws that happen to be most diverse, and
    # their spread understates the sampling error (EPIC v2 myometrium, 30 donors on 10 chips, gave intervals of
    # 0.01 SD); no interval then
    return {
        k: (float(np.quantile(v, 0.025)), float(np.quantile(v, 0.975)), len(v))
        if len(v) >= n / 2
        else (float("nan"), float("nan"), len(v))
        for k, v in draws.items()
    }


def mouse_features(paths, rng: np.random.Generator, spec: dict) -> tuple[list[str], pl.DataFrame, pl.DataFrame]:
    """Mouse promoter CpGs, their annotation and gene-set membership. ``annotation = "mm285"``: the MM285 annotation
    built from Zhou's manifest and UCSC islands (``build_mm285_probes``); ``"mouse"``: the GPL32543 table, whose
    island flags are sparse. ``strip_probe_suffix = true`` for arrays that write ``cg...`` without the ``_BC21``
    suffix (the 40k mammalian array)."""
    from entropy_age.data.mouse import MouseAnnotationLayout, build_mm285_probes, mouse_membership
    from entropy_age.experiments.bxd_liver import probes_without_suffix
    from entropy_age.experiments.mouse_tissues import BACKGROUND, SETS

    lay = MouseAnnotationLayout(paths)
    if spec["annotation"] == "mm285":
        path = lay.probes_mm285 if lay.probes_mm285.exists() else build_mm285_probes(paths)
    else:
        path = lay.probes
    probes = pl.read_parquet(path).filter(pl.col("promoter") & ~pl.col("chr").is_in(["X", "Y", "MT"]))
    if spec.get("strip_probe_suffix"):
        probes = probes_without_suffix(probes)
    membership = mouse_membership(probes, pl.read_parquet(lay.gene_sets)).filter(pl.col("group").is_in(SETS))
    members = set(membership["feature"].to_list())
    rest = sorted(set(probes["cpg"].to_list()) - members)
    background = sorted(rng.choice(rest, size=min(BACKGROUND, len(rest)), replace=False).tolist())
    return sorted(members | set(background)), probes, membership


def run_series(ctx: Context, name: str) -> list[Path]:
    spec = specs()[name]
    gse = spec["accession"]
    matrix = ctx.paths.input / "geo" / gse / spec.get("matrix", f"{gse}_series_matrix.txt.gz")
    if not (reprocessed.available(ctx.paths, gse) and matrix.exists()):
        print(f"series {name}: no SeSAMe rebuild or series matrix for {gse}; nothing run", file=sys.stderr)
        return []
    out = ctx.out_dir(f"series_{name}")
    dropped: set[str] = set()
    if "exclude_cpgs" in spec:  # a table with a ``cpg`` column, relative to the repository (e.g. tumour-changed CpGs)
        dropped = set(pl.read_csv(CONFIG.parents[1] / spec["exclude_cpgs"])["cpg"].to_list())
    if spec.get("annotation", "human") in ("mouse", "mm285"):
        features, probes, membership = mouse_features(ctx.paths, ctx.rng(120_000), spec)
        probes = probes.filter(pl.col("cpg").is_in(features))
    else:
        features, _, membership = _features(ctx.paths, ctx.rng(120_000))
        probes = pl.read_parquet(AnnotationLayout(ctx.paths).probes_450k).filter(pl.col("cpg").is_in(features))
    if dropped:
        features = [f for f in features if f not in dropped]
        membership = membership.filter(~pl.col("feature").is_in(list(dropped)))
        print(f"series {name}: {len(dropped)} listed CpGs left out", file=sys.stderr)
    island_of = dict(probes.select("cpg", "cpg_island").iter_rows())
    qc = reprocessed.read_qc(ctx.paths, gse).filter(pl.col("passed"))
    chip_of = dict(qc.select("sample", "chip").iter_rows())
    exclude: set[str] = set()
    if spec.get("exclude_jaffe_donors"):
        jaffe = pl.read_parquet(cohort_layout("jaffe", ctx.paths).samples)
        exclude = set(jaffe.filter((pl.col("group") == "Control") & (pl.col("age") >= 20))["donor"].to_list())
    chars = _characteristics(matrix).with_columns(
        source_name=pl.Series(titles(matrix, "!Sample_source_name_ch1"))  # for series that keep tissue only there
    )
    table = sample_table(chars, titles(matrix), spec)
    hold = spec.get("hold", [])
    rows, sizes, overlap = [], [], []
    platform_of = dict(qc.select("sample", "platform").iter_rows()) if "platform" in qc.columns else {}
    for index, (label, s) in enumerate(groups_of(table, spec, qc["sample"].to_list(), exclude).items()):
        platforms = {platform_of[v] for v in s["sample"].to_list() if v in platform_of}
        if len(platforms) > 1:
            raise ValueError(f"series {name}: {label} mixes array platforms {sorted(platforms)}; they are never pooled")
        if s.height < MIN_SAMPLES:
            print(f"series {name}: {label}: {s.height} arrays, skipped", file=sys.stderr)
            continue
        betas = reprocessed.by_accession(ctx.paths, gse, {v: v for v in s["sample"].to_list()})
        s = s.filter(pl.col("sample").is_in(betas.columns)).sort("age", "sample")
        x, feats = complete_matrix(betas, features, s["sample"].to_list())
        if not feats:  # a probe-id mismatch between array and annotation must not pass as an empty result
            raise ValueError(
                f"series {name}: {label}: no annotated promoter CpG is measured in every array "
                f"({len(set(features) & set(betas['cpg'].to_list()))} of {len(features)} annotated CpGs on the array)"
            )
        ages = s["age"].to_numpy()
        chips = np.array([str(chip_of.get(v) or "NA") for v in s["sample"].to_list()])
        held = np.column_stack(
            [dummies(s["sex"].fill_null("NA").to_numpy()), *(dummies(s[c].fill_null("NA").to_numpy()) for c in hold)]
        )
        island = np.array([bool(island_of.get(f, False)) for f in feats])
        sets = _sets(feats, membership)
        lines = s[hold[0]].fill_null("NA").to_numpy() if hold else np.repeat("all", s.height)
        overlap += [{"group": label, **r} for r in chip_overlap(ages, lines, chips).to_dicts()]
        sizes.append(
            {
                "group": label,
                "donors": s.height,
                "age_min": float(ages.min()),
                "age_max": float(ages.max()),
                "cpgs": len(feats),
                "chips": int(np.unique(chips).size),
            }
        )
        designs = {"held": held, "held + chip": np.column_stack([held, dummies(chips)])}
        for k, (design, cov) in enumerate(designs.items()):
            rng = ctx.rng(350_000 + 10 * index + k)
            cov = cov[:, np.ptp(cov, axis=0) > 1e-12] if cov.shape[1] else cov
            if s.height - cov.shape[1] - 2 < 10:
                continue
            by_donor = donor_intervals(x, ages, cov, feats, island, sets, N_DONOR_BOOT[str(ctx.scale.name)], rng)
            for r in group_gaps(x, ages, cov, feats, island_of, membership, N_RESAMPLES[str(ctx.scale.name)], ctx.seed):
                lo, hi, kept = by_donor.get((r["band"], r["metric"], r["group"]), (float("nan"), float("nan"), 0))
                rows.append(
                    {
                        "cohort": label,
                        "covariates": design,
                        **r,
                        "donor_ci_low": lo,
                        "donor_ci_high": hi,
                        "donor_resamples": kept,
                    }
                )
        print(f"series {name}: {label}: {s.height} donors, {len(feats)} CpGs", file=sys.stderr, flush=True)
    for table_name, data in (("sizes", sizes), ("gaps", rows), ("chip_overlap", overlap)):
        if data:
            write_table(pl.DataFrame(data), out / table_name)
    return []
