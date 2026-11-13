"""Do the headline gaps survive reprocessing from raw IDATs? (roadmap RM18 (3))

Every cohort test runs on the betas the authors deposited. Here Jaffe's prefrontal cortex series (GSE74193) is
reprocessed from its IDATs (``scripts/reprocess_idats.py``: methylprep with SeSAMe-like steps: noob, non-linear
dye-bias correction, channel-switch inference, quality mask) into ``data/interim/cohorts/jaffe_betas_sesame.parquet``.
The adult controls of ``figures tissues`` (same samples, plate, sex, race and neuron fraction held) are then
measured on both: the age slopes of entropy and dispersion, each gene set against look-alikes of the same island
status and level decile, in SD, on the CpGs present in both. A gap that is a preprocessing artefact changes when the
preprocessing does; one in the biology does not.

``entropy-age figures reprocessed_jaffe`` as a Slurm job, after R61a wrote the reprocessed betas.
"""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import polars as pl

from entropy_age import plots
from entropy_age.data.matrix import FeatureMatrix
from entropy_age.experiments.common import Context, dummies, write_table
from entropy_age.experiments.tissue_contrast import (
    LIFESPAN_METRICS,
    MAX_MISSING,
    Prepared,
    _features,
    _lifespan_metrics,
    _matched_sd,
    _with_strata,
)

REPROCESSED = "cohorts/jaffe_betas_sesame.parquet"
DEPOSITED = "cohorts/jaffe_betas.parquet"
SAMPLES = "cohorts/jaffe_samples.parquet"
SOURCES = ("deposited", "reprocessed from IDATs")


def adult_controls(samples: pl.DataFrame) -> pl.DataFrame:
    """The samples ``figures tissues`` uses: adult controls, one per donor (best QC first), with age and neurons."""
    return (
        samples.filter((pl.col("group") == "Control") & ~pl.col("drop_sample") & (pl.col("age") >= 20))
        .sort("best_qc", descending=True)
        .unique(subset="donor", keep="first", maintain_order=True)
        .drop_nulls(["age", "neuron_fraction"])
        .sort("age", "sample")
    )


def _prepared(path: Path, s: pl.DataFrame, features: list[str]) -> Prepared:
    frame = pl.scan_parquet(path).select("cpg", *s["sample"].to_list()).filter(pl.col("cpg").is_in(features)).collect()
    m = FeatureMatrix.from_frame(frame, feature_col="cpg", samples=s["sample"].to_list(), max_missing=MAX_MISSING)
    cov = np.column_stack(
        [*(dummies(s[c].to_numpy()) for c in ("plate", "sex", "race")), s["neuron_fraction"].to_numpy()]
    )
    return Prepared("cortex (Jaffe)", m.values, s["age"].to_numpy(), cov, m.features, s["donor"].to_numpy(), cov)


def run(ctx: Context) -> list[Path]:
    interim = ctx.paths.interim
    if not (interim / REPROCESSED).exists():
        print(f"reprocessed_jaffe: {interim / REPROCESSED} missing (run R61a); nothing run", file=sys.stderr)
        return []
    plots.apply_style()
    out = ctx.out_dir("reprocessed_jaffe")
    features, probes, membership = _features(ctx.paths, ctx.rng(120_000))
    s = adult_controls(pl.read_parquet(interim / SAMPLES))
    have = set(pl.read_parquet(interim / REPROCESSED, n_rows=0).columns)
    missing = [x for x in s["sample"].to_list() if x not in have]
    if missing:
        print(f"reprocessed_jaffe: {len(missing)} samples not reprocessed, left out", file=sys.stderr)
        s = s.filter(pl.col("sample").is_in(list(have)))
    preps = {
        SOURCES[0]: _prepared(interim / DEPOSITED, s, features),
        SOURCES[1]: _prepared(interim / REPROCESSED, s, features),
    }
    common = sorted(set(preps[SOURCES[0]].features) & set(preps[SOURCES[1]].features))
    print(f"reprocessed_jaffe: {s.height} donors, {len(common)} CpGs in both", file=sys.stderr, flush=True)
    rows, agree = [], []
    metrics = {}
    for source, p in preps.items():
        idx = [p.features.index(c) for c in common]
        q = Prepared(p.label, p.x[:, idx], p.ages, p.covariates, common, p.donors, p.technical)
        metrics[source] = _lifespan_metrics(q)
        table = _matched_sd(_with_strata(metrics[source], probes), membership, LIFESPAN_METRICS, 1000, ctx.seed)
        rows.append(table.with_columns(source=pl.lit(source), donors=pl.lit(s.height), cpgs=pl.lit(len(common))))
    for metric in LIFESPAN_METRICS:
        a = metrics[SOURCES[0]][metric].to_numpy()
        b = metrics[SOURCES[1]][metric].to_numpy()
        agree.append({"metric": metric, "per_cpg_correlation": float(np.corrcoef(a, b)[0, 1])})
    table = pl.concat(rows)
    write_table(table, out / "gaps")
    write_table(pl.DataFrame(agree), out / "per_cpg_agreement")
    return [_figure(ctx, table)]


def _figure(ctx: Context, table: pl.DataFrame) -> Path:
    fig, axes = plots.panels(1, ncols=1, size=5.0)
    ax = axes[0]
    keys = table.filter(pl.col("source") == SOURCES[0]).select("metric", "group")
    y = np.arange(keys.height)
    for i, source in enumerate(SOURCES):
        t = keys.join(table.filter(pl.col("source") == source), on=["metric", "group"], how="left")
        ax.barh(y + 0.38 * i, t["matched_diff_sd"], height=0.36, color=plots.PALETTE[i], label=source)
    ax.set_yticks(y + 0.19, [f"{m.replace('_slope', '')} · {g.replace('_', ' ')}" for m, g in keys.iter_rows()])
    ax.axvline(0, color=plots.MUTED, lw=0.8)
    ax.invert_yaxis()
    ax.set(xlabel="age gap, set minus look-alikes (SD)", title="Jaffe cortex: deposited against IDAT-reprocessed")
    ax.legend(fontsize=7, loc="lower right")
    return plots.save(fig, ctx.figure_path("reprocessed_jaffe"))
