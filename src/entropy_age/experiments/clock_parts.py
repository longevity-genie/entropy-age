"""Shielded and unshielded parts of a clock, tissue by tissue (Anton, 2026-11-24: "check how shielded and non-shielded
components of the clock contribute to correlation, unless shielded components are too small").

The share check came first (clock_weights, R64): shielded CpGs (promoters of core-essential or housekeeping genes)
carry 1-8 % of the age signal of Hannum, Zhang, the cortical and the Ying clocks, too little to split, and 9-26 % of
Horvath 2013, StocH, StocP, Lin and PhenoAge. Those five are split here.

Per clock and rebuilt series group (the groups of ``config/series``, one array per donor, as the gap runs use them):
the clock's linear predictor over its shielded CpGs and over its other CpGs (``sum of weight x beta``, before any
transform; a CpG missing in an array takes the group's mean for that CpG, and the measured share is reported). For each
part: its correlation with age, and its share of the age signal in that group (|cov(part, age)| over the sum of both
parts'). Also the correlation between the two parts' accelerations (residual on age): whether the shielded part calls
the same donors old. Writes ``clock_parts/parts`` and a figure. ``entropy-age figures clock_parts`` (needs the clock
atlas inputs, data/interim/clocks/series/<GSE>_cpgs.parquet).
"""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import polars as pl

from entropy_age import plots
from entropy_age.data import reprocessed
from entropy_age.data.annotation import AnnotationLayout, link_membership
from entropy_age.data.cohorts import _characteristics
from entropy_age.experiments.clock_cpgs import clock_table
from entropy_age.experiments.clock_weights import PROTECTED
from entropy_age.experiments.common import Context, write_table
from entropy_age.experiments.series import groups_of, sample_table, specs, titles

CLOCKS = ("Horvath 2013", "StocH", "StocP", "Lin", "PhenoAge")
BETAS = Path("data/interim/clocks/series")
MIN_ARRAYS = 15


def residual(y: np.ndarray, x: np.ndarray) -> np.ndarray:
    return y - np.polyval(np.polyfit(x, y, 1), x)


def run(ctx: Context) -> list[Path]:
    ann = AnnotationLayout(ctx.paths)
    mem = link_membership(pl.read_parquet(ann.links_450k), pl.read_parquet(ann.gene_sets))
    shielded = set(mem.filter(pl.col("group").is_in(PROTECTED))["feature"].to_list())
    weights = clock_table(ctx.paths).filter(pl.col("clock").is_in(CLOCKS) & pl.col("cpg").str.starts_with("cg"))
    rows = []
    for name, spec in specs().items():
        if spec.get("annotation", "human") != "human":
            continue
        gse = spec["accession"]
        path = ctx.paths.root / BETAS / f"{gse}_cpgs.parquet"
        matrix = ctx.paths.input / "geo" / gse / spec.get("matrix", f"{gse}_series_matrix.txt.gz")
        if not (path.exists() and matrix.exists() and reprocessed.available(ctx.paths, gse)):
            continue
        cpgs = sorted(set(weights["cpg"].to_list()))
        betas = pl.read_parquet(path).filter(pl.col("cpg").is_in(cpgs))
        chars = _characteristics(matrix).with_columns(source_name=pl.Series(titles(matrix, "!Sample_source_name_ch1")))
        table = sample_table(chars, titles(matrix), spec)
        usable = reprocessed.read_qc(ctx.paths, gse).filter(pl.col("passed"))["sample"].to_list()
        for label, s in groups_of(table, spec, usable, set()).items():
            s = s.drop_nulls("age").filter(pl.col("sample").is_in(betas.columns))
            if s.height < MIN_ARRAYS:
                continue
            age = s["age"].to_numpy().astype(float)
            sub = betas.select("cpg", *s["sample"].to_list())
            for clock in CLOCKS:
                w = weights.filter(pl.col("clock") == clock).select("cpg", "weight")
                m = w.join(sub, on="cpg", how="left")
                x = m.select(s["sample"].to_list()).to_numpy().astype(float)  # CpGs x arrays
                present = ~np.isnan(x)
                row_mean = np.where(present.any(axis=1), np.nanmean(np.where(present, x, np.nan), axis=1), np.nan)
                keep = ~np.isnan(row_mean)  # a CpG absent from the whole group drops out of both parts
                x = np.where(present, x, row_mean[:, None])[keep]
                wt = m["weight"].to_numpy()[keep]
                prot = np.array([c in shielded for c in m["cpg"].to_list()])[keep]
                if prot.sum() < 3 or (~prot).sum() < 3:
                    continue
                part_s = wt[prot] @ x[prot]
                part_o = wt[~prot] @ x[~prot]
                cov_s = abs(np.cov(part_s, age)[0, 1])
                cov_o = abs(np.cov(part_o, age)[0, 1])
                total_w = np.abs(m["weight"].to_numpy()).sum()
                rows.append(
                    {
                        "series": name,
                        "group": label,
                        "clock": clock,
                        "arrays": s.height,
                        "shielded_cpgs": int(prot.sum()),
                        "other_cpgs": int((~prot).sum()),
                        "weight_measured": float(np.abs(wt).sum() / total_w),
                        "values_measured": float(present[keep].mean()),
                        "r_age_shielded": float(np.corrcoef(part_s, age)[0, 1]),
                        "r_age_other": float(np.corrcoef(part_o, age)[0, 1]),
                        "r_age_whole": float(np.corrcoef(part_s + part_o, age)[0, 1]),
                        "signal_share_shielded": float(cov_s / (cov_s + cov_o)) if cov_s + cov_o > 0 else float("nan"),
                        "r_acc_parts": float(np.corrcoef(residual(part_s, age), residual(part_o, age))[0, 1]),
                    }
                )
    if not rows:
        print("clock_parts: no series with exported betas", file=sys.stderr)
        return []
    out = ctx.out_dir("clock_parts")
    parts = pl.DataFrame(rows)
    write_table(parts, out / "parts")
    return [draw(parts, ctx.figure_path("clock_parts"))]


def draw(parts: pl.DataFrame, path: Path) -> Path:
    plots.apply_style()
    fig, axes = plots.plt.subplots(1, 2, figsize=(11, 4.2))
    colours = {c: plots.PALETTE[i] for i, c in enumerate(CLOCKS)}
    ax = axes[0]
    for clock in CLOCKS:
        d = parts.filter(pl.col("clock") == clock)
        ax.scatter(d["r_age_other"], d["r_age_shielded"], s=18, color=colours[clock], label=clock, alpha=0.8)
    ax.plot([-0.2, 1], [-0.2, 1], color=plots.MUTED, ls="--", lw=0.8)
    ax.set(
        xlabel="r with age, clock's other CpGs",
        ylabel="r with age, clock's shielded CpGs",
        title="a  Each part's correlation with age, per tissue group",
    )
    ax.legend(fontsize=7, frameon=False)
    ax = axes[1]
    for i, clock in enumerate(CLOCKS):
        d = parts.filter(pl.col("clock") == clock)
        ax.scatter(
            np.full(d.height, i) + np.random.default_rng(i).normal(0, 0.06, d.height),
            d["signal_share_shielded"],
            s=16,
            color=colours[clock],
            alpha=0.8,
        )
    ax.set_xticks(range(len(CLOCKS)), CLOCKS, rotation=20, fontsize=8)
    ax.set(
        ylabel="shielded CpGs' share of the age signal",
        title="b  Share of the age signal from shielded CpGs",
        ylim=(0, 1),
    )
    fig.tight_layout()
    return plots.save(fig, path)
