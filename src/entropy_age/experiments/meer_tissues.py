"""Four mouse tissues across four ages, RRBS: is the shield missing where cells do not divide? (roadmap RM2)

Meer et al. 2018 (GSE121141) profiled heart, liver, lung and cerebral cortex of male C57BL/6 mice at 6, 12, 20
and 30 months, about 20 per tissue (RRBS, Bismark coverage files). On the mammalian array (GSE279437) the shield
was missing in cortex and hippocampus but kept in heart. This is the RRBS replication with exact ages.

Per tissue: promoter CpGs (RefSeq TSS −1500 to +500, mm10) covered ``MIN_COVERAGE`` times in at least 90 % of the
tissue's mice (gaps filled by the CpG's median); the age slopes of entropy and dispersion, each set against
look-alikes of the same island status and 6-month level decile (``interventions.aging_gaps``), in SD, with an
interval from resampling mice; and the per-mouse shield share ``F`` against the 6-month mice, its change per 10
months.

``entropy-age figures meer_tissues`` as a Slurm job, after GSE121141's RAW.tar is unpacked into
``data/input/geo/GSE121141/cov``.
"""

from __future__ import annotations

import gzip
import re
import sys
from pathlib import Path

import numpy as np
import polars as pl

from entropy_age import plots
from entropy_age.data.matrix import FeatureMatrix
from entropy_age.data.rrbs import TRACKS, RrbsLayout, in_intervals, islands, promoters
from entropy_age.experiments.bxd_liver import ols
from entropy_age.experiments.common import Context, scaled_intervals, write_table
from entropy_age.experiments.interventions import _sets, aging_gaps, membership
from entropy_age.experiments.per_person import score_against
from entropy_age.experiments.reprogramming_rrbs import read_cov
from entropy_age.protection import fit_reference

GSE = "GSE121141"
MIN_COVERAGE = 10
MIN_PRESENT = 0.9
YOUNG_MONTHS = 6.0
N_BOOT = {"paper": 200, "quick": 30}


def read_mice(matrix: Path) -> pl.DataFrame:
    """``sample`` (GSM), ``tissue`` and ``age_months`` from the series matrix."""
    rows: dict[str, list[str]] = {}
    with gzip.open(matrix, "rt") as handle:
        for line in handle:
            cells = [c.strip().strip('"') for c in line.rstrip("\n").split("\t")]
            if cells[0] == "!Sample_geo_accession":
                rows["sample"] = cells[1:]
            elif cells[0] == "!Sample_source_name_ch1":
                rows["tissue"] = cells[1:]
            elif cells[0] == "!Sample_characteristics_ch1" and cells[1].startswith("age"):
                rows["age_months"] = [c.split(":", 1)[1].strip() for c in cells[1:]]
    return pl.DataFrame(rows).with_columns(pl.col("age_months").cast(pl.Float64))


def _tissue_matrix(files: dict[str, Path], samples: list[str], prom: pl.DataFrame) -> pl.DataFrame:
    long = []
    for s in samples:
        t = read_cov(files[s]).filter(pl.col("coverage") >= MIN_COVERAGE)
        chrom = t["cpg"].str.split(":").list.get(0).to_numpy()
        pos = t["cpg"].str.split(":").list.get(1).cast(pl.Int64).to_numpy() - 1
        long.append(
            t.filter(pl.Series(in_intervals(chrom, pos, prom))).select("cpg", "beta").with_columns(sample=pl.lit(s))
        )
    long = pl.concat(long)
    kept = long.group_by("cpg").len().filter(pl.col("len") >= MIN_PRESENT * len(samples)).select("cpg")
    return long.join(kept, on="cpg", how="semi").pivot(on="sample", index="cpg", values="beta").sort("cpg")


def run(ctx: Context) -> list[Path]:
    base = ctx.paths.input / "geo" / GSE
    if not (base / "cov").exists():
        print(f"meer_tissues: {base / 'cov'} missing; nothing run", file=sys.stderr)
        return []
    plots.apply_style()
    out = ctx.out_dir("meer_tissues")
    lay = RrbsLayout(ctx.paths)
    prom = promoters(lay.track(TRACKS[0]), lay.track(TRACKS[1]))
    isl_tracks = islands(lay.track(TRACKS[2]))
    mice = read_mice(base / f"{GSE}_series_matrix.txt.gz")
    files = {}
    for f in (base / "cov").glob("GSM*cov*"):
        m = re.match(r"(GSM\d+)", f.name)
        if m:
            files[m.group(1)] = f
    mice = mice.filter(pl.col("sample").is_in(list(files)))
    rng = ctx.rng(290_000)
    gap_rows, f_rows = [], []
    for (tissue,), sel in mice.group_by("tissue", maintain_order=True):
        sel = sel.sort("age_months", "sample")
        wide = _tissue_matrix(files, sel["sample"].to_list(), prom)
        m = FeatureMatrix.from_frame(
            wide, feature_col="cpg", samples=sel["sample"].to_list(), max_missing=1 - MIN_PRESENT
        )
        features = m.features
        cpgs = pl.DataFrame({"cpg": features}).with_columns(
            chr=pl.col("cpg").str.split(":").list.get(0), pos=pl.col("cpg").str.split(":").list.get(1).cast(pl.Int64)
        )
        isl = in_intervals(cpgs["chr"].to_numpy(), cpgs["pos"].to_numpy() - 1, isl_tracks)
        members = membership(ctx.paths, features, cpgs)
        x, ages = m.values, sel["age_months"].to_numpy()
        # aging_gaps reads its young level from 3-6 month mice: the 6-month mice here
        point = aging_gaps(x, ages, features, isl, members)
        boot = []
        for b in range(N_BOOT[str(ctx.scale.name)]):
            rows = np.sort(rng.choice(x.shape[0], size=x.shape[0], replace=True))
            boot += [{"boot": b, **r} for r in aging_gaps(x, ages, features, isl, members, rows)]
        gaps = scaled_intervals(pl.DataFrame(point), pl.DataFrame(boot), ["metric", "group"])
        gap_rows.append(gaps.with_columns(tissue=pl.lit(tissue), mice=pl.lit(x.shape[0]), cpgs=pl.lit(len(features))))
        sets = _sets(features, members)
        member = np.any([s for s, _ in sets.values()], axis=0)
        young = np.flatnonzero(ages <= YOUNG_MONTHS)
        older = np.flatnonzero(ages > YOUNG_MONTHS)
        sc = score_against(x[older], fit_reference(x[young], isl, ~member), sets)
        design = np.column_stack([np.ones(older.size), ages[older] / 10.0])
        for group in sets:
            coef, se, p = ols(sc[f"F_{group}"], design, 1)
            f_rows.append({"tissue": tissue, "group": group, "F_per_10_months": coef, "se": se, "p": p})
        print(f"meer_tissues: {tissue}: {len(features)} CpGs, {x.shape[0]} mice", file=sys.stderr, flush=True)
    gaps = pl.concat(gap_rows)
    shares = pl.DataFrame(f_rows)
    write_table(gaps, out / "aging_gaps")
    write_table(shares, out / "share_trend")
    return [_figure(ctx, gaps)]


def _figure(ctx: Context, gaps: pl.DataFrame) -> Path:
    fig, axes = plots.panels(2, ncols=2, size=4.4)
    groups = ["core_essential", "housekeeping", "non_essential"]
    tissues = gaps["tissue"].unique(maintain_order=True).to_list()
    for ax, metric, letter in zip(axes, ("entropy_slope", "dispersion_slope"), "ab", strict=True):
        for i, group in enumerate(groups):
            t = gaps.filter((pl.col("metric") == metric) & (pl.col("group") == group))
            t = pl.DataFrame({"tissue": tissues}).join(t, on="tissue", how="left")
            y = np.arange(len(tissues)) + 0.25 * (i - 1)
            ax.errorbar(
                t["gap_sd"],
                y,
                xerr=[t["gap_sd"] - t["ci_low"], t["ci_high"] - t["gap_sd"]],
                fmt="o",
                ms=4,
                lw=1,
                color=plots.PALETTE[i],
                label=group.replace("_", " "),
            )
        ax.set_yticks(range(len(tissues)), tissues if letter == "a" else [""] * len(tissues))
        ax.axvline(0, color=plots.MUTED, lw=0.8)
        ax.invert_yaxis()
        ax.set(xlabel="set minus look-alikes (SD)", title=metric.replace("_slope", " slope"))
        if letter == "a":
            ax.legend(fontsize=7)
        plots.label_panel(ax, letter)
    fig.suptitle("Four mouse tissues, 6-30 months (Meer 2018, RRBS)", x=0.02, y=1.03, ha="left", fontsize=11)
    fig.subplots_adjust(wspace=0.15)
    return plots.save(fig, ctx.figure_path("meer_tissues"))
