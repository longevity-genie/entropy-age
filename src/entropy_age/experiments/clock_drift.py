"""Which part of drift each clock generation reads, tissue by tissue (RESEARCH_PLAN idea 2; Figure 10b in EPIC-Italy
blood, tested here in every rebuilt series group with ages).

Per series group (the groups of ``config/series``; one array per donor): each person's drift at housekeeping look-alike
promoters (B) and shield share (F), scored against a reference fitted on everyone else in the group as
``incident_cancer`` does (``loo_scores``), on the CpGs measured in every array. Each scored clock's acceleration is its
prediction's residual on age within the group. Per group and clock: the correlation of the acceleration with B and
with F, and the clock's generation. Clocks outside their tissue are kept, labelled by ``clock_atlas.valid``.

Writes ``clock_drift/correlations`` and a figure (per generation, the spread of r with F and with B across tissue
groups). ``entropy-age figures clock_drift`` after the clock panel per series (data/interim/clocks/series/).
"""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import polars as pl

from entropy_age import plots
from entropy_age.data import reprocessed
from entropy_age.data.annotation import AnnotationLayout
from entropy_age.experiments.clock_atlas import valid
from entropy_age.experiments.clock_increment import GENERATION
from entropy_age.experiments.common import Context, write_table
from entropy_age.experiments.health_states import loo_scores
from entropy_age.experiments.per_person import _islands, _sets
from entropy_age.experiments.series import groups_of, series_table, specs
from entropy_age.experiments.tissue_contrast import _features

SERIES = Path("data/interim/clocks/series")
GROUP = "housekeeping"
MIN_ARRAYS = 25
GENS = ("1st: chronological", "2nd: phenotypic and mortality", "3rd: pace of ageing", "mitotic and stochastic")


def residual(y: np.ndarray, x: np.ndarray) -> np.ndarray:
    return y - np.polyval(np.polyfit(x, y, 1), x)


def run(ctx: Context) -> list[Path]:
    features, _, membership = _features(ctx.paths, ctx.rng(120_000))
    probes = pl.read_parquet(AnnotationLayout(ctx.paths).probes_450k).filter(pl.col("cpg").is_in(features))
    rows = []
    for name, spec in specs().items():
        if spec.get("annotation", "human") != "human":
            continue
        gse = spec["accession"]
        betas_path = ctx.paths.root / SERIES / f"{gse}_cpgs.parquet"
        clocks_path = ctx.paths.root / SERIES / f"{gse}_clocks.csv"
        table = series_table(ctx.paths, spec)
        if not (betas_path.exists() and clocks_path.exists() and table is not None):
            continue
        if not reprocessed.available(ctx.paths, gse):
            continue
        betas = pl.read_parquet(betas_path).filter(pl.col("cpg").is_in(features))
        scores = pl.read_csv(clocks_path).pivot(on="clock", index="sample", values="prediction")
        usable = reprocessed.read_qc(ctx.paths, gse).filter(pl.col("passed"))["sample"].to_list()
        for label, s in groups_of(table, spec, usable, set()).items():
            s = s.drop_nulls("age").filter(
                pl.col("sample").is_in(betas.columns) & pl.col("sample").is_in(scores["sample"])
            )
            if s.height < MIN_ARRAYS:
                continue
            samples = s["sample"].to_list()
            sub = betas.select("cpg", *samples)
            x = sub.select(samples).to_numpy().T.astype(float)
            ok = ~np.isnan(x).any(axis=0)
            kept = [c for c, k in zip(sub["cpg"].to_list(), ok, strict=True) if k]
            sets = _sets(kept, membership)
            if GROUP not in sets:
                continue
            print(
                f"clock_drift: {name} / {label}: {len(samples)} arrays, {len(kept)} CpGs", file=sys.stderr, flush=True
            )
            per = loo_scores(x[:, ok], _islands(kept, probes), sets)
            b, f = per[f"B_{GROUP}"], per[f"F_{GROUP}"]
            age = s["age"].to_numpy().astype(float)
            clk = s.select("sample").join(scores, on="sample", how="left")
            for clock in [c for c in clk.columns if c != "sample"]:
                v = clk[clock].to_numpy().astype(float)
                if not np.isfinite(v).all() or np.std(v) == 0:
                    continue
                acc = residual(v, age)
                rows.append(
                    {
                        "series": name,
                        "group": label,
                        "arrays": len(samples),
                        "clock": clock,
                        "generation": GENERATION.get(clock, "other"),
                        "valid": valid(clock, label),
                        "r_B": float(np.corrcoef(acc, b)[0, 1]) if np.std(b) > 0 else float("nan"),
                        "r_F": float(np.corrcoef(acc, f)[0, 1]) if np.std(f) > 0 else float("nan"),
                    }
                )
    if not rows:
        print("clock_drift: nothing to score", file=sys.stderr)
        return []
    out = ctx.out_dir("clock_drift")
    table = pl.DataFrame(rows)
    write_table(table, out / "correlations")
    return [draw(table, ctx.figure_path("clock_drift"))]


def draw(table: pl.DataFrame, path: Path) -> Path:
    plots.apply_style()
    fig, axes = plots.plt.subplots(1, 2, figsize=(11, 4), sharey=True)
    d = table.filter(pl.col("generation").is_in(GENS))
    for ax, col, title in (
        (axes[0], "r_F", "a  with the shield share F"),
        (axes[1], "r_B", "b  with look-alike drift B"),
    ):
        for i, g in enumerate(GENS):
            v = d.filter(pl.col("generation") == g)[col].drop_nans().to_numpy()
            if v.size:
                ax.scatter(
                    np.full(v.size, i) + np.random.default_rng(i).normal(0, 0.07, v.size),
                    v,
                    s=10,
                    alpha=0.5,
                    color=plots.PALETTE[i],
                )
                ax.plot([i - 0.25, i + 0.25], [np.median(v)] * 2, color=plots.INK, lw=2)
        ax.axhline(0, color=plots.MUTED, lw=0.8, ls="--")
        ax.set_xticks(range(len(GENS)), [g.split(":")[0] + "\n" + g.split(": ")[1] for g in GENS], fontsize=7)
        ax.set_title(title, loc="left", fontsize=9)
    axes[0].set_ylabel("r, clock acceleration with the drift measure\n(each dot: one clock in one tissue group)")
    fig.tight_layout()
    return plots.save(fig, path)
