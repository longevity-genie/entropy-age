"""Which part of drift each clock generation reads, tissue by tissue (RESEARCH_PLAN idea 2; Figure 10b in EPIC-Italy
blood, tested here in every rebuilt series group with ages).

Per series group (the groups of ``config/series``; one array per donor): each person's drift at housekeeping look-alike
promoters (B) and shield share (F), scored against a reference fitted on everyone else in the group as
``incident_cancer`` does (``loo_scores``), on the CpGs measured in every array. Each scored clock's acceleration is its
prediction's residual on age within the group. Per group and clock: the correlation of the acceleration with B and
with F, and the clock's generation. Clocks outside their tissue are kept, labelled by ``clock_atlas.valid``.

Arrays are chosen by ``series.selection`` as the gap runner chooses them (QC, configured donor and CpG exclusions,
one array per donor), among arrays present in every input. The mitotic clocks of ``mitotic_clocks`` are scored on
each group's own arrays from the clock track, with coverage judged there; their CpGs and EpiTOC1's are kept out of B
and F. Rows from sensitivity configs (``series.SENSITIVITY``) are flagged and left out of the figure.

Writes ``clock_drift/correlations`` and a figure (per generation, the spread of r with F and with B across tissue
groups). ``entropy-age figures clock_drift`` after the clock panel per series (data/interim/clocks/series/).
"""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import polars as pl

from entropy_age import plots
from entropy_age.data.annotation import AnnotationLayout
from entropy_age.experiments import mitotic_clocks
from entropy_age.experiments.clock_atlas import valid
from entropy_age.experiments.clock_increment import GENERATION
from entropy_age.experiments.common import Context, write_table
from entropy_age.experiments.health_states import loo_scores
from entropy_age.experiments.per_person import _islands, _sets
from entropy_age.experiments.series import SENSITIVITY, selection, specs
from entropy_age.experiments.tissue_contrast import _features

#: drift (B, F) on the standard rebuild; clock scores on the clock track (protocol change 4)
SERIES = Path("data/interim/clocks/series")
CLOCK_SCORES = Path("data/interim/clocks/series_unmasked")
GROUP = "housekeeping"
MIN_ARRAYS = 25
GENS = ("1st: chronological", "2nd: phenotypic and mortality", "3rd: pace of ageing", "mitotic and stochastic")


def residual(y: np.ndarray, x: np.ndarray) -> np.ndarray:
    return y - np.polyval(np.polyfit(x, y, 1), x)


def run(ctx: Context) -> list[Path]:
    features, _, membership = _features(ctx.paths, ctx.rng(120_000))
    probes = pl.read_parquet(AnnotationLayout(ctx.paths).probes_450k).filter(pl.col("cpg").is_in(features))
    mitotic = mitotic_clocks.load(ctx.paths.root)  # epiTOC2/3, HypoClock, RepliTali, epiCMIT, scored per group here
    # every CpG a mitotic clock reads (EpiTOC1's included) stays out of B and F, so a clock and the drift measures
    # never share a probe (Codex review 2026-11-25: epiTOC and epiCMIT CpGs sat in the matched background)
    clock_cpgs = mitotic_clocks.cpgs(mitotic) if mitotic is not None else set()
    rows, coverage = [], []
    for name, spec in specs().items():
        if spec.get("annotation", "human") != "human":
            continue
        gse = spec["accession"]
        betas_path = ctx.paths.root / SERIES / f"{gse}_cpgs.parquet"
        clocks_path = ctx.paths.root / CLOCK_SCORES / f"{gse}_clocks.csv"
        unmasked = ctx.paths.root / CLOCK_SCORES / f"{gse}_cpgs.parquet"
        if not (betas_path.exists() and clocks_path.exists()):
            continue
        long = pl.read_csv(clocks_path, schema_overrides={"sample": pl.String}).select("sample", "clock", "prediction")
        betas = pl.read_parquet(betas_path).filter(pl.col("cpg").is_in(features) & ~pl.col("cpg").is_in(clock_cpgs))
        mb = None
        available = set(betas.columns) & set(long["sample"].to_list())
        if mitotic is not None and unmasked.exists():
            mb = pl.scan_parquet(unmasked).filter(pl.col("cpg").is_in(clock_cpgs)).collect()
            available &= set(mb.columns)
        chosen = selection(ctx.paths, name, spec, available)
        if chosen is None:
            continue
        groups, dropped = chosen
        if dropped:
            betas = betas.filter(~pl.col("cpg").is_in(dropped))
        base = long.cast({"prediction": pl.Float64}).pivot(on="clock", index="sample", values="prediction")
        if base.height != long["sample"].n_unique():
            raise ValueError(f"clock_drift: {name}: duplicate sample ids in {clocks_path.name}")
        for label, s in groups.items():
            s = s.drop_nulls("age")
            if s.height < MIN_ARRAYS:
                continue
            samples = s["sample"].to_list()
            scores = base
            if mb is not None:
                extra, cover = mitotic_clocks.score(mb, samples, mitotic)
                coverage += [{"series": name, "group": label, **r} for r in cover.to_dicts()]
                if extra.height:
                    scores = base.join(
                        extra.pivot(on="clock", index="sample", values="prediction"), on="sample", how="left"
                    )
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
            clk = s.select("sample").join(scores, on="sample", how="left", validate="1:1")
            if clk["sample"].to_list() != samples:
                raise ValueError(f"clock_drift: {name} / {label}: clock scores out of the arrays' order")
            for clock in [c for c in clk.columns if c != "sample"]:
                v = clk[clock].to_numpy().astype(float)
                if not np.isfinite(v).all() or np.std(v) == 0:
                    continue
                acc = residual(v, age)
                rows.append(
                    {
                        "series": name,
                        "group": label,
                        "sensitivity": bool(SENSITIVITY.search(name)),
                        "arrays": len(samples),
                        "cpgs": len(kept),
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
    write_table(table, out / "correlations", csv_max_rows=10_000)
    if coverage:
        write_table(pl.DataFrame(coverage), out / "mitotic_coverage", csv_max_rows=10_000)
    return [draw(table.filter(~pl.col("sensitivity")), ctx.figure_path("clock_drift"))]


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
        ax.set_xticks(range(len(GENS)), [g.replace(": ", "\n") for g in GENS], fontsize=7)
        ax.set_title(title, loc="left", fontsize=9)
    axes[0].set_ylabel("r, clock acceleration with the drift measure\n(each dot: one clock in one tissue group)")
    fig.tight_layout()
    return plots.save(fig, path)
