"""Do the headline gaps depend on probe design or unreliable probes? (roadmap RM18)

The 450k array carries two probe designs. Infinium I probes measure both alleles with two bead types and have a
wider dynamic range; Infinium II probes use one bead and two colours and compress toward the middle. CpG-island
probes, where essential promoters sit, are more often type I. Matching on island status and starting level
(``tissue_contrast._with_strata``) does not separate the designs, so part of a gap could be design.

Some probes are unreliable whatever the design: they cross-hybridise, carry a SNP at the CpG or the extension
base, or map poorly. Zhou, Laird and Shen 2017 list them; ``M_general`` is their recommended mask.

This recomputes the point estimates of ``figures tissues`` (Hannum blood and Jaffe cortex age slopes, Lunnon's
same-donor tissues) four ways (``VARIANTS``): as published here; with probe type added to the strata; with masked
probes left out; with both. A gap that is design or bad probes moves; a gap that is biology stays.

Run as a job after ``entropy-age fetch annotation`` (which now also fetches the mask) and the cohorts of
``figures tissues``: ``entropy-age figures probe_robustness``.
"""

from __future__ import annotations

import sys
from pathlib import Path

import polars as pl
from matplotlib.lines import Line2D

from entropy_age import plots
from entropy_age.data.annotation import ZHOU_MASK_450K, AnnotationLayout, read_probe_mask
from entropy_age.experiments.common import Context, write_table
from entropy_age.experiments.tissue_contrast import (
    LIFESPAN_METRICS,
    SET_LABELS,
    SETS,
    _dispersion_metrics,
    _inputs_missing,
    _lifespan_metrics,
    _load,
    _matched_sd,
)
from entropy_age.stochasticity import baseline_strata

VARIANTS = {
    "as published": (False, False),
    "probe type in strata": (True, False),
    "masked probes out": (False, True),
    "both": (True, True),
}
PERMUTATIONS = 1000


def strata(metrics: pl.DataFrame, probes: pl.DataFrame, *, by_type: bool) -> pl.DataFrame:
    """``metrics`` with ``stratum``: island status x baseline decile, and the probe design when ``by_type``."""
    cols = ["cpg_island", "probe_type"] if by_type else ["cpg_island"]
    out = metrics.join(probes.select(pl.col("cpg").alias("feature"), *cols), on="feature", how="left")
    if by_type and out["probe_type"].null_count():
        raise ValueError("probe_robustness: CpGs without a probe type; re-run `entropy-age fetch annotation`")
    decile = pl.Series(baseline_strata(out["baseline"].to_numpy(), 10))
    parts = [pl.col("cpg_island").fill_null(False).cast(pl.String), decile.cast(pl.String)]
    if by_type:
        parts.append(pl.col("probe_type"))
    return out.with_columns(stratum=pl.concat_str(parts, separator="|"))


def keep_unmasked(metrics: pl.DataFrame, masked: set[str]) -> pl.DataFrame:
    return metrics.filter(~pl.col("feature").is_in(list(masked)))


def run(ctx: Context) -> list[Path]:
    mask_path = AnnotationLayout(ctx.paths).input / ZHOU_MASK_450K
    if _inputs_missing(ctx, "probe_robustness"):
        return []
    if not mask_path.exists():
        print(
            f"probe_robustness: {mask_path} missing; run `entropy-age fetch annotation`. Nothing run.", file=sys.stderr
        )
        return []
    plots.apply_style()
    out = ctx.out_dir("probe_robustness")
    lifespan, lunnon, probes, membership = _load(ctx)
    if "probe_type" not in probes.columns:
        raise ValueError(
            "probe_robustness: the 450k annotation has no probe_type; re-run `entropy-age fetch annotation`"
        )
    mask = read_probe_mask(mask_path)
    masked = set(mask.filter(pl.col("masked"))["cpg"].to_list())
    tables = []
    jobs = [("lifespan", k, _lifespan_metrics(p), LIFESPAN_METRICS) for k, p in lifespan.items()]
    jobs += [("same donors", k, _dispersion_metrics(p), ("dispersion",)) for k, p in lunnon.items()]
    for part, tissue, metrics, names in jobs:
        for variant, (by_type, drop_masked) in VARIANTS.items():
            m = keep_unmasked(metrics, masked) if drop_masked else metrics
            t = _matched_sd(strata(m, probes, by_type=by_type), membership, names, PERMUTATIONS, 1)
            tables.append(
                t.with_columns(part=pl.lit(part), tissue=pl.lit(tissue), variant=pl.lit(variant), cpgs=pl.lit(m.height))
            )
            print(f"probe_robustness: {tissue}, {variant}: done", file=sys.stderr, flush=True)
    table = pl.concat(tables, how="diagonal")
    write_table(table, out / "gaps")
    write_table(_design_share(probes, membership, masked), out / "design_share")
    return [_figure(ctx, table)]


def _design_share(probes: pl.DataFrame, membership: pl.DataFrame, masked: set[str]) -> pl.DataFrame:
    """Share of type I and of masked probes among each set's CpGs and among all promoter CpGs."""
    p = probes.with_columns(masked=pl.col("cpg").is_in(list(masked)))
    rows = [
        {
            "group": "all promoter CpGs",
            "cpgs": p.height,
            "type_I": float((p["probe_type"] == "I").mean()),
            "masked": float(p["masked"].mean()),
        }
    ]
    for group in SETS:
        cpgs = membership.filter(pl.col("group") == group)["feature"].unique()
        s = p.filter(pl.col("cpg").is_in(cpgs.implode()))
        rows.append(
            {
                "group": group,
                "cpgs": s.height,
                "type_I": float((s["probe_type"] == "I").mean()),
                "masked": float(s["masked"].mean()),
            }
        )
    return pl.DataFrame(rows)


def _figure(ctx: Context, table: pl.DataFrame) -> Path:
    parts = [("lifespan", "entropy_slope"), ("lifespan", "dispersion_slope"), ("same donors", "dispersion")]
    fig, axes = plots.panels(3, ncols=3, size=4.2)
    for i, (part, metric) in enumerate(parts):
        ax = axes[i]
        sel = table.filter((pl.col("part") == part) & (pl.col("metric") == metric))
        labels, y = [], 0.0
        for tissue in sel["tissue"].unique(maintain_order=True).to_list():
            for g, group in enumerate(SETS):
                for v, variant in enumerate(VARIANTS):
                    r = sel.filter(
                        (pl.col("tissue") == tissue) & (pl.col("group") == group) & (pl.col("variant") == variant)
                    )
                    if r.is_empty():
                        continue
                    ax.scatter(
                        r["matched_diff_sd"][0],
                        y + 0.18 * v,
                        s=14,
                        color=plots.PALETTE[g],
                        marker="osD^"[v],
                        alpha=1.0 if v else 0.6,
                    )
                labels.append((y + 0.27, f"{tissue} · {SET_LABELS[group]}"))
                y += 1
            y += 0.5
        ax.axvline(0, color=plots.MUTED, lw=0.8)
        ax.set_yticks([p for p, _ in labels], [t for _, t in labels], fontsize=6)
        ax.invert_yaxis()
        ax.set(xlabel="Set minus matched CpGs (SD)", title=f"{part}: {metric.replace('_', ' ')}")
        plots.label_panel(ax, "abc"[i])
    handles = [
        Line2D([], [], marker="osD^"[v], ls="", color=plots.MUTED, label=variant) for v, variant in enumerate(VARIANTS)
    ]
    fig.legend(handles=handles, loc="lower center", ncol=4, fontsize=8, frameon=False, bbox_to_anchor=(0.5, -0.04))
    fig.suptitle(
        "Gaps with probe design in the strata and unreliable probes masked (Zhou 2017)",
        x=0.02,
        y=1.04,
        ha="left",
        fontsize=11,
        fontweight="bold",
    )
    return plots.save(fig, ctx.figure_path("probe_robustness"))
