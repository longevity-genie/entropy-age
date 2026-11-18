"""One false-discovery correction over every p-value the project reports (roadmap RM34 item 3).

Most runs correct within themselves (BH over their own table), and R46 over its 168 tests. The omics re-analysis
practice asks for the correction over the full set tested, not the set reported. This reads every table under
``data/output/*/`` that carries a p-value column (``p``, ``p_permutation``, ``p_bootstrap``, ``p_t``, ``p_perm``,
``p_boot`` and the like; adjusted columns such as ``p_bh`` are skipped), pools them into one ledger and applies one
Benjamini-Hochberg correction across all of them. Per output folder it reports how many tests were nominally
significant and how many survive the global correction. Reads tables only; no data are reprocessed.

``entropy-age figures global_fdr`` (light; seconds).
"""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import polars as pl
from scipy.stats import false_discovery_control

from entropy_age import plots
from entropy_age.experiments.common import Context, write_table

ADJUSTED = ("bh", "fdr", "adj", "holm", "q", "bonferroni")
FOLDER = "global_fdr"


def p_columns(columns: list[str]) -> list[str]:
    """The raw p-value columns of a table, adjusted ones left out.

    Raw: ``p``, ``p_*`` (``p_permutation``, ``p_boot``) and ``*_p`` (``directed_p``, ``B_p_welch`` reads as
    ``*_p_*``). Adjusted: any name with a part ``bh``, ``fdr``, ``adj``, ``holm``, ``q`` or ``bonferroni``.
    (Until 2026-11-18 only ``p`` and ``p_*`` were read; an outside audit found 366 p-values left out.)
    """
    out = []
    for c in columns:
        parts = c.split("_")  # case-sensitive: ``P_core_essential`` is the per-person shield score, not a p
        if "p" not in parts or any(a in (q.lower() for q in parts) for a in ADJUSTED):
            continue
        out.append(c)
    return out


def ledger(output: Path) -> pl.DataFrame:
    """One row per reported p: folder, table, row, column, p and the row's other string fields as a label."""
    rows = []
    for path in sorted(output.glob("*/*.csv")):
        if path.parent.name == FOLDER:
            continue
        try:
            t = pl.read_csv(path, infer_schema_length=10_000)
        except Exception as err:  # a malformed or empty table is not a test; say so
            print(f"global_fdr: skipped {path.parent.name}/{path.name}: {type(err).__name__}", file=sys.stderr)
            continue
        cols = p_columns(t.columns)
        if not cols:
            continue
        label_cols = [c for c, d in t.schema.items() if d == pl.String][:4]
        for c in cols:
            if not t.schema[c].is_numeric():
                continue
            for i, (p, *labels) in enumerate(t.select(c, *label_cols).iter_rows()):
                if p is None or not np.isfinite(p):
                    continue
                if not 0.0 <= p <= 1.0:  # a sentinel or an error, never evidence (audit 2026-11-18)
                    raise ValueError(f"global_fdr: {path.parent.name}/{path.name} row {i} {c} = {p} is not a p-value")
                rows.append(
                    {
                        "folder": path.parent.name,
                        "table": path.stem,
                        "row": i,
                        "column": c,
                        "p": float(p),
                        "label": " | ".join(str(v) for v in labels),
                    }
                )
    return pl.DataFrame(rows)


def run(ctx: Context) -> list[Path]:
    plots.apply_style()
    out = ctx.out_dir(FOLDER)
    led = ledger(ctx.paths.output)
    q = false_discovery_control(led["p"].to_numpy(), method="bh")
    led = led.with_columns(q_global=pl.Series(q))
    write_table(led, out / "ledger")
    summary = (
        led.group_by("folder")
        .agg(
            pl.len().alias("tests"),
            (pl.col("p") < 0.05).sum().alias("nominal"),
            (pl.col("q_global") < 0.05).sum().alias("survive"),
        )
        .with_columns(lost=pl.col("nominal") - pl.col("survive"))
        .sort("folder")
    )
    write_table(summary, out / "by_folder")
    cut = float(led.filter(pl.col("q_global") < 0.05)["p"].max() or 0.0)
    print(
        f"global_fdr: {led.height} tests in {summary.height} folders; {int(summary['nominal'].sum())} at p < 0.05, "
        f"{int(summary['survive'].sum())} survive (largest surviving p {cut:.4f})",
        file=sys.stderr,
    )
    return [_figure(ctx, summary)]


def _figure(ctx: Context, summary: pl.DataFrame) -> Path:
    s = summary.filter(pl.col("nominal") > 0).sort("nominal")
    fig, axes = plots.panels(1, ncols=1, size=max(4.0, 0.16 * s.height))
    ax = axes[0]
    y = np.arange(s.height)
    ax.barh(y, s["nominal"], color=plots.MUTED, label="p < 0.05")
    ax.barh(y, s["survive"], color=plots.PALETTE[0], label="q < 0.05, one correction over all")
    ax.set_yticks(y, s["folder"].to_list(), fontsize=6)
    ax.set(xlabel="tests", title="Every reported p, corrected together")
    ax.legend(fontsize=7, loc="lower right")
    return plots.save(fig, ctx.figure_path("global_fdr"))
