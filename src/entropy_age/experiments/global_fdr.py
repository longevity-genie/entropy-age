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

SKIP_PREFIXES = ("p_bh", "p_fdr", "p_adj", "p_holm")
FOLDER = "global_fdr"


def p_columns(columns: list[str]) -> list[str]:
    """The raw p-value columns of a table: ``p`` or ``p_*``, adjusted ones left out."""
    return [c for c in columns if (c == "p" or c.startswith("p_")) and not c.startswith(SKIP_PREFIXES)]


def ledger(output: Path) -> pl.DataFrame:
    """One row per reported p: folder, table, row, column, p and the row's other string fields as a label."""
    rows = []
    for path in sorted(output.glob("*/*.csv")):
        if path.parent.name == FOLDER:
            continue
        try:
            t = pl.read_csv(path, infer_schema_length=10_000)
        except Exception:  # a malformed or empty table is not a test
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
                rows.append(
                    {
                        "folder": path.parent.name,
                        "table": path.stem,
                        "row": i,
                        "column": c,
                        "p": float(min(max(p, 0.0), 1.0)),
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
