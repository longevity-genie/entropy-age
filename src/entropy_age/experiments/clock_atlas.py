"""The clock atlas: every scored clock in every rebuilt series group (RESEARCH_PLAN, ideas 1 and 6; Anton, 2026-11-24:
"having so many samples and not computing all possible and legit methylation clocks on them is a crime").

Input: ``data/interim/clocks/series/<GSE>_clocks.csv`` (scripts/clocks_biolearn.py on each rebuilt series; clocks
failing the coverage rule are absent) and the groups of ``config/series`` (ages, donors, one array per donor, as the
gap runs use them). Per group and clock: arrays, the correlation of the prediction with age, the median absolute
error where the clock's output is years, and whether the clock is valid in that tissue (trained on it, or
pan-tissue); a clock outside its tissue is kept as a negative control, never dropped. Per group: the median
correlation between the accelerations (residual on age) of its valid adult clocks, the agreement.

Writes ``clock_atlas/clocks`` and ``clock_atlas/agreement`` and a figure: groups by clock, coloured by the
correlation with age, clocks ordered by generation, invalid cells hatched. ``entropy-age figures clock_atlas``.
"""

from __future__ import annotations

import re
import sys
from pathlib import Path

import numpy as np
import polars as pl

from entropy_age import plots
from entropy_age.data import reprocessed
from entropy_age.experiments.clock_increment import GENERATION
from entropy_age.experiments.common import Context, write_table
from entropy_age.experiments.series import groups_of, series_table, specs

#: clock scores from the clock track (protocol change 4: SeSAMe without the design-flag mask)
CLOCKS = Path("data/interim/clocks/series_unmasked")
MIN_ARRAYS = 15
#: clocks trained on blood (or saliva) only; any other kept clock is treated as valid everywhere except where
#: TISSUE_ONLY names its tissue
BLOOD = {
    "Hannum",
    "Lin",
    "PhenoAge",
    "HRSInCHPhenoAge",
    "GrimAgeV1",
    "GrimAgeV2",
    "DunedinPACE",
    "Zhang_10",
    "Weidner",
    "Garagnani",
    "VidalBralo",
    "Bocklandt",
    "YingCausAge",
    "YingDamAge",
    "YingAdaptAge",
}
TISSUE_ONLY = {"DNAmClockCortical": "cortex|brain|amygdala|hippocamp|neuron", "PEDBE": "buccal|cheek"}
PAN = {"Horvathv1", "Horvathv2", "EpiTOC1", "StocH"}
BLOOD_WORDS = re.compile(r"blood|pbmc|leukocyte|monocyte|lymphocyte|cd4|cd8|neutrophil", re.I)
YEARS = "Age (Years)"


def valid(clock: str, group: str) -> bool:
    if clock in PAN:
        return True
    if clock in BLOOD:
        return bool(BLOOD_WORDS.search(group))
    if clock in TISSUE_ONLY:
        return bool(re.search(TISSUE_ONLY[clock], group, re.I))
    return False  # trait, gestational and other predictors: negative controls in adult tissue


def residual(y: np.ndarray, x: np.ndarray) -> np.ndarray:
    coef = np.polyfit(x, y, 1)
    return y - np.polyval(coef, x)


def run(ctx: Context) -> list[Path]:
    rows, agreement = [], []
    for name, spec in specs().items():
        if spec.get("annotation", "human") != "human":
            continue
        gse = spec["accession"]
        path = ctx.paths.root / CLOCKS / f"{gse}_clocks.csv"
        table = series_table(ctx.paths, spec)
        if not (path.exists() and table is not None and reprocessed.available(ctx.paths, gse)):
            print(f"clock_atlas: {name}: no clock scores or matrix for {gse}", file=sys.stderr)
            continue
        scores = pl.read_csv(path)
        usable = reprocessed.read_qc(ctx.paths, gse).filter(pl.col("passed"))["sample"].to_list()
        for label, s in groups_of(table, spec, usable, set()).items():
            wide = (
                s.select("sample", "age")
                .join(scores.pivot(on="clock", index="sample", values="prediction"), on="sample", how="inner")
                .drop_nulls("age")
            )
            if wide.height < MIN_ARRAYS:
                continue
            age = wide["age"].to_numpy().astype(float)
            output = dict(scores.select("clock", "output").unique().iter_rows())
            accs = {}
            for clock in [c for c in wide.columns if c not in ("sample", "age")]:
                pred = wide[clock].to_numpy().astype(float)
                if not np.all(np.isfinite(pred)) or np.std(pred) == 0:
                    continue
                ok = valid(clock, label)
                rows.append(
                    {
                        "series": name,
                        "accession": gse,
                        "group": label,
                        "clock": clock,
                        "generation": GENERATION.get(clock, "other"),
                        "valid": ok,
                        "arrays": wide.height,
                        "age_min": float(age.min()),
                        "age_max": float(age.max()),
                        "r_age": float(np.corrcoef(pred, age)[0, 1]),
                        "median_abs_error": float(np.median(np.abs(pred - age)))
                        if output.get(clock) == YEARS
                        else None,
                    }
                )
                if ok and output.get(clock) == YEARS:
                    accs[clock] = residual(pred, age)
            if len(accs) >= 2:
                r = np.corrcoef(np.column_stack(list(accs.values())), rowvar=False)
                pairs = r[np.triu_indices_from(r, k=1)]
                agreement.append(
                    {
                        "series": name,
                        "group": label,
                        "arrays": wide.height,
                        "clocks": len(accs),
                        "median_pair_r": float(np.median(pairs)),
                        "min_pair_r": float(pairs.min()),
                        "clock_names": ", ".join(accs),
                    }
                )
    if not rows:
        print("clock_atlas: no series with clock scores yet", file=sys.stderr)
        return []
    out = ctx.out_dir("clock_atlas")
    table = pl.DataFrame(rows)
    write_table(table, out / "clocks", csv_max_rows=10_000)
    write_table(pl.DataFrame(agreement), out / "agreement")
    return [draw(table, ctx.figure_path("clock_atlas"))]


def draw(table: pl.DataFrame, path: Path) -> Path:
    gens = [
        "1st: chronological",
        "2nd: phenotypic and mortality",
        "3rd: pace of ageing",
        "mitotic and stochastic",
        "other",
    ]
    clocks = (
        table.group_by("clock", "generation")
        .agg(pl.col("valid").sum().alias("n_valid"))
        .with_columns(rank=pl.col("generation").replace_strict({g: i for i, g in enumerate(gens)}, default=9))
        .sort("rank", "n_valid", "clock", descending=[False, True, False])["clock"]
        .to_list()
    )
    groups = table.select("series", "group").unique(maintain_order=True).sort("series", "group")
    labels = [f"{g} ({n})" for n, g in groups.iter_rows()]
    grid = np.full((len(labels), len(clocks)), np.nan)
    invalid = np.zeros_like(grid, dtype=bool)
    index = {(s, g): i for i, (s, g) in enumerate(groups.iter_rows())}
    for row in table.iter_rows(named=True):
        i, j = index[(row["series"], row["group"])], clocks.index(row["clock"])
        grid[i, j] = row["r_age"]
        invalid[i, j] = not row["valid"]
    plots.apply_style()
    fig, ax = plots.plt.subplots(figsize=(0.24 * len(clocks) + 4, 0.26 * len(labels) + 2.2))
    im = ax.imshow(grid, cmap="RdBu_r", vmin=-1, vmax=1, aspect="auto")
    yy, xx = np.nonzero(invalid & np.isfinite(grid))
    ax.scatter(xx, yy, marker="x", s=8, color=plots.MUTED, lw=0.6, label="outside the clock's tissue")
    ax.set_xticks(range(len(clocks)), clocks, rotation=90, fontsize=6)
    ax.set_yticks(range(len(labels)), labels, fontsize=6)
    ax.set_title("Each clock's correlation with age in each tissue group (rebuilt raw arrays)", fontsize=9, loc="left")
    fig.colorbar(im, ax=ax, fraction=0.02, pad=0.01, label="r with age")
    ax.legend(fontsize=6, frameon=False, loc="upper left", bbox_to_anchor=(0, -0.25))
    fig.tight_layout()
    return plots.save(fig, path)
