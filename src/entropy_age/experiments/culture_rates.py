"""Drift per day or per division, across every culture series and agent (roadmap RM39 item 5).

R1 and the RM5 runs (``passaged_cells``) scored every cultured sample against its experiment's first day:
look-alike drift ``B`` matched to each gene set, with each sample's population doublings
and days in culture. Within one standard series doublings and days rise together; the series' other arms pull
them apart: mitomycin arrest and serum withdrawal (days without divisions), reduced oxygen and TERT (Endicott,
GSE197512), dexamethasone, galactose, 2-deoxyglucose, beta-hydroxybutyrate and SURF1 mutants (Sturm,
GSE179847). Over all samples of a series, each score is regressed jointly on doublings and days, each counted
from its experiment's start, with one intercept per experiment (line x sub-experiment): the drift each division
adds with time held, and each day adds with divisions held. Reads the saved per-sample tables; seconds.

``entropy-age figures culture_rates`` (after ``figures passage`` and ``figures passage_sturm``).
"""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import polars as pl
from scipy import stats

from entropy_age.experiments.common import Context, write_table

SERIES = {"Endicott (GSE197512)": "passaged_cells", "Sturm (GSE179847)": "passaged_cells_sturm"}
SCORES = (
    "B_housekeeping",
    "B_core_essential",
    "B_non_essential",
)  # not F: its fall with PD is the array-noise artefact (extension.md, RM5)


def joint_rates(s: pl.DataFrame, y: str) -> dict | None:
    """Score on doublings and days (per 10 of each) with one intercept per experiment; OLS SEs."""
    s = s.drop_nulls(["doublings", "days", y]).filter(pl.col("doublings").is_finite() & pl.col(y).is_finite())
    if s.height < 10:
        return None
    exp = (s["subexperiment"] + "|" + s["donor"]).to_numpy()
    levels = sorted(set(exp))
    design = np.column_stack(
        [s["doublings"].to_numpy() / 10, s["days"].to_numpy() / 10, *((exp == e).astype(float) for e in levels)]
    )
    v = s[y].to_numpy()
    coef, *_ = np.linalg.lstsq(design, v, rcond=None)
    resid = v - design @ coef
    dof = v.size - np.linalg.matrix_rank(design)
    se = np.sqrt(np.diag((resid @ resid / dof) * np.linalg.pinv(design.T @ design)))
    out = {"score": y, "samples": v.size, "experiments": len(levels)}
    out["r_doublings_days"] = float(np.corrcoef(s["doublings"], s["days"])[0, 1])
    for i, name in ((0, "per_10_doublings"), (1, "per_10_days")):
        out[name] = float(coef[i])
        out[f"{name}_se"] = float(se[i])
        out[f"{name}_p"] = float(2 * stats.t.sf(abs(coef[i] / se[i]), dof))
    return out


def run(ctx: Context) -> list[Path]:
    out = ctx.out_dir("culture_rates")
    rows = []
    for series, folder in SERIES.items():
        path = ctx.paths.output / folder / "sample_scores.csv"
        if not path.exists():
            print(f"culture_rates: {path} missing; run the passage figures first", file=sys.stderr)
            continue
        s = pl.read_csv(path)
        for y in SCORES:
            if y in s.columns and (r := joint_rates(s, y)):
                rows.append({"series": series, **r})
    table = pl.DataFrame(rows)
    write_table(table, out / "rates")
    return []
