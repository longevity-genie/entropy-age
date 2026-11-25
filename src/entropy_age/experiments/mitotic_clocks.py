"""Mitotic clocks beyond EpiTOC1 (RESEARCH_PLAN idea 3: is the shield a mitotic clock?).

EpiTOC1 tracks the per-person shield in every dividing tissue (clock_drift). Three other published mitotic clocks test
whether that is a property of cell-division counting or of EpiTOC1's CpGs:

- **epiTOC2** (Teschendorff 2020): the estimated cumulative number of stem-cell divisions,
  ``tnsc = 2 * mean_i (beta_i - beta0_i) / (delta_i * (1 - beta0_i))`` over 163 promoter CpGs, with each CpG's
  per-division gain ``delta`` and ground state ``beta0``.
- **RepliTali** (Endicott 2022): replicative age from 87 solo-WCGW CpGs in partially methylated domains (losing
  methylation with division), a linear model fitted on the EPIC array; most of its CpGs are not on the 450k.
- **epiCMIT** (Duran-Ferrer 2020): the larger of the mean beta over its hypermethylating CpGs and one minus the mean
  over its hypomethylating ones.

Coefficients come from ``scripts/fetch_mitotic_clocks.py`` (data/input/mitotic). A clock is scored on a group of
arrays only when the CpGs it needs (weighted by |coefficient| for RepliTali, counted for the others) are at least 90 %
present (the clock-coverage rule); a missing value in one array takes the group mean for that CpG.
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import polars as pl

DIR = Path("data/input/mitotic")
MIN_COVERAGE = 0.9
NAMES = ("epiTOC2", "RepliTali", "epiCMIT")


def load(root: Path) -> dict[str, pl.DataFrame] | None:
    d = root / DIR
    if not (d / "epitoc2.csv").exists():
        return None
    return {
        "epiTOC2": pl.read_csv(d / "epitoc2.csv"),
        "RepliTali": pl.read_csv(d / "replitali.csv"),
        "epiCMIT": pl.read_csv(d / "epicmit.csv"),
    }


def cpgs(coef: dict[str, pl.DataFrame]) -> set[str]:
    return (
        set(coef["epiTOC2"]["cpg"])
        | set(coef["RepliTali"].filter(pl.col("term") != "(Intercept)")["term"])
        | set(coef["epiCMIT"]["cpg"])
    )


def _matrix(betas: pl.DataFrame, wanted: list[str], samples: list[str]) -> tuple[np.ndarray, np.ndarray]:
    """CpGs x arrays for ``wanted`` (rows in that order, NaN where the CpG is absent), the group mean filling gaps,
    and a mask of rows present in the group."""
    m = pl.DataFrame({"cpg": wanted}).join(betas.select("cpg", *samples), on="cpg", how="left")
    x = m.select(samples).to_numpy().astype(float)
    present = ~np.isnan(x)
    have = present.any(axis=1)
    mean = np.zeros(len(wanted))
    mean[have] = np.nanmean(x[have], axis=1)
    return np.where(present, x, mean[:, None]), have


def score(betas: pl.DataFrame, samples: list[str], coef: dict[str, pl.DataFrame]) -> pl.DataFrame:
    """Long frame (sample, clock, prediction) for every mitotic clock that passes coverage on these arrays."""
    out = []
    e = coef["epiTOC2"]
    x, have = _matrix(betas, e["cpg"].to_list(), samples)
    if have.mean() >= MIN_COVERAGE:
        d, b0 = e["delta"].to_numpy()[have], e["beta0"].to_numpy()[have]
        tnsc = 2 * ((x[have] - b0[:, None]) / (d * (1 - b0))[:, None]).mean(axis=0)
        out.append(("epiTOC2", tnsc))
    r = coef["RepliTali"]
    intercept = float(r.filter(pl.col("term") == "(Intercept)")["coefficient"][0])
    r = r.filter(pl.col("term") != "(Intercept)")
    x, have = _matrix(betas, r["term"].to_list(), samples)
    w = r["coefficient"].to_numpy()
    if np.abs(w[have]).sum() / np.abs(w).sum() >= MIN_COVERAGE:
        out.append(("RepliTali", intercept + w[have] @ x[have]))
    c = coef["epiCMIT"]
    x, have = _matrix(betas, c["cpg"].to_list(), samples)
    if have.mean() >= MIN_COVERAGE:
        hyper = (c["class"] == "hyper").to_numpy()
        h1 = x[have & hyper].mean(axis=0)
        h2 = 1 - x[have & ~hyper].mean(axis=0)
        out.append(("epiCMIT", np.maximum(h1, h2)))
    return pl.DataFrame(
        [
            {"sample": s, "clock": name, "prediction": float(v)}
            for name, vals in out
            for s, v in zip(samples, vals, strict=True)
        ],
        schema={"sample": pl.String, "clock": pl.String, "prediction": pl.Float64},
    )
