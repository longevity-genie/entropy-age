"""Mitotic clocks beyond EpiTOC1 (RESEARCH_PLAN idea 3: is the shield a mitotic clock?).

EpiTOC1 tracks the per-person shield in every dividing tissue (clock_drift). Other published mitotic clocks test
whether that is a property of cell-division counting or of EpiTOC1's own CpGs:

- **epiTOC2** (Teschendorff 2020 Genome Med): the estimated cumulative number of stem-cell divisions,
  ``tnsc = 2 * mean_i (beta_i - beta0_i) / (delta_i * (1 - beta0_i))`` over 163 promoter CpGs that gain methylation,
  with each CpG's per-division gain ``delta`` and ground state ``beta0`` (the table the package's ``epiTOC2()`` uses).
  **epiTOC3** is the package's later version of the same model on 170 CpGs.
- **HypoClock** (Zhou et al. 2018 Nat Genet): mean beta over 678 solo-WCGW CpGs in partially methylated domains,
  which lose methylation with division (lower = more divisions). Measured on the 450k and EPIC alike.
- **RepliTali** (Endicott et al. 2022 Nat Commun): replicative age from 87 CpGs, 18 of them PMD solo-WCGWs, a linear
  model fitted on the EPIC array; most of its CpGs are not on the 450k.
- **epiCMIT** (Duran-Ferrer et al. 2020 Nat Cancer): the larger of the mean beta over its 184 hypermethylating CpGs
  and one minus the mean over its 1,164 hypomethylating ones.

Coefficients come from ``scripts/fetch_mitotic_clocks.py`` (data/input/mitotic, pinned commit in SOURCE).

**Coverage** is judged on the arrays a group analyses (Codex review 2026-11-25: judging it on a whole export let
RepliTali pass on 450k arrays that measure 46 % of its weight): each array's measured share of the clock's weight
(|coefficient| for RepliTali, the per-CpG multiplier ``1 / (delta (1 - beta0))`` for epiTOC2/3, equal weights within
each epiCMIT arm and for HypoClock, each epiCMIT arm on its own) must be at least 90 % in at least 90 % of the arrays,
the rule used for the biolearn clocks. A CpG measured in no array of the group drops out of the model (the score is
then a reduced model; its weight share is reported); a CpG missing in some arrays takes the group mean there.
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import polars as pl

DIR = Path("data/input/mitotic")
MIN_SHARE = 0.9  # of the clock's weight measured in one array
MIN_ARRAYS = 0.9  # share of arrays that must reach MIN_SHARE
NAMES = ("epiTOC2", "epiTOC3", "HypoClock", "RepliTali", "epiCMIT")


def load(root: Path) -> dict[str, pl.DataFrame] | None:
    d = root / DIR
    if not (d / "epitoc2.csv").exists() or not (d / "hypoclock.csv").exists():
        return None
    return {
        "epiTOC2": pl.read_csv(d / "epitoc2.csv"),
        "epiTOC3": pl.read_csv(d / "epitoc3.csv"),
        "HypoClock": pl.read_csv(d / "hypoclock.csv"),
        "RepliTali": pl.read_csv(d / "replitali.csv"),
        "epiCMIT": pl.read_csv(d / "epicmit.csv"),
        "EpiTOC1": pl.read_csv(d / "epitoc1_cpgs.csv"),
    }


def cpgs(coef: dict[str, pl.DataFrame]) -> set[str]:
    """Every CpG any mitotic clock reads, EpiTOC1's included (kept out of the drift measures B and F)."""
    out: set[str] = set()
    for name, frame in coef.items():
        col = "term" if name == "RepliTali" else "cpg"
        out |= set(frame[col].to_list())
    return out - {"(Intercept)"}


def _matrix(betas: pl.DataFrame, wanted: list[str], samples: list[str]) -> tuple[np.ndarray, np.ndarray]:
    """CpGs x arrays for ``wanted`` (rows in that order), NaN where not measured, and the measured mask."""
    m = pl.DataFrame({"cpg": wanted}).join(betas.select("cpg", *samples), on="cpg", how="left")
    if m.height != len(wanted):
        raise ValueError("duplicate CpG rows in the betas")
    x = m.select(samples).to_numpy().astype(float)
    return x, ~np.isnan(x)


def _covered(present: np.ndarray, w: np.ndarray) -> tuple[bool, float]:
    """Whether the weights ``w`` (one per row) are measured well enough, and the median per-array measured share."""
    share = (np.abs(w)[:, None] * present).sum(axis=0) / np.abs(w).sum()
    return bool((share >= MIN_SHARE).mean() >= MIN_ARRAYS), float(np.median(share))


def _filled(x: np.ndarray, present: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """Rows measured in at least one array (mask) and the matrix with the group mean in the gaps."""
    have = present.any(axis=1)
    mean = np.zeros(x.shape[0])
    mean[have] = np.nanmean(x[have], axis=1)
    return have, np.where(present, x, mean[:, None])


def score(betas: pl.DataFrame, samples: list[str], coef: dict[str, pl.DataFrame]) -> tuple[pl.DataFrame, pl.DataFrame]:
    """Score every mitotic clock on these arrays (one analysed group). Returns (long frame sample, clock, prediction
    for the clocks that pass coverage; coverage frame clock, passed, median measured weight share, CpGs used)."""
    out, cover = [], []

    def note(name: str, ok: bool, share: float, used: int) -> None:
        cover.append({"clock": name, "passed": ok, "median_weight_share": share, "cpgs_used": used})

    for name in ("epiTOC2", "epiTOC3"):
        e = coef[name]
        x, present = _matrix(betas, e["cpg"].to_list(), samples)
        d, b0 = e["delta"].to_numpy(), e["beta0"].to_numpy()
        ok, share = _covered(present, 1 / (d * (1 - b0)))
        have, x = _filled(x, present)
        note(name, ok, share, int(have.sum()))
        if ok:
            out.append((name, 2 * ((x[have] - b0[have, None]) / (d[have] * (1 - b0[have]))[:, None]).mean(axis=0)))
    h = coef["HypoClock"]
    x, present = _matrix(betas, h["cpg"].to_list(), samples)
    ok, share = _covered(present, np.ones(h.height))
    have, x = _filled(x, present)
    note("HypoClock", ok, share, int(have.sum()))
    if ok:
        out.append(("HypoClock", x[have].mean(axis=0)))
    r = coef["RepliTali"]
    intercept = float(r.filter(pl.col("term") == "(Intercept)")["coefficient"][0])
    r = r.filter(pl.col("term") != "(Intercept)")
    x, present = _matrix(betas, r["term"].to_list(), samples)
    w = r["coefficient"].to_numpy()
    ok, share = _covered(present, w)
    have, x = _filled(x, present)
    note("RepliTali", ok, share, int(have.sum()))
    if ok:
        out.append(("RepliTali", intercept + w[have] @ x[have]))
    c = coef["epiCMIT"]
    hyper = (c["class"] == "hyper").to_numpy()
    x, present = _matrix(betas, c["cpg"].to_list(), samples)
    ok_hyper, share_hyper = _covered(present[hyper], np.ones(int(hyper.sum())))
    ok_hypo, share_hypo = _covered(present[~hyper], np.ones(int((~hyper).sum())))
    have, x = _filled(x, present)
    note("epiCMIT", ok_hyper and ok_hypo, min(share_hyper, share_hypo), int(have.sum()))
    if ok_hyper and ok_hypo:
        out.append(("epiCMIT", np.maximum(x[have & hyper].mean(axis=0), 1 - x[have & ~hyper].mean(axis=0))))
    scores = pl.DataFrame(
        [
            {"sample": s, "clock": name, "prediction": float(v)}
            for name, vals in out
            for s, v in zip(samples, vals, strict=True)
        ],
        schema={"sample": pl.String, "clock": pl.String, "prediction": pl.Float64},
    )
    return scores, pl.DataFrame(cover)
