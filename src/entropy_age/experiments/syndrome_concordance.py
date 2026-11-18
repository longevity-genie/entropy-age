"""Do chromatin-writer variants age the methylome where aging does? (roadmap RM39 item 6)

R17-R18 (``syndromes``) found more look-alike drift in blood of people with NSD1 (Sotos), EZH2 (Weaver) and KMT2D
or CHD7 variants, and a weaker shield with NSD1 and EZH2. Progeroid mice and NSD1 loss speed up clocks (Sarni
2026; Martin-Herranz 2019); the literature check (docs/LITERATURE_COVERAGE.md, idea 9) left open whether such
variants act at the CpGs where aging is noisy, or at those where it is directed. From Hannum blood (the shield's
features and covariates), each CpG's age slope of the mean and of dispersion; **noisy-aging CpGs** are the top
``TOP`` by dispersion slope, **directed-aging CpGs** the top ``TOP`` by absolute mean slope among CpGs whose
dispersion slope is at or below the median. Then per syndrome arm against its series' controls (``syndromes``'
blood samples):

* **noise concordance:** the mean log variance ratio (carriers over controls) at noisy-aging CpGs, minus that at
  other CpGs of the same island status and young-level decile (weighted as the noisy CpGs fall into strata);
* **directed concordance:** the Spearman correlation, over directed-aging CpGs, of carriers' mean shift with the
  age slope (above 0: carriers look older).

Each is read against ``N_PERM`` permutations of the carrier labels within the arm and controls.

``entropy-age figures syndrome_concordance`` as a Slurm job (after ``fetch syndromes``).
"""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import polars as pl
from scipy import stats

from entropy_age.experiments.common import Context, write_table
from entropy_age.experiments.syndromes import CONTROL, MIN_ARM, SERIES, matrix_path, read_betas, read_people, usable
from entropy_age.experiments.tissue_contrast import _features, _hannum, _lifespan_metrics
from entropy_age.stochasticity import baseline_strata

TOP = 0.05
N_PERM = {"paper": 500, "quick": 50}


def site_classes(m: pl.DataFrame, top: float = TOP) -> pl.DataFrame:
    """``feature``, ``slope`` and ``cls`` (noisy, directed, other) from per-CpG age metrics."""
    disp = m["dispersion_slope"].to_numpy()
    noisy = disp >= np.quantile(disp, 1 - top)
    calm = disp <= np.median(disp)
    absl = m["abs_slope"].to_numpy()
    directed = calm & (absl >= np.quantile(absl[calm], 1 - top * len(absl) / max(calm.sum(), 1)))
    cls = np.where(noisy, "noisy", np.where(directed, "directed", "other"))
    return m.select("feature", "slope", "baseline").with_columns(cls=pl.Series(cls))


def noise_concordance(log_vr: np.ndarray, noisy: np.ndarray, strata: np.ndarray) -> float:
    """Mean log variance ratio at noisy CpGs minus the stratum-matched mean at the others."""
    other = ~noisy
    total, weight = 0.0, 0
    for s in np.unique(strata[noisy]):
        a, b = noisy & (strata == s), other & (strata == s)
        if b.any():
            total += a.sum() * (log_vr[a].mean() - log_vr[b].mean())
            weight += a.sum()
    return total / weight if weight else float("nan")


def _stats(x: np.ndarray, carrier: np.ndarray, noisy, directed, strata, slope_dir) -> tuple[float, float]:
    case, ctrl = x[carrier], x[~carrier]
    log_vr = np.log((case.var(axis=0, ddof=1) + 1e-8) / (ctrl.var(axis=0, ddof=1) + 1e-8))
    shift = case.mean(axis=0) - ctrl.mean(axis=0)
    return noise_concordance(log_vr, noisy, strata), float(stats.spearmanr(shift[directed], slope_dir).statistic)


def run(ctx: Context) -> list[Path]:
    if any(not matrix_path(ctx.paths, g).exists() for g in SERIES):
        print("syndrome_concordance: run `entropy-age fetch syndromes` first", file=sys.stderr)
        return []
    out = ctx.out_dir("syndrome_concordance")
    features, probes, _ = _features(ctx.paths, ctx.rng(140_000))
    classes = site_classes(_lifespan_metrics(_hannum(ctx.paths, features)))
    island = dict(probes.select("cpg", "cpg_island").iter_rows())
    print(f"syndrome_concordance: {classes.group_by('cls').len().to_dicts()}", file=sys.stderr, flush=True)
    rng = ctx.rng(390_006)
    rows = []
    for gse, arms in SERIES.items():
        people = read_people(matrix_path(ctx.paths, gse)).filter(pl.col("tissue") == "whole blood")
        feats = classes["feature"].to_list()
        people = people.filter(pl.col("sample").is_in(usable(matrix_path(ctx.paths, gse), people["sample"].to_list())))
        x_all = read_betas(matrix_path(ctx.paths, gse), feats, people["sample"].to_list())
        ok = ~np.isnan(x_all).any(axis=0)
        x_all = x_all[:, ok]
        c = classes.filter(pl.Series(ok))
        is_control = (people["disease"] == CONTROL).to_numpy()
        base = x_all[is_control].mean(axis=0)
        isl = np.array([bool(island.get(f, False)) for f in c["feature"].to_list()])
        strata = np.char.add(isl.astype(str), baseline_strata(base, 10).astype(str))
        noisy = (c["cls"] == "noisy").to_numpy()
        directed = (c["cls"] == "directed").to_numpy()
        slope_dir = c["slope"].to_numpy()[directed]
        for arm, gene in arms.items():
            is_arm = (people["disease"] == arm).to_numpy()
            if is_arm.sum() < MIN_ARM:
                continue
            keep = is_control | is_arm
            x = x_all[keep]
            carrier = is_arm[keep]
            obs = _stats(x, carrier, noisy, directed, strata, slope_dir)
            null = np.array(
                [
                    _stats(x, rng.permutation(carrier), noisy, directed, strata, slope_dir)
                    for _ in range(N_PERM[str(ctx.scale.name)])
                ]
            )
            rows.append(
                {
                    "series": gse,
                    "arm": arm,
                    "gene": gene,
                    "carriers": int(carrier.sum()),
                    "controls": int((~carrier).sum()),
                    "noise_concordance": obs[0],
                    "noise_null_mean": float(null[:, 0].mean()),
                    "noise_p": float(
                        (1 + (np.abs(null[:, 0] - null[:, 0].mean()) >= abs(obs[0] - null[:, 0].mean())).sum())
                        / (1 + len(null))
                    ),
                    "directed_rho": obs[1],
                    "directed_null_mean": float(null[:, 1].mean()),
                    "directed_p": float((1 + (np.abs(null[:, 1]) >= abs(obs[1])).sum()) / (1 + len(null))),
                }
            )
            print(f"syndrome_concordance: {gse} {arm}: {obs}", file=sys.stderr, flush=True)
    write_table(pl.DataFrame(rows), out / "concordance")
    return []
