"""Drift toward each CpG's own equilibrium: is the shield slower gain, faster repair, or another resting level? (RM32)

Asked by antonkulaga (2026-11-13): methylation repair is not symmetric (DNMT1 maintenance, DNMT3A/B de novo gain,
TET-driven loss), so the binary-entropy slope, which assumes every site drifts toward 0.5, may misread the shield.
In the two-state Markov model of ``docs/model.md`` §2 a site gains methylation at rate ``k_gain`` (per unmethylated
copy and year) and loses it at rate ``k_loss``; its mean then moves as ``dm/dt = k_gain (1 - m) - k_loss m`` toward
its own equilibrium ``z = k_gain / (k_gain + k_loss)``, at speed ``λ = k_gain + k_loss``.

Two readings, on Hannum blood and Jaffe's adult cortex controls (the samples, covariates and CpGs of
``figures tissues``; covariates other than age removed first):

* **Pooled rates per gene set.** Over a few decades most sites are in the linear regime, where a site's slope is
  ``k_gain - (k_gain + k_loss) m``. Regressing the age slope (donors outside the youngest fifth) on the young level
  (the youngest fifth, other donors, so the two errors are independent) within a set and island status gives the
  set's mean ``k_gain`` (intercept) and ``k_loss`` (minus slope, minus intercept). Each set is compared with its
  look-alikes (promoter CpGs in no set, same island status). 95 % intervals resample donors.
* **Per-CpG relaxation.** ``m(t) = z + (m0 - z) e^{-λ t}`` fitted per CpG with ``z`` held in [0, 1]; where the curve
  bends enough to fix ``λ`` (not at the grid's ends), the equilibrium's distance from the young level, the
  fraction of that distance travelled over the cohort's ages, ``k_gain = zλ`` and ``k_loss = (1 - z)λ``; each set
  against look-alikes of the same island status and young-level decile, in SD.

``entropy-age figures equilibrium_drift`` as a Slurm job.
"""

from __future__ import annotations

import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import polars as pl

from entropy_age import plots
from entropy_age.data.chromatin import ChromatinLayout, cpg_states, read_segmentation
from entropy_age.experiments.common import Context, write_table
from entropy_age.experiments.tissue_contrast import SETS, Prepared, _features, _hannum, _jaffe, _matched_sd
from entropy_age.stochasticity import baseline_strata

YOUNG = 0.2
LOOKALIKE = "look-alikes"
N_BOOT = {"paper": 200, "quick": 20}
PERMUTATIONS = {"paper": 1000, "quick": 100}
RATES = np.geomspace(1e-4, 1.0, 64)  # per year
OWN_STATE = {"blood (Hannum)": "E062", "cortex (Jaffe)": "E073"}
ACTIVE = "1_TssA"
PER_CPG = ("headroom", "progress", "k_gain", "k_loss", "equilibrium")


def adjust(x: np.ndarray, ages: np.ndarray, covariates: np.ndarray | None) -> np.ndarray:
    """``x`` with the covariates' part removed and the age trend kept (fit on ``1 + age + covariates``)."""
    if covariates is None or covariates.size == 0:
        return x
    cov = covariates[:, covariates.std(axis=0) > 0]
    design = np.column_stack([np.ones(ages.size), ages, cov])
    coef, *_ = np.linalg.lstsq(design, x, rcond=None)
    centred = cov - cov.mean(axis=0)
    return x - centred @ coef[2:]


def pooled_rates(slope: np.ndarray, level: np.ndarray) -> dict[str, float]:
    """``k_gain``, ``k_loss`` and ``z`` from ``slope = k_gain - (k_gain + k_loss) level`` (OLS over CpGs)."""
    design = np.column_stack([np.ones(level.size), level])
    (a, b), *_ = np.linalg.lstsq(design, slope, rcond=None)
    speed = -b
    return {"k_gain": float(a), "k_loss": float(speed - a), "z": float(a / speed) if speed > 0 else float("nan")}


def bounded_relaxation(x: np.ndarray, ages: np.ndarray, rates: np.ndarray = RATES) -> dict[str, np.ndarray]:
    """Per column: ``m(t) = z + (m0 - z) exp(-λ (t - t_min))``, ``z`` held in [0, 1], λ profiled on ``rates``."""
    t = ages - ages.min()
    n = x.shape[1]
    best_sse = np.full(n, np.inf)
    best = np.zeros((3, n))
    for rate in rates:
        d = np.exp(-rate * t)
        design = np.column_stack([1 - d, d])
        coef, *_ = np.linalg.lstsq(design, x, rcond=None)
        z = np.clip(coef[0], 0, 1)
        m0 = np.where(z == coef[0], coef[1], (d @ (x - np.outer(1 - d, z))) / (d @ d))
        sse = ((x - np.outer(1 - d, z) - np.outer(d, m0)) ** 2).sum(axis=0)
        better = sse < best_sse
        best_sse = np.where(better, sse, best_sse)
        best[:, better] = np.vstack([np.full(n, rate), z, m0])[:, better]
    rate, z, m0 = best
    reach = 1 - np.exp(-rate * t.max())
    return {
        "rate": rate,
        "equilibrium": z,
        "start": m0,
        "headroom": np.abs(z - m0),
        "progress": reach,  # the share of the way from m0 to z travelled between the youngest and oldest donor
        "k_gain": z * rate * 10,
        "k_loss": (1 - z) * rate * 10,
        "identified": (rate > rates[0]) & (rate < rates[-1]),
    }


def _split(ages: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    young = ages <= np.quantile(ages, YOUNG)
    return young, ~young


def _slopes(x: np.ndarray, ages: np.ndarray) -> np.ndarray:
    design = np.column_stack([np.ones(ages.size), ages / 10])
    coef, *_ = np.linalg.lstsq(design, x, rcond=None)
    return coef[1]


def _labels(features: list[str], membership: pl.DataFrame) -> dict[str, np.ndarray]:
    member = set(membership["feature"].to_list())
    out = {LOOKALIKE: np.array([f not in member for f in features])}
    for g in SETS:
        s = set(membership.filter(pl.col("group") == g)["feature"].to_list())
        out[g] = np.array([f in s for f in features])
    return out


def _pooled_table(x, ages, island, labels) -> list[dict]:
    young, rest = _split(ages)
    level, slope = x[young].mean(axis=0), _slopes(x[rest], ages[rest])
    rows = []
    for g, lab in labels.items():
        for isl in (True, False):
            sel = lab & (island == isl)
            if sel.sum() >= 20:
                rows.append({"group": g, "island": isl, "cpgs": int(sel.sum()), **pooled_rates(slope[sel], level[sel])})
    return rows


def _boot_rows(x, ages, island, labels, n_boot, rng) -> pl.DataFrame:
    young, rest = np.flatnonzero(_split(ages)[0]), np.flatnonzero(_split(ages)[1])
    rows = []
    for b in range(n_boot):
        idx = np.concatenate([rng.choice(young, young.size), rng.choice(rest, rest.size)])
        xb, ab = x[idx], ages[idx]
        lb = np.zeros(ab.size, dtype=bool)
        lb[: young.size] = True
        level, slope = xb[lb].mean(axis=0), _slopes(xb[~lb], ab[~lb])
        for g, lab in labels.items():
            for isl in (True, False):
                sel = lab & (island == isl)
                if sel.sum() >= 20:
                    rows.append({"boot": b, "group": g, "island": isl, **pooled_rates(slope[sel], level[sel])})
    return pl.DataFrame(rows)


def _with_intervals(point: pl.DataFrame, boot: pl.DataFrame) -> pl.DataFrame:
    """Each set's rates, its difference from look-alikes of the same island status, and 95 % intervals."""
    look = boot.filter(pl.col("group") == LOOKALIKE).select(
        "boot", "island", *(pl.col(k).alias(f"{k}_la") for k in ("k_gain", "k_loss", "z"))
    )
    b = boot.join(look, on=["boot", "island"])
    b = b.with_columns(**{f"{k}_vs": pl.col(k) - pl.col(f"{k}_la") for k in ("k_gain", "k_loss", "z")})
    agg = []
    for k in ("k_gain", "k_loss", "z"):
        agg += [
            pl.col(k).quantile(0.025).alias(f"{k}_low"),
            pl.col(k).quantile(0.975).alias(f"{k}_high"),
            pl.col(f"{k}_vs").quantile(0.025).alias(f"{k}_vs_low"),
            pl.col(f"{k}_vs").quantile(0.975).alias(f"{k}_vs_high"),
        ]
    ci = b.group_by("group", "island").agg(agg)
    pla = point.filter(pl.col("group") == LOOKALIKE).select(
        "island", *(pl.col(k).alias(f"{k}_la") for k in ("k_gain", "k_loss", "z"))
    )
    out = point.join(pla, on="island").with_columns(
        **{f"{k}_vs": pl.col(k) - pl.col(f"{k}_la") for k in ("k_gain", "k_loss", "z")}
    )
    return out.drop([f"{k}_la" for k in ("k_gain", "k_loss", "z")]).join(ci, on=["group", "island"])


def _per_cpg(p: Prepared, x: np.ndarray, probes: pl.DataFrame, membership: pl.DataFrame, perms: int, seed: int, state):
    fit = bounded_relaxation(x, p.ages)
    young, _ = _split(p.ages)
    frame = pl.DataFrame({"feature": p.features, "baseline": x[young].mean(axis=0), "state": state, **fit})
    island = frame.join(probes.select(pl.col("cpg").alias("feature"), "cpg_island"), on="feature", how="left")
    frame = island.with_columns(
        stratum=pl.format(
            "{}|{}",
            pl.col("cpg_island").fill_null(False),
            pl.Series(baseline_strata(island["baseline"].to_numpy(), 10)),
        )
    )
    share = (
        frame.join(membership, on="feature", how="left")
        .with_columns(pl.col("group").fill_null(LOOKALIKE))
        .group_by("group")
        .agg(pl.col("identified").mean().alias("share_identified"), pl.len().alias("cpgs"))
    )
    ok = frame.filter(pl.col("identified"))
    gaps = _matched_sd(ok, membership, PER_CPG, perms, seed).with_columns(matching=pl.lit("island x level"))
    with_state = ok.with_columns(stratum=pl.format("{}|{}", pl.col("stratum"), pl.col("state")))
    gaps_state = _matched_sd(with_state, membership, PER_CPG, perms, seed).with_columns(
        matching=pl.lit("+ own chromatin state")
    )
    return frame, pl.concat([gaps, gaps_state]), share


def run(ctx: Context) -> list[Path]:
    plots.apply_style()
    out = ctx.out_dir("equilibrium_drift")
    features, probes, membership = _features(ctx.paths, ctx.rng(120_000))
    scale = str(ctx.scale.name)
    rng = ctx.rng(660_000)
    pooled, gaps, shares, per_cpg = [], [], [], []
    for p in (_hannum(ctx.paths, features), _jaffe(ctx.paths, features)):
        x = np.clip(adjust(p.x, p.ages.astype(float), p.covariates), 0, 1)
        ages = p.ages.astype(float)
        isl = dict(zip(probes["cpg"].to_list(), probes["cpg_island"].to_list(), strict=True))
        island = np.array([bool(isl.get(f, False)) for f in p.features])
        labels = _labels(p.features, membership)
        state = _states(ctx, p.features, probes, OWN_STATE[p.label])
        for matching, sel in (
            ("island x level", np.ones(island.size, bool)),
            (f"active TSS only ({ACTIVE})", state == ACTIVE),
        ):
            lab = {k: v[sel] for k, v in labels.items()}
            point = pl.DataFrame(_pooled_table(x[:, sel], ages, island[sel], lab))
            boot = _boot_rows(x[:, sel], ages, island[sel], lab, N_BOOT[scale], rng)
            pooled.append(
                _with_intervals(point, boot).with_columns(
                    tissue=pl.lit(p.label), donors=pl.lit(ages.size), matching=pl.lit(matching)
                )
            )
        frame, g, share = _per_cpg(p, x, probes, membership, PERMUTATIONS[scale], ctx.seed, state)
        gaps.append(g.with_columns(tissue=pl.lit(p.label)))
        shares.append(share.with_columns(tissue=pl.lit(p.label)))
        per_cpg.append(frame.drop("stratum").with_columns(tissue=pl.lit(p.label)))
        print(f"equilibrium_drift: {p.label}: {ages.size} donors, {len(p.features)} CpGs", file=sys.stderr, flush=True)
    pooled_t = pl.concat(pooled, how="diagonal")
    write_table(pooled_t, out / "pooled_rates")
    write_table(pl.concat(gaps), out / "per_cpg_gaps")
    write_table(pl.concat(shares), out / "identified")
    pl.concat(per_cpg).write_parquet(out / "per_cpg.parquet")
    return [_figure(ctx, pooled_t.filter(pl.col("matching") == "island x level"))]


def _states(ctx: Context, features: list[str], probes: pl.DataFrame, epigenome: str) -> np.ndarray:
    """Each CpG's 15-state chromatin segmentation in the tissue's own reference epigenome."""
    frame = pl.DataFrame({"cpg": features}).join(probes.select("cpg", "chr", "pos"), on="cpg", how="left")
    seg = read_segmentation(ChromatinLayout(ctx.paths).segmentation(epigenome))
    return cpg_states(frame, seg).fill_null("none").to_numpy().astype(str)


def _figure(ctx: Context, table: pl.DataFrame) -> Path:
    groups = [*SETS, LOOKALIKE]
    tissues = table["tissue"].unique(maintain_order=True).to_list()
    fig, axes = plots.panels(2 * len(tissues), ncols=4, size=3.0)
    k = 0
    for tissue in tissues:
        for rate, title in (("k_gain", "gain"), ("k_loss", "loss")):
            ax = axes[k]
            for j, isl in enumerate((True, False)):
                t = pl.DataFrame({"group": groups}).join(
                    table.filter((pl.col("tissue") == tissue) & (pl.col("island") == isl)), on="group", how="left"
                )
                y = np.arange(len(groups)) + 0.2 * (j - 0.5)
                v = t[rate].to_numpy()
                ax.errorbar(
                    v,
                    y,
                    xerr=[v - t[f"{rate}_low"].to_numpy(), t[f"{rate}_high"].to_numpy() - v],
                    fmt="o",
                    ms=4,
                    lw=1,
                    color=plots.PALETTE[j],
                    label="CpG island" if isl else "outside islands",
                )
            ax.axvline(0, color=plots.MUTED, lw=0.8)
            ax.xaxis.set_major_locator(plt.MaxNLocator(3))
            ax.set_yticks(
                range(len(groups)),
                [g.replace("_", " ") for g in groups] if k % 4 == 0 else [""] * len(groups),
                fontsize=7,
            )
            ax.invert_yaxis()
            ax.set_title(f"{tissue}: {title} rate", fontsize=8)
            ax.set_xlabel("per decade", fontsize=7)
            if k == 0:
                ax.legend(fontsize=6)
            plots.label_panel(ax, "abcdefgh"[k])
            k += 1
    fig.suptitle(
        "Methylation gain and loss rates by gene set (two-state Markov model)", x=0.02, y=1.05, ha="left", fontsize=11
    )
    fig.subplots_adjust(wspace=0.35)
    return plots.save(fig, ctx.figure_path("equilibrium_drift"))
