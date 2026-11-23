"""Does the shield add to a validated clock? EPIC-Italy breast cancer (protocol P4 and change 2).

Women of EPIC-Italy (GSE51032, rebuilt), breast cancer diagnosed more than 2 years after the draw against women
without a cancer diagnosis, matched on age within 3 years (same chip first) as in ``cancer_matched``. Clocks are
scored with biolearn at the pinned commit (scripts/clocks_biolearn.py) and kept when at least 90 % of their weight is
measured in 90 % of arrays. Per clock, conditional logistic models on identical matched sets, every predictor per SD
among all women, adjusted for age, mean beta and 11 cell fractions:

* M0: clock acceleration (the clock's residual on age among all women);
* M1: M0 + B (drift at housekeeping look-alike promoters);
* M2: M1 + F (the housekeeping shield share).

Reported per clock: the clock's odds ratio in M0, B's in M1 and F's in M2 (bootstrap over matched sets), and the
likelihood-ratio tests M1 vs M0 and M2 vs M1. PhenoAge is the primary clock; every other kept clock is exploratory.
Disagreement: per woman, the SD across the kept age-type clocks of their standardised accelerations; it enters a
model with the panel's mean acceleration (exploratory), and its correlation with B and F is reported.
``entropy-age figures clock_increment`` (needs data/interim/clocks/GSE51032_clocks.csv).
"""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import polars as pl
from scipy import stats

from entropy_age.experiments.cancer_matched import N_BOOT, _score, clogit, match
from entropy_age.experiments.common import Context, write_table

PRIMARY = "PhenoAge"
LAG = 2.0
CLOCKS = Path("data/interim/clocks/GSE51032_clocks.csv")
B, F = "B housekeeping", "F housekeeping"


def z(v: np.ndarray) -> np.ndarray:
    return (v - np.nanmean(v)) / np.nanstd(v)


def acceleration(pred: np.ndarray, age: np.ndarray) -> np.ndarray:
    design = np.column_stack([np.ones(age.size), age])
    coef, *_ = np.linalg.lstsq(design, pred, rcond=None)
    return pred - design @ coef


def run(ctx: Context) -> list[Path]:
    src = ctx.paths.output / "incident_cancer_rsesame" / "per_person.parquet"
    clocks_path = ctx.paths.root / CLOCKS
    if not (src.exists() and clocks_path.exists()):
        print("clock_increment: needs incident_cancer_rsesame per_person and the clock panel", file=sys.stderr)
        return []
    out = ctx.out_dir("clock_increment")
    people = pl.read_parquet(src)
    panel = pl.read_csv(clocks_path).pivot(on="clock", index="sample", values="prediction")
    clocks = [c for c in panel.columns if c != "sample"]
    women = people.filter(pl.col("female")).join(panel, on="sample", how="inner")
    age_all = women["age"].to_numpy()
    for c in clocks:  # accelerations among all women, per SD
        women = women.with_columns(pl.Series(f"acc {c}", z(acceleration(women[c].to_numpy(), age_all))))
    age_type = [
        c for c in clocks if pl.read_csv(clocks_path).filter(pl.col("clock") == c)["output"][0].startswith("Age")
    ]
    acc = np.column_stack([women[f"acc {c}"].to_numpy() for c in age_type])
    women = women.with_columns(disagreement=pl.Series(z(acc.std(axis=1))), panel_mean=pl.Series(z(acc.mean(axis=1))))
    for name in (B, F):
        women = women.with_columns(pl.Series(name, z(women[name].to_numpy())))
    breast = pl.col("icd").fill_null("").str.starts_with("C50")
    cases = women.filter(breast & (pl.col("years_to_diagnosis") > LAG))
    controls = women.filter(pl.col("icd").is_null())
    pool = pl.concat([cases, controls], how="vertical")
    sets = match(cases, controls)
    adj_cols = [c for c in pool.columns if c == "mean_beta" or c.startswith("fraction_")]
    base = np.column_stack([z(pool["age"].to_numpy()), *(z(pool[c].to_numpy()) for c in adj_cols)])
    rng = ctx.rng(380_010)
    n_boot = N_BOOT[str(ctx.scale.name)]
    rows = []
    for c in [PRIMARY, *[k for k in clocks if k != PRIMARY]] + ["panel mean + disagreement"]:
        if c == "panel mean + disagreement":
            first = np.column_stack([pool["panel_mean"].to_numpy(), pool["disagreement"].to_numpy()])
        else:
            first = pool[f"acc {c}"].to_numpy()[:, None]
        m0 = np.column_stack([first, base])
        m1 = np.column_stack([first, pool[B].to_numpy(), base])
        m2 = np.column_stack([first, pool[B].to_numpy(), pool[F].to_numpy(), base])
        k = first.shape[1]
        b0, s0 = clogit(m0, sets)
        b1, s1 = clogit(m1, sets)
        b2, s2 = clogit(m2, sets)
        ll0, ll1, ll2 = (float(_score(m, sets, b)[0]) for m, b in ((m0, b0), (m1, b1), (m2, b2)))
        boot = {"clock": [], "B": [], "F": []}
        for _ in range(n_boot):
            pick = [sets[i] for i in rng.integers(0, len(sets), len(sets))]
            boot["clock"].append(clogit(m0, pick, steps=20)[0][0])
            boot["B"].append(clogit(m1, pick, steps=20)[0][k])
            boot["F"].append(clogit(m2, pick, steps=20)[0][k + 1])

        def ci(v):
            return (float(np.exp(np.quantile(v, 0.025))), float(np.exp(np.quantile(v, 0.975))))

        rows.append(
            {
                "clock": c,
                "role": "primary" if c == PRIMARY else "exploratory",
                "sets": len(sets),
                "clock_or": float(np.exp(b0[0])),
                "clock_low": ci(boot["clock"])[0],
                "clock_high": ci(boot["clock"])[1],
                "disagreement_or": float(np.exp(b0[1])) if k == 2 else float("nan"),
                "B_or_given_clock": float(np.exp(b1[k])),
                "B_low": ci(boot["B"])[0],
                "B_high": ci(boot["B"])[1],
                "lr_p_B": float(stats.chi2.sf(2 * (ll1 - ll0), 1)),
                "F_or_given_B": float(np.exp(b2[k + 1])),
                "F_low": ci(boot["F"])[0],
                "F_high": ci(boot["F"])[1],
                "lr_p_F": float(stats.chi2.sf(2 * (ll2 - ll1), 1)),
            }
        )
        print(f"clock_increment: {c} done", file=sys.stderr, flush=True)
    write_table(pl.DataFrame(rows), out / "nested_models")
    corr = [
        {
            "pair": f"disagreement ~ {m}",
            "r": float(np.corrcoef(women["disagreement"].to_numpy(), women[m].to_numpy())[0, 1]),
        }
        for m in (B, F)
    ]
    write_table(pl.DataFrame(corr), out / "disagreement_correlations")
    return []
