"""The breast-cancer lead as a matched case-control (roadmap RM38 item 4; a simulated referee's concern on round v0.1).

R46 found that drift at look-alike promoters (``B``) predicts breast cancer in EPIC-Italy (GSE51032) by logistic
regression over the whole cohort, with age and sex held. A referee asked for the design the cohort was sampled for:
cases set against controls like them. Here, among women only, each breast cancer case (ICD C50) gets up to
``K`` controls who stayed free of cancer, within ``AGE_CALIPER`` years of her age at the draw, from the same array
chip where possible (chips stand in for the recruitment centre and batch, which GEO does not give), each control
used once. Conditional logistic regression within the matched sets gives the odds ratio per SD of each measure,
with age as a covariate inside the sets, without and with the mean beta and 11 cell fractions of R46b; the 95 %
interval resamples matched sets. Run for all cases and for cases diagnosed over ``LAG`` years after the draw.

Reads ``incident_cancer``'s per-person table (with chip and adjustments); seconds. ``entropy-age figures
cancer_matched``.
"""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import polars as pl
from scipy import stats

from entropy_age.experiments.common import Context, source_suffix, write_table

K = 3
AGE_CALIPER = 3.0
LAG = 2.0
N_BOOT = {"paper": 1000, "quick": 100}
MEASURES = (
    "B housekeeping",
    "B core_essential",
    "F housekeeping",
    "S housekeeping",
    "DunedinPACE",
    "PhenoAge",
    "Horvath 2013",
)


def match(cases: pl.DataFrame, controls: pl.DataFrame, k: int = K, caliper: float = AGE_CALIPER) -> list[list[int]]:
    """Matched sets as row lists ``[case, control, ...]`` into ``cases.vstack(controls)``; controls used once.

    Controls are dealt in rounds, one per case per round (oldest case first, older controls being scarcer), so a
    short supply of controls spreads over all cases instead of filling the first ones; each takes the unused
    control within ``caliper`` years, same chip first, then nearest in age. Cases left without a control drop out.
    """
    n_case = cases.height
    c_age, c_chip = controls["age"].to_numpy(), controls["chip"].to_list()
    used = np.zeros(controls.height, bool)
    chosen: dict[int, list[int]] = {}
    order = [int(i) for i in np.argsort(-cases["age"].to_numpy(), kind="stable")]
    for _ in range(k):
        for i in order:
            age, chip = cases["age"][i], cases["chip"][i]
            d = np.abs(c_age - age)
            ok = np.flatnonzero(~used & (d <= caliper))
            if ok.size == 0:
                continue
            j = min(ok.tolist(), key=lambda j: (c_chip[j] != chip, d[j]))
            used[j] = True
            chosen.setdefault(i, []).append(j)
    return [[i, *(n_case + j for j in chosen[i])] for i in order if i in chosen]


def clogit(x: np.ndarray, sets: list[list[int]], *, steps: int = 50) -> tuple[np.ndarray, np.ndarray]:
    """Conditional logistic regression, one case first in each set; coefficients and SEs by Newton steps, each
    halved until the log-likelihood does not fall (an undamped step overshoots and diverges)."""
    beta = np.zeros(x.shape[1])
    ll, grad, info = _score(x, sets, beta)
    for _ in range(steps):
        step = np.linalg.lstsq(info, grad, rcond=None)[0]
        for _ in range(30):
            new_ll, new_grad, new_info = _score(x, sets, beta + step)
            if new_ll >= ll - 1e-12:
                break
            step = step / 2
        beta += step
        ll, grad, info = new_ll, new_grad, new_info
        if np.abs(step).max() < 1e-10:
            break
    return beta, np.sqrt(np.diag(np.linalg.pinv(info)))


def _score(x, sets, beta):
    ll = 0.0
    grad = np.zeros(x.shape[1])
    info = np.zeros((x.shape[1], x.shape[1]))
    for s in sets:
        xs = x[s]
        e = xs @ beta
        top = e.max()
        w = np.exp(e - top)
        total = w.sum()
        ll += e[0] - top - np.log(total)
        w /= total
        mean = w @ xs
        grad += xs[0] - mean
        centred = xs - mean
        info += centred.T @ (centred * w[:, None])
    return ll, grad, info


def _z(v: np.ndarray) -> np.ndarray:
    return (v - v.mean()) / v.std()


def effects(women: pl.DataFrame, lag: float | None, n_boot: int, rng) -> list[dict]:
    breast = pl.col("icd").fill_null("").str.starts_with("C50")
    cases = women.filter(breast)
    if lag is not None:
        cases = cases.filter(pl.col("years_to_diagnosis") > lag)
    controls = women.filter(pl.col("icd").is_null())
    pool = pl.concat([cases, controls], how="vertical")
    sets = match(cases, controls)
    same_chip = np.mean([pool["chip"][s[0]] == pool["chip"][j] for s in sets for j in s[1:]])
    adj_cols = [c for c in women.columns if c == "mean_beta" or c.startswith("fraction_")]
    age = _z(pool["age"].to_numpy())
    adj = np.column_stack([_z(pool[c].to_numpy()) for c in adj_cols]) if adj_cols else np.empty((pool.height, 0))
    rows = []
    for name in MEASURES:
        if name not in pool.columns:
            continue
        m = (pool[name].to_numpy() - women[name].mean()) / women[name].std()  # per SD among all women
        for adjusted, design in ((False, np.column_stack([m, age])), (True, np.column_stack([m, age, adj]))):
            if adjusted and not adj_cols:
                continue
            beta, se = clogit(design, sets)
            boot = []
            for _ in range(n_boot):
                pick = rng.integers(0, len(sets), len(sets))
                boot.append(clogit(design, [sets[i] for i in pick], steps=20)[0][0])
            lo, hi = np.quantile(boot, [0.025, 0.975]) if boot else (np.nan, np.nan)
            rows.append(
                {
                    "measure": name,
                    "cases": "all" if lag is None else f"diagnosed > {lag:g} y after draw",
                    "adjusted": "mean beta + 11 cell fractions" if adjusted else "age only",
                    "sets": len(sets),
                    "controls": sum(len(s) - 1 for s in sets),
                    "same_chip_share": float(same_chip),
                    "odds_ratio": float(np.exp(beta[0])),
                    "wald_low": float(np.exp(beta[0] - 1.96 * se[0])),
                    "wald_high": float(np.exp(beta[0] + 1.96 * se[0])),
                    "boot_low": float(np.exp(lo)),
                    "boot_high": float(np.exp(hi)),
                    "p": float(2 * stats.norm.sf(abs(beta[0] / se[0]))),
                }
            )
    return rows


def run(ctx: Context) -> list[Path]:
    src = ctx.paths.output / f"incident_cancer{source_suffix()}" / "per_person.parquet"
    if not src.exists() or "chip" not in pl.read_parquet_schema(src):
        print("cancer_matched: run `figures cancer` first (its per-person table with chip)", file=sys.stderr)
        return []
    out = ctx.out_dir("cancer_matched")
    women = pl.read_parquet(src).filter(pl.col("female"))
    rng = ctx.rng(380_004)
    n_boot = N_BOOT[str(ctx.scale.name)]
    rows = effects(women, None, n_boot, rng) + effects(women, LAG, n_boot, rng)
    table = pl.DataFrame(rows)
    write_table(table, out / "odds_ratios")
    print(f"cancer_matched: {table['sets'].max()} matched sets", file=sys.stderr)
    return []
