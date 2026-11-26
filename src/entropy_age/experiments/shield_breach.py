"""Does damage start when drift breaches the shield at the genes that define the cell? (asked by antonkulaga)

Entropy accumulates with age and housekeeping and identity promoters are shielded from it. Two tests on series
already rebuilt from IDATs (``config/series``, arrays chosen by ``series.selection``; betas from
``data/interim/clocks/series``; promoter CpGs, gene sets and matched background as ``tissue_contrast._features`` and
``per_person._sets``):

1. **Threshold.** Each person's background drift ``B`` and shield share ``F`` at housekeeping look-alike promoters
   are scored against a reference fitted on everyone else in the group (``health_states.loo_scores``, as
   ``clock_drift`` uses them). In every group of at least ``MIN_DONORS`` donors, ``F`` is fitted on ``B`` with one
   breakpoint (``piecewise_fit``: continuous, two slopes, the breakpoint chosen over ``GRID`` quantiles of ``B`` by
   least squares, which is the Gaussian profile likelihood) against a straight line, with and without age held. The
   line is the null of a residual bootstrap of the likelihood ratio (``line_null_p``; the breakpoint is not
   identified under the null, so the chi-square reference does not hold); the slope change gets a donor-bootstrap
   interval with the breakpoint re-chosen in each resample. In the six disease series (``DISEASE``) the fit is made on
   disease and control donors together, and the share of each beyond the breakpoint is reported; ``F`` and ``B`` are
   compared between disease and control at a given age by OLS with age, sex and the series' other held columns, and
   again with chip intercepts where the chip is recorded and disease stays estimable beside it.
2. **Identity genes.** Each tissue's identity CpGs are the promoter CpGs (TSS1500, TSS200, 5'UTR, first exon;
   autosomes) of a short curated list of marker genes highly expressed in and specific to the tissue in GTEx and the
   Human Protein Atlas (``IDENTITY``). In disease and control donors separately, per-CpG age slopes of relative
   entropy from the group's own young reference (cross-fitted, ``muscle_identity._kl_slopes``) with sex and held
   columns, and the identity CpGs' gap against promoter background CpGs matched on island status and young-level
   decile, each gene weighted once, in SD of the slopes (``muscle_identity.weighted_matched_gap``). The disease
   minus control difference gets an interval from independent donor resamples of the two groups.

Limits: cross-sectional; a disease may shift ``F`` or ``B`` as a consequence, and treatment or post-mortem state
cannot be separated from it; the leave-one-out reference mixes cases with controls, which pulls contrasts toward
zero; each group's young reference is its own, so a disease group with fewer young donors has a noisier reference;
liver fibrosis is split at stage 2 (0-1 against 2-3); bulk tissues carry composition shifts that move both tests.

Writes ``shield_breach/{threshold,placement,disease_scores,identity_gap}`` and one figure.
Run as ``shield_breach.run(ctx)``.
"""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import polars as pl

from entropy_age import plots
from entropy_age.data import reprocessed
from entropy_age.data.annotation import AnnotationLayout
from entropy_age.experiments.bxd_liver import ols
from entropy_age.experiments.common import Context, dummies, write_table
from entropy_age.experiments.health_states import loo_scores
from entropy_age.experiments.muscle_identity import PROMOTER_REGIONS, _kl_slopes, strata_of, weighted_matched_gap
from entropy_age.experiments.per_person import _islands, _sets
from entropy_age.experiments.series import SENSITIVITY, age_identifiable, selection, specs
from entropy_age.experiments.tissue_contrast import _features

SERIES = Path("data/interim/clocks/series")
GROUP = "housekeeping"
MIN_DONORS = 80
MIN_SUBGROUP = 20
GRID = np.linspace(0.1, 0.9, 33)  # quantiles of B tried as the breakpoint
N_BOOT = {"paper": 300, "quick": 20}
LEVEL_BINS = 10

#: curated tissue identity genes (highly expressed and tissue-specific in GTEx / Human Protein Atlas)
IDENTITY = {
    "lung": ("SFTPC", "SFTPB", "SFTPA1", "NKX2-1", "AGER", "HOPX", "SCGB1A1", "FOXJ1"),
    "retina": ("RHO", "RCVRN", "NRL", "CRX", "PDE6B", "RPE65"),
    "liver": ("ALB", "APOA1", "TF", "CYP3A4", "HNF4A", "SERPINA1"),
    "airway": ("SCGB1A1", "FOXJ1", "KRT5", "MUC5AC"),
    "colon": ("CDX2", "VIL1", "CDH17", "LGALS4"),
    "brain": ("RBFOX3", "SNAP25", "SYT1", "GAD1", "SLC17A7"),
}

#: series -> the combined group (filters), the column telling disease from control, its values, the tissue's genes
DISEASE = {
    "lung_smoking": {
        "label": "lung, healthy and pulmonary fibrosis",
        "filters": {"tissue": "lung", "dx": ["normal", "ipf"]},
        "column": "dx",
        "control": ["normal"],
        "case": ["ipf"],
        "disease": "pulmonary fibrosis",
        "tissue": "lung",
    },
    "retina_amd": {
        "label": "retina, every donor",
        "column": "mgs grade",
        "control": ["1"],
        "case": ["3", "4"],
        "disease": "AMD (MGS 3-4)",
        "tissue": "retina",
    },
    "liver_aatd": {
        "label": "liver (AATD), every donor",
        "column": "fib",
        "control": ["0", "1"],
        "case": ["2", "3"],
        "disease": "fibrosis stage 2-3 (vs 0-1)",
        "tissue": "liver",
    },
    "airway_asthma": {
        "label": "airway, every donor",
        "column": "disease status",
        "control": ["Control"],
        "case": ["Asthma"],
        "disease": "asthma",
        "tissue": "airway",
    },
    "brain_stg": {
        "label": "superior temporal gyrus, every donor",
        "column": "group2",
        "control": ["Control"],
        "case": ["Schizophrenia"],
        "disease": "schizophrenia",
        "tissue": "brain",
    },
    "neurons_parietal": {
        "label": "neurons (parietal), every donor",
        "column": "disease state",
        "control": ["Control"],
        "case": ["PD"],
        "disease": "Parkinson's disease",
        "tissue": "brain",
    },
}


# ---------------------------------------------------------------- pure helpers


def _design(b: np.ndarray, cov: np.ndarray | None, tau: float | None) -> np.ndarray:
    cols = [np.ones(b.size), b]
    if tau is not None:
        cols.append(np.maximum(0.0, b - tau))
    if cov is not None and cov.shape[1]:
        cols.append(cov)
    return np.column_stack(cols)


def _sse(y: np.ndarray, design: np.ndarray) -> float:
    coef, *_ = np.linalg.lstsq(design, y, rcond=None)
    r = y - design @ coef
    return float(r @ r)


def piecewise_fit(b: np.ndarray, f: np.ndarray, cov: np.ndarray | None = None, grid=GRID) -> dict:
    """Continuous one-breakpoint fit ``f ~ b + max(0, b - tau) [+ cov]`` with ``tau`` chosen over ``grid`` quantiles of
    ``b`` by least squares (the profile likelihood), against the straight line ``f ~ b [+ cov]``: the breakpoint, the
    slopes before and after, the slope change, the likelihood ratio ``n log(SSE_line / SSE)`` and the BIC gain of the
    bend (two extra parameters)."""
    b, f = np.asarray(b, float), np.asarray(f, float)
    n = b.size
    line = _design(b, cov, None)
    coef_line, *_ = np.linalg.lstsq(line, f, rcond=None)
    sse_line = _sse(f, line)
    best = None
    for tau in np.unique(np.quantile(b, grid)):
        design = _design(b, cov, tau)
        coef, *_ = np.linalg.lstsq(design, f, rcond=None)
        r = f - design @ coef
        sse = float(r @ r)
        if best is None or sse < best["sse"]:
            best = {
                "tau": float(tau),
                "slope_before": float(coef[1]),
                "slope_change": float(coef[2]),
                "slope_after": float(coef[1] + coef[2]),
                "sse": sse,
            }
    best["slope_line"] = float(coef_line[1])
    best["sse_line"] = sse_line
    best["lr"] = float(n * np.log(sse_line / best["sse"])) if best["sse"] > 0 else float("inf")
    best["delta_bic"] = best["lr"] - 2 * float(np.log(n))
    return best


def line_null_p(b: np.ndarray, f: np.ndarray, cov: np.ndarray | None, n: int, rng: np.random.Generator) -> float:
    """p of the bend's likelihood ratio under the straight line: residual bootstrap from the line fit, the breakpoint
    re-chosen in every draw (add-one estimate)."""
    observed = piecewise_fit(b, f, cov)["lr"]
    line = _design(b, cov, None)
    coef, *_ = np.linalg.lstsq(line, f, rcond=None)
    fitted, resid = line @ coef, f - line @ coef
    hits = sum(piecewise_fit(b, fitted + rng.choice(resid, resid.size), cov)["lr"] >= observed for _ in range(n))
    return (hits + 1) / (n + 1)


def donor_bootstrap(
    b: np.ndarray, f: np.ndarray, cov: np.ndarray | None, n: int, rng: np.random.Generator
) -> dict[str, np.ndarray]:
    """Breakpoint and slope change over ``n`` resamples of donors (rows), the breakpoint re-chosen in each."""
    taus, changes = [], []
    for _ in range(n):
        rows = rng.integers(0, b.size, b.size)
        c = cov[rows] if cov is not None else None
        if c is not None:
            c = c[:, np.ptp(c, axis=0) > 1e-12]
        fit = piecewise_fit(b[rows], f[rows], c)
        taus.append(fit["tau"])
        changes.append(fit["slope_change"])
    return {"tau": np.array(taus), "slope_change": np.array(changes)}


def interval(draws: np.ndarray) -> tuple[float, float]:
    """2.5th and 97.5th percentiles of the finite draws (NaN when fewer than 10)."""
    d = np.asarray(draws, float)
    d = d[np.isfinite(d)]
    if d.size < 10:
        return float("nan"), float("nan")
    return float(np.quantile(d, 0.025)), float(np.quantile(d, 0.975))


def identity_gap(
    x: np.ndarray, ages: np.ndarray, cov: np.ndarray, member: np.ndarray, weights: np.ndarray, island: np.ndarray
) -> float:
    """Identity CpGs' age-slope gap of cross-fitted relative entropy against CpGs matched on island status and
    young-level decile, each gene weighted once, in SD of all slopes."""
    young = ages <= np.quantile(ages, 0.2)
    strata = strata_of(island, x[young].mean(axis=0), LEVEL_BINS)
    slopes = _kl_slopes(x, ages, cov)
    sd = float(np.nanstd(slopes))
    return weighted_matched_gap(slopes, member, np.where(member, weights, 0.0), strata) / sd


def identity_gap_draws(x, ages, cov, member, weights, island, n: int, rng) -> np.ndarray:
    """``identity_gap`` over donor resamples, with the strata and SD of the full sample (as
    ``series.donor_intervals``); resamples leaving age unidentifiable are NaN."""
    young = ages <= np.quantile(ages, 0.2)
    strata = strata_of(island, x[young].mean(axis=0), LEVEL_BINS)
    sd = float(np.nanstd(_kl_slopes(x, ages, cov)))
    w = np.where(member, weights, 0.0)
    out = np.full(n, np.nan)
    for k in range(n):
        rows = rng.integers(0, ages.size, ages.size)
        c = cov[rows]
        c = c[:, np.ptp(c, axis=0) > 1e-12] if c.shape[1] else c
        if not age_identifiable(ages[rows], c, distinct=np.unique(rows).size):
            continue
        out[k] = weighted_matched_gap(_kl_slopes(x[rows], ages[rows], c, rows), member, w, strata) / sd
    return out


# ---------------------------------------------------------------- data


def _held(s: pl.DataFrame, columns: list[str]) -> np.ndarray:
    parts = [dummies(s["sex"].cast(pl.String).fill_null("NA").to_numpy())]
    parts += [dummies(s[c].cast(pl.String).fill_null("NA").to_numpy()) for c in columns]
    m = np.column_stack(parts) if parts else np.zeros((s.height, 0))
    return m[:, np.ptp(m, axis=0) > 1e-12] if m.shape[1] else m


def _jobs() -> list[tuple[str, dict, dict | None]]:
    """(series, spec, disease entry or None); disease series get their combined group added."""
    out = []
    for name, spec in specs().items():
        if spec.get("annotation", "human") != "human" or SENSITIVITY.search(name):
            continue
        d = DISEASE.get(name)
        if d is not None:
            spec = dict(spec)
            groups = dict(spec["groups"])
            groups[d["label"]] = d.get("filters", groups.get(d["label"]))
            spec["groups"] = groups
            spec["hold"] = sorted({*spec.get("hold", []), d["column"]})
        out.append((name, spec, d))
    return out


def run(ctx: Context) -> list[Path]:
    n_boot = N_BOOT[str(ctx.scale.name)]
    rng = ctx.rng(190_000)
    features, _, membership = _features(ctx.paths, ctx.rng(120_000))  # the stream clock_drift uses
    annotation = AnnotationLayout(ctx.paths)
    probes = pl.read_parquet(annotation.probes_450k).filter(~pl.col("chr").is_in(["X", "Y"]))
    links = pl.read_parquet(annotation.links_450k).filter(
        pl.col("region").is_in(list(PROMOTER_REGIONS)) & pl.col("cpg").is_in(probes["cpg"].implode())
    )
    identity = {
        t: links.filter(pl.col("symbol").is_in(list(g))).select("cpg", "symbol").unique(subset="cpg")
        for t, g in IDENTITY.items()
    }
    id_cpgs = set(pl.concat(identity.values())["cpg"].to_list())
    feature_set = set(features)
    wanted = feature_set | id_cpgs
    threshold, placement, disease_rows, gaps, curves = [], [], [], [], {}
    for name, spec, dis in _jobs():
        gse = spec["accession"]
        path = ctx.paths.root / SERIES / f"{gse}_cpgs.parquet"
        if not path.exists():
            continue
        columns = set(pl.read_parquet_schema(path).names()) - {"cpg"}
        chosen = selection(ctx.paths, name, spec, columns)
        if chosen is None:
            continue
        groups, dropped = chosen
        big = {k: s.drop_nulls("age") for k, s in groups.items()}
        big = {k: s for k, s in big.items() if s.height >= MIN_DONORS or (dis and k == dis["label"])}
        if not big:
            continue
        samples = sorted({v for s in big.values() for v in s["sample"].to_list()})
        betas = (
            pl.scan_parquet(path).select("cpg", *samples).filter(pl.col("cpg").is_in(list(wanted - dropped))).collect()
        )
        chip_of = dict(reprocessed.read_qc(ctx.paths, gse).select("sample", "chip").iter_rows())
        for label, s in big.items():
            sub = betas.select("cpg", *s["sample"].to_list())
            x = sub.drop("cpg").to_numpy().T.astype(float)
            ok = ~np.isnan(x).any(axis=0)
            cpgs = [c for c, k in zip(sub["cpg"].to_list(), ok, strict=True) if k]
            x = x[:, ok]
            in_feat = np.array([c in feature_set for c in cpgs])
            fcpgs = [c for c, k in zip(cpgs, in_feat, strict=True) if k]
            sets = _sets(fcpgs, membership)
            if GROUP not in sets:
                continue
            print(f"shield_breach: {name} / {label}: {s.height} donors, {len(fcpgs)} CpGs", file=sys.stderr, flush=True)
            per = loo_scores(x[:, in_feat], _islands(fcpgs, probes), sets)
            b, f = per[f"B_{GROUP}"], per[f"F_{GROUP}"]
            ages = s["age"].to_numpy().astype(float)
            if s.height >= MIN_DONORS:
                for held, cov in (("none", None), ("age", ages[:, None])):
                    fit = piecewise_fit(b, f, cov)
                    boot = donor_bootstrap(b, f, cov, n_boot, rng)
                    lo, hi = interval(boot["slope_change"])
                    tlo, thi = interval(boot["tau"])
                    threshold.append(
                        {
                            "series": name,
                            "group": label,
                            "donors": s.height,
                            "held": held,
                            **{k: fit[k] for k in ("tau", "slope_line", "slope_before", "slope_after")},
                            "slope_change": fit["slope_change"],
                            "slope_change_low": lo,
                            "slope_change_high": hi,
                            "tau_low": tlo,
                            "tau_high": thi,
                            "lr": fit["lr"],
                            "delta_bic": fit["delta_bic"],
                            "p_line": line_null_p(b, f, cov, n_boot, rng),
                        }
                    )
            if not (dis and label == dis["label"]):
                continue
            status = s[dis["column"]].cast(pl.String).to_numpy()
            case = np.isin(status, dis["case"])
            keep = case | np.isin(status, dis["control"])
            fit = piecewise_fit(b[keep], f[keep])
            curves[name] = (b[keep], f[keep], case[keep], fit, dis)
            for who, m in (("control", keep & ~case), ("disease", case)):
                placement.append(
                    {
                        "series": name,
                        "disease": dis["disease"],
                        "status": who,
                        "donors": int(m.sum()),
                        "tau": fit["tau"],
                        "slope_change": fit["slope_change"],
                        "share_beyond_tau": float((b[m] > fit["tau"]).mean()) if m.any() else float("nan"),
                        "median_B": float(np.median(b[m])) if m.any() else float("nan"),
                        "median_F": float(np.median(f[m])) if m.any() else float("nan"),
                    }
                )
            # disease against control at a given age
            sk = s.filter(pl.Series(keep))
            other = [c for c in spec.get("hold", []) if c != dis["column"]]
            held = _held(sk, other)
            chips = np.array([str(chip_of.get(v) or "NA") for v in sk["sample"].to_list()])
            y = case[keep].astype(float)
            designs = {"age, sex, held": held}
            if (chips != "NA").all():
                designs["age, sex, held, chip"] = np.column_stack([held, dummies(chips)])
            for dname, c in designs.items():
                design = np.column_stack([np.ones(y.size), y, ages[keep], c])
                if np.linalg.matrix_rank(design) < design.shape[1] or y.size - design.shape[1] < 10:
                    continue
                for measure, v in (("F", f[keep]), ("B", b[keep])):
                    coef, se, p = ols(v, design, 1)
                    disease_rows.append(
                        {
                            "series": name,
                            "disease": dis["disease"],
                            "measure": measure,
                            "covariates": dname,
                            "controls": int((~case[keep]).sum()),
                            "cases": int(case[keep].sum()),
                            "effect": coef,
                            "se": se,
                            "p": p,
                            "effect_sd": coef / float(np.std(v, ddof=1)),
                        }
                    )
            # identity genes: disease and control separately
            table = identity[dis["tissue"]]
            gene_of = dict(table.iter_rows())
            member = np.array([c in gene_of for c in cpgs])
            per_gene = table.filter(pl.col("cpg").is_in(cpgs))["symbol"].value_counts()
            n_of = dict(per_gene.iter_rows())
            weights = np.array([1.0 / n_of[gene_of[c]] if c in gene_of else 0.0 for c in cpgs])
            # background: promoter CpGs of the matched feature set, never identity CpGs
            use = member | in_feat
            island = _islands(cpgs, probes)
            draws, point = {}, {}
            for who, m in (("control", keep & ~case), ("disease", case)):
                if m.sum() < MIN_SUBGROUP:
                    continue
                xs, a = x[m][:, use], ages[m]
                cov = _held(s.filter(pl.Series(m)), other)
                point[who] = identity_gap(xs, a, cov, member[use], weights[use], island[use])
                draws[who] = identity_gap_draws(xs, a, cov, member[use], weights[use], island[use], n_boot, rng)
                lo, hi = interval(draws[who])
                gaps.append(
                    {
                        "series": name,
                        "disease": dis["disease"],
                        "tissue": dis["tissue"],
                        "status": who,
                        "donors": int(m.sum()),
                        "identity_genes": int(len(n_of)),
                        "identity_cpgs": int(member.sum()),
                        "gap_sd": point[who],
                        "ci_low": lo,
                        "ci_high": hi,
                    }
                )
            if len(point) == 2:
                lo, hi = interval(draws["disease"] - draws["control"])
                gaps.append(
                    {
                        "series": name,
                        "disease": dis["disease"],
                        "tissue": dis["tissue"],
                        "status": "disease minus control",
                        "donors": int(keep.sum()),
                        "identity_genes": int(len(n_of)),
                        "identity_cpgs": int(member.sum()),
                        "gap_sd": point["disease"] - point["control"],
                        "ci_low": lo,
                        "ci_high": hi,
                    }
                )
    out = ctx.out_dir("shield_breach")
    for table_name, rows in (
        ("threshold", threshold),
        ("placement", placement),
        ("disease_scores", disease_rows),
        ("identity_gap", gaps),
    ):
        if rows:
            write_table(pl.DataFrame(rows), out / table_name)
    if not curves:
        return []
    return [_figure(ctx, curves, pl.DataFrame(gaps) if gaps else None)]


def _figure(ctx: Context, curves: dict, gaps: pl.DataFrame | None) -> Path:
    plots.apply_style()
    fig, axes = plots.plt.subplots(2, 4, figsize=(15, 7))
    axes = axes.ravel()
    for ax, (name, (b, f, case, fit, dis)) in zip(axes, curves.items(), strict=False):
        ax.scatter(b[~case], f[~case], s=8, alpha=0.6, color=plots.PALETTE[0], label="control")
        ax.scatter(b[case], f[case], s=8, alpha=0.6, color=plots.PALETTE[1], label=dis["disease"])
        grid = np.linspace(b.min(), b.max(), 100)
        design = _design(b, None, fit["tau"])
        coef, *_ = np.linalg.lstsq(design, f, rcond=None)
        ax.plot(grid, _design(grid, None, fit["tau"]) @ coef, color=plots.INK, lw=1.5)
        ax.axvline(fit["tau"], color=plots.MUTED, ls="--", lw=0.8)
        ax.set_title(f"{name}\nslope change {fit['slope_change']:.2g}, dBIC {fit['delta_bic']:.1f}", fontsize=8)
        ax.set_xlabel("background drift B", fontsize=8)
        ax.set_ylabel("shield share F", fontsize=8)
        ax.legend(fontsize=6, frameon=False)
    ax = axes[-1]
    for extra in axes[len(curves) : -1]:
        extra.axis("off")
    if gaps is not None:
        d = gaps.filter(pl.col("status") == "disease minus control")
        y = np.arange(d.height)
        g, lo, hi = (d[c].to_numpy() for c in ("gap_sd", "ci_low", "ci_high"))
        ax.errorbar(
            g,
            y,
            xerr=np.nan_to_num(np.abs(np.vstack([g - lo, hi - g]))),
            fmt="o",
            color=plots.INK,
        )
        ax.set_yticks(y, d["series"].to_list(), fontsize=7)
        ax.axvline(0, color=plots.MUTED, ls="--", lw=0.8)
        ax.set_xlabel("identity-gene KL slope gap,\ndisease minus control (SD)", fontsize=8)
        ax.set_title("identity genes", fontsize=8)
    fig.tight_layout()
    return plots.save(fig, ctx.figure_path("shield_breach"))
