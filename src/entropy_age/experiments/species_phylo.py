"""Shield against lifespan across mammals, with relatedness held fixed (roadmap RM10b).

Run R12 (``figures species``) regressed each species-tissue group's shield (``gap_sd``) on log maximum lifespan,
log adult weight and tissue, treating species as independent. Related species are not: whales share a shield
because they share ancestors as much as because they share long lives. This refits the same regression by
generalised least squares with the covariance

    V(λ) = λ C + (1 - λ) I,

where ``C`` is the shared NCBI Taxonomy lineage depth of the two groups' species (``entropy_age.data.taxonomy``)
and λ (Pagel's lambda) is chosen by maximum likelihood on a grid. λ = 0 is the independent fit of R12; λ near 1
says the gaps follow relatedness. Fits are made on all groups (tissue as covariate) and within each tissue with
at least ``MIN_SPECIES`` species.

Inputs: R12's ``species_gaps`` table, GSE223748's series matrix (each sample's organism) and the lineages
(``entropy-age fetch taxonomy``). Light: no methylation data is read.
"""

from __future__ import annotations

import gzip
import sys
from pathlib import Path

import numpy as np
import polars as pl
from scipy import stats

from entropy_age import plots
from entropy_age.data.taxonomy import lineage_covariance
from entropy_age.experiments.common import Context, dummies, write_table
from entropy_age.experiments.species_shield import layout, read_samples
from entropy_age.experiments.tissue_contrast import SET_LABELS, SETS

# λ stops short of 1: one species in several tissues shares its whole lineage, and V(1) would be singular
LAMBDAS = np.append(np.linspace(0.0, 0.975, 40), 0.99)
MIN_SPECIES = 10
TERMS = {"log max lifespan": 1, "log adult weight": 2}


def lineages_path(paths) -> Path:
    return paths.input / "taxonomy" / "lineages.csv"


def organisms(matrix: Path) -> pl.DataFrame:
    """Each consortium species' organism name as GEO gives it (``!Sample_organism_ch1``)."""
    accessions = names = None
    with gzip.open(matrix, "rt") as handle:
        for line in handle:
            cells = [c.strip().strip('"') for c in line.rstrip("\n").split("\t")[1:]]
            if line.startswith("!Sample_geo_accession"):
                accessions = cells
            elif line.startswith("!Sample_organism_ch1"):
                names = cells
            if accessions and names:
                break
    per_sample = read_samples(matrix).join(pl.DataFrame({"sample": accessions, "organism": names}), on="sample")
    out = per_sample.select("species", "organism").unique()
    clash = out.group_by("species").len().filter(pl.col("len") > 1)
    if clash.height:
        raise ValueError(f"species_phylo: species with several organisms: {clash['species'].to_list()}")
    return out


def gls(y: np.ndarray, x: np.ndarray, c: np.ndarray, lam: float) -> tuple[np.ndarray, np.ndarray, float, int]:
    """Coefficients, their covariance, the profile log-likelihood and the residual degrees of freedom."""
    v = lam * c + (1.0 - lam) * np.eye(y.size)
    chol = np.linalg.cholesky(v)
    xw, yw = np.linalg.solve(chol, x), np.linalg.solve(chol, y)
    coef, *_ = np.linalg.lstsq(xw, yw, rcond=None)
    resid = yw - xw @ coef
    n, rank = y.size, np.linalg.matrix_rank(xw)
    rss = float(resid @ resid)
    loglik = -0.5 * n * np.log(rss / n) - np.log(np.diag(chol)).sum()
    cov = rss / (n - rank) * np.linalg.pinv(xw.T @ xw)
    return coef, cov, float(loglik), n - rank


def fit_lambda(y: np.ndarray, x: np.ndarray, c: np.ndarray) -> float:
    """Pagel's λ on the grid ``LAMBDAS`` by maximum profile likelihood."""
    return float(LAMBDAS[int(np.argmax([gls(y, x, c, lam)[2] for lam in LAMBDAS]))])


def term_rows(y, x, c, lam: float, label: dict) -> list[dict]:
    coef, cov, _, dof = gls(y, x, c, lam)
    rows = []
    for term, col in TERMS.items():
        se = float(np.sqrt(cov[col, col]))
        t = coef[col] / se
        rows.append(
            {**label, "term": term, "coef": float(coef[col]), "se": se, "p": float(2 * stats.t.sf(abs(t), dof))}
        )
    return rows


def jackknife(y: np.ndarray, x: np.ndarray, species: list[str], col: int = 1) -> dict:
    """The coefficient with each species left out in turn: the species whose removal moves it most, and to what."""
    names = sorted(set(species))
    if len(names) < 3:
        raise ValueError("species_phylo: a jackknife needs at least three species")
    full = float(np.linalg.lstsq(x, y, rcond=None)[0][col])
    moved = {}
    for name in names:
        keep = np.array([s != name for s in species])
        moved[name] = float(np.linalg.lstsq(x[keep], y[keep], rcond=None)[0][col])
    worst = max(moved, key=lambda n: abs(moved[n] - full))
    return {
        "influential_species": worst,
        "coef_without": moved[worst],
        "coef_min": min(moved.values()),
        "coef_max": max(moved.values()),
    }


def run(ctx: Context) -> list[Path]:
    matrix, _ = layout(ctx.paths)
    gaps_path = ctx.paths.output / "species_shield" / "species_gaps.csv"
    missing = [str(p) for p in (matrix, gaps_path, lineages_path(ctx.paths)) if not p.exists()]
    if missing:
        print(
            f"species_phylo: inputs missing ({', '.join(missing)}); run `entropy-age fetch mammals`, "
            "`entropy-age figures species` and `entropy-age fetch taxonomy` first. Nothing run.",
            file=sys.stderr,
        )
        return []
    plots.apply_style()
    out = ctx.out_dir("species_phylo")
    names = organisms(matrix)
    lineage = {
        r["organism"]: r["lineage"].split(";") for r in pl.read_csv(lineages_path(ctx.paths)).iter_rows(named=True)
    }
    gaps = (
        pl.read_csv(gaps_path)
        .join(names, on="species", how="left")
        .drop_nulls(["adult_weight"])
        .filter(pl.col("adult_weight") > 0)
    )
    unknown = sorted(set(gaps["organism"].to_list()) - set(lineage))
    if unknown:
        raise ValueError(f"species_phylo: no lineage for {unknown}; run `entropy-age fetch taxonomy`")
    rows = []
    for metric in ("entropy_slope", "dispersion_slope"):
        for group in SETS:
            sel = gaps.filter((pl.col("metric") == metric) & (pl.col("group") == group))
            scopes = {"all tissues": sel} | {
                t: s for (t,), s in sel.group_by("tissue") if s["species"].n_unique() >= MIN_SPECIES
            }
            for scope, s in scopes.items():
                y = s["gap_sd"].to_numpy()
                cols = [np.ones(y.size), np.log(s["max_lifespan"].to_numpy()), np.log(s["adult_weight"].to_numpy())]
                if scope == "all tissues":
                    cols.append(dummies(s["tissue"].to_numpy()))
                x = np.column_stack(cols)
                c = lineage_covariance([lineage[o] for o in s["organism"].to_list()])
                lam = fit_lambda(y, x, c)
                label = {
                    "metric": metric,
                    "group": group,
                    "scope": scope,
                    "groups": s.height,
                    "species": s["species"].n_unique(),
                }
                jack = jackknife(y, x, s["species"].to_list())
                rows += term_rows(y, x, c, 0.0, label | {"model": "independent", "lambda": 0.0} | jack)
                rows += term_rows(y, x, c, lam, label | {"model": "phylogenetic", "lambda": lam})
    fits = pl.DataFrame(rows)
    write_table(fits, out / "lifespan_fits")
    return [_figure(ctx, fits)]


def _figure(ctx: Context, fits: pl.DataFrame) -> Path:
    sel = fits.filter(pl.col("term") == "log max lifespan")
    fig, axes = plots.panels(2, ncols=2, size=4.6)
    for i, metric in enumerate(("entropy_slope", "dispersion_slope")):
        ax = axes[i]
        sub = sel.filter(pl.col("metric") == metric)
        scopes = ["all tissues", *sorted(set(sub["scope"].to_list()) - {"all tissues"})]
        labels, y = [], 0.0
        for scope in scopes:
            for j, group in enumerate(SETS):
                for k, model in enumerate(("independent", "phylogenetic")):
                    r = sub.filter((pl.col("scope") == scope) & (pl.col("group") == group) & (pl.col("model") == model))
                    if r.is_empty():
                        continue
                    r = r.row(0, named=True)
                    yy = y + 0.3 * k
                    ax.errorbar(
                        r["coef"],
                        yy,
                        xerr=1.96 * r["se"],
                        fmt="o" if k == 0 else "s",
                        ms=4,
                        color=plots.PALETTE[j],
                        alpha=0.5 if k == 0 else 1.0,
                    )
                    if k == 1:
                        ax.annotate(
                            f"λ={r['lambda']:.2f}",
                            (1.0, yy),
                            xycoords=("axes fraction", "data"),
                            fontsize=6,
                            va="center",
                            ha="right",
                        )
                labels.append((y + 0.15, f"{scope} ({r['species']} spp) · {SET_LABELS[group]}"))
                y += 1
            y += 0.6
        ax.axvline(0, color=plots.MUTED, lw=0.8)
        ax.set_yticks([p for p, _ in labels], [t for _, t in labels] if i == 0 else [""] * len(labels), fontsize=6.5)
        ax.invert_yaxis()
        ax.set(
            xlabel="Gap per log maximum lifespan (SD, 95 % CI)",
            title={"entropy_slope": "Entropy gap", "dispersion_slope": "Dispersion gap"}[metric],
        )
        plots.label_panel(ax, "ab"[i])
    fig.suptitle(
        "Shield against lifespan across mammals (> 0: weaker in long-lived species): circles species independent, "
        "squares with relatedness (λ)",
        x=0.02,
        y=1.04,
        ha="left",
        fontsize=11,
        fontweight="bold",
    )
    return plots.save(fig, ctx.figure_path("species_phylo"))
