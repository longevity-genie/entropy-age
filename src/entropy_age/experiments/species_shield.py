"""Do long-lived species shield essential genes more? (roadmap RM10)

The Mammalian Methylation Consortium (GSE223748; ``entropy-age fetch mammals``) measured the same conserved
CpGs, on the mammalian methylation array, in 15,043 samples of 349 species and some 70 tissues, with each
animal's age and each species' maximum lifespan and adult weight. For every species and tissue with at least
``MIN_ANIMALS`` animals spanning a range of ages, the aging gap of each gene set is computed as in
``figures bxd``: entropy and dispersion slopes against age, set promoter probes against probes of the same
CpG-island status and young-reference level decile, in SDs of the metric, with an interval from resampling
the set's genes. The gap in SD does not depend on the unit of age, so species are comparable.

The gaps are then related to longevity across species: ``gap ~ log(max lifespan) + log(adult weight) +
tissue``, so that lifespan is read beyond what body size predicts, and a tissue's general level is held
fixed. Selection, chromatin or anything that protects essential genes more in long-lived species predicts a
negative lifespan coefficient (a stronger shield). Species are not independent (closely related ones share
both lifespan and biology); without a phylogeny this is a first look, and the result is reported per order
as well.

Run as a job: ``entropy-age fetch mammals`` (about 4.5 GB), then ``entropy-age figures species``.
"""

from __future__ import annotations

import gzip
import sys
from pathlib import Path

import numpy as np
import polars as pl

from entropy_age import plots
from entropy_age.data.annotation import AnnotationLayout, link_membership
from entropy_age.data.cohorts import _characteristics
from entropy_age.data.mouse import MouseAnnotationLayout, mouse_membership
from entropy_age.experiments.bxd_liver import ols, probes_without_suffix
from entropy_age.experiments.common import Context, chip_terms, dummies, human_annotation, write_table
from entropy_age.experiments.interventions import _sets
from entropy_age.experiments.tissue_contrast import SET_LABELS, SETS
from entropy_age.protection import matched_scores, resampled_gaps
from entropy_age.stochasticity import baseline_strata, drift_metrics

GSE = "GSE223748"
BETAS = f"{GSE}_datBetaNormalized.csv.gz"
MIN_ANIMALS = 15
MIN_AGE_SPAN = 0.15  # of the species' maximum lifespan
N_RESAMPLES = {"paper": 200, "quick": 40}


def layout(paths) -> tuple[Path, Path]:
    base = paths.input / "geo" / GSE
    return base / f"{GSE}_series_matrix.txt.gz", base / BETAS


def read_samples(matrix: Path) -> pl.DataFrame:
    """One row per sample: barcode, species, tissue, age, maximum lifespan and adult weight (numbers or null)."""
    chars = _characteristics(matrix)
    titles, sources, accessions = [], [], []
    with gzip.open(matrix, "rt") as handle:
        for line in handle:
            cells = [c.strip().strip('"') for c in line.rstrip("\n").split("\t")[1:]]
            if line.startswith("!Sample_title"):
                titles = cells
            elif line.startswith("!Sample_source_name_ch1"):
                sources = cells
            elif line.startswith("!Sample_geo_accession"):
                accessions = cells
    meta = pl.DataFrame({"sample": accessions, "title": titles, "tissue": sources})

    def num(col: str) -> pl.Expr:
        return pl.col(col).cast(pl.Float64, strict=False) if col in chars.columns else pl.lit(None, dtype=pl.Float64)

    species = pl.coalesce(*(pl.col(c) for c in ("speciescommonname",) if c in chars.columns))
    frame = chars.select(
        "sample",
        species.alias("species"),
        num("age").alias("age"),
        num("maximumlifespanyears").alias("max_lifespan"),
        num("average_adultweight").alias("adult_weight"),
        (pl.col("intervention") if "intervention" in chars.columns else pl.lit(None, dtype=pl.String)).alias(
            "intervention"
        ),
        (pl.col("experiment") if "experiment" in chars.columns else pl.lit(None, dtype=pl.String)).alias("experiment"),
    )
    return meta.join(frame, on="sample").with_columns(barcode=pl.col("title").str.extract(r"(\d+_R\d+C\d+)\s*$"))


MIN_RESIDUAL_DF = 10  # animals minus chips minus 2, below which a within-chip age slope is not fitted


def group_gap(
    x: np.ndarray, ages: np.ndarray, features: list[str], island: np.ndarray, sets: dict, n: int, rng, covariates=None
):
    """Per set and metric: gap in SD and its gene-resampling interval, strata from the youngest quarter.

    ``gap`` is the unscaled set-minus-matched difference of the metric's age slope (per year), and
    ``background_rate`` the matched look-alikes' own slope: how fast unprotected CpGs drift in this species
    (roadmap RM20: does a long life lower it toward the protected CpGs' floor?).
    """
    metrics = drift_metrics(x, ages, features, covariates=covariates, spearman=False)
    young = ages <= np.quantile(ages, 0.25)
    strata = np.char.add(island.astype(str), baseline_strata(x[young].mean(axis=0), 10).astype(str))
    out = []
    for metric in ("entropy_slope", "dispersion_slope", "kl_slope"):
        v = metrics[metric].to_numpy()
        sd = float(v.std())
        for group, (in_set, genes) in sets.items():
            gap, background = (float(z[0]) for z in matched_scores(v[None, :], in_set, strata))
            draws = resampled_gaps(v[None, :], genes, strata, n, rng)[:, 0]
            out.append(
                {
                    "metric": metric,
                    "group": group,
                    "gap": gap,
                    "background_rate": background,
                    "gap_sd": gap / sd,
                    "ci_low": float(np.quantile(draws, 0.025)) / sd,
                    "ci_high": float(np.quantile(draws, 0.975)) / sd,
                }
            )
    return out


def outcomes(sel: pl.DataFrame) -> dict[str, np.ndarray]:
    """What is regressed on lifespan: the gap in SD, the look-alikes' drift rate (log), and the gap relative to it.

    Catch-up (roadmap RM20) predicts a lower look-alike rate and a shrinking relative gap in long-lived species; a
    shield that removes a fixed share of drift predicts a lower rate and a constant relative gap.
    """
    rate = sel["background_rate"].to_numpy()
    with np.errstate(divide="ignore", invalid="ignore"):
        return {
            "gap_sd": sel["gap_sd"].to_numpy(),
            "log background rate": np.where(rate > 0, np.log(rate), np.nan),
            "relative gap": np.where(rate > 0, sel["gap"].to_numpy() / rate, np.nan),
        }


def _lifespan_fit(sel: pl.DataFrame, y: np.ndarray, metric: str, group: str, outcome: str) -> list[dict]:
    design = np.column_stack(
        [
            np.ones(y.size),
            np.log(sel["max_lifespan"].to_numpy()),
            np.log(sel["adult_weight"].to_numpy()),
            dummies(sel["tissue"].to_numpy()),
        ]
    )
    rows = []
    for term, col in (("log max lifespan", 1), ("log adult weight", 2)):
        coef, se, p = ols(y, design, col)
        rows.append(
            {
                "metric": metric,
                "group": group,
                "outcome": outcome,
                "term": term,
                "coef": coef,
                "se": se,
                "p": p,
                "groups": sel.height,
                "species": sel["species"].n_unique(),
            }
        )
    return rows


def run(ctx: Context) -> list[Path]:
    matrix, betas_path = layout(ctx.paths)
    if not betas_path.exists():
        print(f"species: {betas_path} missing; run `entropy-age fetch mammals` first. Nothing run.", file=sys.stderr)
        return []
    plots.apply_style()
    out = ctx.out_dir("species_shield")
    samples = read_samples(matrix).filter(
        pl.col("age").is_not_null()
        & pl.col("max_lifespan").is_not_null()
        & pl.col("barcode").is_not_null()
        & (pl.col("intervention").is_null() | pl.col("intervention").is_in(["WT", "Control", ""]))
    )
    lay = MouseAnnotationLayout(ctx.paths)
    probes = probes_without_suffix(
        pl.read_parquet(lay.probes).filter(pl.col("promoter") & ~pl.col("chr").is_in(["X", "Y", "MT"]))
    )
    human = human_annotation()
    if human:  # RM43 item 5: the human groups with promoters, islands and gene sets from the human genome
        samples = samples.filter(pl.col("species") == "Human")
        annotation = AnnotationLayout(ctx.paths)
        probes = pl.read_parquet(annotation.probes_450k).filter(pl.col("promoter") & ~pl.col("chr").is_in(["X", "Y"]))
    header = pl.read_csv(betas_path, n_rows=0).columns
    id_col = header[0]
    groups = (
        samples.filter(pl.col("barcode").is_in(header))
        .group_by("species", "tissue")
        .agg(
            pl.len().alias("n"),
            pl.col("age").min().alias("lo"),
            pl.col("age").max().alias("hi"),
            pl.col("max_lifespan").first(),
        )
        .filter((pl.col("n") >= MIN_ANIMALS) & ((pl.col("hi") - pl.col("lo")) >= MIN_AGE_SPAN * pl.col("max_lifespan")))
    )
    keep = samples.join(groups.select("species", "tissue"), on=["species", "tissue"]).filter(
        pl.col("barcode").is_in(header)
    )
    print(f"species: {groups.height} species-tissue groups, {keep.height} samples", file=sys.stderr, flush=True)
    frame = (
        pl.scan_csv(betas_path, null_values=["NA", ""])
        .select([id_col, *keep["barcode"].unique().to_list()])
        .filter(pl.col(id_col).is_in(probes["cpg"].implode()))
        .collect()
        .rename({id_col: "cpg"})
    )
    features = frame["cpg"].to_list()
    if human:
        members = link_membership(pl.read_parquet(annotation.links_450k), pl.read_parquet(annotation.gene_sets))
        members = members.filter(pl.col("group").is_in(SETS) & pl.col("feature").is_in(features))
    else:
        members = mouse_membership(probes.filter(pl.col("cpg").is_in(features)), pl.read_parquet(lay.gene_sets)).filter(
            pl.col("group").is_in(SETS)
        )
    island = dict(probes.select("cpg", "cpg_island").iter_rows())
    isl = np.array([bool(island.get(f, False)) for f in features])
    sets = _sets(features, members)
    print(
        f"species: {len(features)} promoter probes, set CpGs { ({g: int(s.sum()) for g, (s, _) in sets.items()}) }",
        file=sys.stderr,
        flush=True,
    )
    rng = ctx.rng(190_000)
    rows = []
    for (species, tissue), sel in keep.group_by("species", "tissue", maintain_order=True):
        x = frame.select(sel["barcode"].to_list()).to_numpy().T.astype(float)
        ok = ~np.isnan(x).any(axis=0)
        if ok.sum() < 0.9 * x.shape[1]:
            print(f"species: {species}, {tissue}: {int((~ok).sum())} probes missing; left out", file=sys.stderr)
            continue
        cols = np.flatnonzero(ok)
        sub_sets = {
            g: (s[cols], [np.searchsorted(cols, gi[np.isin(gi, cols)]) for gi in genes if np.isin(gi, cols).any()])
            for g, (s, genes) in sets.items()
        }
        cov = None
        if chip_terms():  # age within chips only (RM43 item 5): chip intercepts, where enough freedom is left
            chips = sel["barcode"].str.extract(r"^(\d+)_").to_numpy()
            cov = dummies(chips)
            cov = cov[:, np.ptp(cov, axis=0) > 0]
            if sel.height - cov.shape[1] - 2 < MIN_RESIDUAL_DF:
                print(
                    f"species: {species}, {tissue}: too few animals per chip for chip terms; left out", file=sys.stderr
                )
                continue
        gaps = group_gap(
            x[:, cols],
            sel["age"].to_numpy(),
            [features[c] for c in cols],
            isl[cols],
            sub_sets,
            N_RESAMPLES[str(ctx.scale.name)],
            rng,
            cov,
        )
        first = sel.row(0, named=True)
        rows += [
            {
                "species": species,
                "tissue": tissue,
                "animals": sel.height,
                "max_lifespan": first["max_lifespan"],
                "adult_weight": first["adult_weight"],
                "age_min": float(sel["age"].min()),
                "age_max": float(sel["age"].max()),
                **g,
            }
            for g in gaps
        ]
    gaps = pl.DataFrame(rows)
    write_table(gaps, out / "species_gaps")
    if human:  # one species: no lifespan fit and no figure
        return []

    fits = []
    for metric in ("entropy_slope", "dispersion_slope", "kl_slope"):
        for group in SETS:
            sel = gaps.filter((pl.col("metric") == metric) & (pl.col("group") == group)).drop_nulls(["adult_weight"])
            sel = sel.filter(pl.col("adult_weight") > 0)
            if sel.height < 10:
                continue
            for outcome, y_all in outcomes(sel).items():
                ok = np.isfinite(y_all)
                if ok.sum() < 10:
                    continue
                fits += _lifespan_fit(sel.filter(pl.Series(ok)), y_all[ok], metric, group, outcome)
    fits = pl.DataFrame(fits)
    write_table(fits, out / "lifespan_fits")
    return [_figure(ctx, gaps, fits)]


METRIC_NAMES = {"kl_slope": "relative entropy", "entropy_slope": "binary entropy", "dispersion_slope": "dispersion"}
#: species marked in panels a and b, each by its own marker rather than a text label per point
MARKED = {"Human": "*", "Naked mole rat": "D", "Mouse": "s"}


def _figure(ctx: Context, gaps: pl.DataFrame, fits: pl.DataFrame) -> Path:
    fig, axes = plots.panels(3, ncols=3, size=4.2)
    for i, metric in enumerate(("kl_slope", "dispersion_slope")):
        ax = axes[i]
        sel = gaps.filter((pl.col("metric") == metric) & (pl.col("group") == "housekeeping"))
        tissues = sel.group_by("tissue").len().sort("len", descending=True)["tissue"].to_list()[:6]
        for j, t in enumerate(tissues):
            sub = sel.filter(pl.col("tissue") == t)
            ax.scatter(
                np.log10(sub["max_lifespan"]),
                sub["gap_sd"],
                s=14,
                color=plots.PALETTE[j % len(plots.PALETTE)],
                label=t,
                alpha=0.8,
            )
        other = sel.filter(~pl.col("tissue").is_in(tissues))
        ax.scatter(
            np.log10(other["max_lifespan"]), other["gap_sd"], s=8, color=plots.MUTED, label="other tissues", alpha=0.5
        )
        for name, marker in MARKED.items():
            hit = sel.filter(pl.col("species").str.to_lowercase() == name.lower())
            ax.scatter(
                np.log10(hit["max_lifespan"]),
                hit["gap_sd"],
                s=40,
                marker=marker,
                facecolors="none",
                edgecolors="black",
                linewidths=0.8,
                label=name,
            )
        ax.axhline(0, color=plots.MUTED, lw=0.8)
        ax.set(
            xlabel="Maximum lifespan (years, log10)",
            ylabel="Housekeeping minus matched CpGs (SD)",
            title=f"{METRIC_NAMES[metric].capitalize()} slope",
        )
        ax.legend(fontsize=5, loc="lower left", ncol=2)
        plots.label_panel(ax, "ab"[i])
    ax = axes[2]
    sel = fits.filter((pl.col("term") == "log max lifespan") & (pl.col("outcome") == "gap_sd"))
    labels = []
    for k, r in enumerate(sel.iter_rows(named=True)):
        ax.errorbar(r["coef"], k, xerr=1.96 * r["se"], fmt="o", color=plots.PALETTE[SETS.index(r["group"])])
        labels.append(f"{SET_LABELS[r['group']]} · {METRIC_NAMES[r['metric']]}")
    ax.axvline(0, color=plots.MUTED, lw=0.8)
    ax.set_yticks(range(len(labels)), labels, fontsize=7)
    ax.invert_yaxis()
    ax.set(
        xlabel="Gap change per ln(max lifespan) (SD), 95 % CI",
        title="Gap against lifespan",
    )
    plots.label_panel(ax, "c")
    fig.suptitle(
        "Shielding across mammals (GSE223748, mammalian array)", x=0.02, ha="left", fontsize=11, fontweight="bold"
    )
    fig.tight_layout()
    return plots.save(fig, ctx.figure_path("species_shield"))
