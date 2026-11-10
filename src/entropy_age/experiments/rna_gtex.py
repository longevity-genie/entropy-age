"""Prediction 2 in RNA: GTEx v10, 47 tissues, along a turnover gradient (roadmap RM4).

The methylation tests found the essential-gene protection in every mouse tissue whose cells divide,
even slowly, and missing in cortex and hippocampus; heart and white blood cells were protected about
equally. GTEx RNA has thirteen brain regions next to blood, gut, skin, heart and muscle in the same
donors, so it asks the same question in human transcripts: is between-donor spread of essential and
housekeeping genes growing more slowly with age than that of genes of the same expression level, and
is that protection missing in brain?

**Inputs.** GTEx v10 counts read from an external parsed copy (``data.gtex_rna``; settings
``ENTROPY_AGE_GTEX_RNA`` and ``ENTROPY_AGE_GTEX_ATTRIBUTES``). Ages are decade bins; the midpoint is
used (25 ... 75). Tissues with fewer than ``MIN_DONORS`` donors (bladder, cervix, fallopian tube,
kidney medulla and the laser-captured sub-samples) and cultured fibroblasts are left out; the
``tissues`` table lists every tissue with the reason.

**Measure.** Per tissue: protein-coding genes (HGNC) with at least 1 CPM in 80 % of the tissue's
samples, as log2(CPM + 1). For each gene the age slope of its absolute residual after
``log CPM ~ age + sex + Hardy death scale + RIN + ischemic time + log2 library size``
(``stochasticity.expression_drift``): how fast donors spread apart at that gene. A second version adds
five composition components of the genes outside every set, made orthogonal to age (the
``tissues_cells`` check), which absorb the cell mix. Each gene set is compared with background genes of
the same expression decile (``compare_groups``), in SDs of the metric. Gene sets: the Hart core-essential
and non-essential sets and the Eisenberg-Levanon housekeeping genes on Ensembl IDs through HGNC
(``annotation.gene_set_ensembl``). Non-essential genes are by definition silent in dividing cell lines,
so few of them pass the expression filter outside testis and brain; their numbers are in the output.

**Size matching.** Fewer donors make per-gene slopes noisier, which inflates the metric's SD and pulls
every SD-scaled difference toward zero. The intervals therefore come from subsamples of ``SUBSAMPLE``
donors per tissue, drawn without replacement, and the class contrasts are taken over those resamples.

**Turnover classes**, decided from the tissue list:

* *renewing* - whole blood, spleen (lymphocytes), transverse colon and terminal ileum (both sampled
  with mucosa), stomach (full thickness, with mucosa), oesophagus mucosa, both skin sites, vagina
  (squamous epithelium), testis (spermatogonia);
* *slow* - liver, kidney cortex, lung, thyroid, pancreas, adrenal, pituitary, both adipose depots,
  breast, prostate, ovary, uterus, minor salivary gland, three arteries, tibial nerve (Schwann cells,
  fibroblasts), sigmoid colon, oesophagus muscularis and gastro-oesophageal junction (all three mostly
  smooth muscle, whose cells divide rarely but do divide);
* *heart & muscle* - both heart samples and skeletal muscle: post-mitotic, not neural;
* *brain* - the thirteen brain regions, reported apart from heart and muscle because the mouse result
  pointed at neurons. Spinal cord and substantia nigra are mostly glia and white matter, cerebellum
  mostly granule neurons; all are post-mitotic for their neurons, not for their glia.
"""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import polars as pl
from scipy.stats import false_discovery_control

from entropy_age import plots
from entropy_age.data.annotation import AnnotationLayout, gene_set_ensembl
from entropy_age.data.gtex_rna import GtexRnaMissing, GtexRnaSource, read_counts, read_genes, read_samples
from entropy_age.experiments.common import Context, dummies, write_table
from entropy_age.experiments.human_tissues import _subsample
from entropy_age.experiments.tissue_contrast import (
    _STATE,
    SET_LABELS,
    SETS,
    Prepared,
    _adjusted,
    _interval,
    _matched_sd,
    _run_parallel,
    _settings,
)
from entropy_age.stochasticity import baseline_strata, expression_drift

RENEWING = (
    "Whole Blood",
    "Spleen",
    "Colon - Transverse",
    "Small Intestine - Terminal Ileum",
    "Stomach",
    "Esophagus - Mucosa",
    "Skin - Not Sun Exposed (Suprapubic)",
    "Skin - Sun Exposed (Lower leg)",
    "Vagina",
    "Testis",
)
SLOW = (
    "Liver",
    "Kidney - Cortex",
    "Lung",
    "Thyroid",
    "Pancreas",
    "Adrenal Gland",
    "Pituitary",
    "Adipose - Subcutaneous",
    "Adipose - Visceral (Omentum)",
    "Breast - Mammary Tissue",
    "Prostate",
    "Ovary",
    "Uterus",
    "Minor Salivary Gland",
    "Artery - Aorta",
    "Artery - Coronary",
    "Artery - Tibial",
    "Nerve - Tibial",
    "Colon - Sigmoid",
    "Esophagus - Muscularis",
    "Esophagus - Gastroesophageal Junction",
)
HEART_MUSCLE = ("Heart - Atrial Appendage", "Heart - Left Ventricle", "Muscle - Skeletal")
BRAIN = (
    "Brain - Amygdala",
    "Brain - Anterior cingulate cortex (BA24)",
    "Brain - Caudate (basal ganglia)",
    "Brain - Cerebellar Hemisphere",
    "Brain - Cerebellum",
    "Brain - Cortex",
    "Brain - Frontal Cortex (BA9)",
    "Brain - Hippocampus",
    "Brain - Hypothalamus",
    "Brain - Nucleus accumbens (basal ganglia)",
    "Brain - Putamen (basal ganglia)",
    "Brain - Spinal cord (cervical c-1)",
    "Brain - Substantia nigra",
)
#: tissue -> turnover class, in plotting order
TURNOVER: dict[str, str] = (
    dict.fromkeys(RENEWING, "renewing")
    | dict.fromkeys(SLOW, "slow")
    | dict.fromkeys(HEART_MUSCLE, "heart & muscle")
    | dict.fromkeys(BRAIN, "brain")
)
CLASSES = ("renewing", "slow", "heart & muscle", "brain")
#: class contrasts reported, first minus second (negative = more protection in the first)
CONTRASTS = (
    ("renewing", "brain"),
    ("slow", "brain"),
    ("heart & muscle", "brain"),
    ("renewing", "heart & muscle"),
    ("slow", "heart & muscle"),
    ("renewing", "slow"),
)
ADJUSTMENTS: dict[str, tuple[int, bool] | None] = {"covariates": None, "pc5_age_free": (5, True)}
METRICS = ("dispersion_slope", "abs_slope")
MIN_DONORS = 100
SUBSAMPLE = 70  # donors per tissue and resample: 70 % of the smallest tissue kept (kidney cortex, 104)
MIN_CPM = 1.0
MIN_FRACTION = 0.8
STRATA = 10


# ----------------------------------------------------------------------------- inputs
def log_cpm(counts: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """log2(CPM + 1) of a samples x genes count matrix, and each sample's library size."""
    counts = np.asarray(counts, dtype=float)
    library = counts.sum(axis=1)
    if np.any(library <= 0):
        raise ValueError("a sample has no counts")
    return np.log2(counts / library[:, None] * 1e6 + 1.0), library


def expressed(log_values: np.ndarray, min_cpm: float = MIN_CPM, fraction: float = MIN_FRACTION) -> np.ndarray:
    """Genes at ``min_cpm`` or more in at least ``fraction`` of samples."""
    return (log_values >= np.log2(min_cpm + 1.0)).mean(axis=0) >= fraction


def technical_covariates(samples: pl.DataFrame, library: np.ndarray) -> tuple[np.ndarray, dict]:
    """Sex and Hardy-scale dummies, RIN, ischemic time and log2 library size, with what was imputed.

    A missing RIN or ischemic time is set to the tissue's median and counted; a missing Hardy scale is
    its own level. The counts go into the ``tissues`` table.
    """
    numeric, imputed = [], {}
    for column in ("rin", "ischemic_minutes"):
        values = samples[column].to_numpy().astype(float)
        missing = np.isnan(values)
        imputed[f"imputed_{column}"] = int(missing.sum())
        if missing.all():
            raise ValueError(f"no {column} for any sample of {samples['tissue'][0]}")
        values[missing] = np.nanmedian(values)
        numeric.append(values)
    hardy = samples["hardy"].fill_null("NA").to_numpy()
    imputed["hardy_missing"] = int((hardy == "NA").sum())
    columns = [dummies(samples["sex"].fill_null("NA").to_numpy()), dummies(hardy), *(v[:, None] for v in numeric)]
    columns.append(np.log2(library)[:, None])
    return np.column_stack(columns), imputed


def prepare_tissue(
    tissue: str, samples: pl.DataFrame, counts: np.ndarray, genes: list[str], protein_coding: set[str]
) -> tuple[Prepared, dict]:
    """One tissue, its expressed protein-coding genes as float32 log CPM, ready for ``_tissue_rows``."""
    values, library = log_cpm(counts)
    coding = np.array([g in protein_coding for g in genes])
    keep = coding & expressed(values)
    technical, info = technical_covariates(samples, library)
    features = [g for g, k in zip(genes, keep, strict=True) if k]
    p = Prepared(
        tissue,
        values[:, keep].astype(np.float32),
        samples["age"].to_numpy().astype(float),
        technical,
        features,
        samples["donor"].to_numpy(),
        technical,
    )
    return p, info | {"genes": len(features)}


def _load(ctx: Context, source: GtexRnaSource, members: set[str], protein_coding: set[str]):
    samples = read_samples(source)
    genes = read_genes(source)
    donors = samples.group_by("tissue").agg(pl.col("donor").n_unique().alias("donors"))
    tissues, table = {}, []
    for tissue, n in donors.sort("tissue").iter_rows():
        row = {"tissue": tissue, "turnover": TURNOVER.get(tissue), "donors": n}
        if tissue not in TURNOVER:
            reason = "cultured cells" if tissue.startswith("Cells") else f"fewer than {MIN_DONORS} donors"
            table.append(row | {"used": False, "reason": reason})
            continue
        if n < MIN_DONORS:
            raise ValueError(f"{tissue} has {n} donors, fewer than MIN_DONORS={MIN_DONORS}: reclassify it")
        t = samples.filter(pl.col("tissue") == tissue).unique(subset="donor", keep="first").sort("donor")
        p, info = prepare_tissue(tissue, t, read_counts(source, t["sample"].to_list()), genes, protein_coding)
        for name, adjustment in ADJUSTMENTS.items():
            tissues[(tissue, name)] = p if adjustment is None else _adjusted(p, members, *adjustment)[0]
        ages = p.ages
        table.append(
            row
            | {"used": True, "reason": None, "age_min_bin_mid": float(ages.min()), "age_max_bin_mid": float(ages.max())}
            | {"age_mean": float(ages.mean())}
            | info
        )
        print(f"[rna_gtex] {tissue}: {n} donors, {info['genes']} genes", file=sys.stderr, flush=True)
    return tissues, pl.DataFrame(table, infer_schema_length=None)


# ----------------------------------------------------------------------------- metrics
def tissue_metrics(p: Prepared, rows: np.ndarray) -> pl.DataFrame:
    cov = p.covariates[rows]
    cov = cov[:, cov.std(axis=0) > 0]  # a subsample can lose a level (or a single-sex tissue has none)
    metrics = expression_drift(p.x[rows].astype(float), p.ages[rows], p.features, covariates=cov)
    return metrics.with_columns(stratum=pl.Series(baseline_strata(metrics["baseline"].to_numpy(), STRATA)))


def _tissue_rows(tissue, adjustment, p, rows, membership, permutations, seed) -> pl.DataFrame:
    sets = membership.filter(pl.col("feature").is_in(p.features))
    table = _matched_sd(tissue_metrics(p, rows), sets, METRICS, permutations, seed)
    return table.with_columns(tissue=pl.lit(tissue), turnover=pl.lit(TURNOVER[tissue]), adjustment=pl.lit(adjustment))


def _boot(seed: int) -> list[dict]:
    """One resample: ``SUBSAMPLE`` donors of every tissue, without replacement, the same in both adjustments."""
    tissues, membership, size = _STATE["tissues"], _STATE["membership"], _STATE["size"]
    rng = np.random.default_rng(seed)
    out = []
    for tissue in TURNOVER:
        if (tissue, "covariates") not in tissues:
            continue
        rows = _subsample(rng, tissues[(tissue, "covariates")].x.shape[0], size)
        for adjustment in ADJUSTMENTS:
            table = _tissue_rows(tissue, adjustment, tissues[(tissue, adjustment)], rows, membership, 0, seed)
            out += [
                {"boot": seed, **r}
                for r in table.select("tissue", "adjustment", "metric", "group", "matched_diff_sd").to_dicts()
            ]
    return out


# ----------------------------------------------------------------------------- summaries (shared with rna_rodents)
def class_contrasts(boot: pl.DataFrame, point: pl.DataFrame, turnover: dict[str, str], contrasts) -> pl.DataFrame:
    """Class contrasts over size-matched resamples (median, 95 % interval, full-data value for reference)."""
    keys = ("adjustment", "metric", "group")
    boot = boot.with_columns(turnover=pl.col("tissue").replace_strict(turnover))
    point = point.with_columns(turnover=pl.col("tissue").replace_strict(turnover))
    rows = []
    for key, sel in boot.group_by(*keys):
        wide = (
            sel.group_by("boot", "turnover")
            .agg(pl.col("matched_diff_sd").mean())
            .pivot(on="turnover", index="boot", values="matched_diff_sd")
            .sort("boot")
        )
        full = point.filter(*(pl.col(k) == v for k, v in zip(keys, key, strict=True)))
        full_class = dict(full.group_by("turnover").agg(pl.col("matched_diff_sd").mean()).iter_rows())
        for first, second in contrasts:
            if first not in wide.columns or second not in wide.columns:
                continue
            diff = (wide[first] - wide[second]).drop_nulls().drop_nans().to_numpy()
            if diff.size < 10:
                continue
            lo, hi, p = _interval(diff)
            rows.append(
                dict(zip(keys, key, strict=True))
                | {
                    "contrast": f"{first} - {second}",
                    "difference": float(np.median(diff)),
                    "ci_low": lo,
                    "ci_high": hi,
                    "p_bootstrap": p,
                    "full_data": full_class[first] - full_class[second],
                }
            )
    return pl.DataFrame(rows).sort("adjustment", "metric", "group", "contrast") if rows else pl.DataFrame()


def tissue_intervals(boot: pl.DataFrame, turnover: dict[str, str]) -> pl.DataFrame:
    keys = ("adjustment", "metric", "group")
    return (
        boot.with_columns(turnover=pl.col("tissue").replace_strict(turnover))
        .group_by("tissue", "turnover", *keys)
        .agg(
            pl.col("matched_diff_sd").median().alias("median"),
            pl.col("matched_diff_sd").quantile(0.025).alias("ci_low"),
            pl.col("matched_diff_sd").quantile(0.975).alias("ci_high"),
            pl.len().alias("resamples"),
        )
        .sort("adjustment", "metric", "group", "tissue")
    )


def summarise(boot: pl.DataFrame, point: pl.DataFrame) -> tuple[pl.DataFrame, pl.DataFrame]:
    """Class contrasts over the size-matched resamples, and each tissue's resample interval.

    A class value is the mean of its tissues' matched differences in one resample; ``difference`` is
    the median over resamples (negative = more protection in the first class), ``full_data`` the same
    contrast of the full-data estimates, for reference only.
    """
    return class_contrasts(boot, point, TURNOVER, CONTRASTS), tissue_intervals(boot, TURNOVER)


# ----------------------------------------------------------------------------- run
def _inputs(ctx: Context) -> tuple[GtexRnaSource, AnnotationLayout] | None:
    annotation = AnnotationLayout(ctx.paths)
    try:
        source = GtexRnaSource.from_environment()
    except GtexRnaMissing as error:
        print(f"rna_gtex: {error}. Nothing was run.", file=sys.stderr)
        return None
    missing = [str(p) for p in (annotation.gene_sets, annotation.hgnc) if not p.exists()]
    if missing:
        print(
            f"rna_gtex: inputs missing ({', '.join(missing)}); run `entropy-age fetch annotation` first. "
            "Nothing was run.",
            file=sys.stderr,
        )
        return None
    return source, annotation


def run(ctx: Context) -> list[Path]:
    inputs = _inputs(ctx)
    if inputs is None:
        return []
    source, annotation = inputs
    plots.apply_style()
    out = ctx.out_dir("rna_gtex")
    workers, n_boot = _settings(ctx)
    hgnc = pl.read_parquet(annotation.hgnc)
    membership = gene_set_ensembl(pl.read_parquet(annotation.gene_sets), hgnc).filter(pl.col("group").is_in(SETS))
    protein_coding = set(hgnc.filter(pl.col("protein_coding") & pl.col("ensembl").is_not_null())["ensembl"].to_list())
    members = set(membership["feature"].to_list())
    tissues, table = _load(ctx, source, members, protein_coding)
    sizes = [
        {
            "tissue": t,
            "group": g,
            "genes": int(membership.filter(pl.col("group") == g)["feature"].is_in(p.features).sum()),
        }
        for (t, adjustment), p in tissues.items()
        if adjustment == "covariates"
        for g in SETS
    ]
    write_table(
        table.join(pl.DataFrame(sizes).pivot(on="group", index="tissue", values="genes"), on="tissue", how="left"),
        out / "tissues",
    )
    write_table(
        membership.group_by("group", "matched_by").agg(pl.len().alias("genes")).sort("group", "matched_by"),
        out / "set_mapping",
    )

    point = pl.concat(
        [
            _tissue_rows(tissue, adjustment, p, np.arange(p.x.shape[0]), membership, 1000, 3)
            for (tissue, adjustment), p in tissues.items()
        ],
        how="diagonal",
    )
    point = point.with_columns(p_bh=pl.Series(false_discovery_control(point["p_permutation"].to_numpy(), method="bh")))
    write_table(point, out / "point_estimates")
    print(f"[rna_gtex] point estimates done; {n_boot} resamples on {workers} workers", file=sys.stderr, flush=True)

    state = {"tissues": tissues, "membership": membership, "size": SUBSAMPLE}
    boot = pl.DataFrame(_run_parallel(_boot, [9000 + b for b in range(n_boot)], workers, state))
    boot.write_parquet(out / "bootstrap.parquet")
    contrasts, per_tissue = summarise(boot, point)
    write_table(contrasts, out / "contrasts")
    write_table(per_tissue, out / "tissue_intervals")
    return [_figure(ctx, per_tissue, contrasts)]


# ----------------------------------------------------------------------------- figure
SHOWN = ("core_essential", "housekeeping")  # non-essential genes are often a handful per tissue: tables only
_ADJ = {"covariates": "cov", "pc5_age_free": "+PCs"}
_SET_SHORT = {"core_essential": "core", "housekeeping": "HK", "non_essential": "non-ess."}


def _short(tissue: str) -> str:
    return tissue.replace("Brain - ", "").replace(" (basal ganglia)", "").replace("Esophagus", "Oesoph.")


def tissue_panel(ax, sel: pl.DataFrame, turnover: dict[str, str], classes, title: str, labels: bool) -> None:
    """Per-tissue resample intervals of the shown sets, tissues grouped by class from top to bottom."""
    present = set(sel["tissue"].to_list())
    tissues = [t for t in turnover if t in present]
    for j, group in enumerate(SHOWN):
        lookup = {r["tissue"]: r for r in sel.filter(pl.col("group") == group).iter_rows(named=True)}
        for i, tissue in enumerate(tissues):
            if tissue in lookup:
                r, y = lookup[tissue], i + (j - 0.5) * 0.32
                ax.plot([r["ci_low"], r["ci_high"]], [y, y], color=plots.PALETTE[j], lw=0.9)
                ax.scatter([r["median"]], [y], color=plots.PALETTE[j], s=9, label=SET_LABELS[group] if i == 0 else None)
    start = 0
    for cls in classes:
        n = sum(turnover[t] == cls for t in tissues)
        if n == 0:
            continue
        if start:
            ax.axhline(start - 0.5, color=plots.MUTED, lw=0.6, ls=":")
        ax.text(
            0.98,
            start - 0.3,
            cls,
            transform=ax.get_yaxis_transform(),
            ha="right",
            va="top",
            fontsize=7,
            color=plots.MUTED,
            style="italic",
        )
        start += n
    ax.axvline(0, color=plots.MUTED, lw=0.8)
    ax.set_yticks(range(len(tissues)), [_short(t) for t in tissues] if labels else [], fontsize=6.5)
    ax.set_ylim(len(tissues) - 0.5, -0.5)
    ax.set(xlabel="Set minus matched genes (SD)", title=title)


def contrast_panel(ax, frame: pl.DataFrame, title: str) -> None:
    """Class contrasts of the shown sets; accent where the interval excludes zero."""
    frame = (
        frame.filter(pl.col("group").is_in(SHOWN))
        .with_columns(
            label=pl.format(
                "{} · {} · {}",
                pl.col("group").replace_strict(_SET_SHORT),
                pl.col("contrast").str.replace(" - ", " − "),
                pl.col("adjustment").replace_strict(_ADJ),
            )
        )
        .sort("group", "contrast", "adjustment")
    )
    for y, r in enumerate(frame.iter_rows(named=True)):
        color = plots.ACCENT if r["ci_high"] < 0 or r["ci_low"] > 0 else plots.MUTED
        ax.plot([r["ci_low"], r["ci_high"]], [y, y], color=color, lw=1.1)
        ax.scatter([r["difference"]], [y], color=color, s=12)
    ax.axvline(0, color=plots.MUTED, lw=0.8)
    ax.set_yticks(range(frame.height), frame["label"].to_list(), fontsize=6)
    ax.set_ylim(max(frame.height - 0.5, 0.5), -0.5)
    ax.set(xlabel="Difference (SD); < 0 = more protection in the first", title=title)


def _figure(ctx: Context, per_tissue: pl.DataFrame, contrasts: pl.DataFrame) -> Path:
    import matplotlib.pyplot as plt

    spread = per_tissue.filter(pl.col("metric") == "dispersion_slope")
    n = spread["tissue"].n_unique()
    fig = plt.figure(figsize=(16, 0.2 * n + 3))
    grid = fig.add_gridspec(1, 3, width_ratios=(1, 0.8, 1), wspace=0.5)
    a, b, c = (fig.add_subplot(grid[i]) for i in range(3))
    tissue_panel(
        a,
        spread.filter(pl.col("adjustment") == "covariates"),
        TURNOVER,
        CLASSES,
        f"Spread gain with age, {SUBSAMPLE} donors/tissue",
        True,
    )
    tissue_panel(
        b,
        spread.filter(pl.col("adjustment") == "pc5_age_free"),
        TURNOVER,
        CLASSES,
        "Same, plus 5 age-free components",
        False,
    )
    contrast_panel(c, contrasts.filter(pl.col("metric") == "dispersion_slope"), "Turnover class contrasts")
    for ax, letter in zip((a, b, c), "abc", strict=True):
        plots.label_panel(ax, letter)
    a.legend(fontsize=7, loc="lower left")
    fig.suptitle(
        f"Do essential genes resist transcriptional drift less in brain? GTEx v10 RNA, {n} tissues",
        x=0.02,
        ha="left",
        fontsize=12,
        fontweight="bold",
    )
    fig.subplots_adjust(left=0.12, right=0.99, top=0.94, bottom=0.06)
    return plots.save(fig, ctx.figure_path("rna_gtex"))
