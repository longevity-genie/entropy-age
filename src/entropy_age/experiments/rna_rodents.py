"""The RNA turnover test in mouse and rat: GSE314817 and GSE314970 (roadmap RM4).

A transcriptomic atlas of aging (C57BL/6J mice, 24 tissues; F344 rats, 28 tissues), male and female
animals at six ages from 6 to 27 months, about 100 mice and 47-85 rats per tissue. Unlike GTEx the
ages are exact and the animals inbred, so between-animal spread has no genetic part, and a test in
rodents repeats the human RNA result (``rna_gtex``) without post-mortem confounders.

**Measure**, as in ``rna_gtex``: per tissue, protein-coding genes with at least 1 CPM in 80 % of the
tissue's samples as log2(CPM + 1); the age slope of the absolute residual after
``log CPM ~ age + sex + log2 library size`` (``expression_drift``), and a second version with five
age-free composition components of the genes outside every set; gene sets against genes of the same
expression decile, in SDs of the metric; size-matched subsamples (``SUBSAMPLE`` animals per tissue,
70 % of the species' smallest tissue, without replacement).

**Genes.** The counts are named by rodent gene symbol. A rodent symbol is put on a human gene by its
upper-case form, as an HGNC approved symbol first, then a previous symbol, then a unique alias
(``symbol_ensembl``); this is a symbol-level orthology, right for the large majority of one-to-one
orthologs and silent for genes named differently in the two species, which drop out (the
``genes`` table counts them). The human gene sets then reach rodent genes through the Ensembl ID.

**Turnover classes:**

* *renewing* - bone marrow, spleen, thymus, lymph node, Peyer's patch, colon, duodenum, jejunum,
  ileum, skin;
* *slow* - liver, kidney, adrenal and pituitary glands, thyroid, gonadal white adipose (gWAT) and
  femur (bone; the series does not say whether its marrow was flushed, and marrow is profiled as a
  tissue of its own, so femur is counted as bone, which turns over slowly);
* *heart & muscle* - the four heart chambers (the skeletal muscles of the atlas are not in these
  count files);
* *brain* - cortex, hippocampus, hypothalamus, cerebellum, spinal cord, and in rat the retina,
  which is neural tissue.

Left out: rat eye tissues of mixed or unclear turnover (lens, cornea, iris, trabecular meshwork,
retinal pigment epithelium with choroid).
"""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import polars as pl
from scipy.stats import false_discovery_control

from entropy_age import plots
from entropy_age.data.annotation import AnnotationLayout, gene_set_ensembl, symbol_ensembl
from entropy_age.data.geo import read_series_matrix
from entropy_age.experiments.common import Context, dummies, write_table
from entropy_age.experiments.human_tissues import _subsample
from entropy_age.experiments.rna_gtex import ADJUSTMENTS, METRICS, expressed, log_cpm, tissue_metrics
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

RENEWING = (
    "bone marrow",
    "spleen",
    "thymus",
    "lymph node",
    "Peyer's patch",
    "colon",
    "duodenum",
    "jejunum",
    "ileum",
    "skin",
)
SLOW = ("liver", "kidney", "adrenal gland", "pituitary gland", "thyroid gland", "gWAT", "femur")
HEART_MUSCLE = ("left atrium", "right atrium", "left ventricle", "right ventricle")
BRAIN = ("cortex", "hippocampus", "hypothalamus", "cerebellum", "spinal cord", "retina")
LEFT_OUT = ("lens", "cornea", "iris", "trabecular meshwork", "RPE and choroid")
TURNOVER: dict[str, str] = (
    dict.fromkeys(RENEWING, "renewing")
    | dict.fromkeys(SLOW, "slow")
    | dict.fromkeys(HEART_MUSCLE, "heart & muscle")
    | dict.fromkeys(BRAIN, "brain")
)
CLASSES = ("renewing", "slow", "heart & muscle", "brain")
CONTRASTS = (
    ("renewing", "brain"),
    ("slow", "brain"),
    ("heart & muscle", "brain"),
    ("renewing", "heart & muscle"),
    ("slow", "heart & muscle"),
    ("renewing", "slow"),
)
#: species -> (GEO series, counts file)
SPECIES: dict[str, tuple[str, str]] = {
    "mouse": ("GSE314817", "GSE314817_mouse_raw_counts_wo_muscle.txt.gz"),
    "rat": ("GSE314970", "GSE314970_rat_raw_counts_wo_muscle.txt.gz"),
}
SUBSAMPLE_FRACTION = 0.7


# ----------------------------------------------------------------------------- inputs
def read_atlas(matrix: Path, counts: Path) -> tuple[pl.DataFrame, list[str], np.ndarray]:
    """Samples (``sample``, ``tissue``, ``sex``, ``age`` in months), gene symbols and a samples x genes count matrix.

    The count columns are the numeric suffix of each GEO sample title (``male_cortex_27 months_1377531``).
    """
    meta = read_series_matrix(matrix).samples
    samples = meta.select(
        pl.col("title").str.extract(r"_(\d+)$").alias("sample"),
        pl.col("sample").alias("gsm"),
        "tissue",
        pl.col("sex").str.to_lowercase(),
        pl.col("age").str.extract(r"([\d.]+)").cast(pl.Float64).alias("age"),
    )
    if samples["sample"].null_count() or samples["age"].null_count():
        raise ValueError(f"{matrix.name}: sample titles or ages that do not parse")
    frame = pl.read_csv(counts, separator="\t", infer_schema_length=0)
    genes = frame[frame.columns[0]].to_list()
    columns = set(frame.columns[1:])
    absent = [s for s in samples["sample"].to_list() if s not in columns]
    if absent:
        print(f"{counts.name}: {len(absent)} GEO samples have no count column, e.g. {absent[:3]}", file=sys.stderr)
    samples = samples.filter(pl.col("sample").is_in(list(columns)))
    x = frame.select(samples["sample"].to_list()).cast(pl.Float64).to_numpy().T
    return samples, genes, x


def _prepare(tissue: str, samples: pl.DataFrame, counts: np.ndarray, features: list[str], coding: np.ndarray):
    values, library = log_cpm(counts)
    keep = coding & expressed(values)
    technical = np.column_stack([dummies(samples["sex"].to_numpy()), np.log2(library)[:, None]])
    names = [f for f, k in zip(features, keep, strict=True) if k]
    animals = samples["sample"].to_numpy()
    p = Prepared(
        tissue, values[:, keep].astype(np.float32), samples["age"].to_numpy(), technical, names, animals, technical
    )
    return p


# ----------------------------------------------------------------------------- shared with rna_gtex
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


def _rows(species, tissue, adjustment, p, rows, membership, permutations, seed) -> pl.DataFrame:
    sets = membership.filter(pl.col("feature").is_in(p.features))
    table = _matched_sd(tissue_metrics(p, rows), sets, METRICS, permutations, seed)
    return table.with_columns(
        species=pl.lit(species), tissue=pl.lit(tissue), turnover=pl.lit(TURNOVER[tissue]), adjustment=pl.lit(adjustment)
    )


def _boot(seed: int) -> list[dict]:
    tissues, membership, size, species = _STATE["tissues"], _STATE["membership"], _STATE["size"], _STATE["species"]
    rng = np.random.default_rng(seed)
    out = []
    for tissue in TURNOVER:
        if (tissue, "covariates") not in tissues:
            continue
        rows = _subsample(rng, tissues[(tissue, "covariates")].x.shape[0], size)
        for adjustment in ADJUSTMENTS:
            table = _rows(species, tissue, adjustment, tissues[(tissue, adjustment)], rows, membership, 0, seed)
            out += [
                {"boot": seed, **r}
                for r in table.select(
                    "species", "tissue", "adjustment", "metric", "group", "matched_diff_sd"
                ).to_dicts()
            ]
    return out


# ----------------------------------------------------------------------------- run
def _species_inputs(ctx: Context) -> dict[str, tuple[Path, Path]]:
    out = {}
    for species, (gse, counts) in SPECIES.items():
        folder = ctx.paths.input / "geo" / gse
        out[species] = (folder / f"{gse}_series_matrix.txt.gz", folder / counts)
    return out


def analyse_species(species, samples, genes, counts, hgnc, gene_sets, n_boot, workers):
    """Every table of one species: tissues, gene mapping, point estimates, resamples, contrasts."""
    mapping = symbol_ensembl(genes, hgnc)
    protein = set(hgnc.filter(pl.col("protein_coding") & pl.col("ensembl").is_not_null())["ensembl"].to_list())
    ensembl = dict(mapping.select("symbol", "ensembl").iter_rows())
    coding = np.array([ensembl.get(g) in protein for g in genes])
    membership = (
        gene_set_ensembl(gene_sets, hgnc)
        .filter(pl.col("group").is_in(SETS))
        .join(mapping.select(pl.col("ensembl").alias("feature"), pl.col("symbol").alias("rodent")), on="feature")
        .select(pl.col("rodent").alias("feature"), "group", "gene")
    )
    members = set(membership["feature"].to_list())
    tissues, table = {}, []
    for tissue, sel in samples.with_row_index("row").group_by("tissue", maintain_order=True):
        name = tissue[0]
        if name not in TURNOVER:
            table.append(
                {
                    "tissue": name,
                    "turnover": None,
                    "animals": sel.height,
                    "used": False,
                    "reason": "eye tissue of unclear turnover" if name in LEFT_OUT else "not classified",
                }
            )
            continue
        p = _prepare(name, sel, counts[sel["row"].to_numpy()], genes, coding)
        for adjustment_name, adjustment in ADJUSTMENTS.items():
            tissues[(name, adjustment_name)] = p if adjustment is None else _adjusted(p, members, *adjustment)[0]
        table.append(
            {
                "tissue": name,
                "turnover": TURNOVER[name],
                "animals": sel.height,
                "used": True,
                "reason": None,
                "age_min": float(p.ages.min()),
                "age_max": float(p.ages.max()),
                "genes": len(p.features),
            }
            | {g: int(membership.filter(pl.col("group") == g)["feature"].is_in(p.features).sum()) for g in SETS}
        )
        print(
            f"[rna_rodents] {species} {name}: {sel.height} animals, {len(p.features)} genes",
            file=sys.stderr,
            flush=True,
        )
    used = [t for t in table if t["used"]]
    size = int(SUBSAMPLE_FRACTION * min(t["animals"] for t in used))
    point = pl.concat(
        [_rows(species, t, a, p, np.arange(p.x.shape[0]), membership, 1000, 3) for (t, a), p in tissues.items()],
        how="diagonal",
    )
    point = point.with_columns(p_bh=pl.Series(false_discovery_control(point["p_permutation"].to_numpy(), method="bh")))
    state = {"tissues": tissues, "membership": membership, "size": size, "species": species}
    boot = pl.DataFrame(_run_parallel(_boot, [11000 + b for b in range(n_boot)], workers, state))
    contrasts = class_contrasts(boot, point, TURNOVER, CONTRASTS).with_columns(species=pl.lit(species))
    intervals = tissue_intervals(boot, TURNOVER).with_columns(species=pl.lit(species))
    genes_table = {
        "species": species,
        "symbols": len(genes),
        "mapped": mapping.height,
        "protein_coding": int(coding.sum()),
        "subsample": size,
    } | {f"matched_by_{k}": v for k, v in mapping.group_by("matched_by").len().iter_rows()}
    tissues_table = pl.DataFrame(table, infer_schema_length=None).with_columns(
        species=pl.lit(species), reason=pl.col("reason").cast(pl.String), turnover=pl.col("turnover").cast(pl.String)
    )
    return tissues_table, genes_table, point, boot, contrasts, intervals


def run(ctx: Context) -> list[Path]:
    inputs = _species_inputs(ctx)
    annotation = AnnotationLayout(ctx.paths)
    needed = [p for pair in inputs.values() for p in pair] + [annotation.gene_sets, annotation.hgnc]
    missing = [str(p) for p in needed if not p.exists()]
    if missing:
        print(
            f"rna_rodents: inputs missing ({', '.join(missing)}); run `entropy-age fetch annotation` and "
            "`entropy-age fetch rodent-rna` first. Nothing was run.",
            file=sys.stderr,
        )
        return []
    plots.apply_style()
    out = ctx.out_dir("rna_rodents")
    workers, n_boot = _settings(ctx)
    hgnc, gene_sets = pl.read_parquet(annotation.hgnc), pl.read_parquet(annotation.gene_sets)
    parts = []
    for species, (matrix, counts) in inputs.items():
        samples, genes, x = read_atlas(matrix, counts)
        parts.append(analyse_species(species, samples, genes, x, hgnc, gene_sets, n_boot, workers))
    tissues, genes_rows, point, boot, contrasts, intervals = zip(*parts, strict=True)
    write_table(pl.concat(tissues, how="diagonal"), out / "tissues")
    write_table(pl.DataFrame(list(genes_rows), infer_schema_length=None), out / "genes")
    write_table(pl.concat(point, how="diagonal"), out / "point_estimates")
    pl.concat(boot, how="diagonal").write_parquet(out / "bootstrap.parquet")
    contrasts = pl.concat(contrasts, how="diagonal")
    intervals = pl.concat(intervals, how="diagonal")
    write_table(contrasts, out / "contrasts")
    write_table(intervals, out / "tissue_intervals")
    return [_figure(ctx, intervals, contrasts)]


def _figure(ctx: Context, intervals: pl.DataFrame, contrasts: pl.DataFrame) -> Path:
    import matplotlib.pyplot as plt

    species = intervals["species"].unique(maintain_order=True).to_list()
    n = max(
        sum(t in set(intervals.filter(pl.col("species") == s)["tissue"].to_list()) for t in TURNOVER) for s in species
    )
    fig = plt.figure(figsize=(15, 0.22 * n * len(species) + 3))
    grid = fig.add_gridspec(len(species), 2, width_ratios=(1, 1.1), wspace=0.55, hspace=0.35)
    colors = {g: plots.PALETTE[i] for i, g in enumerate(SETS)}
    letters = iter("abcdef")
    for row, sp in enumerate(species):
        sel = intervals.filter(
            (pl.col("species") == sp)
            & (pl.col("adjustment") == "covariates")
            & (pl.col("metric") == "dispersion_slope")
        )
        tissues = [t for t in TURNOVER if t in set(sel["tissue"].to_list())]
        ax = fig.add_subplot(grid[row, 0])
        for j, group in enumerate(SETS):
            lookup = {r["tissue"]: r for r in sel.filter(pl.col("group") == group).iter_rows(named=True)}
            for i, tissue in enumerate(tissues):
                if tissue in lookup:
                    r, y = lookup[tissue], i + (j - 1) * 0.25
                    ax.plot([r["ci_low"], r["ci_high"]], [y, y], color=colors[group], lw=0.9)
                    ax.scatter(
                        [r["median"]], [y], color=colors[group], s=9, label=SET_LABELS[group] if i == 0 else None
                    )
        start = 0
        for cls in CLASSES:
            k = sum(TURNOVER[t] == cls for t in tissues)
            if start:
                ax.axhline(start - 0.5, color=plots.MUTED, lw=0.6, ls=":")
            ax.text(
                1.0,
                start + k / 2 - 0.5,
                cls,
                transform=ax.get_yaxis_transform(),
                rotation=-90,
                va="center",
                ha="left",
                fontsize=7,
                color=plots.MUTED,
            )
            start += k
        ax.axvline(0, color=plots.MUTED, lw=0.8)
        ax.set_yticks(range(len(tissues)), tissues, fontsize=6.5)
        ax.set_ylim(len(tissues) - 0.5, -0.5)
        ax.set(xlabel="Set minus matched genes (SD)", title=f"{sp}: spread gain with age, per tissue")
        if row == 0:
            ax.legend(fontsize=7, loc="lower left")
        plots.label_panel(ax, next(letters))
        ax = fig.add_subplot(grid[row, 1])
        frame = (
            contrasts.filter((pl.col("species") == sp) & (pl.col("metric") == "dispersion_slope"))
            .with_columns(
                label=pl.format(
                    "{} · {} · {}",
                    pl.col("group").replace_strict(SET_LABELS),
                    pl.col("contrast"),
                    pl.col("adjustment").replace({"covariates": "cov", "pc5_age_free": "+PCs"}),
                )
            )
            .sort("group", "contrast", "adjustment")
        )
        for y, r in enumerate(frame.iter_rows(named=True)):
            color = plots.ACCENT if r["ci_high"] < 0 or r["ci_low"] > 0 else plots.MUTED
            ax.plot([r["ci_low"], r["ci_high"]], [y, y], color=color, lw=1.1)
            ax.scatter([r["difference"]], [y], color=color, s=12)
        ax.axvline(0, color=plots.MUTED, lw=0.8)
        ax.set_yticks(range(frame.height), frame["label"].to_list(), fontsize=5.5)
        ax.set_ylim(max(frame.height - 0.5, 0.5), -0.5)
        ax.set(xlabel="Difference (SD); < 0 = more protection in the first", title=f"{sp}: class contrasts")
        plots.label_panel(ax, next(letters))
    fig.suptitle(
        "Transcriptional drift by turnover class in mouse and rat (GSE314817, GSE314970)",
        x=0.02,
        ha="left",
        fontsize=12,
        fontweight="bold",
    )
    fig.subplots_adjust(left=0.12, right=0.98, top=0.95, bottom=0.05)
    return plots.save(fig, ctx.figure_path("rna_rodents"))
