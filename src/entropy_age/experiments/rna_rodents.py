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
from entropy_age.experiments.rna_gtex import (
    ADJUSTMENTS,
    METRICS,
    class_contrasts,
    contrast_panel,
    expressed,
    log_cpm,
    tissue_intervals,
    tissue_metrics,
    tissue_panel,
)
from entropy_age.experiments.tissue_contrast import (
    _STATE,
    SETS,
    Prepared,
    _adjusted,
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
    spread = intervals.filter(pl.col("metric") == "dispersion_slope")
    n = max(spread.filter(pl.col("species") == s)["tissue"].n_unique() for s in species)
    fig = plt.figure(figsize=(16, (0.24 * n + 1.5) * len(species) + 1))
    grid = fig.add_gridspec(len(species), 3, width_ratios=(1, 0.8, 1), wspace=0.45, hspace=0.3)
    letters = iter("abcdef")
    for row, sp in enumerate(species):
        sel = spread.filter(pl.col("species") == sp)
        a, b, c = (fig.add_subplot(grid[row, i]) for i in range(3))
        tissue_panel(
            a, sel.filter(pl.col("adjustment") == "covariates"), TURNOVER, CLASSES, f"{sp}: spread gain with age", True
        )
        tissue_panel(
            b,
            sel.filter(pl.col("adjustment") == "pc5_age_free"),
            TURNOVER,
            CLASSES,
            f"{sp}: plus 5 age-free components",
            False,
        )
        contrast_panel(
            c,
            contrasts.filter((pl.col("species") == sp) & (pl.col("metric") == "dispersion_slope")),
            f"{sp}: class contrasts",
        )
        for ax in (a, b, c):
            plots.label_panel(ax, next(letters))
        if row == 0:
            a.legend(fontsize=7, loc="lower left")
    fig.suptitle(
        "Transcriptional drift by turnover class in mouse and rat (GSE314817, GSE314970)",
        x=0.02,
        ha="left",
        fontsize=12,
        fontweight="bold",
    )
    fig.subplots_adjust(left=0.1, right=0.99, top=0.93, bottom=0.04)
    return plots.save(fig, ctx.figure_path("rna_rodents"))
