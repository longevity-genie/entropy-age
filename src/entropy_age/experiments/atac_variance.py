"""Does chromatin accessibility become more variable between people with age, and less so at essential genes?
(roadmap RM19)

Methylation drifts apart between people with age, and essential promoters drift less (the shield). Chromatin
opening could do the same. Morandini et al. 2023 (GSE193140) measured bulk ATAC-seq of blood mononuclear cells
(PBMC) in people aged 20-74, and deposited raw counts over 42,557 peaks (hg38).

Per peak: log2 counts per million, the residual after age and sex (and ``COMPONENTS`` leading principal
components of all peaks, for cell mix), and the age slope of the absolute residual over the residual's SD. That
slope is how fast people spread apart at the peak, in units of the peak's own spread. Peaks are tied to genes
through the 450k promoter CpGs inside them (hg38 positions from Zhou et al.'s manifest), and a gene set's peaks
are compared with the other promoter peaks of the same CpG-island status, mean-accessibility decile and width
tertile (``stochasticity.compare_groups``, labels permuted within strata).

A donor sequenced twice keeps its first library. Light: runs in a minute.

``entropy-age fetch atac`` then ``entropy-age figures atac_variance``.
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
from entropy_age.experiments.common import Context, dummies, write_table
from entropy_age.experiments.tissue_contrast import SET_LABELS, SETS
from entropy_age.stochasticity import baseline_strata, compare_groups

GSE = "GSE193140"
COUNTS = f"{GSE}_atac_counts_qval_0.05.txt.gz"
HG38_MANIFEST = "zhou/HM450.hg38.manifest.tsv.gz"
HG38_MANIFEST_URL = "https://github.com/zhou-lab/InfiniumAnnotationV1/raw/main/Anno/HM450/HM450.hg38.manifest.tsv.gz"
COMPONENTS = (0, 5)
PERMUTATIONS = 2000


def layout(paths) -> tuple[Path, Path, Path]:
    base = paths.input / "geo" / GSE
    return base / COUNTS, base / f"{GSE}_series_matrix.txt.gz", AnnotationLayout(paths).input / HG38_MANIFEST


def read_donors(matrix: Path) -> pl.DataFrame:
    """``donor`` (the title without ``ATAC_``), ``age`` and ``sex`` per sample."""
    chars = _characteristics(matrix)
    titles = accessions = None
    with gzip.open(matrix, "rt") as handle:
        for line in handle:
            cells = [c.strip().strip('"') for c in line.rstrip("\n").split("\t")[1:]]
            if line.startswith("!Sample_title"):
                titles = cells
            elif line.startswith("!Sample_geo_accession"):
                accessions = cells
    frame = pl.DataFrame({"sample": accessions, "donor": [t.removeprefix("ATAC_") for t in titles]})
    return frame.join(
        chars.select("sample", pl.col("age").cast(pl.Float64), pl.col("gender").alias("sex")), on="sample"
    ).select("donor", "age", "sex")


def first_library(columns: list[str]) -> dict[str, str]:
    """Donor -> its first library column (``CR_027_1`` before ``CR_027_2``)."""
    out: dict[str, str] = {}
    for c in sorted(columns, key=lambda c: (c.rsplit("_", 1)[0], int(c.rsplit("_", 1)[1]))):
        out.setdefault(c.rsplit("_", 1)[0], c)
    return out


def log_cpm(counts: np.ndarray) -> np.ndarray:
    """log2 counts per million, samples in rows."""
    return np.log2(counts / counts.sum(axis=1, keepdims=True) * 1e6 + 1.0)


def spread_slopes(y: np.ndarray, age: np.ndarray, covariates: np.ndarray) -> np.ndarray:
    """Per column of ``y``: the age slope of |residual| / residual SD, residuals after age and ``covariates``."""
    design = np.column_stack([np.ones(age.size), age, covariates])
    coef, *_ = np.linalg.lstsq(design, y, rcond=None)
    resid = y - design @ coef
    sd = resid.std(axis=0, ddof=design.shape[1])
    a = np.abs(resid) / np.where(sd > 0, sd, np.nan)
    d2 = np.column_stack([np.ones(age.size), age])
    return np.linalg.lstsq(d2, a, rcond=None)[0][1]


def peak_cpgs(peaks: pl.DataFrame, manifest: pl.DataFrame) -> pl.DataFrame:
    """``peak`` (row index) and ``cpg`` for every 450k CpG inside a peak (half-open [start, end))."""
    rows = []
    for (chrom,), probes in manifest.group_by("chr"):
        p = peaks.with_row_index("peak").filter(pl.col("chr") == chrom).sort("start")
        if p.is_empty():
            continue
        start, end, idx = p["start"].to_numpy(), p["end"].to_numpy(), p["peak"].to_numpy()
        pos = probes["pos"].to_numpy()
        i = np.searchsorted(start, pos, side="right") - 1
        inside = (i >= 0) & (pos < end[np.clip(i, 0, None)])
        rows.append(pl.DataFrame({"peak": idx[i[inside]], "cpg": probes["cpg"].filter(pl.Series(inside))}))
    return pl.concat(rows) if rows else pl.DataFrame(schema={"peak": pl.UInt32, "cpg": pl.String})


def run(ctx: Context) -> list[Path]:
    counts_path, matrix, manifest_path = layout(ctx.paths)
    ann = AnnotationLayout(ctx.paths)
    missing = [str(p) for p in (counts_path, matrix, manifest_path, ann.links_450k, ann.probes_450k) if not p.exists()]
    if missing:
        print(
            f"atac_variance: run `entropy-age fetch atac` first ({', '.join(missing)}). Nothing run.", file=sys.stderr
        )
        return []
    plots.apply_style()
    out = ctx.out_dir("atac_variance")
    table = pl.read_csv(counts_path, separator=" ")
    peaks = table.select("chr", "start", "end")
    libraries = first_library(table.columns[3:])
    donors = read_donors(matrix).filter(pl.col("donor").is_in(list(libraries))).sort("donor")
    y = log_cpm(table.select([libraries[d] for d in donors["donor"]]).to_numpy().T.astype(float))
    age = donors["age"].to_numpy()

    manifest = (
        pl.read_csv(manifest_path, separator="\t", columns=["CpG_chrm", "CpG_beg", "Probe_ID"], null_values=["NA"])
        .drop_nulls()
        .select(pl.col("CpG_chrm").alias("chr"), pl.col("CpG_beg").alias("pos"), pl.col("Probe_ID").alias("cpg"))
    )
    probes = pl.read_parquet(ann.probes_450k).select("cpg", "cpg_island", "promoter")
    inside = peak_cpgs(peaks, manifest).join(probes, on="cpg")
    promoter_peaks = inside.filter(pl.col("promoter"))
    gene_sets = pl.read_parquet(ann.gene_sets)
    cpg_sets = link_membership(pl.read_parquet(ann.links_450k), gene_sets).filter(pl.col("group").is_in(SETS))
    membership = (
        promoter_peaks.join(cpg_sets.rename({"feature": "cpg"}), on="cpg")
        .select(pl.col("peak").cast(pl.String).alias("feature"), "group", "gene")
        .unique()
    )
    info = (
        promoter_peaks.group_by("peak")
        .agg(pl.col("cpg_island").any())
        .with_columns(feature=pl.col("peak").cast(pl.String))
    )
    keep = info["peak"].to_numpy()
    width = (peaks["end"] - peaks["start"]).to_numpy()[keep]
    mean = y[:, keep].mean(axis=0)
    stratum = pl.Series(
        [
            f"{i}|{m}|{w}"
            for i, m, w in zip(info["cpg_island"], baseline_strata(mean, 10), baseline_strata(width, 3), strict=True)
        ]
    )
    rows = []
    sex = dummies(donors["sex"].to_numpy())
    for k in COMPONENTS:
        cov = sex
        if k:
            centred = y - y.mean(axis=0)
            u, s, _ = np.linalg.svd(centred, full_matrices=False)
            cov = np.column_stack([sex, u[:, :k] * s[:k]])
        slopes = spread_slopes(y[:, keep], age, cov)
        metrics = info.select("feature").with_columns(spread_slope=pl.Series(slopes), stratum=stratum)
        res = compare_groups(
            metrics,
            membership,
            metric="spread_slope",
            stratify_by="stratum",
            n_permutations=PERMUTATIONS,
            seed=ctx.seed,
        )
        sd = float(np.nanstd(slopes))
        rows.append(
            res.with_columns(
                components=pl.lit(k),
                gap_sd=pl.col("matched_diff") / sd,
                peaks=pl.lit(int(keep.size)),
                donors=pl.lit(int(age.size)),
            )
        )
        print(f"atac_variance: {k} components done", file=sys.stderr, flush=True)
    result = pl.concat(rows, how="diagonal")
    write_table(result, out / "gaps")
    return [_figure(ctx, result, age)]


def _figure(ctx: Context, result: pl.DataFrame, age: np.ndarray) -> Path:
    fig, axes = plots.panels(2, ncols=2, size=4.2)
    ax = axes[0]
    ax.hist(age, bins=12, color=plots.PALETTE[0])
    ax.set(xlabel="Age (years)", ylabel="Donors", title=f"PBMC ATAC donors (n = {age.size})")
    plots.label_panel(ax, "a")
    ax = axes[1]
    for j, k in enumerate(COMPONENTS):
        sub = result.filter(pl.col("components") == k)
        for g, group in enumerate(SETS):
            r = sub.filter(pl.col("group") == group)
            if r.is_empty():
                continue
            ax.scatter(
                r["gap_sd"][0],
                g + 0.25 * j,
                color=plots.PALETTE[g],
                marker="os"[j],
                s=30,
                label=f"{k} components" if g == 0 else None,
            )
            ax.annotate(
                f"p={r['p_permutation'][0]:.2g}",
                (r["gap_sd"][0], g + 0.25 * j),
                fontsize=6,
                xytext=(4, 0),
                textcoords="offset points",
                va="center",
            )
    ax.axvline(0, color=plots.MUTED, lw=0.8)
    ax.set_yticks(range(len(SETS)), [SET_LABELS[s] for s in SETS])
    ax.invert_yaxis()
    ax.legend(fontsize=7, frameon=False)
    ax.set(
        xlabel="Set minus matched promoter peaks (SD); < 0: spreads less with age",
        title="Between-person accessibility spread with age",
    )
    plots.label_panel(ax, "b")
    return plots.save(fig, ctx.figure_path("atac_variance"))
