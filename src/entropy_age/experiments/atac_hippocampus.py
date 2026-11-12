"""Replication of the cortical accessibility result in hippocampus (roadmap RM19).

In 357 prefrontal cortices (``atac_cortex``) essential promoters spread apart between people with age more than
matched promoters, in nearly every cell type. Zemke et al. 2024 (GSE278576) profiled single-nucleus multiome
(ATAC and RNA) in the hippocampus of 40 donors aged 20-95, with per-donor 10x matrices. Cell Ranger called peaks
per donor, so peaks differ between donors; features here are genes instead: a donor's count for a gene is the
sum, over the donor's peaks that hold one of the gene's 450k promoter CpGs (hg38), of the reads of nuclei of one
subclass. Gene-set genes are compared with the other promoter genes of the same CpG-island status, mean decile
and promoter-CpG-count tertile, as in ``atac_cortex``, with the donor's promoter read fraction and log nucleus
count as covariates (``ADJUSTMENTS``).

Run as a job: ``entropy-age figures atac_hippocampus``.
"""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import polars as pl

from entropy_age import plots
from entropy_age.data.annotation import AnnotationLayout, link_membership
from entropy_age.experiments.atac_cortex import parse_peaks, promoter_fraction
from entropy_age.experiments.atac_variance import HG38_MANIFEST, peak_cpgs, spread_slopes
from entropy_age.experiments.common import Context, dummies, write_table
from entropy_age.experiments.tissue_contrast import SET_LABELS, SETS
from entropy_age.stochasticity import baseline_strata, compare_groups

GSE = "GSE278576"
METADATA = f"{GSE}_hippocampus_RNA_seurat_object_filtered_cells_metadata.tsv.gz"
MIN_NUCLEI = 50
MIN_DONORS = 20
#: name -> (leading components, promoter read fraction as covariate); 40 donors carry at most a few components
ADJUSTMENTS = {"none": (0, False), "promoter fraction": (0, True), "both": (3, True)}
PERMUTATIONS = 2000


def base_dir(paths) -> Path:
    return paths.input / "geo" / GSE


def subclass_peak_sums(indptr, indices, data, columns_by_group: dict[str, np.ndarray], n_features: int):
    """Per group: counts summed over its barcode columns of a CSC (features x barcodes) matrix."""
    out = {}
    for group, cols in columns_by_group.items():
        starts, ends = indptr[cols], indptr[cols + 1]
        idx = (
            np.concatenate([np.arange(a, b) for a, b in zip(starts, ends, strict=True)])
            if cols.size
            else np.array([], int)
        )
        out[group] = np.bincount(indices[idx], weights=data[idx].astype(float), minlength=n_features)
    return out


def gene_counts(peak_sums: np.ndarray, peak_gene: pl.DataFrame, genes: list[str]) -> np.ndarray:
    """Counts per gene: the sum over the distinct peaks that hold any of the gene's promoter CpGs."""
    pairs = peak_gene.unique(["peak", "gene"])
    index = {g: i for i, g in enumerate(genes)}
    out = np.zeros(len(genes))
    np.add.at(out, np.array([index[g] for g in pairs["gene"]], dtype=int), peak_sums[pairs["peak"].to_numpy()])
    return out


def run(ctx: Context) -> list[Path]:
    import h5py

    base = base_dir(ctx.paths)
    ann = AnnotationLayout(ctx.paths)
    manifest_path = ann.input / HG38_MANIFEST
    matrices = sorted(base.glob(f"{GSE}_*_raw_feature_bc_matrix.h5"))
    if not matrices or not (base / METADATA).exists() or not manifest_path.exists():
        print(f"atac_hippocampus: inputs missing under {base}. Nothing run.", file=sys.stderr)
        return []
    plots.apply_style()
    out = ctx.out_dir("atac_hippocampus")
    meta = pl.read_csv(base / METADATA, separator="\t", infer_schema_length=10000).select(
        pl.col("bacrode").alias("cell"),
        pl.col("orig.ident").alias("donor"),
        pl.col("Age").cast(pl.Float64).alias("age"),
        pl.col("Gender").alias("sex"),
        "subclass",
    )
    manifest = (
        pl.read_csv(manifest_path, separator="\t", columns=["CpG_chrm", "CpG_beg", "Probe_ID"], null_values=["NA"])
        .drop_nulls()
        .select(pl.col("CpG_chrm").alias("chr"), pl.col("CpG_beg").alias("pos"), pl.col("Probe_ID").alias("cpg"))
    )
    links = (
        pl.read_parquet(ann.links_450k)
        .filter(pl.col("promoter"))
        .select("cpg", pl.col("symbol").alias("gene"))
        .unique()
    )
    island = pl.read_parquet(ann.probes_450k).select("cpg", "cpg_island")
    genes = sorted(links["gene"].unique().to_list())
    cpg_sets = link_membership(pl.read_parquet(ann.links_450k), pl.read_parquet(ann.gene_sets)).filter(
        pl.col("group").is_in(SETS)
    )
    membership = (
        cpg_sets.select(pl.col("feature").alias("cpg"), "group")
        .join(links, on="cpg")
        .select(pl.col("gene").alias("feature"), "group", "gene")
        .unique()
    )
    gene_info = (
        links.join(island, on="cpg")
        .group_by("gene")
        .agg(pl.col("cpg_island").any(), pl.len().alias("cpgs"))
        .sort("gene")
    )
    rows = []
    for path in matrices:
        donor = path.name.removeprefix(f"{GSE}_").removesuffix("_raw_feature_bc_matrix.h5")
        cells = meta.filter(pl.col("donor") == donor)
        if cells.is_empty():
            print(f"atac_hippocampus: {donor}: no cells in the metadata; left out", file=sys.stderr)
            continue
        with h5py.File(path, "r") as f:
            m = f["matrix"]
            kind = m["features"]["feature_type"][:]
            peak_rows = np.flatnonzero(kind == b"Peaks")
            peaks = parse_peaks([x.decode() for x in m["features"]["id"][:][peak_rows]])
            barcodes = [f"{donor}_{b.decode()}" for b in m["barcodes"][:]]
            col = {b: i for i, b in enumerate(barcodes)}
            by_group = {
                s: np.array([col[c] for c in g["cell"] if c in col], dtype=int)
                for (s,), g in cells.group_by("subclass")
            }
            sums = subclass_peak_sums(m["indptr"][:], m["indices"][:], m["data"][:], by_group, kind.size)
        peak_gene = peak_cpgs(peaks, manifest).join(links, on="cpg")
        for subclass, s in sums.items():
            peak_counts = s[peak_rows]
            rows.append(
                {
                    "donor": donor,
                    "subclass": subclass,
                    "nuclei": int(by_group[subclass].size),
                    "total": float(peak_counts.sum()),
                    "counts": gene_counts(peak_counts, peak_gene, genes),
                }
            )
        print(f"atac_hippocampus: {donor} done", file=sys.stderr, flush=True)
    donors = meta.group_by("donor").agg(pl.col("age").first(), pl.col("sex").first())
    results = []
    for subclass in sorted({r["subclass"] for r in rows}):
        sel = [r for r in rows if r["subclass"] == subclass and r["nuclei"] >= MIN_NUCLEI]
        if len(sel) < MIN_DONORS:
            continue
        info = pl.DataFrame([{k: r[k] for k in ("donor", "nuclei", "total")} for r in sel]).join(donors, on="donor")
        x = np.vstack([r["counts"] for r in sel])
        y = np.log2(x / info["total"].to_numpy()[:, None] * 1e6 + 1.0)
        expressed = (x > 0).mean(axis=0) >= 0.9
        y, gi = y[:, expressed], gene_info.filter(pl.Series(expressed))
        age = info["age"].to_numpy()
        base_cov = np.column_stack([dummies(info["sex"].to_numpy()), np.log(info["nuclei"].to_numpy())])
        quality = promoter_fraction(x, info["total"].to_numpy())
        mean = y.mean(axis=0)
        stratum = pl.Series(
            [
                f"{i}|{m}|{c}"
                for i, m, c in zip(
                    gi["cpg_island"],
                    baseline_strata(mean, 10),
                    baseline_strata(gi["cpgs"].to_numpy().astype(float), 3),
                    strict=True,
                )
            ]
        )
        for adjustment, (k, by_quality) in ADJUSTMENTS.items():
            cov = np.column_stack([base_cov, quality]) if by_quality else base_cov
            if k:
                u, sv, _ = np.linalg.svd(y - mean, full_matrices=False)
                cov = np.column_stack([cov, u[:, :k] * sv[:k]])
            slopes = spread_slopes(y, age, cov)
            metrics = pl.DataFrame({"feature": gi["gene"], "spread_slope": slopes, "stratum": stratum})
            res = compare_groups(
                metrics,
                membership,
                metric="spread_slope",
                stratify_by="stratum",
                n_permutations=PERMUTATIONS,
                seed=ctx.seed,
            )
            sd = float(np.nanstd(slopes))
            results.append(
                res.with_columns(
                    subclass=pl.lit(subclass),
                    adjustment=pl.lit(adjustment),
                    gap_sd=pl.col("matched_diff") / sd,
                    donors=pl.lit(len(sel)),
                    genes=pl.lit(int(expressed.sum())),
                )
            )
        print(f"atac_hippocampus: {subclass}: {len(sel)} donors done", file=sys.stderr, flush=True)
    table = pl.concat(results, how="diagonal")
    write_table(table, out / "gaps")
    return [_figure(ctx, table)]


def _figure(ctx: Context, table: pl.DataFrame) -> Path:
    fig, axes = plots.panels(len(ADJUSTMENTS), ncols=len(ADJUSTMENTS), size=4.0)
    for i, adjustment in enumerate(ADJUSTMENTS):
        ax = axes[i]
        sub = table.filter(pl.col("adjustment") == adjustment)
        cells = sorted(sub["subclass"].unique().to_list())
        for g, group in enumerate(SETS):
            r = sub.filter(pl.col("group") == group).sort("subclass")
            ys = [cells.index(c) + 0.2 * g for c in r["subclass"].to_list()]
            ax.scatter(r["gap_sd"], ys, color=plots.PALETTE[g], s=20, label=SET_LABELS[group])
            for xv, yv, p in zip(r["gap_sd"], ys, r["p_permutation"], strict=True):
                if p < 0.05:
                    ax.annotate("*", (xv, yv), xytext=(4, -3), textcoords="offset points", fontsize=9)
        donors = dict(zip(sub["subclass"], sub["donors"], strict=False))
        ax.set_yticks(
            [j + 0.2 for j in range(len(cells))],
            [f"{c} ({donors[c]})" for c in cells] if i == 0 else [""] * len(cells),
            fontsize=7,
        )
        ax.axvline(0, color=plots.MUTED, lw=0.8)
        ax.invert_yaxis()
        if i == 0:
            ax.legend(fontsize=7, frameon=False, loc="lower left")
        ax.set(xlabel="Set minus matched genes (SD)", title=f"Adjusted for: {adjustment}")
        plots.label_panel(ax, "abc"[i])
    fig.suptitle(
        "Accessibility spread with age, hippocampus subclasses (GSE278576; > 0: essential promoters spread more)",
        x=0.02,
        y=1.04,
        ha="left",
        fontsize=11,
        fontweight="bold",
    )
    return plots.save(fig, ctx.figure_path("atac_hippocampus"))
