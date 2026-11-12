"""Accessibility spread with age per cell type, in 357 human prefrontal cortices (roadmap RM19).

The PBMC test (``atac_variance``) found essential promoters spreading apart between people more slowly with age,
once cell mix was removed. Bulk blood mixes cell types; this repeats the test inside each cell type. Catching et
al. 2026 (Cell Reports; Zenodo 18394349, CC-BY 4.0) profiled single-nucleus ATAC and RNA in the prefrontal cortex
of donors aged 15-100 from four brain banks: 1.5 million nuclei, 521,217 peaks (hg38), seven cell types.

1. Pseudobulk (``pseudobulk``, cached in ``data/interim/atac_cortex``): counts summed per donor and cell type over
   the peaks that hold a 450k promoter CpG, streamed from the 14.5 GB h5ad in row blocks, with each group's total
   counts and nucleus count.
2. Per cell type, donors with at least ``MIN_NUCLEI`` nuclei: log counts per million, residuals after age, sex,
   brain bank, post-mortem interval and log nucleus count (and, per ``ADJUSTMENTS``, leading components and the
   promoter read fraction), and the age slope of the absolute residual in units of the peak's spread
   (``atac_variance.spread_slopes``). A set's peaks
   are compared with the other promoter peaks of the same island status, mean decile and width tertile.

Run as a job: ``entropy-age figures atac_cortex`` (the pseudobulk reads about 67 GB once).
"""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import polars as pl

from entropy_age import plots
from entropy_age.data.annotation import AnnotationLayout, link_membership
from entropy_age.experiments.atac_variance import HG38_MANIFEST, peak_cpgs, spread_slopes
from entropy_age.experiments.common import Context, dummies, write_table
from entropy_age.experiments.tissue_contrast import SET_LABELS, SETS
from entropy_age.stochasticity import baseline_strata, compare_groups

RECORD = "18394349"
H5AD = "final_atac_data.h5ad"
DEMOGRAPHICS = "TableS1_CohortDemographics.xlsx"
MIN_NUCLEI = 50
#: name -> (leading components removed, promoter read fraction as a covariate). The promoter fraction is a
#: donor's share of reads in promoter peaks, a quality score like TSS enrichment: if it varies with age it widens
#: the spread exactly at the strongest promoters, which essential genes own.
ADJUSTMENTS = {"none": (0, False), "5 components": (5, False), "promoter fraction": (0, True), "both": (5, True)}
PERMUTATIONS = 2000
BLOCK = 20_000


def inputs(paths) -> tuple[Path, Path]:
    base = paths.input / "zenodo" / RECORD
    return base / H5AD, base / DEMOGRAPHICS


def parse_peaks(names: list[str]) -> pl.DataFrame:
    """``chr:start-end`` names as ``chr``, ``start``, ``end``."""
    chrom, rest = zip(*(n.split(":") for n in names), strict=True)
    start, end = zip(*(r.split("-") for r in rest), strict=True)
    return pl.DataFrame({"chr": list(chrom), "start": [int(s) for s in start], "end": [int(e) for e in end]})


def donor_id(sample: object, bank: str) -> str:
    """The h5ad's sample id for a demographics row: HBCC donors are bare numbers there (887 -> HBCC-887)."""
    if bank == "HBCC":
        return f"HBCC-{int(float(sample))}"
    return str(sample)


def read_demographics(path: Path) -> pl.DataFrame:
    import openpyxl

    rows = [r for r in openpyxl.load_workbook(path, read_only=True).active.iter_rows(values_only=True)]
    header, body = rows[0], [r for r in rows[1:] if r and r[0] is not None]
    frame = pl.DataFrame([dict(zip(header, r, strict=False)) for r in body], infer_schema_length=None)
    return frame.select(
        pl.struct("SampleID", "BrainBank")
        .map_elements(lambda s: donor_id(s["SampleID"], s["BrainBank"]), return_dtype=pl.String)
        .alias("sample"),
        pl.col("Age").cast(pl.Float64).alias("age"),
        pl.col("Sex").alias("sex"),
        pl.col("PMI").cast(pl.Float64, strict=False).alias("pmi"),
        pl.col("BrainBank").alias("bank"),
    )


def group_sums(indptr, indices, data, groups, keep_col, n_groups: int, n_keep: int):
    """Counts summed per group over the kept columns, and every group's total counts, for one CSR row block."""
    rows = np.repeat(np.arange(indptr.size - 1), np.diff(indptr))
    g = groups[rows]
    totals = np.bincount(g, weights=data, minlength=n_groups)
    col = keep_col[indices]
    hit = col >= 0
    flat = g[hit].astype(np.int64) * n_keep + col[hit]
    sums = np.bincount(flat, weights=data[hit], minlength=n_groups * n_keep).reshape(n_groups, n_keep)
    return sums, totals


def pseudobulk(h5ad: Path, keep: np.ndarray, out: Path) -> Path:
    """Write ``groups.parquet`` (sample, cell_type, nuclei, total) and ``counts.npy`` (groups x kept peaks)."""
    import h5py

    with h5py.File(h5ad, "r") as f:
        obs = f["obs"]
        samples = obs["sample_id"]["codes"][:]
        cells = obs["cell_type"]["codes"][:]
        sample_names = [x.decode() for x in obs["sample_id"]["categories"][:]]
        cell_names = [x.decode() for x in obs["cell_type"]["categories"][:]]
        n_cell = len(cell_names)
        groups = samples.astype(np.int64) * n_cell + cells
        n_groups = len(sample_names) * n_cell
        n_peaks = f["var"]["_index"].shape[0]
        keep_col = np.full(n_peaks, -1, dtype=np.int64)
        keep_col[keep] = np.arange(keep.size)
        indptr = f["X"]["indptr"][:]
        sums = np.zeros((n_groups, keep.size))
        totals = np.zeros(n_groups)
        for lo in range(0, indptr.size - 1, BLOCK):
            hi = min(lo + BLOCK, indptr.size - 1)
            a, b = indptr[lo], indptr[hi]
            s, t = group_sums(
                indptr[lo : hi + 1] - a,
                f["X"]["indices"][a:b],
                f["X"]["data"][a:b].astype(np.float64),
                groups[lo:hi],
                keep_col,
                n_groups,
                keep.size,
            )
            sums += s
            totals += t
            if lo // BLOCK % 10 == 0:
                print(f"atac_cortex: {hi:,} of {indptr.size - 1:,} nuclei", file=sys.stderr, flush=True)
    nuclei = np.bincount(groups, minlength=n_groups)
    table = pl.DataFrame(
        {
            "sample": [s for s in sample_names for _ in cell_names],
            "cell_type": cell_names * len(sample_names),
            "nuclei": nuclei,
            "total": totals,
        }
    )
    out.mkdir(parents=True, exist_ok=True)
    present = nuclei > 0
    table.filter(pl.Series(present)).write_parquet(out / "groups.parquet")
    np.save(out / "counts.npy", sums[present].astype(np.float32))
    return out


def promoter_fraction(counts: np.ndarray, total: np.ndarray) -> np.ndarray:
    """Each group's share of all its reads that fall in the kept promoter peaks."""
    return counts.sum(axis=1) / total


def run(ctx: Context) -> list[Path]:
    h5ad, demographics = inputs(ctx.paths)
    ann = AnnotationLayout(ctx.paths)
    manifest_path = ann.input / HG38_MANIFEST
    missing = [str(p) for p in (h5ad, demographics, manifest_path, ann.links_450k) if not p.exists()]
    if missing:
        print(f"atac_cortex: inputs missing ({', '.join(missing)}). Nothing run.", file=sys.stderr)
        return []
    import h5py

    plots.apply_style()
    out = ctx.out_dir("atac_cortex")
    with h5py.File(h5ad, "r") as f:
        peaks = parse_peaks([x.decode() for x in f["var"]["_index"][:]])
    manifest = (
        pl.read_csv(manifest_path, separator="\t", columns=["CpG_chrm", "CpG_beg", "Probe_ID"], null_values=["NA"])
        .drop_nulls()
        .select(pl.col("CpG_chrm").alias("chr"), pl.col("CpG_beg").alias("pos"), pl.col("Probe_ID").alias("cpg"))
    )
    probes = pl.read_parquet(ann.probes_450k).select("cpg", "cpg_island", "promoter")
    promoter = peak_cpgs(peaks, manifest).join(probes, on="cpg").filter(pl.col("promoter"))
    keep = np.sort(promoter["peak"].unique().to_numpy())
    cache = ctx.paths.interim / "atac_cortex"
    if not (cache / "counts.npy").exists():
        pseudobulk(h5ad, keep, cache)
    groups = pl.read_parquet(cache / "groups.parquet")
    counts = np.load(cache / "counts.npy")
    if counts.shape[1] != keep.size:
        raise ValueError("atac_cortex: the cached pseudobulk has other peaks; delete data/interim/atac_cortex")

    cpg_sets = link_membership(pl.read_parquet(ann.links_450k), pl.read_parquet(ann.gene_sets)).filter(
        pl.col("group").is_in(SETS)
    )
    membership = (
        promoter.join(cpg_sets.rename({"feature": "cpg"}), on="cpg")
        .select(pl.col("peak").cast(pl.String).alias("feature"), "group", "gene")
        .unique()
    )
    island = promoter.group_by("peak").agg(pl.col("cpg_island").any())
    island = dict(zip(island["peak"].to_list(), island["cpg_island"].to_list(), strict=True))
    width = (peaks["end"] - peaks["start"]).to_numpy()[keep]
    demo = read_demographics(demographics)
    rows = []
    for cell_type in sorted(groups["cell_type"].unique().to_list()):
        sel = (
            groups.with_row_index("row")
            .filter((pl.col("cell_type") == cell_type) & (pl.col("nuclei") >= MIN_NUCLEI))
            .join(demo, on="sample")
            .drop_nulls(["age", "sex", "pmi", "bank"])
        )
        if sel.height < 30:
            print(f"atac_cortex: {cell_type}: {sel.height} donors; left out", file=sys.stderr)
            continue
        x = counts[sel["row"].to_numpy()].astype(np.float64)
        y = np.log2(x / sel["total"].to_numpy()[:, None] * 1e6 + 1.0)
        age = sel["age"].to_numpy()
        base_cov = np.column_stack(
            [
                dummies(sel["sex"].to_numpy()),
                dummies(sel["bank"].to_numpy()),
                sel["pmi"].to_numpy(),
                np.log(sel["nuclei"].to_numpy()),
            ]
        )
        mean = y.mean(axis=0)
        stratum = pl.Series(
            [
                f"{island[int(p)]}|{m}|{w}"
                for p, m, w in zip(keep, baseline_strata(mean, 10), baseline_strata(width, 3), strict=True)
            ]
        )
        quality = promoter_fraction(x, sel["total"].to_numpy())
        for adjustment, (k, by_quality) in ADJUSTMENTS.items():
            cov = np.column_stack([base_cov, quality]) if by_quality else base_cov
            if k:
                u, s, _ = np.linalg.svd(y - mean, full_matrices=False)
                cov = np.column_stack([cov, u[:, :k] * s[:k]])
            slopes = spread_slopes(y, age, cov)
            metrics = pl.DataFrame({"feature": [str(p) for p in keep], "spread_slope": slopes, "stratum": stratum})
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
                    cell_type=pl.lit(cell_type),
                    adjustment=pl.lit(adjustment),
                    components=pl.lit(k),
                    gap_sd=pl.col("matched_diff") / sd,
                    donors=pl.lit(sel.height),
                    age_min=pl.lit(float(age.min())),
                    age_max=pl.lit(float(age.max())),
                )
            )
        print(f"atac_cortex: {cell_type}: {sel.height} donors done", file=sys.stderr, flush=True)
    result = pl.concat(rows, how="diagonal")
    write_table(result, out / "gaps")
    return [_figure(ctx, result)]


def _figure(ctx: Context, result: pl.DataFrame) -> Path:
    fig, axes = plots.panels(len(ADJUSTMENTS), ncols=2, size=4.6)
    for i, adjustment in enumerate(ADJUSTMENTS):
        ax = axes[i]
        sub = result.filter(pl.col("adjustment") == adjustment)
        cells = sorted(sub["cell_type"].unique().to_list())
        for g, group in enumerate(SETS):
            r = sub.filter(pl.col("group") == group).sort("cell_type")
            ys = [cells.index(c) + 0.2 * g for c in r["cell_type"].to_list()]
            ax.scatter(r["gap_sd"], ys, color=plots.PALETTE[g], s=22, label=SET_LABELS[group])
            for x, yv, p in zip(r["gap_sd"], ys, r["p_permutation"], strict=True):
                if p < 0.05:
                    ax.annotate("*", (x, yv), xytext=(4, -3), textcoords="offset points", fontsize=9)
        donors = dict(zip(sub["cell_type"], sub["donors"], strict=False))
        labels = [f"{c} ({donors[c]})" for c in cells] if i % 2 == 0 else [""] * len(cells)
        ax.set_yticks([j + 0.2 for j in range(len(cells))], labels)
        ax.axvline(0, color=plots.MUTED, lw=0.8)
        ax.invert_yaxis()
        if i == 0:
            ax.legend(fontsize=7, frameon=False, loc="lower left")
        ax.set(
            xlabel="Set minus matched peaks (SD)",
            title=f"Adjusted for: {adjustment}",
        )
        plots.label_panel(ax, "abcd"[i])
    fig.suptitle(
        "Accessibility spread with age, prefrontal cortex cell types "
        "(> 0: essential promoters spread more; * p < 0.05)",
        x=0.02,
        y=1.03,
        ha="left",
        fontsize=11,
        fontweight="bold",
    )
    return plots.save(fig, ctx.figure_path("atac_cortex"))
