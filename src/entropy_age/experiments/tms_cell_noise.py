"""Per-cell transcriptional noise and the essential-gene shield: Tabula Muris Senis FACS (Smart-seq2).

The methylation shield is that housekeeping and core-essential promoters drift less with age than matched promoters.
The bulk RNA tests (``rna_tms``, ``rna_rodents``, ``rna_gtex``) asked whether between-animal expression spread
grows more slowly at those genes; this module asks the same of cell-to-cell variability inside one cell type of one
mouse, the noise a single cell sees.

**Data.** The official raw object of the TMS FACS release (``tabula-muris-senis-facs-official-raw-obj.h5ad``,
AWS bucket ``czb-tabula-muris-senis``, ``Data-objects/``; Tabula Muris Consortium 2020): full-length Smart-seq2 read
counts of sorted cells from C57BL/6JN mice at 3, 18, 21 and 24 months, with the consortium's cell ontology labels.
Young is 3 months, old is 18-24 months pooled. A tissue x cell type is analysed when each age group has at least
``MIN_CELLS`` cells from at least ``MIN_MICE`` mice (each mouse counted when it gives ``MIN_CELLS_PER_MOUSE``
cells). ERCC spike-ins are left out of genes and depth.

**Genes and sets.** Mouse symbols reach human genes by upper-case symbol through HGNC (``symbol_ensembl``, as in
``rna_rodents``); protein-coding genes only. Sets: Hart core-essential, Eisenberg-Levanon housekeeping and Hart
non-essential, through ``gene_set_ensembl``.

**Per-cell noise.** Per cell type every cell is binomially thinned to one depth, the ``DEPTH_QUANTILE`` quantile of
the cell type's read totals (cells below it are dropped), so that depth differences between cells and between ages
add no counting noise of their own. The larger age group is then subsampled at random to the size of the smaller.
Genes detected in at least ``MIN_DETECTED`` of cells in both age groups are kept. Per age group and gene: the mean
over cells, the variance pooled within mice (sum of squared deviations from each mouse's own mean over
``n_cells - n_mice``), so that differences between mice in mean level do not count as cell-to-cell noise, and
CV2 = variance / mean^2. A trend of log CV2 on log mean is fitted per cell type and age group as running medians in
``TREND_BINS`` quantile bins of log mean, interpolated linearly; a gene's noise is its residual from that trend
(the "distance to median" of Newman et al. 2006, with depth matched first). The age change is the old residual
minus the young one.

**Gap.** Each set's mean age change minus that of background genes (outside every set) in the same stratum, the
strata being deciles of mean expression crossed with tertiles of dropout rate (fraction of cells with zero reads,
within the decile), weighted by the set's genes per stratum. The gap is given in SDs of the age change over all
kept genes of the cell type (the full sample's SD for the point estimate and every resample, ``scaled_intervals``).

**Intervals.** Mice are resampled with replacement within each age group (one draw of the cohort's mice per
resample, shared by every cell type, since the same mice gave most tissues); thinning and strata are held fixed. With
3-6 mice per age the percentile interval is coarse. A gene bootstrap (genes resampled within strata) is given beside
it for every cell type. The across-cell-type mean gap takes its interval from the same mouse resamples.

**Between-mouse version (pseudobulk).** For each mouse with at least ``MIN_CELLS_PER_MOUSE`` depth-matched cells, raw
counts summed over those cells, log2(CPM + 1) over the kept genes. Per gene and age group the mean absolute deviation
from the group mean, times sqrt(n / (n - 1)); the age change is old minus young, the between-animal spread measure
of the bulk modules for two age groups. Gap and strata as above; with 3-6 mice per age a mouse resample would put
duplicate mice into a spread measure, so this version takes its interval from genes resampled within strata.

``uv run python -c`` with a ``Context`` (as ``entropy-age figures``) as a Slurm job; ``ENTROPY_AGE_BOOT`` sets the
resamples, ``ENTROPY_AGE_WORKERS`` the processes.
"""

from __future__ import annotations

import multiprocessing
import os
import sys
from concurrent.futures import ProcessPoolExecutor
from dataclasses import dataclass
from pathlib import Path

import h5py
import numpy as np
import polars as pl
from scipy import sparse

from entropy_age import plots
from entropy_age.data.annotation import AnnotationLayout, gene_set_ensembl, symbol_ensembl
from entropy_age.experiments.common import Context, write_table

H5AD = "tabula-muris-senis-facs-official-raw-obj.h5ad"
URL = f"https://czb-tabula-muris-senis.s3-us-west-2.amazonaws.com/Data-objects/{H5AD}"
SETS = ("core_essential", "housekeeping", "non_essential")
YOUNG = ("3m",)
OLD = ("18m", "21m", "24m")
MIN_CELLS = 50
MIN_MICE = 3
MIN_CELLS_PER_MOUSE = 10
MAX_CELLS = 2000  # per age group, after depth matching, to bound the resampling time
DEPTH_QUANTILE = 0.1
MIN_DETECTED = 0.1
TREND_BINS = 30
MEAN_BINS = 10
DROPOUT_BINS = 3


# ----------------------------------------------------------------------------- pure helpers
def thin_rows(x: sparse.csr_matrix, depth: int, rng: np.random.Generator) -> sparse.csr_matrix:
    """Each row binomially thinned to ``depth`` total counts in expectation; rows must hold at least ``depth``.

    A read kept with probability ``depth / total`` is what a shallower run of the same library would give.
    """
    totals = np.asarray(x.sum(axis=1)).ravel()
    if np.any(totals < depth):
        raise ValueError("thin_rows: a row holds fewer counts than the target depth")
    x = x.tocsr(copy=True)
    p = np.repeat(depth / totals, np.diff(x.indptr))
    x.data = rng.binomial(np.round(x.data).astype(np.int64), p).astype(np.float32)
    x.eliminate_zeros()
    return x


def within_moments(x: np.ndarray, mouse: np.ndarray) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Per column: mean over rows, variance pooled within ``mouse`` groups, and the fraction of zero rows.

    The pooled variance is the sum of squared deviations from each mouse's own mean over ``n - k`` (k mice);
    a mouse with one row adds no deviation and no degree of freedom.
    """
    labels, codes = np.unique(mouse, return_inverse=True)
    n, k = x.shape[0], labels.size
    if n - k < 1:
        raise ValueError("within_moments: no degrees of freedom left within mice")
    ss = np.zeros(x.shape[1])
    for m in range(k):
        rows = x[codes == m]
        ss += ((rows - rows.mean(axis=0)) ** 2).sum(axis=0)
    return x.mean(axis=0), ss / (n - k), (x == 0).mean(axis=0)


def trend_residual(log_mean: np.ndarray, log_cv2: np.ndarray, n_bins: int = TREND_BINS) -> np.ndarray:
    """``log_cv2`` minus a running-median trend on ``log_mean``: medians in quantile bins, interpolated linearly.

    Non-finite inputs give NaN and are left out of the trend.
    """
    ok = np.isfinite(log_mean) & np.isfinite(log_cv2)
    out = np.full(log_mean.shape, np.nan)
    if ok.sum() < 2 * n_bins:
        n_bins = max(1, ok.sum() // 10)
    edges = np.quantile(log_mean[ok], np.linspace(0, 1, n_bins + 1))
    which = np.clip(np.searchsorted(edges, log_mean[ok], side="right") - 1, 0, n_bins - 1)
    centres, medians = [], []
    for b in range(n_bins):
        sel = which == b
        if sel.any():
            centres.append(np.median(log_mean[ok][sel]))
            medians.append(np.median(log_cv2[ok][sel]))
    out[ok] = log_cv2[ok] - np.interp(log_mean[ok], centres, medians)
    return out


def noise_residual(x: np.ndarray, mouse: np.ndarray) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Residual log CV2 (pooled within mice) on its mean trend, with the mean and dropout rate per gene."""
    mean, var, dropout = within_moments(x, mouse)
    with np.errstate(divide="ignore", invalid="ignore"):
        log_mean = np.log(mean)
        log_cv2 = np.log(var / mean**2)
    return trend_residual(log_mean, log_cv2), mean, dropout


def strata(mean: np.ndarray, dropout: np.ndarray, n_mean: int = MEAN_BINS, n_dropout: int = DROPOUT_BINS) -> np.ndarray:
    """Stratum per gene: decile of mean expression crossed with tertile of dropout rate inside the decile."""
    rank = np.argsort(np.argsort(mean, kind="stable"), kind="stable")
    decile = rank * n_mean // mean.size
    out = np.empty(mean.size, dtype=np.int64)
    for d in range(n_mean):
        sel = np.flatnonzero(decile == d)
        r = np.argsort(np.argsort(dropout[sel], kind="stable"), kind="stable")
        out[sel] = d * n_dropout + r * n_dropout // max(sel.size, 1)
    return out


def matched_diff(values: np.ndarray, member: np.ndarray, background: np.ndarray, stratum: np.ndarray) -> float:
    """Mean of ``values`` over members minus the stratum-matched background mean, weighted by members per stratum.

    Members in a stratum without finite background values are left out; NaN values count nowhere.
    """
    ok = np.isfinite(values)
    num, weight = 0.0, 0
    for s in np.unique(stratum[member & ok]):
        m = member & ok & (stratum == s)
        b = background & ok & (stratum == s)
        if not b.any():
            continue
        num += m.sum() * (values[m].mean() - values[b].mean())
        weight += m.sum()
    return num / weight if weight else float("nan")


def gene_boot_diffs(
    values: np.ndarray, member: np.ndarray, background: np.ndarray, stratum: np.ndarray, n: int, seed: int
) -> np.ndarray:
    """``matched_diff`` over ``n`` resamples of genes drawn with replacement within each stratum and group."""
    rng = np.random.default_rng(seed)
    out = np.empty(n)
    idx_m = {s: np.flatnonzero(member & (stratum == s)) for s in np.unique(stratum)}
    idx_b = {s: np.flatnonzero(background & (stratum == s)) for s in np.unique(stratum)}
    for i in range(n):
        picks_v, picks_s, picks_m = [], [], []
        for s in idx_m:
            for idx, flag in ((idx_m[s], True), (idx_b[s], False)):
                if idx.size:
                    draw = rng.choice(idx, size=idx.size, replace=True)
                    picks_v.append(values[draw])
                    picks_s.append(np.full(draw.size, s))
                    picks_m.append(np.full(draw.size, flag))
        mem = np.concatenate(picks_m)
        out[i] = matched_diff(np.concatenate(picks_v), mem, ~mem, np.concatenate(picks_s))
    return out


def between_spread(logcpm: np.ndarray) -> np.ndarray:
    """Per column: mean absolute deviation from the column mean, times sqrt(n / (n - 1)) (rows are mice)."""
    n = logcpm.shape[0]
    if n < 2:
        raise ValueError("between_spread: needs at least two mice")
    return np.abs(logcpm - logcpm.mean(axis=0)).mean(axis=0) * np.sqrt(n / (n - 1))


def mouse_draw(mice: dict[str, list[str]], rng: np.random.Generator) -> dict[str, list[str]]:
    """Mice of each age group drawn with replacement, as many as the group has."""
    return {age: list(rng.choice(ids, size=len(ids), replace=True)) for age, ids in mice.items()}


# ----------------------------------------------------------------------------- reading the h5ad
def _decode(values) -> np.ndarray:
    arr = np.asarray(values)
    return arr.astype(str) if arr.dtype.kind in "SO" else arr


def _column(f: h5py.File, name: str) -> np.ndarray:
    obs = f["obs"]
    if isinstance(obs, h5py.Dataset):  # old anndata: a compound dataset, categories under uns/__categories
        values = obs[name]
        cats = f.get(f"uns/{name}_categories") or f.get(f"obs/__categories/{name}")
        if cats is not None and values.dtype.kind in "iu":
            return np.where(values >= 0, _decode(cats[()])[np.maximum(values, 0)], "nan")
        return _decode(values)
    item = obs[name]
    if isinstance(item, h5py.Group):  # anndata >= 0.8 categorical
        return _decode(item["categories"][()])[item["codes"][()]]
    cats = obs.get(f"__categories/{name}")
    if cats is not None:
        return _decode(cats[()])[item[()]]
    return _decode(item[()])


def _var_names(f: h5py.File) -> np.ndarray:
    var = f["var"]
    if isinstance(var, h5py.Dataset):
        return _decode(var["index"])
    key = var.attrs.get("_index", "index")
    return _decode(var[key][()])


def read_obs(path: Path) -> tuple[pl.DataFrame, np.ndarray]:
    """Cells (``row``, ``age``, ``mouse``, ``tissue``, ``cell_type``, ``sex``) and gene symbols of the h5ad."""
    with h5py.File(path, "r") as f:
        cols = {
            "age": "age",
            "mouse": "mouse.id",
            "tissue": "tissue",
            "cell_type": "cell_ontology_class",
            "sex": "sex",
        }
        frame = pl.DataFrame({k: _column(f, v).astype(str) for k, v in cols.items()}).with_row_index("row")
        genes = _var_names(f)
    return frame, genes


def read_rows(path: Path, rows: np.ndarray, n_genes: int) -> sparse.csr_matrix:
    """The given rows (sorted) of the CSR count matrix ``X``."""
    with h5py.File(path, "r") as f:
        x = f["X"]
        indptr = x["indptr"][()]
        data, indices, ptr = [], [], [0]
        for r in rows:
            a, b = indptr[r], indptr[r + 1]
            data.append(x["data"][a:b])
            indices.append(x["indices"][a:b])
            ptr.append(ptr[-1] + (b - a))
    return sparse.csr_matrix(
        (np.concatenate(data).astype(np.float32), np.concatenate(indices), np.array(ptr)), shape=(rows.size, n_genes)
    )


# ----------------------------------------------------------------------------- one cell type
@dataclass
class Task:
    path: Path
    tissue: str
    cell_type: str
    cells: pl.DataFrame  # row, age_group, mouse
    gene_index: np.ndarray  # columns of X kept (protein-coding, mapped, not ERCC)
    genes: np.ndarray
    sets: dict[str, np.ndarray]  # set -> boolean over genes
    draws: list[dict[str, list[str]]]
    n_gene_boot: int
    seed: int


def _sc_delta(x, cells) -> np.ndarray:
    """Old minus young residual log CV2 per gene (columns of ``x``)."""
    young = (cells["age_group"] == "young").to_numpy()
    r_y, _, _ = noise_residual(x[young], cells["mouse"].to_numpy()[young])
    r_o, _, _ = noise_residual(x[~young], cells["mouse"].to_numpy()[~young])
    return r_o - r_y


def analyse(task: Task) -> tuple[list[dict], list[dict], dict]:
    """Point gaps (both versions), mouse-resample gaps and a cell-type summary row."""
    print(
        f"[tms_cell_noise] start {task.tissue} / {task.cell_type}: {task.cells.height} cells",
        file=sys.stderr,
        flush=True,
    )
    rng = np.random.default_rng(task.seed)
    rows = task.cells["row"].to_numpy()
    order = np.argsort(rows)
    cells = task.cells[order.tolist()]
    raw = read_rows(task.path, rows[order], _n_genes(task.path))
    raw = raw[:, task.gene_index].tocsr()
    totals = np.asarray(raw.sum(axis=1)).ravel()
    depth = int(np.quantile(totals, DEPTH_QUANTILE))
    ok = totals >= depth
    raw, cells, totals = raw[ok], cells.filter(pl.Series(ok)), totals[ok]
    # equal cells per age group
    young = (cells["age_group"] == "young").to_numpy()
    n_eq = min(young.sum(), (~young).sum(), MAX_CELLS)
    pick = np.sort(
        np.concatenate(
            [
                rng.choice(np.flatnonzero(young), n_eq, replace=False),
                rng.choice(np.flatnonzero(~young), n_eq, replace=False),
            ]
        )
    )
    raw, cells = raw[pick], cells[pick.tolist()]
    thinned = thin_rows(raw, depth, rng).toarray()
    young = (cells["age_group"] == "young").to_numpy()
    detected = np.minimum((thinned[young] > 0).mean(axis=0), (thinned[~young] > 0).mean(axis=0))
    keep = detected >= MIN_DETECTED
    x = thinned[:, keep]
    del thinned
    sets = {s: m[keep] for s, m in task.sets.items()}
    background = ~np.any(np.column_stack(list(sets.values())), axis=1)
    stratum = strata(x.mean(axis=0), (x == 0).mean(axis=0))
    delta = _sc_delta(x, cells)
    sd = float(np.nanstd(delta))
    # young level of noise, for context
    r_y, _, _ = noise_residual(x[young], cells["mouse"].to_numpy()[young])
    point = []
    for s, member in sets.items():
        d = matched_diff(delta, member, background, stratum)
        gb = gene_boot_diffs(delta, member, background, stratum, task.n_gene_boot, task.seed + 1)
        point.append(
            {
                "version": "single_cell",
                "set": s,
                "n_set": int(member.sum()),
                "matched_diff": d,
                "sd": sd,
                "gene_low": float(np.nanquantile(gb, 0.025)),
                "gene_high": float(np.nanquantile(gb, 0.975)),
                "young_level_gap_sd": matched_diff(r_y, member, background, stratum) / float(np.nanstd(r_y)),
            }
        )
    # mouse resamples
    mouse = cells["mouse"].to_numpy()
    boot = []
    for b, draw in enumerate(task.draws):
        idx, labels = [], []
        for age, ids in draw.items():
            for j, m in enumerate(ids):
                sel = np.flatnonzero(mouse == m)
                idx.append(sel)
                labels.append(np.full(sel.size, f"{age}:{j}"))
        idx = np.concatenate(idx) if idx else np.array([], int)
        sub = cells[idx.tolist()].with_columns(mouse=pl.Series(np.concatenate(labels) if labels else []))
        yb = (sub["age_group"] == "young").to_numpy()
        n_mice = [sub.filter(pl.Series(g))["mouse"].n_unique() for g in (yb, ~yb)]
        if min(n_mice) < 2 or min(yb.sum(), (~yb).sum()) - max(n_mice) < 5:
            continue
        with np.errstate(invalid="ignore", divide="ignore"):
            dlt = _sc_delta(x[idx], sub)
        for s, member in sets.items():
            boot.append({"boot": b, "set": s, "matched_diff": matched_diff(dlt, member, background, stratum)})
    # pseudobulk
    pb_rows = _pseudobulk(task, raw, cells, keep, sets, background, stratum)
    point += pb_rows
    for r in point:
        r |= {"tissue": task.tissue, "cell_type": task.cell_type}
    for r in boot:
        r |= {"tissue": task.tissue, "cell_type": task.cell_type, "version": "single_cell"}
    cy = cells.filter(pl.Series(young))
    co = cells.filter(pl.Series(~young))
    summary = {
        "tissue": task.tissue,
        "cell_type": task.cell_type,
        "depth": depth,
        "cells_young": cy.height,
        "cells_old": co.height,
        "mice_young": cy["mouse"].n_unique(),
        "mice_old": co["mouse"].n_unique(),
        "genes": int(keep.sum()),
        "mouse_resamples": len({r["boot"] for r in boot}),
    } | {f"n_{s}": int(m.sum()) for s, m in sets.items()}
    print(f"[tms_cell_noise] {task.tissue} / {task.cell_type}: {summary}", file=sys.stderr, flush=True)
    return point, boot, summary


def _n_genes(path: Path) -> int:
    with h5py.File(path, "r") as f:
        x = f["X"]
        shape = x.attrs.get("shape", None)
        if shape is None:
            shape = x.attrs["h5sparse_shape"]
        return int(shape[1])


def _pseudobulk(task, raw, cells, keep, sets, background, stratum) -> list[dict]:
    """Between-mouse spread of pseudobulk log CPM; its age change and matched gaps with a gene bootstrap."""
    mouse = cells["mouse"].to_numpy()
    young = (cells["age_group"] == "young").to_numpy()
    spreads = {}
    for name, sel in (("young", young), ("old", ~young)):
        mats = []
        for m in np.unique(mouse[sel]):
            r = sel & (mouse == m)
            if r.sum() >= MIN_CELLS_PER_MOUSE:
                mats.append(np.asarray(raw[r].sum(axis=0)).ravel())
        if len(mats) < MIN_MICE:
            return []
        counts = np.vstack(mats)
        cpm = counts / counts.sum(axis=1, keepdims=True) * 1e6
        spreads[name] = (between_spread(np.log2(cpm[:, keep] + 1)), len(mats))
    delta = spreads["old"][0] - spreads["young"][0]
    sd = float(np.std(delta))
    out = []
    for s, member in sets.items():
        gb = gene_boot_diffs(delta, member, background, stratum, task.n_gene_boot, task.seed + 2)
        out.append(
            {
                "version": "pseudobulk",
                "set": s,
                "n_set": int(member.sum()),
                "matched_diff": matched_diff(delta, member, background, stratum),
                "sd": sd,
                "gene_low": float(np.nanquantile(gb, 0.025)),
                "gene_high": float(np.nanquantile(gb, 0.975)),
                "pb_mice_young": spreads["young"][1],
                "pb_mice_old": spreads["old"][1],
            }
        )
    return out


# ----------------------------------------------------------------------------- run
def eligible(obs: pl.DataFrame) -> pl.DataFrame:
    """Tissue x cell types with ``MIN_CELLS`` cells per age group from ``MIN_MICE`` mice of ``MIN_CELLS_PER_MOUSE``."""
    per_mouse = obs.group_by("tissue", "cell_type", "age_group", "mouse").len()
    mice = (
        per_mouse.filter(pl.col("len") >= MIN_CELLS_PER_MOUSE)
        .group_by("tissue", "cell_type", "age_group")
        .agg(pl.len().alias("mice"))
    )
    cells = obs.group_by("tissue", "cell_type", "age_group").agg(pl.len().alias("cells"))
    both = cells.join(mice, on=["tissue", "cell_type", "age_group"], how="left").fill_null(0)
    ok = both.filter((pl.col("cells") >= MIN_CELLS) & (pl.col("mice") >= MIN_MICE))
    return (
        ok.group_by("tissue", "cell_type")
        .agg(pl.len().alias("k"))
        .filter(pl.col("k") == 2)
        .drop("k")
        .sort("tissue", "cell_type")
    )


def _settings() -> tuple[int, int]:
    return int(os.environ.get("ENTROPY_AGE_WORKERS", "1")), int(os.environ.get("ENTROPY_AGE_BOOT", "200"))


def run(ctx: Context) -> list[Path]:
    path = ctx.paths.input / "tms_facs" / H5AD
    annotation = AnnotationLayout(ctx.paths)
    needed = [path, annotation.gene_sets, annotation.hgnc]
    if any(not p.exists() for p in needed):
        print(f"tms_cell_noise: inputs missing ({', '.join(str(p) for p in needed if not p.exists())}); {URL}")
        return []
    plots.apply_style()
    out = ctx.out_dir("tms_cell_noise")
    workers, n_boot = _settings()
    obs, symbols = read_obs(path)
    # cells of two pooled mice ("3_10_M/3_11_M") belong to neither mouse and are left out
    obs = obs.filter(
        pl.col("age").is_in(YOUNG + OLD)
        & ~pl.col("cell_type").is_in(["nan", "", "None"])
        & ~pl.col("mouse").str.contains("/", literal=True)
    )
    obs = obs.with_columns(age_group=pl.when(pl.col("age").is_in(YOUNG)).then(pl.lit("young")).otherwise(pl.lit("old")))
    hgnc, gene_sets = pl.read_parquet(annotation.hgnc), pl.read_parquet(annotation.gene_sets)
    mapping = symbol_ensembl([str(s) for s in symbols], hgnc)
    protein = set(hgnc.filter(pl.col("protein_coding") & pl.col("ensembl").is_not_null())["ensembl"].to_list())
    ensembl = dict(mapping.select("symbol", "ensembl").iter_rows())
    coding = np.array([ensembl.get(str(g)) in protein and not str(g).startswith("ERCC-") for g in symbols])
    gene_index = np.flatnonzero(coding)
    genes = symbols[gene_index]
    membership = (
        gene_set_ensembl(gene_sets, hgnc)
        .filter(pl.col("group").is_in(SETS))
        .join(mapping.select(pl.col("ensembl").alias("feature"), pl.col("symbol").alias("mouse")), on="feature")
    )
    sets = {s: np.isin(genes, membership.filter(pl.col("group") == s)["mouse"].to_list()) for s in SETS}
    groups = eligible(obs)
    if limit := int(os.environ.get("ENTROPY_AGE_TMS_CELL_TYPES", "0")):  # smoke tests
        groups = groups.head(limit)
    mice = {g: sorted(obs.filter(pl.col("age_group") == g)["mouse"].unique().to_list()) for g in ("young", "old")}
    rng = ctx.rng(4242)
    draws = [mouse_draw(mice, rng) for _ in range(n_boot)]
    tasks = []
    for i, (tissue, cell_type) in enumerate(groups.iter_rows()):
        cells = obs.filter((pl.col("tissue") == tissue) & (pl.col("cell_type") == cell_type)).select(
            "row", "age_group", "mouse"
        )
        tasks.append(Task(path, tissue, cell_type, cells, gene_index, genes, sets, draws, n_boot, 1000 + i))
    print(f"[tms_cell_noise] {len(tasks)} cell types, {len(genes)} genes, {n_boot} resamples", file=sys.stderr)
    if workers > 1:
        # spawn: a forked child inherits polars' locked thread pool and hangs (smoke job 6657459, 2026-11-26)
        with ProcessPoolExecutor(workers, mp_context=multiprocessing.get_context("spawn")) as pool:
            results = list(pool.map(analyse, tasks))
    else:
        results = [analyse(t) for t in tasks]
    point = pl.DataFrame([r for p, _, _ in results for r in p], infer_schema_length=None)
    boot = pl.DataFrame([r for _, b, _ in results for r in b], infer_schema_length=None)
    cell_types = pl.DataFrame([s for _, _, s in results], infer_schema_length=None)
    gaps = _intervals(point, boot)
    pooled = _pooled(gaps, boot)
    write_table(cell_types, out / "cell_types")
    write_table(gaps, out / "gaps")
    write_table(pooled, out / "pooled")
    boot.write_parquet(out / "mouse_resamples.parquet")
    write_table(
        pl.DataFrame(
            {
                "symbols": [len(symbols)],
                "protein_coding_mapped": [len(genes)],
                **{f"set_{s}": [int(m.sum())] for s, m in sets.items()},
                "mouse_resamples": [n_boot],
                "mice_young": [len(mice["young"])],
                "mice_old": [len(mice["old"])],
            }
        ),
        out / "genes",
    )
    return [_figure(ctx, gaps, pooled)]


def _intervals(point: pl.DataFrame, boot: pl.DataFrame) -> pl.DataFrame:
    keys = ["tissue", "cell_type", "version", "set"]
    q = boot.group_by(keys).agg(
        pl.col("matched_diff").quantile(0.025).alias("mouse_low"),
        pl.col("matched_diff").quantile(0.975).alias("mouse_high"),
    )
    return (
        point.join(q, on=keys, how="left")
        .with_columns(
            gap_sd=pl.col("matched_diff") / pl.col("sd"),
            mouse_ci_low=pl.col("mouse_low") / pl.col("sd"),
            mouse_ci_high=pl.col("mouse_high") / pl.col("sd"),
            gene_ci_low=pl.col("gene_low") / pl.col("sd"),
            gene_ci_high=pl.col("gene_high") / pl.col("sd"),
        )
        .drop("mouse_low", "mouse_high", "gene_low", "gene_high")
        .sort("version", "set", "tissue", "cell_type")
    )


def _pooled(gaps: pl.DataFrame, boot: pl.DataFrame) -> pl.DataFrame:
    """Mean gap over cell types; single-cell interval from the shared mouse resamples (complete resamples only)."""
    rows = []
    for (version, s), sel in gaps.group_by("version", "set", maintain_order=True):
        row = {
            "version": version,
            "set": s,
            "cell_types": sel.height,
            "mean_gap_sd": float(sel["gap_sd"].mean()),
            "cell_types_below_zero": int((sel["gap_sd"] < 0).sum()),
        }
        if version == "single_cell":
            b = (
                boot.filter(pl.col("set") == s)
                .join(sel.select("tissue", "cell_type", "sd"), on=["tissue", "cell_type"])
                .with_columns(g=pl.col("matched_diff") / pl.col("sd"))
                .group_by("boot")
                .agg(pl.col("g").mean(), pl.len().alias("n"))
                .filter(pl.col("n") == sel.height)
            )
            if b.height:
                row |= {
                    "ci_low": float(b["g"].quantile(0.025)),
                    "ci_high": float(b["g"].quantile(0.975)),
                    "complete_resamples": b.height,
                }
        rows.append(row)
    return pl.DataFrame(rows, infer_schema_length=None)


def _figure(ctx: Context, gaps: pl.DataFrame, pooled: pl.DataFrame) -> Path:
    fig, axes = plots.panels(2, ncols=2, size=6.0)
    sel = gaps.filter(pl.col("set").is_in(["housekeeping", "core_essential"]))
    labels = sel.select("tissue", "cell_type").unique().sort("tissue", "cell_type")
    names = [f"{t}: {c}" for t, c in labels.iter_rows()]
    pos = {k: i for i, k in enumerate(labels.iter_rows())}
    for ax, version, title, lo, hi in (
        (axes[0], "single_cell", "cell-to-cell noise (mouse resamples)", "mouse_ci_low", "mouse_ci_high"),
        (axes[1], "pseudobulk", "between-mouse spread (gene resamples)", "gene_ci_low", "gene_ci_high"),
    ):
        for k, (s, colour) in enumerate((("housekeeping", plots.PALETTE[0]), ("core_essential", plots.PALETTE[1]))):
            part = sel.filter((pl.col("version") == version) & (pl.col("set") == s))
            y = np.array([pos[(t, c)] for t, c in part.select("tissue", "cell_type").iter_rows()]) + (k - 0.5) * 0.3
            x = part["gap_sd"].to_numpy()
            # a resampling interval need not contain the point estimate (few mice): clip the bar lengths at zero
            err = np.vstack([x - part[lo].fill_null(np.nan).to_numpy(), part[hi].fill_null(np.nan).to_numpy() - x])
            err = np.clip(np.nan_to_num(err, nan=0.0), 0, None)
            ax.errorbar(x, y, xerr=err, fmt="o", ms=3, color=colour, lw=0.8, label=s.replace("_", " "))
        ax.axvline(0, color=plots.MUTED, lw=0.8)
        ax.set_yticks(range(len(names)), names, fontsize=5.5)
        ax.invert_yaxis()
        ax.set_xlabel("gap in age change, SD (set minus matched genes)")
        ax.set_title(title)
    axes[0].legend(loc="lower left")
    axes[1].set_yticklabels([])
    for ax, letter in zip(axes, "ab", strict=True):
        plots.label_panel(ax, letter)
    fig.set_size_inches(12, max(4, 0.2 * len(names) + 1.5))
    fig.suptitle("Tabula Muris Senis FACS: essential-gene shield in per-cell noise (3 m vs 18-24 m)", x=0.02, ha="left")
    return plots.save(fig, ctx.figure_path("tms_cell_noise"))
