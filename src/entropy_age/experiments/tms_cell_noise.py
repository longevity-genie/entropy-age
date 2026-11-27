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

**Sensitivities (Codex review 8; ``run_sensitivity``, tables ``sensitivity_*``, ``lomo``, ``planted``,
``contrasts``, ``orthology``).** Every one reports the single-cell gap pooled over its estimable cell types (finite
point gap), with an interval from shared Bayesian-bootstrap mouse weights (Dirichlet within age group, every mouse
kept, so no cell type drops out of a resample). (a) The main multinomial mouse draws re-gated on distinct original
mouse IDs (a mouse drawn twice no longer counts as two). (b) Leave one mouse out, each of the cohort's mice in turn.
(c) Planted effects: in old cells, housekeeping and core-essential genes get multiplicative Gamma noise raising their
within-mouse log CV2 by 0.05-0.5 SD of the age change, cells and mice as observed; the recovered gap and its
interval calibrate what the estimator can see. (d) 3 vs 18 months, and same-sex contrasts where three mice per age
allow it. (e) Depth matched on all endogenous reads (thinning rate from all genes). (f) Gene sets through MGI
one-to-one orthologs (human Entrez to HGNC), and the subset where both routes name the same human gene. (g) Young
reference: genes kept on young detection and matched on young mean and dropout, ties kept in one bin; matched set
sizes and covariate balance reported. Main estimates keep their definition; pooling now excludes cell types with no
finite gap (no retained non-essential genes) instead of returning NaN.

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
from scipy.stats import rankdata

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
def thin_rows(
    x: sparse.csr_matrix, depth: int, rng: np.random.Generator, totals: np.ndarray | None = None
) -> sparse.csr_matrix:
    """Each row binomially thinned to ``depth`` total counts in expectation; rows must hold at least ``depth``.

    A read kept with probability ``depth / total`` is what a shallower run of the same library would give.
    ``totals`` (per row) defaults to the row sums of ``x``; given, it is the depth over a wider gene set (all
    endogenous reads) of which ``x`` holds some columns, and the columns are thinned at that set's rate.
    """
    if totals is None:
        totals = np.asarray(x.sum(axis=1)).ravel()
    if np.any(totals < depth):
        raise ValueError("thin_rows: a row holds fewer counts than the target depth")
    x = x.tocsr(copy=True)
    p = np.repeat(depth / totals, np.diff(x.indptr))
    x.data = rng.binomial(np.round(x.data).astype(np.int64), p).astype(np.float32)
    x.eliminate_zeros()
    return x


def within_moments(
    x: np.ndarray, mouse: np.ndarray, weights: dict[str, float] | None = None
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Per column: mean over rows, variance pooled within ``mouse`` groups, and the fraction of zero rows.

    The pooled variance is the sum of squared deviations from each mouse's own mean over ``n - k`` (k mice);
    a mouse with one row adds no deviation and no degree of freedom. ``weights`` (mouse -> weight, a shared
    mouse-weight resample) weight each mouse's cells, sums of squares and degrees of freedom.
    """
    labels, codes = np.unique(mouse, return_inverse=True)
    n, k = x.shape[0], labels.size
    if n - k < 1:
        raise ValueError("within_moments: no degrees of freedom left within mice")
    if weights is not None:
        total, ss, zeros = np.zeros(x.shape[1]), np.zeros(x.shape[1]), np.zeros(x.shape[1])
        n_w = df_w = 0.0
        for m, label in enumerate(labels):
            w = float(weights.get(str(label), 0.0))
            rows = x[codes == m]
            total += w * rows.sum(axis=0)
            ss += w * ((rows - rows.mean(axis=0)) ** 2).sum(axis=0)
            zeros += w * (rows == 0).sum(axis=0)
            n_w += w * rows.shape[0]
            df_w += w * (rows.shape[0] - 1)
        if df_w <= 0:
            raise ValueError("within_moments: no weighted degrees of freedom")
        return total / n_w, ss / df_w, zeros / n_w
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


def noise_residual(
    x: np.ndarray, mouse: np.ndarray, weights: dict[str, float] | None = None
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Residual log CV2 (pooled within mice) on its mean trend, with the mean and dropout rate per gene."""
    mean, var, dropout = within_moments(x, mouse, weights)
    with np.errstate(divide="ignore", invalid="ignore"):
        log_mean = np.log(mean)
        log_cv2 = np.log(var / mean**2)
    return trend_residual(log_mean, log_cv2), mean, dropout


def strata(
    mean: np.ndarray,
    dropout: np.ndarray,
    n_mean: int = MEAN_BINS,
    n_dropout: int = DROPOUT_BINS,
    *,
    ties: bool = False,
) -> np.ndarray:
    """Stratum per gene: decile of mean expression crossed with tertile of dropout rate inside the decile.

    By default tied values are ranked in input order; with ``ties`` tied values share their lowest rank, so
    genes with equal mean or dropout always fall into the same bin.
    """

    def ranks(v: np.ndarray) -> np.ndarray:
        if ties:
            return rankdata(v, method="min").astype(np.int64) - 1
        return np.argsort(np.argsort(v, kind="stable"), kind="stable")

    rank = ranks(mean)
    decile = rank * n_mean // mean.size
    out = np.empty(mean.size, dtype=np.int64)
    for d in range(n_mean):
        sel = np.flatnonzero(decile == d)
        if not sel.size:
            continue
        r = ranks(dropout[sel])
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


def mouse_weights(mice: dict[str, list[str]], rng: np.random.Generator) -> dict[str, float]:
    """Bayesian-bootstrap weights per mouse: Dirichlet(1, ..., 1) within each age group, scaled to mean 1.

    Every mouse keeps a positive weight, so no cell type loses its support in a resample; drawn once per resample
    and shared by every cell type, as the same mice gave most tissues.
    """
    out = {}
    for ids in mice.values():
        w = rng.dirichlet(np.ones(len(ids))) * len(ids)
        out |= {str(m): float(v) for m, v in zip(ids, w, strict=True)}
    return out


def distinct_mice_ok(draw_ids: list[str], present: set[str], minimum: int = 2) -> bool:
    """Whether a multinomial draw of one age group holds ``minimum`` distinct original mice present in a cell type.

    A mouse drawn twice counts once: duplicate labels must not stand in for a second mouse.
    """
    return len(set(draw_ids) & present) >= minimum


def plant_noise(
    x: np.ndarray, rows: np.ndarray, cols: np.ndarray, shift: float, mouse: np.ndarray, rng: np.random.Generator
) -> np.ndarray:
    """A copy of ``x`` whose ``rows`` x ``cols`` block carries extra cell-to-cell noise of known size.

    Each entry is multiplied by Gamma noise of mean 1 and variance phi per column, chosen so that the column's
    within-mouse variance grows by ``var * (exp(shift) - 1)``: log CV2 (pooled within mice) rises by ``shift`` in
    expectation while the mean, the zeros and the cells and mice stay as they were.
    """
    out = x.astype(np.float64, copy=True)
    if shift <= 0 or not rows.any() or not cols.any():
        return out
    ri, ci = np.flatnonzero(rows), np.flatnonzero(cols)
    sub = out[np.ix_(ri, ci)]
    _, var, _ = within_moments(sub, mouse[ri])
    mean_sq = (sub**2).mean(axis=0)
    with np.errstate(divide="ignore", invalid="ignore"):
        phi = np.where(mean_sq > 0, var * np.expm1(shift) / mean_sq, 0.0)
    ok = phi > 0
    g = np.ones_like(sub)
    shape = np.broadcast_to(1 / phi[ok], (sub.shape[0], int(ok.sum())))
    g[:, ok] = rng.gamma(shape, 1 / shape)
    out[np.ix_(ri, ci)] = sub * g
    return out


def one_to_one_orthologs(orthologs: pl.DataFrame) -> pl.DataFrame:
    """``mouse_symbol, human_entrez, human_symbol`` for MGI homology classes with one mouse and one human gene."""
    counts = orthologs.group_by("class_key").agg(
        (pl.col("organism") == "mouse").sum().alias("n_mouse"), (pl.col("organism") == "human").sum().alias("n_human")
    )
    single = counts.filter((pl.col("n_mouse") == 1) & (pl.col("n_human") == 1)).select("class_key")
    mouse = orthologs.filter(pl.col("organism") == "mouse").join(single, on="class_key")
    human = orthologs.filter(pl.col("organism") == "human").join(single, on="class_key")
    out = mouse.select("class_key", pl.col("symbol").alias("mouse_symbol")).join(
        human.select("class_key", pl.col("entrez").alias("human_entrez"), pl.col("symbol").alias("human_symbol")),
        on="class_key",
    )
    dup = out.group_by("mouse_symbol").len().filter(pl.col("len") > 1)["mouse_symbol"]
    return out.filter(~pl.col("mouse_symbol").is_in(dup.implode())).drop("class_key").sort("mouse_symbol")


def pool_gaps(gaps: pl.DataFrame, boot: pl.DataFrame | None, keys: list[str]) -> pl.DataFrame:
    """Mean ``gap_sd`` over cell types per ``keys``, with a resampling interval when ``boot`` has rows for the key.

    The cell types that enter are those with a finite point gap (the estimable set of the contrast); the others are
    counted as excluded. A resample enters the interval only when it gives a finite gap for every one of them.
    """
    rows = []
    for key, sel in gaps.group_by(keys, maintain_order=True):
        est = sel.filter(pl.col("gap_sd").is_finite())
        row = dict(zip(keys, key, strict=True)) | {
            "cell_types": est.height,
            "cell_types_excluded": sel.height - est.height,
            "mean_gap_sd": float(est["gap_sd"].mean()) if est.height else float("nan"),
            "cell_types_below_zero": int((est["gap_sd"] < 0).sum()),
        }
        if boot is not None and boot.height:
            mine = boot.join(sel.select(keys).unique(), on=keys)
            if mine.height:
                b = (
                    mine.join(est.select(*keys, "tissue", "cell_type", "sd"), on=[*keys, "tissue", "cell_type"])
                    .with_columns(g=pl.col("matched_diff") / pl.col("sd"))
                    .filter(pl.col("g").is_finite())
                    .group_by("boot")
                    .agg(pl.col("g").mean(), pl.len().alias("n"))
                    .filter(pl.col("n") == est.height)
                )
                row |= {
                    "ci_low": float(b["g"].quantile(0.025)) if b.height else float("nan"),
                    "ci_high": float(b["g"].quantile(0.975)) if b.height else float("nan"),
                    "complete_resamples": b.height,
                    "resamples": mine["boot"].n_unique(),
                }
        rows.append(row)
    return pl.DataFrame(rows, infer_schema_length=None)


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


def _sc_delta(x, cells, weights: dict[str, float] | None = None) -> np.ndarray:
    """Old minus young residual log CV2 per gene (columns of ``x``)."""
    young = (cells["age_group"] == "young").to_numpy()
    r_y, _, _ = noise_residual(x[young], cells["mouse"].to_numpy()[young], weights)
    r_o, _, _ = noise_residual(x[~young], cells["mouse"].to_numpy()[~young], weights)
    return r_o - r_y


def _prepare(
    path: Path, cells: pl.DataFrame, gene_index: np.ndarray, rng: np.random.Generator, depth_index=None
) -> tuple[sparse.csr_matrix, pl.DataFrame, np.ndarray, int]:
    """Raw counts of the kept genes, cells, depth-matched counts (dense) and the depth, for one cell type.

    Depth is the ``DEPTH_QUANTILE`` quantile of read totals over ``gene_index`` (the mapped genes), or over
    ``depth_index`` (e.g. all endogenous genes) when given; shallower cells are dropped, the larger age group is
    subsampled to the smaller (at most ``MAX_CELLS``), and the kept genes are thinned at the depth's rate.
    """
    rows = cells["row"].to_numpy()
    order = np.argsort(rows)
    cells = cells[order.tolist()]
    raw_all = read_rows(path, rows[order], _n_genes(path))
    raw = raw_all[:, gene_index].tocsr()
    totals = np.asarray((raw if depth_index is None else raw_all[:, depth_index]).sum(axis=1)).ravel()
    del raw_all
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
    raw, cells, totals = raw[pick], cells[pick.tolist()], totals[pick]
    thinned = thin_rows(raw, depth, rng, None if depth_index is None else totals).toarray()
    return raw, cells, thinned, depth


def analyse(task: Task) -> tuple[list[dict], list[dict], dict]:
    """Point gaps (both versions), mouse-resample gaps and a cell-type summary row."""
    print(
        f"[tms_cell_noise] start {task.tissue} / {task.cell_type}: {task.cells.height} cells",
        file=sys.stderr,
        flush=True,
    )
    rng = np.random.default_rng(task.seed)
    raw, cells, thinned, depth = _prepare(task.path, task.cells, task.gene_index, rng)
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


@dataclass
class Loaded:
    path: Path
    obs: pl.DataFrame  # every usable cell of every age: row, age, mouse, tissue, cell_type, sex
    symbols: np.ndarray
    hgnc: pl.DataFrame
    set_genes: pl.DataFrame  # feature (human Ensembl), group
    ensembl: dict[str, str]  # mouse symbol -> human Ensembl, upper-case symbol route
    protein: set[str]
    gene_index: np.ndarray
    genes: np.ndarray
    sets: dict[str, np.ndarray]


def _inputs(ctx: Context) -> tuple[Path, list[Path]]:
    path = ctx.paths.input / "tms_facs" / H5AD
    annotation = AnnotationLayout(ctx.paths)
    return path, [path, annotation.gene_sets, annotation.hgnc]


def _load(ctx: Context) -> Loaded:
    path, _ = _inputs(ctx)
    annotation = AnnotationLayout(ctx.paths)
    obs, symbols = read_obs(path)
    # cells of two pooled mice ("3_10_M/3_11_M") belong to neither mouse and are left out
    obs = obs.filter(
        pl.col("age").is_in(YOUNG + OLD)
        & ~pl.col("cell_type").is_in(["nan", "", "None"])
        & ~pl.col("mouse").str.contains("/", literal=True)
    )
    hgnc, gene_sets = pl.read_parquet(annotation.hgnc), pl.read_parquet(annotation.gene_sets)
    mapping = symbol_ensembl([str(s) for s in symbols], hgnc)
    protein = set(hgnc.filter(pl.col("protein_coding") & pl.col("ensembl").is_not_null())["ensembl"].to_list())
    ensembl = dict(mapping.select("symbol", "ensembl").iter_rows())
    coding = np.array([ensembl.get(str(g)) in protein and not str(g).startswith("ERCC-") for g in symbols])
    gene_index = np.flatnonzero(coding)
    genes = symbols[gene_index]
    set_genes = gene_set_ensembl(gene_sets, hgnc).filter(pl.col("group").is_in(SETS))
    membership = set_genes.join(
        mapping.select(pl.col("ensembl").alias("feature"), pl.col("symbol").alias("mouse")), on="feature"
    )
    sets = {s: np.isin(genes, membership.filter(pl.col("group") == s)["mouse"].to_list()) for s in SETS}
    return Loaded(path, obs, symbols, hgnc, set_genes, ensembl, protein, gene_index, genes, sets)


def _age_groups(obs: pl.DataFrame, young: tuple[str, ...], old: tuple[str, ...], sex: str | None = None):
    sel = obs.filter(pl.col("age").is_in(young + old))
    if sex is not None:
        sel = sel.filter(pl.col("sex") == sex)
    return sel.with_columns(
        age_group=pl.when(pl.col("age").is_in(young)).then(pl.lit("young")).otherwise(pl.lit("old"))
    )


def run(ctx: Context) -> list[Path]:
    path, needed = _inputs(ctx)
    if any(not p.exists() for p in needed):
        print(f"tms_cell_noise: inputs missing ({', '.join(str(p) for p in needed if not p.exists())}); {URL}")
        return []
    plots.apply_style()
    out = ctx.out_dir("tms_cell_noise")
    workers, n_boot = _settings()
    data = _load(ctx)
    obs = _age_groups(data.obs, YOUNG, OLD)
    symbols, gene_index, genes, sets = data.symbols, data.gene_index, data.genes, data.sets
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
    """Mean gap over the estimable cell types; single-cell interval from the shared mouse resamples.

    A cell type without a finite gap (no retained non-essential genes) is excluded and counted; a resample enters
    only when it gives a finite gap for every estimable cell type (``pool_gaps``).
    """
    return pool_gaps(gaps, boot, ["version", "set"])


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


# ----------------------------------------------------------------------------- sensitivities (Codex review 8)
#: contrasts beside the main 3 m vs pooled 18-24 m; the 24 m mice are all male and the 21 m mice all female, so
#: same-sex contrasts take exact ages where three mice allow it (male 3 vs 18 m and female 3 vs 18 m have two old
#: mice and are listed as not estimable)
CONTRASTS: tuple[tuple[str, tuple[str, ...], tuple[str, ...], str | None], ...] = (
    ("3m_vs_18m", ("3m",), ("18m",), None),
    ("male_3m_vs_24m", ("3m",), ("24m",), "male"),
    ("female_3m_vs_21m", ("3m",), ("21m",), "female"),
    ("male_3m_vs_18m", ("3m",), ("18m",), "male"),
    ("female_3m_vs_18m", ("3m",), ("18m",), "female"),
    ("male_3m_vs_18-24m", ("3m",), ("18m", "24m"), "male"),
    ("female_3m_vs_18-21m", ("3m",), ("18m", "21m"), "female"),
)
#: planted rises of log CV2 at housekeeping and core-essential genes in old cells, in SDs of the age change
PLANTED = (0.05, 0.1, 0.2, 0.3, 0.5)


@dataclass
class SensTask:
    variant: str
    path: Path
    tissue: str
    cell_type: str
    cells: pl.DataFrame  # row, age_group, mouse
    gene_index: np.ndarray
    depth_index: np.ndarray | None
    sets: dict[str, tuple[np.ndarray, np.ndarray]]  # name -> (member, background) over gene_index
    weights: list[dict[str, float]]
    seed: int
    draws: list[dict[str, list[str]]]  # multinomial mouse draws, gated on distinct original mouse IDs
    mice: list[str]  # leave-one-mouse-out over these; empty for none
    planted: tuple[float, ...]
    plant_sets: tuple[str, ...]
    n_plant_weights: int
    young_reference: bool


def matched_count(values: np.ndarray, member: np.ndarray, background: np.ndarray, stratum: np.ndarray) -> int:
    """Members with a finite value in a stratum that holds finite background values (those ``matched_diff`` uses)."""
    ok = np.isfinite(values)
    have = np.unique(stratum[background & ok])
    return int((member & ok & np.isin(stratum, have)).sum())


def _gate(sub_cells: pl.DataFrame) -> bool:
    yb = (sub_cells["age_group"] == "young").to_numpy()
    n_mice = [sub_cells.filter(pl.Series(g))["mouse"].n_unique() for g in (yb, ~yb)]
    return min(n_mice) >= 2 and min(yb.sum(), (~yb).sum()) - max(n_mice) >= 5


def analyse_sensitivity(task: SensTask) -> tuple[list[dict], list[dict], dict]:
    """Point gaps per matching, resampled gaps (mouse weights, ID-gated mouse draws, left-out mice, planted)."""
    print(f"[tms_cell_noise:{task.variant}] start {task.tissue} / {task.cell_type}", file=sys.stderr, flush=True)
    rng = np.random.default_rng(task.seed)
    _, cells, thinned, depth = _prepare(task.path, task.cells, task.gene_index, rng, task.depth_index)
    young = (cells["age_group"] == "young").to_numpy()
    mouse = cells["mouse"].to_numpy()
    point, boot = [], []
    modes = ["pooled"] + (["young_reference"] if task.young_reference else [])
    for mode in modes:
        if mode == "pooled":
            detected = np.minimum((thinned[young] > 0).mean(axis=0), (thinned[~young] > 0).mean(axis=0))
            keep = detected >= MIN_DETECTED
            x = thinned[:, keep]
            ref = x
            stratum = strata(x.mean(axis=0), (x == 0).mean(axis=0))
        else:  # eligibility and matching on the young cells only, ties kept together
            keep = (thinned[young] > 0).mean(axis=0) >= MIN_DETECTED
            x = thinned[:, keep]
            ref = x[young]
            stratum = strata(ref.mean(axis=0), (ref == 0).mean(axis=0), ties=True)
        sets = {s: (m[keep], b[keep]) for s, (m, b) in task.sets.items()}
        ref_log_mean = np.log(ref.mean(axis=0))
        ref_dropout = (ref == 0).mean(axis=0)
        with np.errstate(invalid="ignore", divide="ignore"):
            delta = _sc_delta(x, cells)
        sd = float(np.nanstd(delta))
        key = {"variant": task.variant, "matching": mode, "tissue": task.tissue, "cell_type": task.cell_type}
        for s, (m, b) in sets.items():
            point.append(
                key
                | {
                    "set": s,
                    "n_set": int(m.sum()),
                    "n_matched": matched_count(delta, m, b, stratum),
                    "matched_diff": matched_diff(delta, m, b, stratum),
                    "sd": sd,
                    "balance_log_mean": matched_diff(ref_log_mean, m, b, stratum),
                    "balance_dropout": matched_diff(ref_dropout, m, b, stratum),
                }
            )

        def emit(kind, b, dlt, extra=None, sets=sets, key=key, stratum=stratum) -> None:
            for s, (m, bg) in sets.items():
                boot.append(
                    key
                    | {"kind": kind, "boot": str(b), "set": s, "matched_diff": matched_diff(dlt, m, bg, stratum)}
                    | (extra or {})
                )

        for b, w in enumerate(task.weights):
            with np.errstate(invalid="ignore", divide="ignore"):
                emit("weights", b, _sc_delta(x, cells, w))
        if mode != "pooled":
            continue
        for b, draw in enumerate(task.draws):
            present = set(mouse.tolist())
            if not all(distinct_mice_ok(ids, present) for ids in draw.values()):
                continue
            idx, labels = [], []
            for age, ids in draw.items():
                for j, m in enumerate(ids):
                    sel = np.flatnonzero(mouse == m)
                    idx.append(sel)
                    labels.append(np.full(sel.size, f"{age}:{j}"))
            idx = np.concatenate(idx)
            sub = cells[idx.tolist()].with_columns(mouse=pl.Series(np.concatenate(labels)))
            if not _gate(sub):
                continue
            with np.errstate(invalid="ignore", divide="ignore"):
                emit("mouse_id", b, _sc_delta(x[idx], sub))
        for m in task.mice:
            if m not in mouse:
                emit("lomo", m, delta)
                continue
            keep_rows = mouse != m
            sub = cells.filter(pl.Series(keep_rows))
            if not _gate(sub):
                emit("lomo", m, np.full(delta.shape, np.nan))
                continue
            with np.errstate(invalid="ignore", divide="ignore"):
                emit("lomo", m, _sc_delta(x[keep_rows], sub))
        if task.planted:
            union = np.any(np.column_stack([sets[s][0] for s in task.plant_sets]), axis=1)
            for k, target in enumerate(task.planted):
                xp = plant_noise(
                    x, ~young, union, target * sd, mouse, np.random.default_rng(task.seed + 7919 * (k + 1))
                )
                with np.errstate(invalid="ignore", divide="ignore"):
                    emit("planted_point", f"{target}", _sc_delta(xp, cells), {"target": target})
                    for b, w in enumerate(task.weights[: task.n_plant_weights]):
                        emit("planted", b, _sc_delta(xp, cells, w), {"target": target})
    cy, co = cells.filter(pl.Series(young)), cells.filter(pl.Series(~young))
    summary = {
        "variant": task.variant,
        "tissue": task.tissue,
        "cell_type": task.cell_type,
        "depth": depth,
        "cells_young": cy.height,
        "cells_old": co.height,
        "mice_young": cy["mouse"].n_unique(),
        "mice_old": co["mouse"].n_unique(),
    }
    return point, boot, summary


def _mgi_sets(ctx: Context, data: Loaded) -> tuple[np.ndarray, dict[str, tuple[np.ndarray, np.ndarray]], pl.DataFrame]:
    """Gene columns and sets by MGI one-to-one orthology; the agreeing subset; a coverage table."""
    from entropy_age.data.mouse import ORTHOLOG_FILE, MouseAnnotationLayout, read_orthologs

    report = MouseAnnotationLayout(ctx.paths).input / ORTHOLOG_FILE
    pairs = one_to_one_orthologs(read_orthologs(report))
    hgnc = data.hgnc.filter(pl.col("entrez").is_not_null() & pl.col("ensembl").is_not_null())
    by_entrez = pairs.join(hgnc.select(pl.col("entrez").alias("human_entrez"), "ensembl"), on="human_entrez")
    mgi = {m: e for m, e in by_entrez.select("mouse_symbol", "ensembl").iter_rows() if e in data.protein}
    symbols = [str(g) for g in data.symbols]
    index = np.array([i for i, g in enumerate(symbols) if g in mgi and not g.startswith("ERCC-")])
    genes = [symbols[i] for i in index]
    ens = np.array([mgi[g] for g in genes])
    sym_route = np.array([data.ensembl.get(g) if data.ensembl.get(g) in data.protein else None for g in genes])
    agree = sym_route == ens
    in_set = {s: np.isin(ens, data.set_genes.filter(pl.col("group") == s)["feature"].to_list()) for s in SETS}
    background = ~np.any(np.column_stack(list(in_set.values())), axis=1)
    sets = {s: (m, background) for s, m in in_set.items()}
    sets |= {f"{s}|agree": (m & agree, background & agree) for s, m in in_set.items()}
    # coverage of each route over the TMS symbols
    sym_member = {s: set(np.asarray(data.genes)[data.sets[s]].tolist()) for s in SETS}
    mgi_member = {s: set(np.asarray(genes)[in_set[s]].tolist()) for s in SETS}
    rows = [
        {
            "set": "all_protein_coding",
            "symbol_route": len(data.genes),
            "mgi_one_to_one": len(genes),
            "both_same_gene": int(agree.sum()),
        }
    ] + [
        {
            "set": s,
            "symbol_route": len(sym_member[s]),
            "mgi_one_to_one": len(mgi_member[s]),
            "both_same_gene": int((in_set[s] & agree).sum()),
            "symbol_route_only": len(sym_member[s] - mgi_member[s]),
            "mgi_only": len(mgi_member[s] - sym_member[s]),
        }
        for s in SETS
    ]
    coverage = pl.DataFrame(rows, infer_schema_length=None).with_columns(
        one_to_one_pairs=pl.lit(pairs.height), tms_symbols=pl.lit(len(symbols))
    )
    return index, sets, coverage


def run_sensitivity(ctx: Context) -> list[Path]:
    """Codex review 8 sensitivities beside the main tables (``sensitivity_*``, ``lomo``, ``planted``, ...)."""
    path, needed = _inputs(ctx)
    if any(not p.exists() for p in needed):
        print("tms_cell_noise sensitivities: inputs missing")
        return []
    out = ctx.out_dir("tms_cell_noise")
    workers, n_boot = _settings()
    n_weights = int(os.environ.get("ENTROPY_AGE_SENS_BOOT", "200"))
    n_plant = int(os.environ.get("ENTROPY_AGE_PLANT_BOOT", "50"))
    limit = int(os.environ.get("ENTROPY_AGE_TMS_CELL_TYPES", "0"))
    data = _load(ctx)
    symbol_sets = {s: (m, ~np.any(np.column_stack(list(data.sets.values())), axis=1)) for s, m in data.sets.items()}
    endogenous = np.array([i for i, g in enumerate(data.symbols) if not str(g).startswith("ERCC-")])
    mgi_index, mgi_sets, coverage = _mgi_sets(ctx, data)
    write_table(coverage, out / "orthology")

    def mice_of(obs: pl.DataFrame) -> dict[str, list[str]]:
        return {g: sorted(obs.filter(pl.col("age_group") == g)["mouse"].unique().to_list()) for g in ("young", "old")}

    tasks, contrast_rows = [], []
    main_obs = _age_groups(data.obs, YOUNG, OLD)
    main_mice = mice_of(main_obs)
    rng = ctx.rng(4242)  # the main run's multinomial draws, re-gated here on distinct original mouse IDs
    draws = [mouse_draw(main_mice, rng) for _ in range(n_boot)]
    wrng = ctx.rng(4343)
    main_weights = [mouse_weights(main_mice, wrng) for _ in range(n_weights)]
    all_mice = main_mice["young"] + main_mice["old"]
    main_groups = eligible(main_obs)
    if limit:
        main_groups = main_groups.head(limit)
    for i, (tissue, cell_type) in enumerate(main_groups.iter_rows()):
        cells = main_obs.filter((pl.col("tissue") == tissue) & (pl.col("cell_type") == cell_type)).select(
            "row", "age_group", "mouse"
        )
        common = {"path": path, "tissue": tissue, "cell_type": cell_type, "cells": cells, "seed": 1000 + i}
        tasks.append(
            SensTask(
                "main",
                gene_index=data.gene_index,
                depth_index=None,
                sets=symbol_sets,
                weights=main_weights,
                draws=draws,
                mice=all_mice,
                planted=PLANTED,
                plant_sets=("housekeeping", "core_essential"),
                n_plant_weights=n_plant,
                young_reference=True,
                **common,
            )
        )
        tasks.append(
            SensTask(
                "all_reads_depth",
                gene_index=data.gene_index,
                depth_index=endogenous,
                sets=symbol_sets,
                weights=main_weights,
                draws=[],
                mice=[],
                planted=(),
                plant_sets=(),
                n_plant_weights=0,
                young_reference=False,
                **common,
            )
        )
        tasks.append(
            SensTask(
                "mgi_one_to_one",
                gene_index=mgi_index,
                depth_index=None,
                sets=mgi_sets,
                weights=main_weights,
                draws=[],
                mice=[],
                planted=(),
                plant_sets=(),
                n_plant_weights=0,
                young_reference=False,
                **common,
            )
        )
    contrast_rows.append(
        {
            "contrast": "main",
            "young": "3m",
            "old": "18m,21m,24m",
            "sex": "both",
            "mice_young": len(main_mice["young"]),
            "mice_old": len(main_mice["old"]),
            "cell_types": main_groups.height,
        }
    )
    for j, (name, y_ages, o_ages, sex) in enumerate(CONTRASTS):
        obs = _age_groups(data.obs, y_ages, o_ages, sex)
        mice = mice_of(obs)
        groups = eligible(obs) if min(len(v) for v in mice.values()) >= MIN_MICE else eligible(obs).head(0)
        if limit:
            groups = groups.head(limit)
        contrast_rows.append(
            {
                "contrast": name,
                "young": ",".join(y_ages),
                "old": ",".join(o_ages),
                "sex": sex or "both",
                "mice_young": len(mice["young"]),
                "mice_old": len(mice["old"]),
                "cell_types": groups.height,
            }
        )
        if not groups.height:
            continue
        crng = ctx.rng(4400 + j)
        weights = [mouse_weights(mice, crng) for _ in range(n_weights)]
        for i, (tissue, cell_type) in enumerate(groups.iter_rows()):
            cells = obs.filter((pl.col("tissue") == tissue) & (pl.col("cell_type") == cell_type)).select(
                "row", "age_group", "mouse"
            )
            tasks.append(
                SensTask(
                    name,
                    path,
                    tissue,
                    cell_type,
                    cells,
                    data.gene_index,
                    None,
                    symbol_sets,
                    weights,
                    2000 + 100 * j + i,
                    [],
                    [],
                    (),
                    (),
                    0,
                    False,
                )
            )
    write_table(pl.DataFrame(contrast_rows), out / "contrasts")
    print(f"[tms_cell_noise] sensitivities: {len(tasks)} tasks, {n_weights} mouse-weight resamples", file=sys.stderr)
    if workers > 1:
        # spawn: a forked child inherits polars' locked thread pool and hangs (smoke job 6657459, 2026-11-26)
        with ProcessPoolExecutor(workers, mp_context=multiprocessing.get_context("spawn")) as pool:
            results = list(pool.map(analyse_sensitivity, tasks))
    else:
        results = [analyse_sensitivity(t) for t in tasks]
    point = pl.DataFrame([r for p, _, _ in results for r in p], infer_schema_length=None)
    boot = pl.DataFrame([r for _, b, _ in results for r in b], infer_schema_length=None)
    summary = pl.DataFrame([s for _, _, s in results], infer_schema_length=None)
    point = point.with_columns(gap_sd=pl.col("matched_diff") / pl.col("sd"))
    boot.write_parquet(out / "sensitivity_resamples.parquet")
    keys = ["variant", "matching", "set"]
    weights_b = boot.filter(pl.col("kind") == "weights").select(*keys, "tissue", "cell_type", "boot", "matched_diff")
    pooled = pool_gaps(point, weights_b.with_columns(pl.col("boot").cast(pl.Int64)), keys)
    id_b = boot.filter(pl.col("kind") == "mouse_id").select(*keys, "tissue", "cell_type", "boot", "matched_diff")
    id_pooled = pool_gaps(
        point.filter((pl.col("variant") == "main") & (pl.col("matching") == "pooled")),
        id_b.with_columns(pl.col("boot").cast(pl.Int64)),
        keys,
    ).select(
        *keys,
        pl.col("ci_low").alias("id_ci_low"),
        pl.col("ci_high").alias("id_ci_high"),
        pl.col("complete_resamples").alias("id_complete_resamples"),
        pl.col("resamples").alias("id_resamples"),
    )
    pooled = pooled.join(id_pooled, on=keys, how="left")
    # leave one mouse out: the pooled gap over the main estimable cell types with each mouse left out
    lomo = (
        boot.filter(pl.col("kind") == "lomo")
        .join(point.select(*keys, "tissue", "cell_type", "sd"), on=[*keys, "tissue", "cell_type"])
        .with_columns(g=pl.col("matched_diff") / pl.col("sd"))
        .filter(pl.col("g").is_finite())
        .group_by("set", pl.col("boot").alias("mouse_left_out"))
        .agg(pl.col("g").mean().alias("mean_gap_sd"), pl.len().alias("cell_types"))
        .sort("set", "mouse_left_out")
    )
    # planted effects: recovered pooled gap with its mouse-weight interval
    pp = (
        boot.filter(pl.col("kind") == "planted_point")
        .join(point.select(*keys, "tissue", "cell_type", "sd"), on=[*keys, "tissue", "cell_type"])
        .with_columns(gap_sd=pl.col("matched_diff") / pl.col("sd"))
    )
    pb = boot.filter(pl.col("kind") == "planted").select(
        "target", *keys, "tissue", "cell_type", pl.col("boot").cast(pl.Int64), "matched_diff"
    )
    planted = pool_gaps(pp.select("target", *keys, "tissue", "cell_type", "sd", "gap_sd"), pb, ["target", *keys])
    write_table(point.sort(*keys, "tissue", "cell_type"), out / "sensitivity_gaps", csv_max_rows=20_000)
    write_table(pooled.sort(*keys), out / "sensitivity_pooled")
    write_table(summary, out / "sensitivity_cell_types", csv_max_rows=5_000)
    write_table(lomo, out / "lomo")
    write_table(planted.sort("set", "target"), out / "planted")
    return [out / "sensitivity_pooled.csv"]
