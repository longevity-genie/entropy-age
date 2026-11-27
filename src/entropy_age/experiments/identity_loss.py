"""Do cells lose their methylation identity with age? Marker fidelity and classifier confidence (GSE186458)

Loyfer et al. 2023 (Nature 613:355; GSE186458) sorted 27 cell-type and tissue groups of donors aged 19-83 and called,
per cell type, the regions unmethylated in that type alone (the UXM atlas markers, ``Markers.U250.hg38.tsv`` of
github.com/nloyfer/UXM_deconv, ``data/input/loyfer_markers/``: up to 250 per type, 39 types, 9,520 regions). Every
published marker is a U marker; the paper's supplementary Table S4 lists U markers only as well, so no M markers exist
to test. Donor ids, sorting purity and depth come from the paper's Table S1 (``loyfer2023_table_s1_samples.csv``).

Each sample's PAT reads are cut to the marker regions. A fragment is informative at a marker when it covers at least
``MIN_CPGS`` = 3 of the region's CpGs; it carries the U pattern when all its CpGs there are unmethylated (strict,
the main measure) or at most ``U_MAX`` = 25 % are (UXM's threshold, a variant). Per sample:

- **own fidelity**: the mean over its own cell type's markers (``MIN_FRAGS`` informative fragments each) of the share
  of fragments with the U pattern; falling with age means fewer molecules keep the type's mark;
- **other-marker U share** (control): the same share at every other type's markers, where a sorted cell should be
  methylated; rising with age means other types' programmes leak into the cell;
- **identity contrast**: own fidelity minus the other-marker U share.

Age slopes: OLS with sex held and each cell-type group as a fixed effect (only age differences within a group count),
intervals by percentile bootstrap over samples within groups (sample-level, as ``read_entropy``) and over donors
(Table S1's PatientID, the real donor of each sample, resampled whole). Per cell type (marker type), the slope over its
samples with a t interval (group effects held where a type spans several groups; 3-10 samples each).

Link to stochastic drift: per sample, Shannon entropy of the 4-CpG epialleles at the promoter windows of
``read_entropy`` (counts from its caches), rarefied to ``RAREFY`` fragments per window so depth does not move it,
averaged over housekeeping windows and over their matched look-alikes (same CpG-island status and pooled-methylation
decile, ``read_entropy.matched_sides``). The association of fidelity with entropy at a given age: the entropy
coefficient (per SD of entropy) in fidelity ~ entropy + age + sex + group effects, and the partial correlation.

Classifier confidence: per-sample mean methylation at every marker region (all calls); a multinomial L2 logistic
regression (C = ``C_REG``) predicts the marker type (Endocrine colon cells, with no marker type, as their own class).
Leave-one-sample-out, and a split that trains on the samples younger than their class's median age and tests the
rest. The log-odds of the true class's probability is then regressed on age and on entropy as above.

Held-out markers (``heldout_*`` tables). The published markers were selected on these very samples (with others), so
their fidelity starts near its ceiling. Markers are chosen again from Loyfer's genome-wide blocks (``BLOCKS``; at least
5 CpGs and at most 1,500 bp, as the published markers) by ``select_u_markers``: every training sample of the type
covers the block, and the type's lower-quartile UXM U share beats the highest other-group mean by ``SEL_DELTA``; the
250 widest margins are kept. Scored on the strict U share of samples left out of the choice:
leave-donor-out (every other donor, all ages), and younger-only (samples below their type's median age choose, the
older rest are scored; donors of scored samples leave the background as well). An in-sample choice under the same rule
shows the ceiling for comparison. Limits: the blocks themselves come from Loyfer's unsupervised segmentation of these
and 106 other samples (no labels or ages used, but the held-out samples shaped the boundaries); only the 26 groups
deposited here form the background (Loyfer had 39 types); bootstrap and leave-one-out intervals hold the scores fixed,
so a donor's part in other samples' marker choice is not resampled.

Purity (``purity``, ``loso``, ``loso_summary``): Table S1's estimated purity held as a covariate, samples below 0.90
and 0.95 dropped, and every slope refitted without each sample in turn. ``entropy_sides``: the fidelity-on-entropy
coefficient and the entropy age slope for housekeeping and matched windows, and their difference from the same
resamples.

Caveats: slopes are within-type and small groups (3-10 samples) hold most types; the classifier is near perfect.

``entropy-age figures identity_loss`` would need registration; run as ``uv run python -c 'from
entropy_age.experiments import identity_loss; ...'`` in a Slurm job (reads about 33 GB of PAT files once, cached).
"""

from __future__ import annotations

import hashlib
import multiprocessing
import os
import sys
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

import numpy as np
import polars as pl
from scipy import stats

from entropy_age import plots
from entropy_age.experiments.common import Context, write_table

GSE = "GSE186458"
MARKERS = "loyfer_markers/Markers.U250.hg38.tsv"
TABLE_S1 = "loyfer_markers/loyfer2023_table_s1_samples.csv"
MIN_CPGS = 3
U_MAX = 0.25
MIN_FRAGS = 5
MIN_CALLS = 10
RAREFY = 10
C_REG = 1.0
N_BOOT = 1000
CACHE = "identity_loss_markers"
#: Loyfer's genome-wide segmentation of these and 106 other samples (GEO supplementary of GSE186458), the candidate
#: pool for markers chosen again on held-out donors
BLOCKS = "geo/GSE186458/blocks/GSE186458_blocks.s207.hg38.bed.gz"
BLOCK_CACHE = "identity_loss_blocks"
BLOCK_MIN_CPGS = 5  # the published U250 markers span 5-77 CpGs and at most 1,500 bp
BLOCK_MAX_BP = 1500
SEL_TOP = 250
SEL_DELTA = 0.4
SEL_TG_QUANT = 0.25
SEL_BG_COVER = 0.8
READ_CHUNK = 5_000_000
PURITY_CUTS = (0.9, 0.95)
#: atlas group (cell type | tissue) -> UXM marker type
TYPE_OF = {
    "Acinar|Pancreas": "Pancreas-Acinar",
    "Adipocytes|Abdominal Subcut.": "Adipocytes",
    "Alpha|Pancreas": "Pancreas-Alpha",
    "Basal epithelial|Breast": "Breast-Basal-Ep",
    "Beta|Pancreas": "Pancreas-Beta",
    "Cardiomyocyte|Heart": "Heart-Cardio",
    "Delta|Pancreas": "Pancreas-Delta",
    "Duct|Pancreas": "Pancreas-Duct",
    "Endothelium|Kidney glomerular": "Endothel",
    "Endothelium|Lung alveolar": "Endothel",
    "Endothelium|Vascular saphenous": "Endothel",
    "Epithelium|Colon": "Colon-Ep",
    "Epithelium|Gastric antrum": "Gastric-Ep",
    "Epithelium|Gastric body": "Gastric-Ep",
    "Epithelium|Gastric fundus": "Gastric-Ep",
    "Epithelium|Lung alveolar": "Lung-Ep-Alveo",
    "Epithelium|Prostate": "Prostate-Ep",
    "Epithelium|Tonsil palatine": "Head-Neck-Ep",
    "Fibroblast|Heart": "Heart-Fibro",
    "Hepatocyte|Liver": "Liver-Hep",
    "Luminal epithelial|Breast": "Breast-Luminal-Ep",
    "Neuronal|Brain": "Neuron",
    "Oligodendrocytes|Brain": "Oligodend",
    "T central memory CD4|Blood": "Blood-T",
    "T effector cell CD8|Blood": "Blood-T",
    "T effector memory CD4|Blood": "Blood-T",
}


# ----------------------------------------------------------------------------- pure helpers
def overlap_pieces(
    read_start: np.ndarray, read_len: np.ndarray, m_start: np.ndarray, m_end: np.ndarray
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """Reads against disjoint, sorted marker regions ``[m_start, m_end)`` (CpG indices): for every (read, marker)
    overlap, the read's row, the marker's row, and the offset and length of the read's calls inside the marker. A read
    spanning two markers yields two pieces."""
    read_start, read_len = np.asarray(read_start, dtype=np.int64), np.asarray(read_len, dtype=np.int64)
    read_end = read_start + read_len
    last = np.searchsorted(m_start, read_end, side="left") - 1  # the last marker starting before the read ends
    rows, marks = [], []
    back = 0
    while True:  # step back one marker at a time; disjoint sorted regions end in order, so a miss ends the read
        j = last - back
        ok = j >= 0
        ok[ok] = m_end[j[ok]] > read_start[ok]
        if not ok.any():
            break
        rows.append(np.flatnonzero(ok))
        marks.append(j[ok])
        back += 1
    if not rows:
        empty = np.zeros(0, dtype=np.int64)
        return empty, empty, empty, empty
    rows, marks = np.concatenate(rows), np.concatenate(marks)
    lo = np.maximum(m_start[marks], read_start[rows])
    hi = np.minimum(m_end[marks], read_end[rows])
    return rows, marks, lo - read_start[rows], hi - lo


def classify_fragments(
    n_meth: np.ndarray, n_unmeth: np.ndarray, min_cpgs: int = MIN_CPGS, u_max: float = U_MAX
) -> dict[str, np.ndarray]:
    """Per fragment piece (its methylated and unmethylated calls inside one marker): ``informative`` (at least
    ``min_cpgs`` called CpGs), ``u_full`` / ``m_full`` (informative and all unmethylated / all methylated), ``u_uxm`` /
    ``m_uxm`` (informative and a methylated share at most ``u_max`` / at least ``1 - u_max``, UXM's rule)."""
    n_meth, n_unmeth = np.asarray(n_meth), np.asarray(n_unmeth)
    called = n_meth + n_unmeth
    informative = called >= min_cpgs
    share = n_meth / np.maximum(called, 1)
    return {
        "informative": informative,
        "u_full": informative & (n_meth == 0),
        "m_full": informative & (n_unmeth == 0),
        "u_uxm": informative & (share <= u_max),
        "m_uxm": informative & (share >= 1 - u_max),
    }


def fidelity(u: np.ndarray, frags: np.ndarray, use: np.ndarray, min_frags: int = MIN_FRAGS) -> np.ndarray:
    """Per sample (row of ``(samples, markers)`` counts), the mean over markers in ``use`` (a boolean matrix of the same
    shape, or a vector) with at least ``min_frags`` informative fragments of the share carrying the pattern."""
    u, frags = np.asarray(u, dtype=float), np.asarray(frags, dtype=float)
    use = np.broadcast_to(use, u.shape) & (frags >= min_frags)
    with np.errstate(divide="ignore", invalid="ignore"):
        share = np.where(use, u / np.maximum(frags, 1), 0.0)
        return share.sum(axis=1) / use.sum(axis=1)


def within(x: np.ndarray, group: np.ndarray) -> np.ndarray:
    """Columns of ``x`` centred on their group means (the fixed-effect projection)."""
    x = np.asarray(x, dtype=float)
    out = x.copy()
    for g in np.unique(group):
        rows = group == g
        out[rows] -= x[rows].mean(axis=0)
    return out


def fe_ols(y: np.ndarray, x: np.ndarray, group: np.ndarray) -> tuple[np.ndarray, np.ndarray, int]:
    """OLS of ``y`` on the columns of ``x`` with ``group`` fixed effects: coefficients, SEs and residual df. Constant
    columns (after centring) get NaN."""
    x = np.asarray(x, dtype=float).reshape(len(y), -1)
    yc, xc = within(np.asarray(y, dtype=float)[:, None], group)[:, 0], within(x, group)
    live = xc.std(axis=0) > 1e-12
    coef, se = np.full(x.shape[1], np.nan), np.full(x.shape[1], np.nan)
    df = len(y) - np.unique(group).size - int(live.sum())
    if df < 1 or not live.any():
        return coef, se, df
    b, *_ = np.linalg.lstsq(xc[:, live], yc, rcond=None)
    resid = yc - xc[:, live] @ b
    cov = np.linalg.pinv(xc[:, live].T @ xc[:, live]) * (resid @ resid) / df
    coef[live], se[live] = b, np.sqrt(np.diag(cov))
    return coef, se, df


def partial_r(y: np.ndarray, x: np.ndarray, controls: np.ndarray, group: np.ndarray) -> float:
    """Correlation of ``y`` and ``x`` after both are freed of ``controls`` and group means."""
    c = within(np.asarray(controls, dtype=float).reshape(len(y), -1), group)
    c = c[:, c.std(axis=0) > 1e-12]

    def resid(v):
        v = within(np.asarray(v, dtype=float)[:, None], group)[:, 0]
        return v - c @ np.linalg.lstsq(c, v, rcond=None)[0] if c.size else v

    a, b = resid(y), resid(x)
    return float(np.corrcoef(a, b)[0, 1]) if a.std() > 0 and b.std() > 0 else float("nan")


def young_old_split(label: np.ndarray, age: np.ndarray) -> np.ndarray:
    """True for training samples: those younger than their class's median age (the rest are tested)."""
    train = np.zeros(len(label), dtype=bool)
    for c in np.unique(label):
        rows = label == c
        train[rows] = age[rows] < np.median(age[rows])
    return train


def u_shares(u: np.ndarray, frags: np.ndarray, min_frags: int = MIN_FRAGS) -> np.ndarray:
    """Per sample and region, the share of informative fragments with the U pattern; NaN below ``min_frags``."""
    u, frags = np.asarray(u, dtype=np.float32), np.asarray(frags, dtype=np.float32)
    with np.errstate(divide="ignore", invalid="ignore"):
        return np.where(frags >= min_frags, u / np.maximum(frags, 1), np.nan).astype(np.float32)


def select_u_markers(
    share: np.ndarray,
    train: np.ndarray,
    type_of: np.ndarray,
    group_of: np.ndarray,
    target: str,
    n_top: int = SEL_TOP,
    min_delta: float = SEL_DELTA,
    tg_quant: float = SEL_TG_QUANT,
    bg_cover: float = SEL_BG_COVER,
) -> np.ndarray:
    """U markers of ``target`` chosen on the ``train`` rows of ``share`` (samples x regions U shares, NaN where thin).

    A region qualifies when every training sample of the target type covers it, at least ``bg_cover`` of the other
    groups (every atlas group whose type differs, unlabelled groups included) cover it, and the target's
    ``tg_quant`` quantile of U share exceeds the highest other-group mean by at least ``min_delta``. The ``n_top``
    regions with the largest margin are returned (indices into the columns, best first). This mirrors the rule of
    Loyfer's markers (target unmethylated, every background group methylated, ranked by the margin) on UXM U shares."""
    type_of, group_of = np.asarray(type_of, dtype=object), np.asarray(group_of, dtype=object)
    tg_rows = np.flatnonzero(train & (type_of == target))
    if tg_rows.size == 0:
        return np.zeros(0, dtype=np.int64)
    tg = share[tg_rows]
    covered = np.isfinite(tg).all(axis=0)
    tg_low = np.quantile(np.where(covered, tg, 0.0), tg_quant, axis=0)
    bg_groups = [g for g in np.unique(group_of[train]) if not np.any((group_of == g) & (type_of == target))]
    means = np.full((len(bg_groups), share.shape[1]), np.nan, dtype=np.float32)
    for k, g in enumerate(bg_groups):
        rows = share[train & (group_of == g)]
        n = np.isfinite(rows).sum(axis=0)
        with np.errstate(invalid="ignore", divide="ignore"):
            means[k] = np.where(n > 0, np.nansum(rows, axis=0) / np.maximum(n, 1), np.nan)
    cover = np.isfinite(means).mean(axis=0) >= bg_cover if bg_groups else np.zeros(share.shape[1], dtype=bool)
    bg_max = np.max(np.where(np.isfinite(means), means, -np.inf), axis=0) if bg_groups else np.zeros(share.shape[1])
    delta = np.where(covered & cover, tg_low - bg_max, -np.inf)
    ok = np.flatnonzero(delta >= min_delta)
    return ok[np.argsort(-delta[ok], kind="stable")][:n_top]


def logit(p: np.ndarray, eps: float = 1e-6) -> np.ndarray:
    p = np.clip(np.asarray(p, dtype=float), eps, 1 - eps)
    return np.log(p / (1 - p))


# ----------------------------------------------------------------------------- data
def read_markers(path: Path) -> pl.DataFrame:
    return pl.read_csv(path, separator="\t").select(
        pl.col("startCpG").alias("start"), pl.col("endCpG").alias("end"), pl.col("target").alias("type")
    )


def _cache_key(regions: np.ndarray, pat: Path) -> str:
    st = pat.stat()
    h = hashlib.sha1(np.ascontiguousarray(regions, dtype=np.int64).tobytes())
    # "walk=all": overlap_pieces follows a read over any number of regions (caches before 2026-11-28 stopped at 3)
    h.update(f"{MIN_CPGS}|{U_MAX}|{pat.name}|{st.st_size}|{int(st.st_mtime)}|walk=all".encode())
    return h.hexdigest()[:16]


def marker_counts(pat: Path, m_start: np.ndarray, m_end: np.ndarray, light: bool = False) -> pl.DataFrame:
    """Per marker region (row of ``m_start``/``m_end``): weighted fragment counts of each class of
    ``classify_fragments`` and the methylated and total calls of all overlapping reads. ``light`` keeps only the
    informative, strict-U and UXM-U counts (as uint32), for the genome-wide blocks. Reads are cut in slices of
    ``READ_CHUNK`` so the pieces of a deep sample against 1.6 M blocks fit in memory."""
    reads = pl.read_csv(
        pat,
        separator="\t",
        has_header=False,
        columns=[1, 2, 3],
        new_columns=["chr", "start", "calls", "n"],
        schema_overrides={"column_2": pl.Int64, "column_3": pl.String, "column_4": pl.Int64},
        quote_char=None,
    )
    k = m_start.size
    names = ("informative", "u_full", "u_uxm") if light else ("informative", "u_full", "m_full", "u_uxm", "m_uxm")
    acc = {name: np.zeros(k) for name in (*names, *(() if light else ("meth_calls", "calls")))}
    for lo in range(0, reads.height, READ_CHUNK):
        chunk = reads.slice(lo, READ_CHUNK)
        rows, marks, offset, length = overlap_pieces(
            chunk["start"].to_numpy(), chunk["calls"].str.len_bytes().to_numpy(), m_start, m_end
        )
        piece = pl.DataFrame(
            {"calls": chunk["calls"].gather(rows), "offset": offset, "length": length, "n": chunk["n"].gather(rows)}
        ).select(pl.col("calls").str.slice(pl.col("offset"), pl.col("length")).alias("p"), "n")
        n_meth = piece["p"].str.count_matches("C").to_numpy()
        n_unmeth = piece["p"].str.count_matches("T").to_numpy()
        w = piece["n"].to_numpy()
        flags = classify_fragments(n_meth, n_unmeth)
        for name in names:
            acc[name] += np.bincount(marks, weights=w * flags[name], minlength=k)
        if not light:
            acc["meth_calls"] += np.bincount(marks, weights=w * n_meth, minlength=k)
            acc["calls"] += np.bincount(marks, weights=w * (n_meth + n_unmeth), minlength=k)
    del reads
    out = {"marker": np.arange(k)}
    out.update({c: (v.astype(np.uint32) if light else v) for c, v in acc.items()})
    return pl.DataFrame(out)


def _count_job(job: tuple[str, np.ndarray, np.ndarray, str, bool]) -> str:
    pat, m_start, m_end, cache, light = job
    if not Path(cache).exists():
        tmp = Path(cache).with_suffix(".tmp")
        marker_counts(Path(pat), m_start, m_end, light).write_parquet(tmp)
        tmp.rename(cache)
    return cache


def _samples(ctx: Context, folder: Path) -> pl.DataFrame:
    from entropy_age.experiments.read_entropy import read_samples

    s1 = pl.read_csv(ctx.paths.input / TABLE_S1, infer_schema_length=0).select(
        pl.col("Sample name").alias("title"),
        pl.col("PatientID").alias("patient"),
        pl.col("Estimated purity (Table S10)").cast(pl.Float64, strict=False).alias("purity"),
        pl.col("CpG coverage (fragments)").cast(pl.Float64, strict=False).alias("coverage"),
    )
    s = read_samples(folder).sort("age", "sample").join(s1, on="title", how="left")
    return s.with_columns(
        pl.col("group").replace_strict(TYPE_OF, default=None).alias("type"),
        pl.col("patient").fill_null(pl.col("sample")),
    )


def _load_marker_counts(
    ctx: Context, samples: pl.DataFrame, regions: pl.DataFrame, cache: str = CACHE, light: bool = False
) -> dict[str, np.ndarray]:
    m_start, m_end = regions["start"].to_numpy(), regions["end"].to_numpy()
    key_arr = np.concatenate([m_start, m_end])
    cache_dir = ctx.paths.interim / cache
    cache_dir.mkdir(parents=True, exist_ok=True)
    jobs = [
        (p, m_start, m_end, str(cache_dir / f"{s}_{_cache_key(key_arr, Path(p))}.parquet"), light)
        for s, p in samples.select("sample", "pat").iter_rows()
    ]
    missing = [j for j in jobs if not Path(j[3]).exists()]
    if missing:
        print(f"identity_loss: reading {len(missing)} PAT files into {cache}", file=sys.stderr, flush=True)
        os.environ["POLARS_MAX_THREADS"] = "1"  # inherited by the spawned workers only
        workers = int(os.environ.get("SLURM_CPUS_PER_TASK", "2"))
        with ProcessPoolExecutor(workers, mp_context=multiprocessing.get_context("spawn")) as pool:
            for i, _ in enumerate(pool.map(_count_job, missing), 1):
                if i % 10 == 0:
                    print(f"identity_loss: {i}/{len(missing)} PAT files", file=sys.stderr, flush=True)
    if light:
        cols = ("informative", "u_full", "u_uxm")
        out = {c: np.zeros((len(jobs), m_start.size), dtype=np.float32) for c in cols}
        for i, j in enumerate(jobs):
            f = pl.read_parquet(j[3]).sort("marker")
            for c in cols:
                out[c][i] = f[c].to_numpy()
        return out
    frames = [pl.read_parquet(j[3]).sort("marker") for j in jobs]
    cols = ("informative", "u_full", "m_full", "u_uxm", "m_uxm", "meth_calls", "calls")
    return {c: np.vstack([f[c].to_numpy() for f in frames]) for c in cols}


def _entropy(ctx: Context, samples: pl.DataFrame) -> pl.DataFrame:
    """Per sample: rarefied epiallele entropy at housekeeping windows and at their matched look-alikes."""
    from entropy_age.experiments.read_entropy import _load_counts, epiallele_metrics, matched_sides, rarefy
    from entropy_age.experiments.read_level import probe_windows
    from entropy_age.experiments.tissue_contrast import _features
    from entropy_age.stochasticity import baseline_strata

    features, probes, membership = _features(ctx.paths, ctx.rng(140_000))  # read_entropy's windows and caches
    windows = probe_windows(ctx.paths, features)
    counts = _load_counts(ctx, samples, windows["start"].to_numpy())
    counts = rarefy(counts, RAREFY, ctx.rng(186_459))
    n, w = counts.shape[:2]
    ent = epiallele_metrics(counts.reshape(n * w, -1), min_reads=RAREFY)
    entropy, mean = ent["entropy"].reshape(n, w), ent["mean"].reshape(n, w)
    keep = np.isfinite(entropy).mean(axis=0) >= 0.9
    entropy, mean = entropy[:, keep], mean[:, keep]
    feat = windows["cpg"].to_numpy()[keep]
    island = (
        pl.DataFrame({"feature": feat})
        .join(probes.select(pl.col("cpg").alias("feature"), "cpg_island"), on="feature", how="left")["cpg_island"]
        .fill_null(False)
        .to_numpy()
    )
    strata = np.char.add(island.astype(str), baseline_strata(np.nanmean(mean, axis=0), 10).astype(str))
    hk = np.isin(feat, membership.filter(pl.col("group") == "housekeeping")["feature"].to_numpy())
    sides = np.array([matched_sides(entropy[i], hk, strata) for i in range(n)])
    print(f"identity_loss: entropy over {keep.sum()} windows ({hk.sum()} housekeeping)", file=sys.stderr, flush=True)
    return pl.DataFrame(
        {
            "sample": samples["sample"],
            "entropy_housekeeping": sides[:, 0],
            "entropy_matched": sides[:, 1],
            "entropy_all": np.nanmean(entropy, axis=1),
        }
    )


# ----------------------------------------------------------------------------- analysis
def _boot(stat, n: int, draws: list[np.ndarray] | None, rng, n_boot: int, cluster: np.ndarray | None):
    """Percentile 95 % interval of ``stat(rows)`` over resamples within groups (``draws``) or whole clusters."""
    vals = []
    ids = np.unique(cluster) if cluster is not None else None
    members = [np.flatnonzero(cluster == c) for c in ids] if cluster is not None else None
    for _ in range(n_boot):
        if cluster is None:
            rows = np.concatenate([rng.choice(r, r.size, replace=True) for r in draws])
        else:
            rows = np.concatenate([members[i] for i in rng.integers(0, len(members), len(members))])
        with np.errstate(all="ignore"):
            vals.append(stat(rows))
    v = np.asarray(vals, dtype=float)
    v = v[np.isfinite(v)]
    if v.size < 10:
        return float("nan"), float("nan"), int(v.size)
    return float(np.quantile(v, 0.025)), float(np.quantile(v, 0.975)), int(v.size)


ENTROPY_TERMS = ("entropy_housekeeping", "entropy_matched", "entropy_all")


def _finite(frame: pl.DataFrame, cols) -> pl.DataFrame:
    """Rows where every named column present in ``frame`` is finite."""
    cols = [c for c in cols if c in frame.columns]
    return frame.filter(pl.all_horizontal([pl.col(c).cast(pl.Float64).is_finite() for c in cols])) if cols else frame


def _pooled_rows(
    frame: pl.DataFrame, outcomes: list[str], ctx: Context, extra: str = "", controls: tuple[str, ...] = ()
) -> list[dict]:
    """Pooled within-group age slopes and entropy links of each outcome, with both bootstrap intervals; ``controls``
    are further columns held in every model (purity). Each outcome is fitted on its own finite rows; an outcome with
    fewer than 10 such rows is skipped."""
    rows = []
    for y_name in outcomes:
        sub = _finite(frame, [y_name, "age", *ENTROPY_TERMS, *controls])
        if sub.height >= 10:
            rows += _pooled_rows_one(sub, [y_name], ctx, extra, controls)
    return rows


def _pooled_rows_one(
    frame: pl.DataFrame, outcomes: list[str], ctx: Context, extra: str = "", controls: tuple[str, ...] = ()
) -> list[dict]:
    age = frame["age"].to_numpy()
    female = np.column_stack([frame["female"].to_numpy().astype(float)] + [frame[c].to_numpy() for c in controls])
    if controls:
        extra = f" + {' + '.join(controls)}{extra}"
    group, patient = frame["group"].to_numpy(), frame["patient"].to_numpy()
    draws = [np.flatnonzero(group == g) for g in np.unique(group)]
    rows = []
    for y_name in outcomes:
        y = frame[y_name].to_numpy()
        specs = [("age", None)] + [
            (e, e) for e in ("entropy_housekeeping", "entropy_matched", "entropy_all") if e in frame.columns
        ]
        for term, ent in specs:
            if ent is None:

                def stat(r, y=y):
                    return fe_ols(y[r], np.column_stack([age[r], female[r]]), group[r])[0][0]

                coef, se, df = fe_ols(y, np.column_stack([age, female]), group)
                r = float("nan")
                unit = "per decade"
                scale = 10.0
            else:
                z = frame[ent].to_numpy()
                sd = float(np.std(within(z[:, None], group)[:, 0], ddof=1))

                def stat(r, y=y, z=z, sd=sd):
                    return fe_ols(y[r], np.column_stack([z[r] / sd, age[r], female[r]]), group[r])[0][0]

                coef, se, df = fe_ols(y, np.column_stack([z / sd, age, female]), group)
                r = partial_r(y, z, np.column_stack([age, female]), group)
                unit = "per within-group SD of entropy"
                scale = 1.0
            lo_s, hi_s, n_s = _boot(stat, len(y), draws, ctx.rng(7001), N_BOOT, None)
            lo_d, hi_d, n_d = _boot(stat, len(y), None, ctx.rng(7002), N_BOOT, patient)
            p = 2 * stats.t.sf(abs(coef[0] / se[0]), df) if np.isfinite(se[0]) and se[0] > 0 else float("nan")
            rows.append(
                {
                    "outcome": y_name,
                    "term": term,
                    "unit": unit,
                    "estimate": coef[0] * scale,
                    "se_ols": se[0] * scale,
                    "p_ols": float(p),
                    "partial_r": r,
                    "lo_samples": lo_s * scale,
                    "hi_samples": hi_s * scale,
                    "lo_donors": lo_d * scale,
                    "hi_donors": hi_d * scale,
                    "boots_usable": min(n_s, n_d),
                    "samples": len(y),
                    "groups": len(draws),
                    "donors": int(np.unique(patient).size),
                    "model": f"y ~ {term if ent is None else 'entropy + age'} + sex + group FE{extra}",
                    "outcome_mean": float(np.mean(y)),
                }
            )
    return rows


def _per_type(frame: pl.DataFrame, outcome: str) -> list[dict]:
    rows = []
    for t in sorted(frame["type"].drop_nulls().unique().to_list()):
        sub = frame.filter(pl.col("type") == t)
        y, age, group = sub[outcome].to_numpy(), sub["age"].to_numpy(), sub["group"].to_numpy()
        coef, se, df = fe_ols(y, age[:, None], group)
        half = stats.t.ppf(0.975, df) * se[0] if df >= 1 else np.nan
        rows.append(
            {
                "type": t,
                "outcome": outcome,
                "samples": sub.height,
                "donors": sub["patient"].n_unique(),
                "ages": ", ".join(str(int(a)) for a in sorted(age)),
                "mean": float(np.mean(y)),
                "slope_per_decade": 10 * coef[0],
                "lo": 10 * (coef[0] - half),
                "hi": 10 * (coef[0] + half),
                "df": df,
            }
        )
    return rows


def _boot_vec(stat, draws: list[np.ndarray] | None, rng, n_boot: int, cluster: np.ndarray | None):
    """``_boot`` for a statistic returning a vector: per-element 2.5 and 97.5 % quantiles over usable draws."""
    if cluster is not None:
        members = [np.flatnonzero(cluster == c) for c in np.unique(cluster)]
    vals = []
    for _ in range(n_boot):
        if cluster is None:
            rows = np.concatenate([rng.choice(r, r.size, replace=True) for r in draws])
        else:
            rows = np.concatenate([members[i] for i in rng.integers(0, len(members), len(members))])
        with np.errstate(all="ignore"):
            vals.append(np.asarray(stat(rows), dtype=float))
    v = np.asarray(vals)
    v = v[np.isfinite(v).all(axis=1)]
    if v.shape[0] < 10:
        nan = np.full(v.shape[1], np.nan)
        return nan, nan, int(v.shape[0])
    return np.quantile(v, 0.025, axis=0), np.quantile(v, 0.975, axis=0), int(v.shape[0])


def _side_rows(frame: pl.DataFrame, outcome: str, ctx: Context) -> list[dict]:
    """Housekeeping and matched windows side by side, with their difference from the same resamples: (a) the
    fidelity-on-entropy coefficient (per within-group SD of each side's entropy; age, sex and group held), and (b) each
    side's entropy age slope (per decade, in within-group SD of that side's entropy; sex and group held)."""
    frame = _finite(frame, [outcome, "age", "entropy_housekeeping", "entropy_matched"])  # one mask for both sides
    if frame.height < 10:
        return []
    age, female = frame["age"].to_numpy(), frame["female"].to_numpy().astype(float)
    group, patient = frame["group"].to_numpy(), frame["patient"].to_numpy()
    y = frame[outcome].to_numpy()
    z = {side: frame[f"entropy_{side}"].to_numpy() for side in ("housekeeping", "matched")}
    sd = {side: float(np.std(within(v[:, None], group)[:, 0], ddof=1)) for side, v in z.items()}

    def link(r):
        b = [fe_ols(y[r], np.column_stack([z[k][r] / sd[k], age[r], female[r]]), group[r])[0][0] for k in z]
        return [b[0], b[1], b[0] - b[1]]

    def drift(r):
        b = [10 * fe_ols(z[k][r] / sd[k], np.column_stack([age[r], female[r]]), group[r])[0][0] for k in z]
        return [b[0], b[1], b[0] - b[1]]

    draws = [np.flatnonzero(group == g) for g in np.unique(group)]
    rows = []
    for quantity, stat, unit in (
        ("fidelity on entropy", link, "per within-group SD of entropy"),
        ("entropy on age", drift, "within-group SD of entropy per decade"),
    ):
        est = stat(np.arange(len(y)))
        lo_s, hi_s, n_s = _boot_vec(stat, draws, ctx.rng(7011), N_BOOT, None)
        lo_d, hi_d, n_d = _boot_vec(stat, None, ctx.rng(7012), N_BOOT, patient)
        for k, side in enumerate(("housekeeping", "matched", "housekeeping minus matched")):
            rows.append(
                {
                    "outcome": outcome,
                    "quantity": quantity,
                    "windows": side,
                    "unit": unit,
                    "estimate": float(est[k]),
                    "lo_samples": float(lo_s[k]),
                    "hi_samples": float(hi_s[k]),
                    "lo_donors": float(lo_d[k]),
                    "hi_donors": float(hi_d[k]),
                    "boots_usable": min(n_s, n_d),
                    "samples": len(y),
                    "donors": int(np.unique(patient).size),
                }
            )
    return rows


def _loso_rows(frame: pl.DataFrame, outcome: str) -> list[dict]:
    """Leave-one-sample-out: the age slope (per decade) and the entropy coefficients (per within-group SD) without
    each sample in turn. The entropy SD stays that of all samples, so a change is a change of the coefficient alone.
    Scores are held fixed: for held-out fidelity, the dropped sample's part in other samples' marker choice stays."""
    frame = _finite(frame, [outcome, "age", *ENTROPY_TERMS])
    age, female, group = frame["age"].to_numpy(), frame["female"].to_numpy().astype(float), frame["group"].to_numpy()
    y = frame[outcome].to_numpy()
    sds = {
        t: float(np.std(within(frame[t].to_numpy()[:, None], group)[:, 0], ddof=1))
        for t in ENTROPY_TERMS
        if t in frame.columns
    }

    def est(term, r):
        if term == "age":
            return 10 * fe_ols(y[r], np.column_stack([age[r], female[r]]), group[r])[0][0]
        z = frame[term].to_numpy()[r]
        return fe_ols(y[r], np.column_stack([z / sds[term], age[r], female[r]]), group[r])[0][0]

    n, rows = len(y), []
    for term in ("age", *sds):
        full = est(term, np.arange(n))
        for i in range(n):
            without = est(term, np.flatnonzero(np.arange(n) != i))
            rows.append(
                {
                    "outcome": outcome,
                    "term": term,
                    "sample": frame["sample"][i],
                    "title": frame["title"][i],
                    "group": group[i],
                    "age": float(age[i]),
                    "purity": frame["purity"][i],
                    "outcome_value": float(y[i]),
                    "estimate_full": float(full),
                    "estimate_without": float(without),
                    "change": float(without - full),
                }
            )
    return rows


def _loso_summary(loso: pl.DataFrame) -> pl.DataFrame:
    """Per outcome and term: the leave-one-sample-out range and the three samples that move the estimate most."""
    rows = []
    for (outcome, term), sub in loso.group_by(["outcome", "term"], maintain_order=True):
        top = sub.with_columns(pl.col("change").abs().alias("a")).sort("a", descending=True).head(3)
        rows.append(
            {
                "outcome": outcome,
                "term": term,
                "estimate_full": sub["estimate_full"][0],
                "min_without_one": sub["estimate_without"].min(),
                "max_without_one": sub["estimate_without"].max(),
                "most_influential": "; ".join(
                    f"{t} (purity {pu if pu is not None else 'NA'}, {c:+.4f})"
                    for t, pu, c in top.select("title", "purity", "change").iter_rows()
                ),
            }
        )
    return pl.DataFrame(rows)


def _heldout_fidelity(ctx: Context, samples: pl.DataFrame) -> tuple[pl.DataFrame, pl.DataFrame]:
    """Own-marker fidelity on markers chosen without the scored sample (``select_u_markers`` on Loyfer's
    genome-wide blocks):

    - ``fidelity_ldo``: leave-donor-out, markers chosen on every other donor's samples (all ages);
    - ``fidelity_young_sel``: markers chosen on the samples younger than their type's median age, scored on the rest
      (NaN for the training samples), so drift-resistant regions cannot be favoured by old samples;
    - ``fidelity_insample``: the same rule on all samples (the circular version, for comparison).

    Selection uses UXM U shares; scoring uses the strict all-unmethylated share (``_uxm`` columns: UXM shares)."""
    b = pl.read_csv(
        ctx.paths.input / BLOCKS,
        separator="\t",
        has_header=False,
        new_columns=["chr", "bp_start", "bp_end", "start", "end"],
        schema_overrides={"column_4": pl.Int64, "column_5": pl.Int64},
    )
    blocks = (
        b.filter(
            (pl.col("end") - pl.col("start") >= BLOCK_MIN_CPGS)
            & (pl.col("bp_end") - pl.col("bp_start") <= BLOCK_MAX_BP)
        )
        .select("start", "end")
        .unique()
        .sort("start")
    )
    print(f"identity_loss: {blocks.height} candidate blocks", file=sys.stderr, flush=True)
    counts = _load_marker_counts(ctx, samples, blocks, cache=BLOCK_CACHE, light=True)
    frags = counts["informative"]
    share = u_shares(counts["u_uxm"], frags)
    type_of = samples["type"].to_numpy().astype(object)
    group_of, patient = samples["group"].to_numpy().astype(object), samples["patient"].to_numpy()
    age = samples["age"].to_numpy()
    n = samples.height
    typed = np.array([t is not None for t in type_of])

    def score(i, sel, u):
        use = np.zeros(frags.shape[1], dtype=bool)
        use[sel] = True
        return float(fidelity(u[i : i + 1], frags[i : i + 1], use)[0]) if sel.size else float("nan")

    def scored(i, sel):
        return float((frags[i, sel] >= MIN_FRAGS).sum())

    out = {
        c: np.full(n, np.nan)
        for c in (
            "fidelity_ldo",
            "fidelity_ldo_uxm",
            "ldo_markers",
            "ldo_markers_scored",
            "young_sel_markers_scored",
            "fidelity_young_sel",
            "fidelity_young_sel_uxm",
            "fidelity_insample",
        )
    }
    everyone = np.ones(n, dtype=bool)
    insample = {t: select_u_markers(share, everyone, type_of, group_of, t) for t in sorted(set(type_of[typed]))}
    for i in np.flatnonzero(typed):
        out["fidelity_insample"][i] = score(i, insample[type_of[i]], counts["u_full"])
    for d in np.unique(patient[typed]):
        train = patient != d
        for t in sorted(set(type_of[typed & ~train])):
            sel = select_u_markers(share, train, type_of, group_of, t)
            for i in np.flatnonzero(~train & (type_of == t)):
                out["fidelity_ldo"][i] = score(i, sel, counts["u_full"])
                out["fidelity_ldo_uxm"][i] = score(i, sel, counts["u_uxm"])
                out["ldo_markers"][i] = sel.size
                out["ldo_markers_scored"][i] = scored(i, sel)
    label = np.where(typed, type_of, "").astype(str)
    young = young_old_split(label, age) & typed
    young_train = young | ~typed  # unlabelled groups (endocrine colon) stay in the background
    young_sel = {}
    for t in sorted(set(type_of[typed])):
        scored_donors = np.unique(patient[typed & ~young & (type_of == t)])
        # a scored donor's samples of other types stay out of the background too
        young_sel[t] = select_u_markers(share, young_train & ~np.isin(patient, scored_donors), type_of, group_of, t)
    for i in np.flatnonzero(typed & ~young):
        out["fidelity_young_sel"][i] = score(i, young_sel[type_of[i]], counts["u_full"])
        out["fidelity_young_sel_uxm"][i] = score(i, young_sel[type_of[i]], counts["u_uxm"])
        out["young_sel_markers_scored"][i] = scored(i, young_sel[type_of[i]])
    print("identity_loss: held-out marker selection done", file=sys.stderr, flush=True)

    published = read_markers(ctx.paths.input / MARKERS)
    pub_key = {t: set(zip(sub["start"], sub["end"], strict=True)) for (t,), sub in published.group_by(["type"])}
    bs, be = blocks["start"].to_numpy(), blocks["end"].to_numpy()
    summary = []
    for t, sel in insample.items():
        keys = set(zip(bs[sel].tolist(), be[sel].tolist(), strict=True))
        ysel = young_sel[t]
        ykeys = set(zip(bs[ysel].tolist(), be[ysel].tolist(), strict=True))
        rows_t = typed & (type_of == t)
        summary.append(
            {
                "type": t,
                "samples": int(rows_t.sum()),
                "young_training_samples": int((young & rows_t).sum()),
                "insample_markers": len(keys),
                "insample_in_published": len(keys & pub_key.get(t, set())),
                "published_markers": len(pub_key.get(t, set())),
                "young_markers": len(ykeys),
                "young_shared_with_insample": len(ykeys & keys),
                "ldo_markers_median": float(np.nanmedian(out["ldo_markers"][rows_t])) if rows_t.any() else np.nan,
            }
        )
    frame = samples.select("sample").with_columns(
        **{k: pl.Series(v) for k, v in out.items()}, young_training=pl.Series(young)
    )
    return frame, pl.DataFrame(summary)


def _classifier(x: np.ndarray, label: np.ndarray, age: np.ndarray) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """True-class probability, leave-one-sample-out and young-to-old (NaN for training samples), and LOSO accuracy."""
    from sklearn.impute import SimpleImputer
    from sklearn.linear_model import LogisticRegression
    from sklearn.pipeline import make_pipeline

    def model():
        return make_pipeline(SimpleImputer(strategy="mean"), LogisticRegression(C=C_REG, max_iter=5000))

    n = len(label)
    loso, hit = np.full(n, np.nan), np.zeros(n, dtype=bool)
    for i in range(n):
        tr = np.arange(n) != i
        m = model().fit(x[tr], label[tr])
        p = m.predict_proba(x[i : i + 1])[0]
        classes = list(m.classes_)
        loso[i] = p[classes.index(label[i])]
        hit[i] = classes[int(np.argmax(p))] == label[i]
    train = young_old_split(label, age)
    split = np.full(n, np.nan)
    m = model().fit(x[train], label[train])
    classes = list(m.classes_)
    test = np.flatnonzero(~train & np.isin(label, classes))
    p = m.predict_proba(x[test])
    split[test] = [p[k, classes.index(label[i])] for k, i in enumerate(test)]
    return loso, split, hit


def run(ctx: Context) -> list[Path]:
    """Marker fidelity, other-marker leakage and classifier confidence against age and read-level entropy."""
    folder = ctx.paths.input / "geo" / GSE / "atlas"
    if not (ctx.paths.input / MARKERS).exists() or not folder.exists():
        print("identity_loss: markers or GSE186458 atlas missing. Nothing run.", file=sys.stderr)
        return []
    plots.apply_style()
    out = ctx.out_dir("identity_loss")
    markers = read_markers(ctx.paths.input / MARKERS)
    regions = markers.select("start", "end").unique().sort("start")
    samples = _samples(ctx, folder)
    counts = _load_marker_counts(ctx, samples, regions)
    region_types = regions.join(markers, on=["start", "end"]).group_by("start", "end").agg("type").sort("start")
    types_of_region = region_types["type"].to_list()
    own = np.array([[t in ts for ts in types_of_region] for t in samples["type"].to_list()])  # samples x regions
    other = ~own  # every other type's markers
    frags = counts["informative"]
    measures = {
        "fidelity_own": fidelity(counts["u_full"], frags, own),
        "fidelity_own_uxm": fidelity(counts["u_uxm"], frags, own),
        "u_share_other": fidelity(counts["u_full"], frags, other),
        "u_share_other_uxm": fidelity(counts["u_uxm"], frags, other),
        "own_markers_used": (own & (frags >= MIN_FRAGS)).sum(axis=1),
        "own_median_frags": np.array(
            [np.median(frags[i, own[i]]) if own[i].any() else np.nan for i in range(len(own))]
        ),
    }
    measures["identity_contrast"] = measures["fidelity_own"] - measures["u_share_other"]
    with np.errstate(divide="ignore", invalid="ignore"):
        beta = np.where(counts["calls"] >= MIN_CALLS, counts["meth_calls"] / counts["calls"], np.nan)
    measures["own_mean_methylation"] = np.array(
        [np.nanmean(beta[i, own[i]]) if own[i].any() else np.nan for i in range(len(own))]
    )
    frame = samples.drop("pat").with_columns(**{k: pl.Series(v) for k, v in measures.items()})
    frame = frame.join(_entropy(ctx, samples), on="sample", how="left")

    label = frame["type"].fill_null(frame["group"].str.replace(r"\|", "-")).to_numpy()
    loso, split, hit = _classifier(beta, label, frame["age"].to_numpy())
    frame = frame.with_columns(
        class_label=pl.Series(label),
        p_true_loso=pl.Series(loso),
        logit_p_true_loso=pl.Series(logit(loso)),
        loso_correct=pl.Series(hit),
        p_true_young_to_old=pl.Series(split),
        logit_p_true_young_to_old=pl.Series(logit(split)),
    )
    write_table(frame, out / "samples")

    typed = frame.filter(pl.col("type").is_not_null() & pl.col("fidelity_own").is_finite())
    slope_rows = _pooled_rows(
        typed, ["fidelity_own", "fidelity_own_uxm", "u_share_other", "identity_contrast"], ctx
    ) + _pooled_rows(frame, ["logit_p_true_loso"], ctx)
    split_frame = frame.filter(pl.col("p_true_young_to_old").is_finite())
    if split_frame.height >= 10:
        slope_rows += _pooled_rows(split_frame, ["logit_p_true_young_to_old"], ctx)
    slopes = pl.DataFrame(slope_rows, infer_schema_length=None)
    write_table(slopes, out / "slopes")

    # held-out markers (selection without the scored donor, or on younger donors only)
    held, selection = _heldout_fidelity(ctx, samples)
    write_table(selection, out / "heldout_markers")
    typed = typed.join(held, on="sample", how="left")
    write_table(
        typed.select("sample", "title", "group", "type", "age", "patient", "purity", *held.columns[1:]),
        out / "heldout_samples",
    )
    held_outcomes = ["fidelity_ldo", "fidelity_ldo_uxm", "fidelity_insample"]
    held_rows = [r | {"rows": "all scored"} for r in _pooled_rows(typed, held_outcomes, ctx)]
    old_rows = typed.filter(pl.col("fidelity_young_sel").is_finite())
    held_rows += [
        r | {"rows": "older than type median (young-selection test rows)"}
        for r in _pooled_rows(old_rows, ["fidelity_young_sel", "fidelity_young_sel_uxm", "fidelity_ldo"], ctx)
    ]
    write_table(pl.DataFrame(held_rows, infer_schema_length=None), out / "heldout_slopes")

    # purity: held as a covariate, and low-purity samples dropped
    main = ["fidelity_own", "fidelity_ldo"]
    pur = typed.filter(pl.col("purity").is_finite())
    pur_rows = [r | {"subset": "purity known"} for r in _pooled_rows(pur, main, ctx)]
    pur_rows += [
        r | {"subset": "purity known, purity held"} for r in _pooled_rows(pur, main, ctx, controls=("purity",))
    ]
    for cut in PURITY_CUTS:
        sub = pur.filter(pl.col("purity") >= cut)
        pur_rows += [r | {"subset": f"purity >= {cut}"} for r in _pooled_rows(sub, main, ctx)]
    write_table(pl.DataFrame(pur_rows, infer_schema_length=None), out / "purity")
    loso = pl.DataFrame(_loso_rows(typed, "fidelity_own") + _loso_rows(typed, "fidelity_ldo"))
    write_table(loso, out / "loso")
    write_table(_loso_summary(loso), out / "loso_summary")

    # housekeeping vs matched windows, and their difference
    side_rows = _side_rows(typed, "fidelity_own", ctx) + _side_rows(typed, "fidelity_ldo", ctx)
    write_table(pl.DataFrame(side_rows), out / "entropy_sides")
    per_type = pl.DataFrame(
        _per_type(typed, "fidelity_own") + _per_type(typed, "u_share_other") + _per_type(frame, "logit_p_true_loso")
    )
    write_table(per_type, out / "per_type")
    summary = (
        markers.group_by("type")
        .agg(pl.len().alias("markers"))
        .join(
            frame.group_by("type").agg(
                pl.len().alias("samples"),
                pl.col("age").sort().alias("ages_list"),
                pl.col("fidelity_own").mean(),
                pl.col("own_markers_used").median(),
            ),
            on="type",
            how="left",
        )
        .with_columns(pl.col("ages_list").list.eval(pl.element().cast(pl.Int64).cast(pl.String)).list.join(", "))
        .rename({"ages_list": "ages"})
        .sort("type")
    )
    write_table(summary, out / "markers")
    print(f"identity_loss: LOSO accuracy {hit.mean():.3f}", file=sys.stderr, flush=True)
    return [_figure(ctx, typed, frame, per_type, slopes)]


def _figure(ctx, typed: pl.DataFrame, frame: pl.DataFrame, per_type: pl.DataFrame, slopes: pl.DataFrame) -> Path:
    fig, axes = plots.panels(4, ncols=2, size=4.2)
    g = typed["group"].to_numpy()
    age_c = within(typed["age"].to_numpy()[:, None], g)[:, 0]

    def pooled(outcome, term):
        r = slopes.filter((pl.col("outcome") == outcome) & (pl.col("term") == term)).row(0, named=True)
        return f"{r['estimate']:+.4f} [{r['lo_donors']:+.4f}, {r['hi_donors']:+.4f}]"

    ax = axes[0]
    for k, (col, name) in enumerate((("fidelity_own", "own markers: U share"), ("u_share_other", "other markers"))):
        y = within(typed[col].to_numpy()[:, None], g)[:, 0]
        ax.scatter(age_c, y, s=12, color=plots.PALETTE[k], alpha=0.75, label=f"{name}, {pooled(col, 'age')}/decade")
    ax.axhline(0, color=plots.MUTED, lw=0.8)
    ax.set(xlabel="age minus group mean (years)", ylabel="fragment share minus group mean", title="Marker fidelity")
    ax.legend(fontsize=6.5, loc="lower left")

    ax = axes[1]
    pt = per_type.filter(pl.col("outcome") == "fidelity_own").sort("slope_per_decade")
    yy = np.arange(pt.height)
    lo, hi = pt["lo"].to_numpy(), pt["hi"].to_numpy()
    est = pt["slope_per_decade"].to_numpy()
    ok = np.isfinite(lo)
    ax.errorbar(est[ok], yy[ok], xerr=[est[ok] - lo[ok], hi[ok] - est[ok]], fmt="o", ms=3, color=plots.PALETTE[0])
    ax.plot(est[~ok], yy[~ok], "o", ms=3, mfc="white", color=plots.PALETTE[0])
    ax.set_yticks(yy, [f"{t} ({n})" for t, n in zip(pt["type"], pt["samples"], strict=True)], fontsize=6)
    ax.axvline(0, color=plots.MUTED, lw=0.8)
    ax.set_xlim(-0.15, 0.15)
    ax.set(xlabel="own-marker fidelity slope per decade (95 % t)", title="Per cell type")

    ax = axes[2]
    z = within(typed["entropy_housekeeping"].to_numpy()[:, None], g)[:, 0]
    y = within(typed["fidelity_own"].to_numpy()[:, None], g)[:, 0]
    ax.scatter(z, y, s=12, c=typed["age"].to_numpy(), cmap="viridis")
    ax.set(
        xlabel="housekeeping read entropy minus group mean (bits/CpG)",
        ylabel="own fidelity minus group mean",
        title=f"Fidelity and drift: {pooled('fidelity_own', 'entropy_housekeeping')} per SD",
    )
    ax.title.set_fontsize(8)

    ax = axes[3]
    g2 = frame["group"].to_numpy()
    ax.scatter(
        within(frame["age"].to_numpy()[:, None], g2)[:, 0],
        within(frame["logit_p_true_loso"].to_numpy()[:, None], g2)[:, 0],
        s=12,
        color=plots.PALETTE[2],
    )
    ax.axhline(0, color=plots.MUTED, lw=0.8)
    ax.set(
        xlabel="age minus group mean (years)",
        ylabel="logit P(true type) minus group mean",
        title=f"Classifier (LOSO): {pooled('logit_p_true_loso', 'age')}/decade",
    )
    ax.title.set_fontsize(8)
    for k, a in enumerate(axes):
        plots.label_panel(a, "abcd"[k])
    fig.suptitle(
        "Cell identity at Loyfer's U markers against age, sorted cells (GSE186458); [donor-bootstrap 95 %]",
        x=0.02,
        y=1.02,
        ha="left",
        fontsize=10,
        fontweight="bold",
    )
    fig.tight_layout()
    return plots.save(fig, ctx.figure_path("identity_loss"))
