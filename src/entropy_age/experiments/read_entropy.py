"""Is there a molecule-level counterpart of the shield? Read-level epiallele measures and age (GSE186458)

The paper's shield is measured on CpG averages over cells. A sequencing read reports several neighbouring CpGs of one
DNA molecule, so the mix of epialleles in a window says how much the cells disagree with each other. Loyfer et al.
2023 (GSE186458; whole-genome bisulfite, wgbstools PAT reads) give three sets with ages: white blood cells and
plasma cfDNA of the same 23 donors (21-75 years) and 101 sorted cell-type samples (27 cell-type and tissue groups,
19-83 years). Windows are the ``WIDTH`` = 4 consecutive CpGs starting at each 450k promoter probe's CpG (hg38),
read off every read that covers all four (``read_level.epiallele_counts``). A window counts in a sample with at
least ``MIN_READS`` such reads. Per sample and window, from the 16 epiallele shares p:

- **mean** methylation over the window's CpGs;
- **entropy**, Shannon entropy of the epiallele distribution in bits divided by the window's CpGs (the methylation
  entropy of Xie et al. 2011, Nucleic Acids Res 39:4099);
- **epipolymorphism**, 1 - sum p^2 (Landan et al. 2012, Nat Genet 44:1207);
- **MHL**, methylation haplotype load (Guo et al. 2017, Nat Genet 49:635): sum over lengths i = 1..4 of
  i x (share of fully methylated stretches of i CpGs among all stretches of i CpGs on the reads), divided by sum i;
- **relative entropy**, KL(p || p_young) in bits per CpG: the read-level form of the paper's headline measure, with
  p_young the pooled epiallele shares of the young reference (ages at or below the 20 % quantile; in the sorted
  cells, of each cell-type group), the sample itself left out, and a pseudo-count of ``KL_PSEUDO`` on every share;
- two depth checks on the Shannon-entropy windows: **log depth** (log(1 + complete fragments)) and **eligible** (1
  when the window reaches ``MIN_READS``), whose age slopes say whether older samples were read more or less deeply.

Per window and measure, the age slope: OLS with sex held (when it varies) and, for the sorted cells, each cell-type
group as a fixed effect, so only age differences within a cell type count. Windows need a value of that measure in
three quarters of the samples (and ``MIN_SAMPLES``); each measure keeps its own windows, so relative entropy, which
loses samples alone at their group's young end, does not narrow the Shannon-entropy panel. The contrast compares a
set's windows (housekeeping, core-essential, non-essential) with all other windows of the same CpG-island status and
young-level decile (the matching of ``tissue_contrast._with_strata``; the decile from the young reference's mean
methylation). Both sides are reported: the set's mean slope and the matched mean slope (stratum means of the
non-members, weighted by the set's windows per stratum), and their difference, in raw units and in SD of the
slopes. Intervals: percentile bootstrap over samples, everything (slopes, young reference, strata, contrast)
recomputed per resample; ``boots_usable`` counts the resamples with a finite contrast and ``p_boot`` is a two-sided
sign-tail summary of them, twice ``(1 + draws on the minority side) / (1 + usable)``, never zero.

Variants (``variant`` column):

- **raw**: the counts as read;
- **rarefied to D**: every scored sample-window cut to exactly D complete fragments by sampling without replacement
  (one fixed draw), so plug-in entropy carries the same downward small-sample bias everywhere; D = 10 on all windows
  and D = 20 on windows with balanced coverage (the share of samples reaching D differs by at most ``BALANCE_TOL``
  between the younger and the older half of each group);
- **leukocyte-adjusted** (white blood cells): leukocyte fractions estimated from each donor's own reads (mean
  methylation at marker windows of six sorted blood cell types of the same atlas, ``blood_refs/``, by non-negative
  least squares), their leading principal components held in every slope;
- **composition-only mixture** (white blood cells): each donor rebuilt as the mixture, at its estimated fractions,
  of the pooled sorted-cell epiallele distributions (the same at every age), sampled at the donor's own depth per
  window; any contrast it shows comes from composition and depth alone;
- **putative-donor clusters** (sorted cells): GEO gives no donor id for the atlas; samples sharing age, sex, lab and
  organ are taken as one donor and resampled together. The main sorted-cell intervals resample samples within groups
  and are sample-level. Leave-one-group-out point estimates (``leave_one_group_out``) show each group's influence.

``young_reference`` gives per source and group the young samples, the putative donors among them and the depth of
the pooled reference. Count caches are keyed by the ordered window starts and each PAT file's name, size and time.

``entropy-age figures read_entropy`` as a Slurm job; reads about 2.5 GB of counts (counted once per panel).
"""

from __future__ import annotations

import hashlib
import os
import re
import sys
from collections.abc import Callable
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

import numpy as np
import polars as pl
from scipy import optimize, stats

from entropy_age import plots
from entropy_age.experiments.common import Context, write_table
from entropy_age.experiments.read_level import GSE, INDEX, SOURCES, WIDTH, _count_one, probe_windows
from entropy_age.experiments.tissue_contrast import _features
from entropy_age.stochasticity import baseline_strata

MIN_READS = 10
MIN_SAMPLES = 18
KL_PSEUDO = 0.01
YOUNG_QUANTILE = 0.2
MEASURES = ("mean", "entropy", "epipolymorphism", "mhl", "relative_entropy")
DEPTH_MEASURES = ("log_depth", "eligible")
SETS = ("housekeeping", "core_essential", "non_essential")
RAREFY = ((10, False), (20, True))  # (fragments, balanced-coverage windows only)
BALANCE_TOL = 0.1
BLOOD_TYPES = ("Granulocytes", "Monocytes", "B", "NK", "T-CD4", "T-CD8")
MARKERS_PER_SIDE = 25  # per blood cell type: the most hyper- and the most hypomethylated windows against the others
FRACTION_VARIANCE = 0.9  # leading fraction PCs held: enough for this share of the variance, at most 3
CACHE = "read_entropy_counts"


# ----------------------------------------------------------------------------- pure helpers
def _ones(width: int) -> np.ndarray:
    """Methylated CpGs of every epiallele code (most significant bit = first CpG)."""
    return np.array([bin(c).count("1") for c in range(2**width)])


def methylated_runs(width: int) -> np.ndarray:
    """(2^width, width): for each epiallele code, how many stretches of i + 1 consecutive CpGs are all methylated."""
    out = np.zeros((2**width, width), dtype=int)
    for code in range(2**width):
        bits = [(code >> (width - 1 - k)) & 1 for k in range(width)]
        for length in range(1, width + 1):
            out[code, length - 1] = sum(all(bits[s : s + length]) for s in range(width - length + 1))
    return out


def epiallele_metrics(counts: np.ndarray, width: int = WIDTH, min_reads: int = MIN_READS) -> dict[str, np.ndarray]:
    """Per window (row of ``(windows, 2^width)`` epiallele counts): reads, mean, entropy (bits per CpG),
    epipolymorphism and MHL; NaN below ``min_reads``."""
    counts = np.asarray(counts, dtype=float)
    reads = counts.sum(axis=1)
    p = counts / np.maximum(reads, 1)[:, None]
    with np.errstate(divide="ignore", invalid="ignore"):
        entropy = -np.where(p > 0, p * np.log2(p), 0.0).sum(axis=1) / width
    lengths = np.arange(1, width + 1)
    stretches = width - lengths + 1  # stretches of each length on one read of the window
    full = p @ methylated_runs(width) / stretches  # share of fully methylated stretches, per length
    out = {
        "reads": reads,
        "mean": p @ _ones(width) / width,
        "entropy": entropy,
        "epipolymorphism": 1 - (p**2).sum(axis=1),
        "mhl": full @ lengths / lengths.sum(),
    }
    low = reads < min_reads
    for key in MEASURES[:-1]:
        out[key] = np.where(low, np.nan, out[key])
    return out


def epiallele_relative_entropy(
    counts: np.ndarray, reference: np.ndarray, width: int = WIDTH, min_reads: int = MIN_READS, a: float = KL_PSEUDO
) -> np.ndarray:
    """KL(p || p_ref) in bits per CpG per window; both share vectors get pseudo-count ``a`` on every epiallele.
    NaN where the sample or the reference has fewer than ``min_reads`` reads."""
    counts, reference = np.asarray(counts, dtype=float), np.asarray(reference, dtype=float)
    n, m = counts.sum(axis=-1), reference.sum(axis=-1)
    k = 2**width
    p = (counts / np.maximum(n, 1)[..., None] + a) / (1 + k * a)
    q = (reference / np.maximum(m, 1)[..., None] + a) / (1 + k * a)
    kl = (p * np.log2(p / q)).sum(axis=-1) / width
    return np.where((n < min_reads) | (m < min_reads), np.nan, kl)


def fe_age_slopes(
    y: np.ndarray,
    age: np.ndarray,
    female: np.ndarray,
    group: np.ndarray,
    min_samples: int = MIN_SAMPLES,
    eligible: int | None = None,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Per column of ``y`` (samples x windows, NaN = missing): OLS age slope with sex and ``group`` fixed effects held
    (Frisch-Waugh within groups, the observed samples of each column only), its SE and residual df. NaN where fewer
    than ``max(min_samples, 3/4 of the samples)`` are observed (3/4 of ``eligible`` samples, when given) or age does
    not vary within groups."""
    y = np.asarray(y, dtype=float)
    ok = np.isfinite(y)
    y0 = np.where(ok, y, 0.0)
    mask = ok.astype(float)
    codes = np.unique(group, return_inverse=True)[1]
    g = np.eye(codes.max() + 1)[codes]  # samples x groups
    n_g = g.T @ mask
    inv_n = np.where(n_g > 0, 1 / np.maximum(n_g, 1), 0.0)
    z = [np.asarray(age, dtype=float), np.asarray(female, dtype=float)]
    sz = [g.T @ (mask * v[:, None]) for v in z]
    sy = g.T @ y0

    def within(a_col, b_col, sa, sb):
        return (a_col * b_col).sum(axis=0) - (sa * sb * inv_n).sum(axis=0)

    zc = [mask * v[:, None] for v in z]
    saa = within(zc[0], z[0][:, None], sz[0], sz[0])
    sss = within(zc[1], z[1][:, None], sz[1], sz[1])
    sas = within(zc[0], z[1][:, None], sz[0], sz[1])
    say = within(zc[0], y0, sz[0], sy)
    ssy = within(zc[1], y0, sz[1], sy)
    syy = within(y0, y0, sy, sy)
    with np.errstate(divide="ignore", invalid="ignore"):
        det = saa * sss - sas**2
        use_sex = (sss > 1e-9) & (det > 1e-9 * np.maximum(saa * sss, 1e-300))
        b_age = np.where(use_sex, (sss * say - sas * ssy) / det, say / saa)
        b_sex = np.where(use_sex, (saa * ssy - sas * say) / det, 0.0)
        rss = syy - b_age * say - b_sex * ssy
        df = ok.sum(axis=0) - (n_g > 0).sum(axis=0) - 1 - use_sex
        var_age = np.where(use_sex, sss / det, 1 / saa)
        se = np.sqrt(np.maximum(rss, 0) / df * var_age)
    need = max(min_samples, int(0.75 * (y.shape[0] if eligible is None else eligible)))
    bad = (ok.sum(axis=0) < need) | ~(saa > 1e-9) | (df < 2)
    return np.where(bad, np.nan, b_age), np.where(bad, np.nan, se), np.where(bad, np.nan, df)


def cov_age_slopes(
    y: np.ndarray, age: np.ndarray, covariates: np.ndarray, min_samples: int = MIN_SAMPLES, eligible: int | None = None
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Per column of ``y`` (samples x windows, NaN = missing): OLS age slope in ``y ~ 1 + age + covariates`` over the
    column's observed samples, its SE and residual df (batched normal equations). NaN where the design is rank
    deficient on the observed samples or fewer than ``max(min_samples, 3/4 of the samples)`` are observed."""
    y = np.asarray(y, dtype=float)
    cov = np.asarray(covariates, dtype=float).reshape(y.shape[0], -1)
    cov = cov[:, np.ptp(cov, axis=0) > 1e-12] if cov.shape[1] else cov
    x = np.column_stack([np.ones(y.shape[0]), np.asarray(age, dtype=float), cov])
    x = x / np.maximum(np.abs(x).max(axis=0), 1e-12)  # columns on one scale: the batched solves stay conditioned
    unit = np.abs(np.asarray(age, dtype=float)).max()
    ok = np.isfinite(y)
    m = ok.astype(float)
    y0 = np.where(ok, y, 0.0)
    xtx = np.einsum("ij,ia,ib->jab", m, x, x)
    xty = np.einsum("ij,ia->ja", y0, x)
    yty = (y0 * y0).sum(axis=0)
    p = x.shape[1]
    rank = np.linalg.matrix_rank(xtx, hermitian=True)
    inv = np.linalg.pinv(xtx, hermitian=True)
    b = np.einsum("jab,jb->ja", inv, xty)
    rss = yty - np.einsum("ja,ja->j", b, xty)
    df = ok.sum(axis=0) - p
    with np.errstate(divide="ignore", invalid="ignore"):
        se = np.sqrt(np.maximum(rss, 0) / df * inv[:, 1, 1]) / unit
    need = max(min_samples, int(0.75 * (y.shape[0] if eligible is None else eligible)))
    bad = (rank < p) | (ok.sum(axis=0) < need) | (df < 2)
    return np.where(bad, np.nan, b[:, 1] / unit), np.where(bad, np.nan, se), np.where(bad, np.nan, df.astype(float))


def matched_sides(values: np.ndarray, member: np.ndarray, strata: np.ndarray) -> tuple[float, float]:
    """The two sides of the stratified contrast of ``stochasticity.compare_groups``: the mean over member windows
    (strata holding both kinds only), and the matched mean (each such stratum's non-member mean, weighted by its
    member windows). NaN with fewer than two matched members."""
    ok = np.isfinite(values)
    values, member, strata = values[ok], member[ok], strata[ok]
    codes = np.unique(strata, return_inverse=True)[1]
    k = codes.max() + 1 if codes.size else 0
    n_in = np.bincount(codes, weights=member.astype(float), minlength=k)
    n_all = np.bincount(codes, minlength=k).astype(float)
    s_in = np.bincount(codes, weights=values * member, minlength=k)
    s_all = np.bincount(codes, weights=values, minlength=k)
    mixed = (n_in > 0) & (n_in < n_all)
    if n_in[mixed].sum() < 2:
        return float("nan"), float("nan")
    rest = (s_all - s_in)[mixed] / (n_all - n_in)[mixed]
    w = n_in[mixed]
    return float(s_in[mixed].sum() / w.sum()), float((rest * w).sum() / w.sum())


def matched_gap(values: np.ndarray, member: np.ndarray, strata: np.ndarray) -> float:
    """Set side minus matched side of ``matched_sides``."""
    inside, matched = matched_sides(values, member, strata)
    return inside - matched


def bh_q(p: np.ndarray) -> np.ndarray:
    """Benjamini-Hochberg q values (NaN kept)."""
    p = np.asarray(p, dtype=float)
    out = np.full(p.shape, np.nan)
    ok = np.flatnonzero(np.isfinite(p))
    if ok.size == 0:
        return out
    order = ok[np.argsort(p[ok])]
    ranked = p[order] * ok.size / np.arange(1, ok.size + 1)
    out[order] = np.minimum(np.minimum.accumulate(ranked[::-1])[::-1], 1.0)
    return out


def sign_tail(draws: np.ndarray) -> tuple[float, int]:
    """Two-sided sign-tail summary of bootstrap draws, and the usable (finite) draws: twice the smaller of
    ``(1 + draws above 0) / (1 + usable)`` and ``(1 + draws at or below 0) / (1 + usable)``, capped at 1. NaN draws are
    left out (never counted as non-positive); the summary is never zero."""
    d = np.asarray(draws, dtype=float)
    d = d[np.isfinite(d)]
    if d.size == 0:
        return float("nan"), 0
    pos = int((d > 0).sum())
    tail = min(pos, d.size - pos)
    return float(min(1.0, 2 * (1 + tail) / (1 + d.size))), int(d.size)


def young_reference(
    counts: np.ndarray, age: np.ndarray, group: np.ndarray, ids: np.ndarray | None = None, q: float = YOUNG_QUANTILE
) -> np.ndarray:
    """(samples, windows, 2^w) pooled epiallele counts of each sample's young reference: the samples of its group at
    or below the group's ``q`` age quantile, leaving out the sample itself (every copy of it in a bootstrap resample,
    by ``ids``). A young sample scored against a reference holding its own reads would sit closer to it by sampling
    alone, and relative entropy would rise with age for that reason only."""
    ids = np.arange(counts.shape[0]) if ids is None else np.asarray(ids)
    ref = np.empty_like(counts)
    for g in np.unique(group):
        rows = np.flatnonzero(group == g)
        young = rows[age[rows] <= np.quantile(age[rows], q)]
        pooled = counts[young].sum(axis=0)
        for i in rows:
            copies = int((ids[young] == ids[i]).sum())
            ref[i] = pooled - copies * counts[i]
    return ref


def rarefy(counts: np.ndarray, depth: int, rng: np.random.Generator) -> np.ndarray:
    """Each row (last axis: epiallele counts) with at least ``depth`` fragments cut to exactly ``depth`` by sampling
    without replacement (multivariate hypergeometric, drawn category by category); rows below ``depth`` become 0."""
    counts = np.asarray(counts, dtype=np.int64)
    flat_in = counts.reshape(-1, counts.shape[-1])
    flat_out = np.zeros_like(flat_in)
    total = flat_in.sum(axis=1)
    left = np.where(total >= depth, depth, 0)
    pool = total.copy()
    for k in range(flat_in.shape[1] - 1):
        c = flat_in[:, k]
        draw = rng.hypergeometric(c, pool - c, left)
        flat_out[:, k] = draw
        left = left - draw
        pool = pool - c
    flat_out[:, -1] = left
    return flat_out.reshape(counts.shape)


def balanced_windows(
    reads: np.ndarray, age: np.ndarray, group: np.ndarray, depth: int, tol: float = BALANCE_TOL
) -> np.ndarray:
    """Windows whose share of samples reaching ``depth`` fragments differs by at most ``tol`` between the younger
    half (at or below its group's median age) and the older half of the groups with an age spread."""
    young = np.zeros(age.size, dtype=bool)
    used = np.zeros(age.size, dtype=bool)
    for g in np.unique(group):
        rows = np.flatnonzero(group == g)
        if np.ptp(age[rows]) == 0:
            continue
        used[rows] = True
        young[rows] = age[rows] <= np.median(age[rows])
    ok = np.asarray(reads) >= depth
    return np.abs(ok[used & young].mean(axis=0) - ok[used & ~young].mean(axis=0)) <= tol


def marker_windows(type_means: np.ndarray, per_side: int = MARKERS_PER_SIDE) -> np.ndarray:
    """Indices of windows that tell cell types apart, from ``(types, windows)`` mean methylation: per type, the
    ``per_side`` windows where it sits furthest above every other type and the ``per_side`` furthest below."""
    ok = np.isfinite(type_means).all(axis=0)
    idx = np.flatnonzero(ok)
    m = type_means[:, ok]
    chosen: set[int] = set()
    for t in range(m.shape[0]):
        others = np.delete(m, t, axis=0)
        for margin in (m[t] - others.max(axis=0), others.min(axis=0) - m[t]):
            order = np.argsort(-margin)[:per_side]
            chosen |= {int(idx[i]) for i in order if margin[i] > 0}
    return np.array(sorted(chosen), dtype=int)


def estimate_fractions(sample_means: np.ndarray, type_means: np.ndarray, markers: np.ndarray) -> np.ndarray:
    """(samples, types) fractions: non-negative least squares of each sample's mean methylation at ``markers`` on
    the type profiles (markers missing in the sample left out), scaled to sum to 1."""
    out = np.full((sample_means.shape[0], type_means.shape[0]), np.nan)
    ref = type_means[:, markers].T
    for i, y in enumerate(sample_means[:, markers]):
        ok = np.isfinite(y)
        if ok.sum() < 2 * type_means.shape[0]:
            continue
        f, _ = optimize.nnls(ref[ok], y[ok])
        if f.sum() > 0:
            out[i] = f / f.sum()
    return out


def fraction_components(fractions: np.ndarray, share: float = FRACTION_VARIANCE, most: int = 3) -> np.ndarray:
    """Leading principal-component scores of the centred fractions: as many as reach ``share`` of the variance,
    at most ``most``."""
    f = fractions - fractions.mean(axis=0)
    u, s, _ = np.linalg.svd(f, full_matrices=False)
    var = s**2 / max(float((s**2).sum()), 1e-300)
    k = int(min(most, np.searchsorted(np.cumsum(var), share) + 1))
    return u[:, :k] * s[:k]


def mixture_counts(
    type_shares: np.ndarray, fractions: np.ndarray, reads: np.ndarray, rng: np.random.Generator
) -> np.ndarray:
    """(samples, windows, 2^w) counts of samples made of cell types only: shares ``fractions @ type_shares`` (types,
    windows, 2^w; the same for every sample), ``reads`` fragments per sample and window drawn multinomially; windows
    with a missing type share get none."""
    shares = np.einsum("st,twk->swk", fractions, type_shares)
    shares = shares / shares.sum(axis=-1, keepdims=True)
    n = np.where(np.isfinite(shares).all(axis=-1), reads, 0).astype(np.int64)
    return rng.multinomial(n, np.nan_to_num(shares, nan=1 / shares.shape[-1]))


def putative_donors(samples: pl.DataFrame) -> np.ndarray:
    """One label per sample for the donor it may share with other samples: GEO gives no donor id for the atlas, so
    samples with the same age, sex, lab and organ (first word of the tissue) are taken as one donor."""
    organ = pl.col("tissue").fill_null("").str.split(" ").list.first().str.to_lowercase()
    key = pl.concat_str(
        [pl.col("age").cast(pl.String), pl.col("female").cast(pl.String), pl.col("lab").fill_null(""), organ],
        separator="|",
    )
    return samples.select(key.alias("key"))["key"].to_numpy().astype(str)


def cache_key(starts: np.ndarray, pat: Path) -> str:
    """Key of one PAT file's counts: the ordered window starts and the file's name, size and modification time."""
    st = Path(pat).stat()
    h = hashlib.sha1(np.ascontiguousarray(starts, dtype=np.int64).tobytes())
    h.update(f"w{WIDTH}|{Path(pat).name}|{st.st_size}|{int(st.st_mtime)}".encode())
    return h.hexdigest()[:16]


# ----------------------------------------------------------------------------- data
def read_samples(folder: Path) -> pl.DataFrame:
    """``sample``, ``title``, ``age``, ``female``, ``group`` (cell type | tissue), ``tissue``, ``lab`` and ``pat`` of
    every sample with an age."""
    rows = []
    for soft in sorted(folder.glob("GSM*.soft")):
        text = soft.read_text()
        age = re.search(r"age: (\d+(?:\.\d+)?)", text)
        sex = re.search(r"(?:gender|sex): (\w+)", text, flags=re.I)
        cell = re.search(r"cell type: (.+)", text)
        tissue = re.search(r"tissue: (.+)", text)
        lab = re.search(r"lab: (.+)", text)
        title = re.search(r"!Sample_title = (.+)", text)
        pats = sorted(folder.glob(f"{soft.stem}_*.hg38.pat.gz"))
        if age and pats:
            rows.append(
                {
                    "sample": soft.stem,
                    "title": title.group(1).strip() if title else "",
                    "age": float(age.group(1)),
                    "female": bool(sex and sex.group(1).lower().startswith("f")),
                    "group": f"{cell.group(1).strip()}|{tissue.group(1).strip() if tissue else ''}" if cell else "",
                    "tissue": tissue.group(1).strip() if tissue else "",
                    "lab": lab.group(1).strip() if lab else "",
                    "pat": str(pats[0]),
                }
            )
    return pl.DataFrame(rows)


def _load_counts(ctx: Context, samples: pl.DataFrame, starts: np.ndarray) -> np.ndarray:
    cache_dir = ctx.paths.interim / CACHE
    cache_dir.mkdir(parents=True, exist_ok=True)
    jobs = [
        (p, starts, str(cache_dir / f"{s}_w{WIDTH}_{starts.size}_{cache_key(starts, Path(p))}.npy"))
        for s, p in samples.select("sample", "pat").iter_rows()
    ]
    missing = [j for j in jobs if not Path(j[2]).exists()]
    if missing:
        print(f"read_entropy: counting {len(missing)} PAT files", file=sys.stderr, flush=True)
        with ProcessPoolExecutor(int(os.environ.get("SLURM_CPUS_PER_TASK", "2"))) as pool:
            list(pool.map(_count_one, missing))
    return np.stack([np.load(j[2]).astype(np.int32) for j in jobs])


# ----------------------------------------------------------------------------- analysis
def _metrics(counts: np.ndarray, min_reads: int) -> dict[str, np.ndarray]:
    m = epiallele_metrics(counts.reshape(-1, counts.shape[-1]), min_reads=min_reads)
    return {k: m[k].reshape(counts.shape[:2]) for k in (*MEASURES[:-1], "reads")}


def _fit(counts, base, age, female, group, cols, ids, min_reads, cov=None):
    """Slopes (and SE, df) of every measure over windows ``cols``, the young-reference mean level and its depth."""
    c = counts[:, cols]
    ref = young_reference(c, age, group, ids)
    kl = epiallele_relative_entropy(c, ref, min_reads=min_reads)
    ref_reads = ref.sum(axis=-1).astype(float)
    with np.errstate(divide="ignore", invalid="ignore"):
        young_mean = np.where(ref_reads >= min_reads, (ref @ _ones(WIDTH)) / WIDTH / ref_reads, np.nan)
        baseline = np.nanmean(young_mean, axis=0)
    reads = base["reads"][:, cols]
    values = {m: base[m][:, cols] for m in MEASURES[:-1]} | {
        "relative_entropy": kl,
        "log_depth": np.log1p(reads),
        "eligible": (reads >= min_reads).astype(float),
    }
    # a sample alone at its group's young end has no reference once it is left out: it cannot carry relative entropy
    scorable = int((ref_reads.sum(axis=1) > 0).sum())
    fits = {}
    for m, v in values.items():
        elig = scorable if m == "relative_entropy" else None
        if cov is None:
            fits[m] = fe_age_slopes(v, age, female, group, eligible=elig)
        else:
            fits[m] = cov_age_slopes(v, age, np.column_stack([female, cov]), eligible=elig)
    return fits, baseline, ref_reads


def _keeps(fits, baseline) -> dict[str, np.ndarray]:
    """Per measure, the windows (of the fitted columns) with a slope and a young level; the depth checks take the
    Shannon-entropy windows."""
    keep = {m: np.isfinite(fits[m][0]) & np.isfinite(baseline) for m in MEASURES}
    for m in DEPTH_MEASURES:
        keep[m] = keep["entropy"]
    return keep


def _contrasts(slopes, baseline, island, members, keep) -> dict[tuple[str, str], tuple[float, float, float, float]]:
    """(measure, set) -> set side, matched side, gap and gap in SD of the measure's slopes, over its own windows."""
    out = {}
    for m, s in slopes.items():
        k = keep[m]
        if k.sum() < 2:
            continue
        sd = float(np.nanstd(s[k], ddof=1))
        strata = np.char.add(island[k].astype(str), baseline_strata(baseline[k], 10).astype(str))
        for name, member in members.items():
            inside, matched = matched_sides(s[k], member[k], strata)
            out[(m, name)] = (inside, matched, inside - matched, (inside - matched) / sd)
    return out


def _variant(
    *,
    tag: str,
    label: str,
    resampling: str,
    counts: np.ndarray,
    min_reads: int,
    age,
    female,
    group,
    island,
    members,
    draw: Callable[[np.random.Generator], np.ndarray],
    n_boot: int,
    rng: np.random.Generator,
    cov: np.ndarray | None = None,
    window_mask: np.ndarray | None = None,
    point_only: bool = False,
) -> dict:
    """Point fit, contrasts and bootstrap of one variant: its gap rows and what the raw tables need."""
    base = _metrics(counts, min_reads)
    cols0 = np.arange(counts.shape[1]) if window_mask is None else np.flatnonzero(window_mask)
    fits, baseline, ref_reads = _fit(counts, base, age, female, group, cols0, None, min_reads, cov)
    keep0 = _keeps(fits, baseline)
    union = np.any([keep0[m] for m in MEASURES], axis=0)
    cols = cols0[union]
    keep = {m: k[union] for m, k in keep0.items()}
    slopes = {m: f[0][union] for m, f in fits.items()}
    isl, mem = island[cols], {k: v[cols] for k, v in members.items()}
    point = _contrasts(slopes, baseline[union], isl, mem, keep)
    boots: dict[tuple[str, str], list] = {k: [] for k in point}
    medians: dict[str, list] = {m: [] for m in MEASURES}
    for b in range(0 if point_only else n_boot):
        rows = np.sort(draw(rng))
        sub = {m: base[m][rows] for m in (*MEASURES[:-1], "reads")}
        c = None if cov is None else cov[rows]
        f, bl, _ = _fit(counts[rows], sub, age[rows], female[rows], group[rows], cols, rows, min_reads, c)
        s = {m: v[0] for m, v in f.items()}
        got = _contrasts(s, bl, isl, mem, keep)
        for k in point:
            boots[k].append(got.get(k, (np.nan,) * 4))
        for m in MEASURES:
            medians[m].append(np.nanmedian(s[m][keep[m]]) if np.isfinite(s[m][keep[m]]).any() else np.nan)
        if (b + 1) % 50 == 0:
            print(f"read_entropy: {tag}: {label}: bootstrap {b + 1}/{n_boot}", file=sys.stderr, flush=True)
    rows_out = []
    for (m, name), (inside, matched, gap, gap_sd) in point.items():
        arr = np.array(boots[(m, name)], dtype=float).reshape(-1, 4)
        usable = arr[np.isfinite(arr).all(axis=1)]
        enough = not point_only and usable.shape[0] >= n_boot / 2
        lo, hi = np.quantile(usable, [0.025, 0.975], axis=0) if enough else (np.full(4, np.nan), np.full(4, np.nan))
        p, n_ok = sign_tail(arr[:, 3])
        rows_out.append(
            {
                "source": tag,
                "variant": label,
                "resampling": resampling,
                "measure": m,
                "group": name,
                "windows": int(mem[name][keep[m]].sum()),
                "all_windows": int(keep[m].sum()),
                "set_slope": inside,
                "set_lo": lo[0],
                "set_hi": hi[0],
                "matched_slope": matched,
                "matched_lo": lo[1],
                "matched_hi": hi[1],
                "gap": gap,
                "gap_lo": lo[2],
                "gap_hi": hi[2],
                "gap_sd": gap_sd,
                "gap_sd_lo": lo[3],
                "gap_sd_hi": hi[3],
                "p_boot": p if enough else np.nan,
                "samples": int(age.size),
                "boots": 0 if point_only else n_boot,
                "boots_usable": n_ok,
            }
        )
    return {
        "rows": rows_out,
        "fits": fits,
        "keep0": keep0,
        "cols0": cols0,
        "baseline": baseline,
        "ref_reads": ref_reads,
        "base": base,
        "union": union,
        "slopes": slopes,
        "keep": keep,
        "medians": medians,
    }


def _young_rows(tag, age, group, donors, ref_reads, keep_kl) -> list[dict]:
    """Per group: the young reference's samples, putative donors and pooled depth over the relative-entropy windows."""
    out = []
    for g in np.unique(group):
        rows = np.flatnonzero(group == g)
        is_young = age[rows] <= np.quantile(age[rows], YOUNG_QUANTILE)
        depth = ref_reads[rows][:, keep_kl]
        out.append(
            {
                "source": tag,
                "group": str(g) if tag == "sorted cells" else "all",
                "samples": int(rows.size),
                "young_samples": int(is_young.sum()),
                "young_donors": int(np.unique(donors[rows[is_young]]).size),
                "young_ages": ", ".join(f"{a:g}" for a in age[rows[is_young]]),
                "scorable_samples": int((ref_reads[rows].sum(axis=1) > 0).sum()),
                "ref_reads_median_older": float(np.median(depth[~is_young])) if (~is_young).any() else np.nan,
                "ref_reads_median_young_left_out": float(np.median(depth[is_young])) if is_young.any() else np.nan,
                "windows": int(keep_kl.sum()),
            }
        )
    return out


def _blood_types(ctx: Context, folder: Path, starts: np.ndarray) -> tuple[list[str], np.ndarray] | None:
    """Pooled epiallele counts (types, windows, 2^w) of the sorted blood references, or None when absent."""
    refs = read_samples(folder) if folder.exists() else pl.DataFrame()
    if refs.is_empty():
        return None
    refs = refs.with_columns(kind=pl.col("title").str.extract(r"^Blood-(.+)-Z[0-9A-Z]+$")).filter(
        pl.col("kind").is_in(list(BLOOD_TYPES))
    )
    names = [t for t in BLOOD_TYPES if t in refs["kind"].to_list()]
    if len(names) < 3:
        return None
    counts = _load_counts(ctx, refs, starts)
    kinds = refs["kind"].to_numpy()
    print(f"read_entropy: blood references: {refs.height} samples, types {names}", file=sys.stderr, flush=True)
    return names, np.stack([counts[kinds == t].sum(axis=0) for t in names]).astype(np.int64)


def _source(ctx: Context, tag: str, folder: Path, windows: pl.DataFrame, probes, membership, out: Path, n_boot: int):
    samples = read_samples(folder).sort("age", "sample") if folder.exists() else pl.DataFrame()
    if samples.height < MIN_SAMPLES:
        return None
    starts = windows["start"].to_numpy()
    counts = _load_counts(ctx, samples, starts)
    age, female = samples["age"].to_numpy(), samples["female"].to_numpy()
    sorted_cells = tag == "sorted cells"
    group = samples["group"].to_numpy() if sorted_cells else np.zeros(samples.height, dtype=int)
    donors = putative_donors(samples) if sorted_cells else samples["sample"].to_numpy()
    print(f"read_entropy: {tag}: {samples.height} samples, {starts.size} windows", file=sys.stderr, flush=True)
    feature_all = windows["cpg"].to_numpy()
    island = (
        pl.DataFrame({"feature": feature_all})
        .join(probes.select(pl.col("cpg").alias("feature"), "cpg_island"), on="feature", how="left")["cpg_island"]
        .fill_null(False)
        .to_numpy()
    )
    members = {
        s: np.isin(feature_all, membership.filter(pl.col("group") == s)["feature"].to_numpy())
        for s in SETS
        if s in membership["group"].unique()
    }
    rows_of = [np.flatnonzero(group == g) for g in np.unique(group)]
    clusters = [np.flatnonzero(donors == d) for d in np.unique(donors)]

    def within_groups(rng):
        return np.concatenate([rng.choice(r, r.size, replace=True) for r in rows_of])

    def whole_clusters(rng):
        return np.concatenate([clusters[i] for i in rng.integers(0, len(clusters), len(clusters))])

    sample_level = "samples within cell-type groups (sample-level)" if sorted_cells else "donors"
    seed = 186_458 + len(tag)
    common = {
        "tag": tag, "age": age, "female": female, "group": group, "island": island, "members": members,
        "n_boot": n_boot,
    }  # fmt: skip
    raw = _variant(
        **common, label="raw", resampling=sample_level, counts=counts, min_reads=MIN_READS, draw=within_groups,
        rng=ctx.rng(seed),
    )  # fmt: skip
    gap_rows = list(raw["rows"])

    # windows table and age-estimability, raw counts
    union_cols = raw["cols0"][raw["union"]]
    frame = {
        "feature": feature_all[union_cols],
        "young_mean": raw["baseline"][raw["union"]],
        "median_reads": np.median(raw["base"]["reads"][:, union_cols], 0),
        "young_ref_reads_median": np.median(raw["ref_reads"][:, raw["union"]], 0),
    }
    for m, (b, se, df) in raw["fits"].items():
        k = raw["keep0"][m][raw["union"]]
        b, se, df = b[raw["union"]], se[raw["union"]], df[raw["union"]]
        with np.errstate(divide="ignore", invalid="ignore"):
            p = np.where(k, 2 * stats.t.sf(np.abs(b / se), df), np.nan)
        frame |= {f"{m}_slope": np.where(k, b, np.nan), f"{m}_q": bh_q(p)}
        if m in raw["base"]:
            frame[f"{m}_level"] = np.nanmean(raw["base"][m][:, union_cols], 0)
    suffix = SOURCES[tag] or "wbc"
    write_table(pl.DataFrame(frame), out / f"windows_{suffix}")
    within_age = np.concatenate([age[r] - age[r].mean() for r in rows_of])
    reads_med = np.median(raw["base"]["reads"][:, union_cols], axis=1)
    est_rows = []
    for m in MEASURES:
        arr = np.array(raw["medians"][m], dtype=float)
        arr = arr[np.isfinite(arr)]
        k = raw["keep"][m]
        q = frame[f"{m}_q"][k]
        est_rows.append(
            {
                "source": tag,
                "measure": m,
                "samples": samples.height,
                "groups_with_age_spread": int(sum(np.ptp(age[r]) > 0 for r in rows_of)),
                "within_group_age_sd": float(within_age.std(ddof=1)),
                "age_min": float(age.min()),
                "age_max": float(age.max()),
                "windows": int(k.sum()),
                "median_slope": float(np.nanmedian(raw["slopes"][m][k])),
                "median_slope_lo": float(np.quantile(arr, 0.025)) if arr.size else np.nan,
                "median_slope_hi": float(np.quantile(arr, 0.975)) if arr.size else np.nan,
                "boots_usable": int(arr.size),
                "q05_windows": int(np.sum(q < 0.05)),
                "share_q05": float(np.mean(q < 0.05)),
                "age_reads_r": float(np.corrcoef(within_age, np.concatenate([reads_med[r] for r in rows_of]))[0, 1]),
            }
        )
    young_rows = _young_rows(tag, age, group, donors, raw["ref_reads"], raw["keep0"]["relative_entropy"])
    del raw

    # rarefied counts
    for depth, balanced in RAREFY:
        rare = rarefy(counts, depth, ctx.rng(seed + 10 * depth)).astype(np.int32)
        mask = balanced_windows(counts.sum(axis=-1), age, group, depth) if balanced else None
        label = f"rarefied to {depth}" + (", balanced coverage" if balanced else "")
        v = _variant(
            **common, label=label, resampling=sample_level, counts=rare, min_reads=depth, draw=within_groups,
            rng=ctx.rng(seed + 10 * depth + 1), window_mask=mask,
        )  # fmt: skip
        gap_rows += v["rows"]
        del rare, v

    extra: dict[str, pl.DataFrame] = {}
    if sorted_cells:
        v = _variant(
            **common, label="raw", resampling="putative-donor clusters (age, sex, lab, organ)", counts=counts,
            min_reads=MIN_READS, draw=whole_clusters, rng=ctx.rng(seed + 3),
        )  # fmt: skip
        gap_rows += v["rows"]
        logo = []
        for g in np.unique(group):
            kept = group != g
            w = _variant(
                **(common | {"age": age[kept], "female": female[kept], "group": group[kept]}),
                label=f"without {g}", resampling="none", counts=counts[kept], min_reads=MIN_READS,
                draw=within_groups, rng=ctx.rng(seed), point_only=True,
            )  # fmt: skip
            logo += [r | {"left_out": str(g)} for r in w["rows"]]
        extra["leave_one_group_out"] = pl.DataFrame(logo)
        extra["putative_donors"] = pl.DataFrame(
            {"sample": samples["sample"], "title": samples["title"], "group": group, "putative_donor": donors}
        )
    if tag == "white blood cells":
        blood = _blood_types(ctx, folder / "blood_refs", starts)
        if blood is None:
            print("read_entropy: no sorted blood references; WBC composition not adjusted", file=sys.stderr)
        else:
            names, pooled = blood
            type_means = _metrics(pooled, MIN_READS)["mean"]
            markers = marker_windows(type_means)
            fractions = estimate_fractions(_metrics(counts, MIN_READS)["mean"], type_means, markers)
            pcs = fraction_components(fractions)
            extra["wbc_fractions"] = pl.DataFrame(
                {"sample": samples["sample"], "age": age, "female": female}
                | {f"fraction_{t}": fractions[:, i] for i, t in enumerate(names)}
                | {f"pc{i + 1}": pcs[:, i] for i in range(pcs.shape[1])}
                | {"markers": np.full(age.size, markers.size)}
            )
            v = _variant(
                **common, label=f"leukocyte-adjusted ({pcs.shape[1]} fraction PCs)", resampling=sample_level,
                counts=counts, min_reads=MIN_READS, draw=within_groups, rng=ctx.rng(seed + 4), cov=pcs,
            )  # fmt: skip
            gap_rows += v["rows"]
            shares = pooled / np.maximum(pooled.sum(axis=-1, keepdims=True), 1)
            covered = (pooled.sum(axis=-1) >= MIN_READS).all(axis=0)
            synthetic = mixture_counts(
                np.where(covered[None, :, None], shares, np.nan), fractions, counts.sum(axis=-1), ctx.rng(seed + 5)
            ).astype(np.int32)
            v = _variant(
                **common, label="composition-only mixture", resampling=sample_level, counts=synthetic,
                min_reads=MIN_READS, draw=within_groups, rng=ctx.rng(seed + 6), window_mask=covered,
            )  # fmt: skip
            gap_rows += v["rows"]
            del synthetic
    print(f"read_entropy: {tag}: done", file=sys.stderr, flush=True)
    return pl.DataFrame(gap_rows), pl.DataFrame(est_rows), pl.DataFrame(young_rows), extra


def run(ctx: Context) -> list[Path]:
    folder = ctx.paths.input / "geo" / GSE
    if not (ctx.paths.interim / INDEX).exists() or not folder.exists():
        print("read_entropy: CpG index or GSE186458 files missing. Nothing run.", file=sys.stderr)
        return []
    plots.apply_style()
    out = ctx.out_dir("read_entropy")
    features, probes, membership = _features(ctx.paths, ctx.rng(140_000))  # the stream read_level uses: same windows
    windows = probe_windows(ctx.paths, features)
    n_boot = 200 if ctx.scale.name == "paper" else 20
    results = [
        r
        for tag, sub in SOURCES.items()
        if (r := _source(ctx, tag, folder / sub if sub else folder, windows, probes, membership, out, n_boot))
        is not None
    ]
    gaps = pl.concat([r[0] for r in results], how="diagonal_relaxed")
    write_table(gaps, out / "gaps")
    write_table(pl.concat([r[1] for r in results]), out / "age_estimable")
    write_table(pl.concat([r[2] for r in results]), out / "young_reference")
    for r in results:
        for name, frame in r[3].items():
            write_table(frame, out / name)
    return [_figure(ctx, gaps)]


def _figure(ctx: Context, table: pl.DataFrame) -> Path:
    sources = table["source"].unique(maintain_order=True).to_list()
    fig, axes = plots.panels(2 * len(sources), ncols=len(sources), size=4.4)
    axes = np.ravel(axes)
    raw = table.filter((pl.col("variant") == "raw") & ~pl.col("resampling").str.starts_with("putative"))
    for k, source in enumerate(sources):
        ax = axes[k]
        sub = raw.filter(pl.col("source") == source)
        for i, name in enumerate(SETS):
            t = sub.filter((pl.col("group") == name) & pl.col("measure").is_in(list(MEASURES)))
            if t.is_empty():
                continue
            order = [MEASURES.index(m) for m in t["measure"]]
            y = np.array(order) + 0.22 * (i - 1)
            gap, lo, hi = (t[c].to_numpy() for c in ("gap_sd", "gap_sd_lo", "gap_sd_hi"))
            ax.errorbar(gap, y, xerr=[gap - lo, hi - gap], fmt="o", ms=4, color=plots.PALETTE[i], capsize=2,
                        label=name.replace("_", "-"))  # fmt: skip
        ax.set_yticks(range(len(MEASURES)), [m.replace("_", " ") for m in MEASURES] if k == 0 else [""] * len(MEASURES))
        ax.axvline(0, color=plots.MUTED, lw=0.8)
        ax.invert_yaxis()
        if k == 0:
            ax.legend(fontsize=7, loc="lower left")
        ax.set(xlabel="age slope, set minus matched windows (SD)", title=f"{source} (n = {sub['samples'][0]})")
        plots.label_panel(ax, "abcdef"[k])
        # Shannon entropy at housekeeping windows: both sides of the contrast, per variant
        ax = axes[len(sources) + k]
        t = table.filter(
            (pl.col("source") == source) & (pl.col("measure") == "entropy") & (pl.col("group") == "housekeeping")
        )
        labels = [
            f"{v}{' (donor clusters)' if r.startswith('putative') else ''}"
            for v, r in zip(t["variant"], t["resampling"], strict=True)
        ]
        for j, (side, color, name) in enumerate(
            (("set", plots.PALETTE[0], "housekeeping"), ("matched", plots.MUTED, "matched"))
        ):
            v, lo, hi = (t[c].to_numpy() for c in (f"{side}_slope", f"{side}_lo", f"{side}_hi"))
            y = np.arange(t.height) + 0.18 * (j - 0.5)
            ax.errorbar(v, y, xerr=[v - lo, hi - v], fmt="o", ms=4, color=color, capsize=2, label=name)
        ax.set_yticks(range(t.height), labels, fontsize=7)
        ax.axvline(0, color=plots.MUTED, lw=0.8)
        ax.invert_yaxis()
        if k == 0:
            ax.legend(fontsize=7, loc="lower left")
        ax.set(xlabel="Shannon-entropy age slope (bits per CpG per year)")
        plots.label_panel(ax, "abcdef"[len(sources) + k])
    fig.suptitle(
        "Read-level epiallele measures: age slopes at essential promoters and at matched windows "
        "(GSE186458, 95 % bootstrap)",
        x=0.02,
        y=1.02,
        ha="left",
        fontsize=11,
        fontweight="bold",
    )
    return plots.save(fig, ctx.figure_path("read_entropy"))
