"""Is there a molecule-level counterpart of the shield? Read-level epiallele measures and age (GSE186458)

The paper's shield is measured on CpG averages over cells. A sequencing read reports several neighbouring CpGs of one
DNA molecule, so the mix of epialleles in a window says how much the cells disagree with each other. Loyfer et al.
2023 (GSE186458; whole-genome bisulfite, wgbstools PAT reads) give three sets with ages: white blood cells and
plasma cfDNA of the same 23 donors (21-75 years) and 101 sorted cell-type samples (27 cell-type and tissue groups,
19-83 years). Windows are the ``WIDTH`` = 4 consecutive CpGs starting at each 450k promoter probe's CpG (hg38), read
off every read that covers all four (``read_level.epiallele_counts``; its per-sample counts are reused). A window
counts in a sample with at least ``MIN_READS`` such reads. Per sample and window, from the 16 epiallele shares p:

- **mean** methylation over the window's CpGs;
- **entropy**, Shannon entropy of the epiallele distribution in bits divided by the window's CpGs (the methylation
  entropy of Xie et al. 2011, Nucleic Acids Res 39:4099);
- **epipolymorphism**, 1 - sum p^2 (Landan et al. 2012, Nat Genet 44:1207);
- **MHL**, methylation haplotype load (Guo et al. 2017, Nat Genet 49:635): sum over lengths i = 1..4 of
  i x (share of fully methylated stretches of i CpGs among all stretches of i CpGs on the reads), divided by sum i;
- **relative entropy**, KL(p || p_young) in bits per CpG: the read-level form of the paper's headline measure, with
  p_young the pooled epiallele shares of the young reference (ages at or below the 20 % quantile; in the sorted
  cells, of each cell-type group), the sample itself left out, and a pseudo-count of ``KL_PSEUDO`` on every share.

Per window and measure, the age slope: OLS with sex held (when it varies) and, for the sorted cells, each cell-type
group as a fixed effect, so only age differences within a cell type count. Windows need a value in three quarters
of the samples (and ``MIN_SAMPLES``). The gap is the stratified difference of the slopes between a set's windows
(housekeeping, core-essential, non-essential) and all other windows of the same CpG-island status and young-level
decile (the matching of ``tissue_contrast._with_strata``; the decile is taken from the young reference's
mean methylation), in raw units and in SD of the slopes. Intervals: percentile bootstrap over samples (donors; in
the sorted cells, samples resampled within their group so the groups stay), slopes, young reference, strata and gap
recomputed in each resample. Whether age is estimable at all: per measure, the median slope over windows with its
interval and the share of windows with an age effect at BH q < 0.05.

``run(ctx)`` as a Slurm job (not yet in ``experiments.FIGURES``); reads about 2.5 GB of cached counts.
"""

from __future__ import annotations

import os
import re
import sys
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

import numpy as np
import polars as pl
from scipy import stats

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
SETS = ("housekeeping", "core_essential", "non_essential")


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


def matched_gap(values: np.ndarray, member: np.ndarray, strata: np.ndarray) -> float:
    """Mean over member windows of (value minus the mean of non-members in its stratum), strata holding both only;
    the stratified contrast of ``stochasticity.compare_groups``."""
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
        return float("nan")
    rest = (s_all - s_in)[mixed] / (n_all - n_in)[mixed]
    return float(((s_in[mixed] / n_in[mixed] - rest) * n_in[mixed]).sum() / n_in[mixed].sum())


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


# ----------------------------------------------------------------------------- data
def read_samples(folder: Path) -> pl.DataFrame:
    """``sample``, ``age``, ``female``, ``group`` (cell type | tissue) and ``pat`` of every sample with an age."""
    rows = []
    for soft in sorted(folder.glob("GSM*.soft")):
        text = soft.read_text()
        age = re.search(r"age: (\d+(?:\.\d+)?)", text)
        sex = re.search(r"(?:gender|sex): (\w+)", text, flags=re.I)
        cell = re.search(r"cell type: (.+)", text)
        tissue = re.search(r"tissue: (.+)", text)
        pats = sorted(folder.glob(f"{soft.stem}_*.hg38.pat.gz"))
        if age and pats:
            rows.append(
                {
                    "sample": soft.stem,
                    "age": float(age.group(1)),
                    "female": bool(sex and sex.group(1).lower().startswith("f")),
                    "group": f"{cell.group(1).strip()}|{tissue.group(1).strip() if tissue else ''}" if cell else "",
                    "pat": str(pats[0]),
                }
            )
    return pl.DataFrame(rows)


def _load_counts(ctx: Context, samples: pl.DataFrame, starts: np.ndarray) -> np.ndarray:
    cache_dir = ctx.paths.interim / "read_level"
    cache_dir.mkdir(parents=True, exist_ok=True)
    jobs = [
        (p, starts, str(cache_dir / f"{s}_w{WIDTH}_{starts.size}.npy"))
        for s, p in samples.select("sample", "pat").iter_rows()
    ]
    missing = [j for j in jobs if not Path(j[2]).exists()]
    if missing:
        print(f"read_entropy: counting {len(missing)} PAT files", file=sys.stderr, flush=True)
        with ProcessPoolExecutor(int(os.environ.get("SLURM_CPUS_PER_TASK", "2"))) as pool:
            list(pool.map(_count_one, missing))
    return np.stack([np.load(j[2]).astype(np.int32) for j in jobs])


# ----------------------------------------------------------------------------- analysis
def _fit(counts, base: dict[str, np.ndarray], age, female, group, cols=None, ids=None):
    """Slopes (and SE, df) of every measure over windows ``cols``, and the young-reference mean level."""
    cols = np.arange(counts.shape[1]) if cols is None else cols
    ref = young_reference(counts[:, cols], age, group, ids)
    kl = epiallele_relative_entropy(counts[:, cols], ref)
    young_mean = epiallele_metrics(ref.reshape(-1, ref.shape[-1]))["mean"].reshape(ref.shape[:2])
    values = {m: base[m][:, cols] for m in MEASURES[:-1]} | {"relative_entropy": kl}
    # a sample alone at its group's young end has no reference once it is left out: it cannot carry relative entropy
    scorable = int((ref.reshape(ref.shape[0], -1).sum(axis=1) > 0).sum())
    fits = {
        m: fe_age_slopes(v, age, female, group, eligible=scorable if m == "relative_entropy" else None)
        for m, v in values.items()
    }
    with np.errstate(invalid="ignore"):
        baseline = np.nanmean(young_mean, axis=0)
    return fits, baseline


def _gaps(slopes: dict[str, np.ndarray], baseline, island, members) -> dict[tuple[str, str], tuple[float, float]]:
    strata = np.char.add(island.astype(str), baseline_strata(baseline, 10).astype(str))
    out = {}
    for m, s in slopes.items():
        sd = float(np.nanstd(s, ddof=1))
        for name, member in members.items():
            gap = matched_gap(s, member, strata)
            out[(m, name)] = (gap, gap / sd)
    return out


def _source(ctx: Context, tag: str, folder: Path, windows: pl.DataFrame, probes, membership, out: Path, n_boot: int):
    samples = read_samples(folder).sort("age", "sample") if folder.exists() else pl.DataFrame()
    if samples.height < MIN_SAMPLES:
        return None
    starts = windows["start"].to_numpy()
    counts = _load_counts(ctx, samples, starts)
    age, female = samples["age"].to_numpy(), samples["female"].to_numpy()
    group = samples["group"].to_numpy() if tag == "sorted cells" else np.zeros(samples.height, dtype=int)
    print(f"read_entropy: {tag}: {samples.height} samples, {starts.size} windows", file=sys.stderr, flush=True)
    base = epiallele_metrics(counts.reshape(-1, counts.shape[-1]))
    base = {m: base[m].reshape(counts.shape[:2]) for m in (*MEASURES[:-1], "reads")}
    fits, baseline = _fit(counts, base, age, female, group)
    keep = np.flatnonzero(np.all([np.isfinite(f[0]) for f in fits.values()], axis=0) & np.isfinite(baseline))
    feature = windows["cpg"].to_numpy()[keep]
    island = (
        pl.DataFrame({"feature": feature})
        .join(probes.select(pl.col("cpg").alias("feature"), "cpg_island"), on="feature", how="left")["cpg_island"]
        .fill_null(False)
        .to_numpy()
    )
    members = {
        s: np.isin(feature, membership.filter(pl.col("group") == s)["feature"].to_numpy())
        for s in SETS
        if s in membership["group"].unique()
    }
    slopes = {m: f[0][keep] for m, f in fits.items()}
    frame = {"feature": feature, "young_mean": baseline[keep], "median_reads": np.median(base["reads"][:, keep], 0)}
    for m, (b, se, df) in fits.items():
        p = 2 * stats.t.sf(np.abs(b[keep] / se[keep]), df[keep])
        frame |= {f"{m}_slope": b[keep], f"{m}_q": bh_q(p)}
        if m in base:
            frame[f"{m}_level"] = np.nanmean(base[m][:, keep], 0)
    suffix = SOURCES[tag] or "wbc"
    write_table(pl.DataFrame(frame), out / f"windows_{suffix}")
    point = _gaps(slopes, baseline[keep], island, members)

    # bootstrap over samples (within groups for the sorted cells)
    rng = ctx.rng(186_458 + len(tag))
    rows_of = [np.flatnonzero(group == g) for g in np.unique(group)]
    boots: dict[tuple[str, str], list] = {k: [] for k in point}
    medians: dict[str, list] = {m: [] for m in MEASURES}
    for b in range(n_boot):
        rows = np.sort(np.concatenate([rng.choice(r, r.size, replace=True) for r in rows_of]))
        sub = {m: base[m][rows] for m in MEASURES[:-1]}
        f, bl = _fit(counts[rows], sub, age[rows], female[rows], group[rows], keep, rows)
        s = {m: v[0] for m, v in f.items()}
        for k, v in _gaps(s, bl, island, members).items():
            boots[k].append(v)
        for m in MEASURES:
            medians[m].append(np.nanmedian(s[m]))
        if (b + 1) % 25 == 0:
            print(f"read_entropy: {tag}: bootstrap {b + 1}/{n_boot}", file=sys.stderr, flush=True)
    gap_rows = []
    for (m, name), (gap, gap_sd) in point.items():
        arr = np.array(boots[(m, name)], dtype=float)
        lo, hi = np.nanquantile(arr, [0.025, 0.975], axis=0) if arr.size else ([np.nan] * 2, [np.nan] * 2)
        share_pos = np.nanmean(arr[:, 1] > 0) if arr.size else np.nan
        gap_rows.append(
            {
                "source": tag,
                "measure": m,
                "group": name,
                "windows": int(members[name].sum()),
                "gap": gap,
                "gap_lo": lo[0],
                "gap_hi": hi[0],
                "gap_sd": gap_sd,
                "gap_sd_lo": lo[1],
                "gap_sd_hi": hi[1],
                "p_boot": float(2 * min(share_pos, 1 - share_pos)),
                "samples": samples.height,
                "boots": n_boot,
            }
        )
    within_age = np.concatenate([age[r] - age[r].mean() for r in rows_of])
    reads_med = np.median(base["reads"][:, keep], axis=1)
    est_rows = []
    for m in MEASURES:
        arr = np.array(medians[m])
        est_rows.append(
            {
                "source": tag,
                "measure": m,
                "samples": samples.height,
                "groups_with_age_spread": int(sum(np.ptp(age[r]) > 0 for r in rows_of)),
                "within_group_age_sd": float(within_age.std(ddof=1)),
                "age_min": float(age.min()),
                "age_max": float(age.max()),
                "windows": int(keep.size),
                "median_slope": float(np.nanmedian(slopes[m])),
                "median_slope_lo": float(np.nanquantile(arr, 0.025)) if arr.size else np.nan,
                "median_slope_hi": float(np.nanquantile(arr, 0.975)) if arr.size else np.nan,
                "share_q05": float(np.mean(frame[f"{m}_q"] < 0.05)),
                "age_reads_r": float(np.corrcoef(within_age, np.concatenate([reads_med[r] for r in rows_of]))[0, 1]),
            }
        )
    print(f"read_entropy: {tag}: {keep.size} windows fitted", file=sys.stderr, flush=True)
    return pl.DataFrame(gap_rows), pl.DataFrame(est_rows)


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
    gaps = pl.concat([r[0] for r in results])
    estimable = pl.concat([r[1] for r in results])
    write_table(gaps, out / "gaps")
    write_table(estimable, out / "age_estimable")
    return [_figure(ctx, gaps)]


def _figure(ctx: Context, table: pl.DataFrame) -> Path:
    sources = table["source"].unique(maintain_order=True).to_list()
    fig, axes = plots.panels(len(sources), ncols=len(sources), size=4.4)
    for k, source in enumerate(sources):
        ax = axes[k]
        sub = table.filter(pl.col("source") == source)
        for i, name in enumerate(SETS):
            t = sub.filter(pl.col("group") == name)
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
        plots.label_panel(ax, "abc"[k])
    fig.suptitle(
        "Read-level epiallele measures: age-slope gap at essential promoters (GSE186458, 95 % bootstrap)",
        x=0.02,
        y=1.04,
        ha="left",
        fontsize=11,
        fontweight="bold",
    )
    return plots.save(fig, ctx.figure_path("read_entropy"))
