"""Do essential promoters stay uniform molecule by molecule with age? Read-level WGBS (roadmap RM22)

Arrays give each CpG's average over cells; a read gives several CpGs of one molecule from one cell. Loyfer et al.
2023 (GSE186458) sequenced white blood cells of 23 donors aged 21-75 (whole-genome bisulfite, 150 bp pairs), in
wgbstools' PAT format: per read, its first CpG's genome-wide index and its calls (``C`` methylated, ``T`` not,
``.`` unknown). For each 450k promoter probe, the window of ``WIDTH`` CpGs starting at the probe's CpG (hg38,
``data.cpg_index``) is read off every read covering it whole, as one of 2^WIDTH epialleles. Per donor and window:

- mean methylation; Shannon **entropy** of the epiallele distribution (bits, divided by ``WIDTH``);
- **PDR**, the share of discordant reads (some CpGs methylated, some not).

Per window, the age slope of entropy and PDR over the donors (sex held), and the set (essential, housekeeping)
against look-alikes of the same island status and mean-level decile (``tissue_contrast``'s matching), in SD. A
shield at the molecule level means a smaller rise of discordance where the array saw a smaller rise of spread.
Caveat: white blood cells are a mixture whose myeloid share grows with age; one donor's cells differing from each
other and cell types differing are not told apart here.

``entropy-age figures read_level`` (as a Slurm job, after the CpG index and the 23 WBC files).
"""

from __future__ import annotations

import gzip
import os
import re
import sys
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

import numpy as np
import polars as pl

from entropy_age import plots
from entropy_age.experiments.common import Context, write_table
from entropy_age.experiments.tissue_contrast import _features, _matched_sd, _with_strata

GSE = "GSE186458"
WIDTH = 4
MIN_READS = 10
MIN_DONORS = 18
METRICS = ("entropy_slope", "pdr_slope", "mean_entropy", "mean_pdr")
MANIFEST = "annotation/zhou/HM450.hg38.manifest.tsv.gz"
INDEX = "annotation/cpg_index_hg38.parquet"


def read_donors(folder: Path) -> pl.DataFrame:
    """``sample``, ``age``, ``female`` and ``pat`` path for every WBC sample with a SOFT record and a PAT file."""
    rows = []
    for soft in sorted(folder.glob("GSM*.soft")):
        text = soft.read_text()
        age = re.search(r"age: (\d+(?:\.\d+)?)", text)
        sex = re.search(r"gender: (\w+)", text)
        pats = sorted(folder.glob(f"{soft.stem}_*.hg38.pat.gz"))
        if age and pats:
            rows.append(
                {
                    "sample": soft.stem,
                    "age": float(age.group(1)),
                    "female": bool(sex and sex.group(1).lower() == "female"),
                    "pat": str(pats[0]),
                }
            )
    return pl.DataFrame(rows)


def epiallele_counts(pat: Path, starts: np.ndarray, width: int = WIDTH) -> np.ndarray:
    """(windows, 2^width) counts of the epialleles of reads covering each window ``[start, start + width)``."""
    row = {int(s): i for i, s in enumerate(starts)}
    counts = np.zeros((starts.size, 2**width), dtype=np.int64)
    opener = gzip.open if str(pat).endswith(".gz") else open
    with opener(pat, "rt") as handle:
        for line in handle:
            _, s, calls, n = line.split("\t")
            s = int(s)
            for a in range(s, s + len(calls) - width + 1):
                i = row.get(a)
                if i is None:
                    continue
                piece = calls[a - s : a - s + width]
                if "." in piece:
                    continue
                code = 0
                for c in piece:
                    code = 2 * code + (c == "C")
                counts[i, code] += int(n)
    return counts


def window_metrics(counts: np.ndarray, width: int = WIDTH) -> dict[str, np.ndarray]:
    """Per window: reads, mean methylation, epiallele entropy (bits per CpG) and the share of discordant reads."""
    reads = counts.sum(axis=1)
    p = counts / np.maximum(reads, 1)[:, None]
    with np.errstate(divide="ignore", invalid="ignore"):
        entropy = -np.nansum(np.where(p > 0, p * np.log2(p), 0.0), axis=1) / width
    ones = np.array([bin(c).count("1") for c in range(2**width)])
    mean = (p * ones).sum(axis=1) / width
    concordant = p[:, 0] + p[:, -1]
    out = {"reads": reads, "mean": mean, "entropy": entropy, "pdr": 1 - concordant}
    for k in ("mean", "entropy", "pdr"):
        out[k] = np.where(reads >= MIN_READS, out[k], np.nan)
    return out


def age_slopes(values: np.ndarray, age: np.ndarray, female: np.ndarray) -> np.ndarray:
    """Per window (column), the OLS age slope with sex held, over donors with a value; NaN with too few donors."""
    out = np.full(values.shape[1], np.nan)
    for j in range(values.shape[1]):
        ok = ~np.isnan(values[:, j])
        if ok.sum() < MIN_DONORS:
            continue
        cols = [np.ones(ok.sum()), age[ok]]
        if np.unique(female[ok]).size > 1:
            cols.append(female[ok].astype(float))
        coef, *_ = np.linalg.lstsq(np.column_stack(cols), values[ok, j], rcond=None)
        out[j] = coef[1]
    return out


def probe_windows(paths, features: list[str]) -> pl.DataFrame:
    """``cpg`` (probe) and ``start`` (its CpG's genome-wide index, hg38) for the probes found in the index."""
    manifest = (
        pl.read_csv(
            paths.input / MANIFEST, separator="\t", columns=["CpG_chrm", "CpG_beg", "Probe_ID"], null_values=["NA"]
        )
        .drop_nulls()
        .select(
            pl.col("Probe_ID").alias("cpg"),
            pl.col("CpG_chrm").alias("chr"),
            (pl.col("CpG_beg").cast(pl.Int64) + 1).alias("pos"),
        )
        .filter(pl.col("cpg").is_in(features))
    )
    index = pl.read_parquet(paths.interim / INDEX)
    return manifest.join(index, on=["chr", "pos"]).select("cpg", pl.col("index").alias("start")).sort("start")


def _count_one(job: tuple[str, np.ndarray, str]) -> str:
    pat, starts, cache = job
    if not Path(cache).exists():
        np.save(cache, epiallele_counts(Path(pat), starts))
    return cache


def run(ctx: Context) -> list[Path]:
    folder = ctx.paths.input / "geo" / GSE
    if not (ctx.paths.interim / INDEX).exists() or not folder.exists():
        print("read_level: CpG index or GSE186458 files missing. Nothing run.", file=sys.stderr)
        return []
    plots.apply_style()
    out = ctx.out_dir("read_level")
    donors = read_donors(folder).sort("age")
    features, probes, membership = _features(ctx.paths, ctx.rng(140_000))
    windows = probe_windows(ctx.paths, features)
    starts = windows["start"].to_numpy()
    cache_dir = ctx.paths.interim / "read_level"
    cache_dir.mkdir(parents=True, exist_ok=True)
    jobs = [
        (p, starts, str(cache_dir / f"{s}_w{WIDTH}_{starts.size}.npy"))
        for s, p in donors.select("sample", "pat").iter_rows()
    ]
    workers = int(os.environ.get("SLURM_CPUS_PER_TASK", "2"))
    print(f"read_level: {donors.height} donors, {starts.size} windows, {workers} workers", file=sys.stderr, flush=True)
    with ProcessPoolExecutor(workers) as pool:
        caches = list(pool.map(_count_one, jobs))
    per = [window_metrics(np.load(c)) for c in caches]
    age, female = donors["age"].to_numpy(), donors["female"].to_numpy()
    stack = {k: np.vstack([m[k] for m in per]) for k in ("mean", "entropy", "pdr", "reads")}
    metrics = (
        pl.DataFrame(
            {
                "feature": windows["cpg"],
                "baseline": np.nanmean(stack["mean"], axis=0),
                "mean_entropy": np.nanmean(stack["entropy"], axis=0),
                "mean_pdr": np.nanmean(stack["pdr"], axis=0),
                "entropy_slope": age_slopes(stack["entropy"], age, female),
                "pdr_slope": age_slopes(stack["pdr"], age, female),
                "median_reads": np.median(stack["reads"], axis=0),
            }
        )
        .drop_nulls()
        .filter(pl.all_horizontal(pl.col(METRICS).is_not_nan()))
    )
    write_table(metrics, out / "windows")
    table = _matched_sd(
        _with_strata(metrics, probes),
        membership.filter(pl.col("feature").is_in(metrics["feature"])),
        METRICS,
        1000,
        ctx.seed,
    )
    write_table(table.with_columns(donors=pl.lit(donors.height)), out / "gaps")
    print(f"read_level: {metrics.height} windows with coverage", file=sys.stderr, flush=True)
    return [_figure(ctx, table)]


def _figure(ctx: Context, table: pl.DataFrame) -> Path:
    fig, axes = plots.panels(1, ncols=1, size=4.6)
    ax = axes[0]
    groups = [g for g in ("core_essential", "housekeeping", "non_essential") if g in table["group"].unique()]
    for i, group in enumerate(groups):
        t = table.filter(pl.col("group") == group)
        ax.barh(
            np.arange(t.height) + 0.27 * i,
            t["matched_diff_sd"],
            height=0.25,
            color=plots.PALETTE[i],
            label=group.replace("_", " "),
        )
    ax.set_yticks(
        np.arange(len(METRICS)) + 0.27,
        [m.replace("_", " ") for m in table.filter(pl.col("group") == groups[0])["metric"].to_list()],
    )
    ax.axvline(0, color=plots.MUTED, lw=0.8)
    ax.invert_yaxis()
    ax.legend(fontsize=7)
    ax.set(xlabel="set minus look-alikes (SD)", title="Molecule-level disorder, white blood cells (Loyfer 2023)")
    return plots.save(fig, ctx.figure_path("read_level"))
