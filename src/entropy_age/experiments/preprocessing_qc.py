"""QC, batch structure and agreement with deposited betas for every series rebuilt from IDATs (roadmap RM41).

For each series in ``data/interim/sesame`` (``scripts/sesame_preprocess.R``):

* **Sample QC:** SeSAMe's detection rate (pOOBAH), mean intensity and bisulfite-conversion index per sample;
  samples detected below ``reprocessed.MIN_DETECTED`` fail and are listed.
* **Batch:** the 10 leading principal components of the passing samples' betas (the 20,000 most variable complete
  probes); for chip, array row and position, the share of each component's variance it explains (eta squared),
  summed over components weighted by their variance: how much of the main structure is technical.
* **Design confounding:** for every characteristic GEO records for the series (age, group, tissue, sex ...), the
  share of its variance chip explains: 1 means the variable cannot be separated from chip in this series (the
  2026-11-18 audit's finding, made routine).
* **Agreement:** where the project also holds the authors' deposited betas, per-sample Pearson r between deposited
  and rebuilt betas over shared probes, and per-probe r across samples (median), matched by GSM or by chip_position.

``entropy-age figures preprocessing_qc``.
"""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import polars as pl

from entropy_age.data import reprocessed
from entropy_age.data.cohorts import _characteristics
from entropy_age.experiments.common import Context, write_table

N_PCS = 10
TOP_PROBES = 20_000
DEPOSITED = {  # GSE -> the project's table of the authors' betas
    "GSE74193": "jaffe",
    "GSE279437": "huang",
    "GSE197512": "endicott",
    "GSE179847": "sturm",
}


def eta_squared(values: np.ndarray, labels: np.ndarray) -> float:
    """Share of ``values``' variance explained by the categorical ``labels`` (0 when labels do not vary)."""
    total = ((values - values.mean()) ** 2).sum()
    if total == 0 or np.unique(labels).size < 2:
        return 0.0
    within = sum(((values[labels == k] - values[labels == k].mean()) ** 2).sum() for k in np.unique(labels))
    return float(1 - within / total)


def pcs(x: np.ndarray, k: int = N_PCS) -> tuple[np.ndarray, np.ndarray]:
    """Sample scores on the ``k`` leading components of samples x probes ``x`` and each component's variance share."""
    c = x - x.mean(axis=0)
    u, s, _ = np.linalg.svd(c, full_matrices=False)
    var = s**2 / (s**2).sum()
    k = min(k, u.shape[1])
    return u[:, :k] * s[:k], var[:k]


def batch_shares(scores: np.ndarray, var: np.ndarray, factors: dict[str, np.ndarray]) -> dict[str, float]:
    """Per factor: variance-weighted eta squared over the components (share of the leading structure)."""
    return {
        name: float(sum(v * eta_squared(scores[:, i], lab) for i, v in enumerate(var)) / var.sum())
        for name, lab in factors.items()
    }


def confounding(design: pl.DataFrame, chips: np.ndarray) -> list[dict]:
    """For each design column: the share of its variance chip explains (numeric as is, categories as indicators)."""
    rows = []
    for col in design.columns:
        v = design[col]
        if v.null_count() == v.len() or v.n_unique() < 2:
            continue
        num = v.cast(pl.Float64, strict=False)
        if num.null_count() <= 0.1 * v.len() and num.n_unique() > 5:
            ok = num.is_not_null().to_numpy()
            r2 = eta_squared(num.to_numpy()[ok], chips[ok])
            kind = "numeric"
        else:
            labs = v.cast(pl.String).fill_null("NA").to_numpy()
            r2 = float(np.mean([eta_squared((labs == k).astype(float), chips) for k in np.unique(labs)]))
            kind = f"{np.unique(labs).size} levels"
        rows.append({"variable": col, "kind": kind, "r2_on_chip": r2})
    return rows


def agreement(rebuilt: pl.DataFrame, qc: pl.DataFrame, deposited: pl.DataFrame) -> dict:
    """Per-sample and per-probe Pearson r between deposited and rebuilt betas, samples matched by GSM or position."""
    key = dict(zip(qc["sample"], (qc["chip"] + "_" + qc["position"]).fill_null(""), strict=True))
    pairs = []
    for col in rebuilt.columns[1:]:
        if col in deposited.columns:
            pairs.append((col, col))
        elif key.get(col) in deposited.columns:
            pairs.append((col, key[col]))
    if len(pairs) < 3:
        return {"matched_samples": len(pairs)}
    dep = deposited.select(
        pl.col(deposited.columns[0]).alias("probe"), *(pl.col(b).alias(f"dep::{a}") for a, b in pairs)
    )
    j = rebuilt.select("probe", *(a for a, _ in pairs)).join(dep, on="probe")
    a = j.select(a for a, _ in pairs).to_numpy().astype(float)
    b = j.select(f"dep::{a}" for a, _ in pairs).to_numpy().astype(float)
    ok = ~(np.isnan(a).any(axis=1) | np.isnan(b).any(axis=1))
    a, b = a[ok], b[ok]
    per_sample = [np.corrcoef(a[:, i], b[:, i])[0, 1] for i in range(a.shape[1])]
    ac, bc = a - a.mean(axis=1, keepdims=True), b - b.mean(axis=1, keepdims=True)
    den = np.sqrt((ac**2).sum(axis=1) * (bc**2).sum(axis=1))
    per_probe = np.where(den > 0, (ac * bc).sum(axis=1) / np.where(den > 0, den, 1), np.nan)
    return {
        "matched_samples": len(pairs),
        "shared_probes": int(ok.sum()),
        "sample_r_median": float(np.median(per_sample)),
        "sample_r_min": float(np.min(per_sample)),
        "probe_r_median": float(np.nanmedian(per_probe)),
        "mean_abs_diff": float(np.abs(a - b).mean()),
    }


def run(ctx: Context) -> list[Path]:
    out = ctx.out_dir("preprocessing_qc")
    root = ctx.paths.interim / "sesame"
    series = sorted(p.name for p in root.iterdir()) if root.exists() else []
    qc_rows, batch_rows, conf_rows, agree_rows = [], [], [], []
    for gse in series:
        if not reprocessed.available(ctx.paths, gse):
            continue
        qc = reprocessed.read_qc(ctx.paths, gse)
        qc_rows.append(qc.with_columns(series=pl.lit(gse)))
        betas = reprocessed.read_betas(ctx.paths, gse)
        samples = betas.columns[1:]
        q = qc.filter(pl.col("sample").is_in(samples)).with_columns(pl.col("sample").cast(pl.String))
        q = pl.DataFrame({"sample": samples}).join(q, on="sample", how="left")
        x = betas.drop("probe").to_numpy().T.astype(float)
        x = x[:, ~np.isnan(x).any(axis=0)]
        if x.shape[1] > TOP_PROBES:
            x = x[:, np.argsort(x.var(axis=0))[-TOP_PROBES:]]
        chips = q["chip"].fill_null("NA").to_numpy()
        pos = q["position"].fill_null("NA").to_numpy()
        row = np.array([p[:3] for p in pos])
        scores, var = pcs(x)
        shares = batch_shares(scores, var, {"chip": chips, "row": row, "position": pos})
        batch_rows.append({"series": gse, "samples": len(samples), "probes_in_pca": x.shape[1], **shares})
        matrix = ctx.paths.input / "geo" / gse / f"{gse}_series_matrix.txt.gz"
        if matrix.exists():
            design = pl.DataFrame({"sample": samples}).join(_characteristics(matrix), on="sample", how="left")
            conf_rows += [{"series": gse, **r} for r in confounding(design.drop("sample"), chips)]
        if gse in DEPOSITED:
            dep_path = ctx.paths.interim / "cohorts" / f"{DEPOSITED[gse]}_betas.parquet"
            if dep_path.exists():
                agree_rows.append({"series": gse, **agreement(betas, qc, pl.read_parquet(dep_path))})
        print(f"preprocessing_qc: {gse}: {len(samples)} passing samples", file=sys.stderr, flush=True)
    if qc_rows:
        write_table(pl.concat(qc_rows, how="diagonal"), out / "sample_qc")
    write_table(pl.DataFrame(batch_rows), out / "batch")
    write_table(pl.DataFrame(conf_rows), out / "design_confounding")
    write_table(pl.DataFrame(agree_rows), out / "agreement")
    return []
