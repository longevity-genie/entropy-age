# /// script
# requires-python = ">=3.10,<3.13"
# dependencies = ["biolearn @ git+https://github.com/bio-learn/biolearn@0d714f5a0c0ad5743b6315f3edab1f86112d394f", "torch", "seaborn", "matplotlib", "requests", "pyyaml", "appdirs", "scikit-learn", "pyarrow"]
# ///
"""Score published clocks with biolearn itself (pinned commit, its documented missing-CpG imputation): protocol P4.

Input: a parquet of betas, ``cpg`` plus one column per sample (written by scripts/export_clock_cpgs.py). Output: a
CSV with one row per sample and clock: the prediction and the share of the clock's absolute weight that was measured
(not imputed) in that sample. Runs in its own environment: ``uv run --script scripts/clocks_biolearn.py IN OUT [META]``.
META, a CSV of ``sample, age, female``, lets the GrimAge models run: they take age and sex as inputs beside the CpGs.
biolearn's GrimAge sums each component over the CpGs present and skips missing values (a missing CpG counts as 0), so
here a missing value is filled with that CpG's cohort mean first, and the coverage rule applies as for every clock.
biolearn's gallery wraps PhenoAge in its imputation; clocks without a wrapper refuse missing CpGs and are skipped for
samples that lack any.
"""

import os
import sys

import numpy as np
import pandas as pd
from biolearn.data_library import GeoData
from biolearn.model_gallery import ModelGallery

MIN_WEIGHT_SHARE = 0.9  # a clock is kept when at least this share of its absolute weight is measured ...
MIN_ARRAY_SHARE = 0.9  # ... in at least this share of the arrays (protocol change 2)
HYBRID_ROW_SHARE = 0.8  # biolearn's hybrid imputation replaces rows measured in fewer arrays than this


def site_weights(model) -> pd.Series:
    """Absolute weight per required CpG; equal weights where the model has no coefficient table."""
    clock = getattr(model, "clock", model)
    coef = getattr(clock, "coefficients", None)
    if isinstance(coef, pd.DataFrame) and coef.shape[1] >= 1:
        w = coef.iloc[:, 0].drop(index=[i for i in coef.index if not str(i).startswith(("cg", "ch"))], errors="ignore")
        return w.abs()
    sites = list(clock.methylation_sites())
    return pd.Series(1.0, index=sites)


def grimage_weights(model) -> pd.Series:
    coef = model.coefficients
    cpg = coef[coef["var"].astype(str).str.startswith(("cg", "ch"))]
    return cpg.groupby("var")["beta"].apply(lambda b: b.abs().sum())


def main(src: str, dst: str, meta_path: str | None = None) -> None:
    betas = pd.read_parquet(src).set_index("cpg")
    meta = None
    if meta_path:
        meta = pd.read_csv(meta_path).set_index("sample").reindex(betas.columns)
        meta = pd.DataFrame({"age": meta["age"].astype(float), "sex": np.where(meta["female"].astype(bool), 0, 1)},
                            index=betas.columns)  # biolearn: sex 0 is female
    gallery = ModelGallery()
    data = GeoData(pd.DataFrame(index=betas.columns), dnam=betas)
    rows, summary = [], []
    from biolearn.model import model_definitions

    for name, entry in model_definitions.items():
        only = os.environ.get("CLOCKS_ONLY")  # comma-separated name prefixes, for a quick check
        if only and not name.startswith(tuple(only.split(","))):
            continue
        if entry.get("species") != "Human":
            continue
        grim = entry["model"]["type"] == "GrimageModel"
        if grim and meta is None:
            summary.append({"clock": name, "status": "GrimageModel needs META (age, sex)"})
            continue
        if not grim and "Methylation" not in entry["model"]["type"]:  # listed, not silently skipped
            summary.append({"clock": name, "status": f"model type {entry['model']['type']} not scored"})
            continue
        try:
            model = gallery.get(name)
            weight = grimage_weights(getattr(model, "clock", model)) if grim else site_weights(model)
        except Exception as err:  # a model the gallery cannot build here is listed, not scored
            summary.append({"clock": name, "status": f"not built: {type(err).__name__}"})
            continue
        if weight.sum() == 0:
            summary.append({"clock": name, "status": "no weights"})
            continue
        present = betas.reindex(weight.index).notna()
        # biolearn's hybrid imputation replaces a CpG row wholesale where it is measured in under 80 % of arrays,
        # observed entries included: such rows do not count as measured
        kept_rows = present.mean(axis=1) >= (0.0 if grim else HYBRID_ROW_SHARE)  # GrimAge: entry-wise fill below
        measured = present.loc[kept_rows].mul(weight[kept_rows], axis=0).sum(axis=0) / weight.sum()
        array_share = float((measured >= MIN_WEIGHT_SHARE).mean())
        if array_share < MIN_ARRAY_SHARE:
            summary.append({"clock": name, "status": "coverage", "arrays_covered": array_share,
                            "median_weight_share": float(measured.median())})
            continue
        try:
            if grim:
                sub = betas.reindex(weight.index).dropna(how="all")  # a CpG off the array stays out (coverage)
                sub = sub.apply(lambda r: r.fillna(r.mean()), axis=1)  # cohort mean; GrimAge would count a gap as 0
                pred = model.predict(GeoData(meta.copy(), dnam=sub.copy()))["DNAmGrimAge"]
            else:
                pred = model.predict(data)
        except Exception as err:
            summary.append({"clock": name, "status": f"predict failed: {type(err).__name__}: {str(err)[:80]}"})
            continue
        values = pred.iloc[:, 0] if isinstance(pred, pd.DataFrame) else pred
        for sample in betas.columns:
            rows.append({"sample": sample, "clock": name, "output": entry.get("output", ""),
                         "prediction": float(values.loc[sample]),
                         "measured_weight_share": float(measured.loc[sample])})
        summary.append({"clock": name, "status": "kept", "arrays_covered": array_share,
                        "median_weight_share": float(measured.median()), "output": entry.get("output", "")})
    pd.DataFrame(rows).to_csv(dst, index=False)
    pd.DataFrame(summary).to_csv(dst.removesuffix(".tmp").replace(".csv", "_panel.csv"), index=False)
    kept = [s["clock"] for s in summary if s["status"] == "kept"]
    print(f"{len(kept)} clocks kept of {len(summary)}: {', '.join(kept)}")


if __name__ == "__main__":
    main(*sys.argv[1:4])
