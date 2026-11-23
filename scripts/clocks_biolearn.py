# /// script
# requires-python = ">=3.10,<3.13"
# dependencies = ["biolearn @ git+https://github.com/bio-learn/biolearn@0d714f5a0c0ad5743b6315f3edab1f86112d394f", "torch", "seaborn", "matplotlib", "requests", "pyyaml", "appdirs", "scikit-learn", "pyarrow"]
# ///
"""Score published clocks with biolearn itself (pinned commit, its documented missing-CpG imputation): protocol P4.

Input: a parquet of betas, ``cpg`` plus one column per sample (written by scripts/export_clock_cpgs.py). Output: a
CSV with one row per sample and clock: the prediction and the share of the clock's absolute weight that was measured
(not imputed) in that sample. Runs in its own environment: ``uv run --script scripts/clocks_biolearn.py IN OUT``.
biolearn's gallery wraps PhenoAge in its imputation; clocks without a wrapper refuse missing CpGs and are skipped for
samples that lack any.
"""

import sys

import numpy as np
import pandas as pd
from biolearn.data_library import GeoData
from biolearn.model_gallery import ModelGallery

CLOCKS = ("PhenoAge",)


def main(src: str, dst: str) -> None:
    betas = pd.read_parquet(src).set_index("cpg")
    gallery = ModelGallery()
    rows = []
    for name in CLOCKS:
        model = gallery.get(name)
        clock = getattr(model, "clock", model)
        coef = clock.coefficients
        coef = coef.drop(index=[i for i in coef.index if i == "intercept"])
        weight = coef.iloc[:, 0].abs()
        present = betas.reindex(weight.index)
        measured = present.notna().mul(weight, axis=0).sum(axis=0) / weight.sum()
        data = GeoData(pd.DataFrame(index=betas.columns), dnam=betas)
        pred = model.predict(data)
        col = pred.columns[0] if isinstance(pred, pd.DataFrame) else None
        values = pred[col] if col else pred
        for sample in betas.columns:
            rows.append({"sample": sample, "clock": name, "prediction": float(values.loc[sample]),
                         "measured_weight_share": float(measured.loc[sample])})
    pd.DataFrame(rows).to_csv(dst, index=False)
    print(f"{len(rows)} predictions; median measured weight share "
          f"{np.median([r['measured_weight_share'] for r in rows]):.3f}")


if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2])
