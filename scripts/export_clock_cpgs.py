"""Export a rebuilt series' betas (every CpG probe) for the clock panel of scripts/clocks_biolearn.py (protocol P4).

``uv run python scripts/export_clock_cpgs.py GSE51032 OUT.parquet``: the passing arrays of the SeSAMe rebuild (one
array per sample, detection-masked values left missing so the clock's own imputation fills them and the measured
share is reported). Beside it, ``<OUT stem>.manifest.csv`` records each exported sample's physical array
(``sample``, ``array``: chip and position, the file prefix where either is missing), so the clock analyses can
check at their join that the scores come from the array the current rebuild keeps (Codex review 5, finding 7).
"""

import sys
from pathlib import Path

import polars as pl

from entropy_age.data import reprocessed
from entropy_age.experiments.clock_drift_controls import array_ids, manifest_path
from entropy_age.paths import ProjectPaths


def main(gse: str, out: str) -> None:
    paths = ProjectPaths.default()
    betas = reprocessed.read_betas(paths, gse).rename({"probe": "cpg"}).filter(pl.col("cpg").str.contains(r"^c[gh]"))
    betas.write_parquet(out)
    qc = reprocessed.read_qc(paths, gse)  # the same rebuild (ENTROPY_AGE_SESAME_DIR) the betas come from
    arrays = array_ids(qc).filter(pl.col("sample").is_in(betas.columns[1:]))
    arrays.write_csv(manifest_path(Path(out)))
    print(f"{gse}: {betas.height} CpGs, {betas.width - 1} arrays (every clock's CpGs are read from here)")


if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2])
