"""Export a rebuilt series' betas at the CpGs of the clocks scored by scripts/clocks_biolearn.py (protocol P4).

``uv run python scripts/export_clock_cpgs.py GSE51032 OUT.parquet``: the passing arrays of the SeSAMe rebuild (one
array per sample, detection-masked values left missing so the clock's own imputation fills them and the measured
share is reported).
"""

import sys

import polars as pl

from entropy_age.data import reprocessed
from entropy_age.data.clocks import ClockLayout, read_clock
from entropy_age.paths import ProjectPaths


def main(gse: str, out: str) -> None:
    paths = ProjectPaths.default()
    cpgs = set(read_clock(ClockLayout(paths).file("PhenoAge"))["cpg"].to_list())
    betas = reprocessed.read_betas(paths, gse).rename({"probe": "cpg"}).filter(pl.col("cpg").is_in(list(cpgs)))
    betas.write_parquet(out)
    print(f"{gse}: {betas.height} of {len(cpgs)} PhenoAge CpGs on the array, {betas.width - 1} arrays")


if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2])
