# /// script
# requires-python = ">=3.11"
# dependencies = ["methylprep==1.7.1", "pandas", "pyarrow"]
# ///
"""Reprocess a GEO series from its raw IDATs into betas, SeSAMe-style (roadmap RM18 (3)).

methylprep's pipeline with ``sesame=True`` applies noob background correction, the non-linear dye-bias correction,
channel-switch inference and SeSAMe's quality mask, then writes betas. The result is saved as a parquet of
``cpg`` plus one column per sample, named by the GSM accession that starts every GEO IDAT file name.

    uv run --script scripts/reprocess_idats.py IDAT_DIR OUT.parquet [--batch 100]

Run as a Slurm job (it downloads the manifest once, through the proxy).
"""

from __future__ import annotations

import argparse
import re
import sys
from pathlib import Path

import pandas as pd


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("idats", type=Path)
    parser.add_argument("out", type=Path)
    parser.add_argument("--batch", type=int, default=100)
    args = parser.parse_args()
    import methylprep

    betas = methylprep.run_pipeline(
        args.idats,
        betas=True,
        make_sample_sheet=True,
        batch_size=args.batch,
        sesame=True,
        export=False,
        low_memory=True,
    )
    if isinstance(betas, list):
        betas = pd.concat(betas, axis=1)
    # GEO names IDATs GSMnnn_<sentrix id>_<position>_Grn.idat(.gz); methylprep names samples by the last two
    by_position = {}
    for f in args.idats.rglob("*_Grn.idat*"):
        m = re.match(r"(GSM\d+)_(\d+_R\d+C\d+)_Grn", f.name)
        if m:
            by_position[m.group(2)] = m.group(1)
    rename = {}
    for c in betas.columns:
        gsm = re.match(r"(GSM\d+)", str(c))
        pos = re.search(r"(\d+_R\d+C\d+)", str(c))
        rename[c] = gsm.group(1) if gsm else by_position.get(pos.group(1), str(c)) if pos else str(c)
    betas = betas.rename(columns=rename)
    betas.index.name = "cpg"
    betas = betas.reset_index()
    args.out.parent.mkdir(parents=True, exist_ok=True)
    betas.to_parquet(args.out, index=False)
    print(f"reprocess_idats: {betas.shape[0]} CpGs, {betas.shape[1] - 1} samples -> {args.out}", file=sys.stderr)


if __name__ == "__main__":
    main()
