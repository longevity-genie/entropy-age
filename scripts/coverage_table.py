#!/usr/bin/env -S uv run --script
# /// script
# requires-python = ">=3.11"
# dependencies = ["polars"]
# ///
"""Render docs/analysis_coverage.csv as Supplementary Table S5 between its markers in manuscript/supplement.md.

One row per methylation series, one column per analysis: the run that did it, or "not run". It shows how far each
claim reaches: a statement about what explains the gap rests only on the rows whose cell for that analysis is filled.
"""

from pathlib import Path

import polars as pl

ROOT = Path(__file__).resolve().parents[1]
START, END = "<!-- coverage:start -->", "<!-- coverage:end -->"
COLUMNS = {
    "series": "series",
    "species": "species",
    "tissue_or_cells": "tissue or cells",
    "accession": "accession",
    "platform": "platform",
    "processing": "processing",
    "units": "samples",
    "gap": "gap",
    "off_floor": "off the 0 % floor",
    "chip_or_batch_held": "chip or batch held",
    "chromatin_state": "chromatin state matched",
    "expression_same_samples": "expression from the same samples",
    "accessibility": "accessibility",
    "gain_loss_rates": "gain and loss rates",
    "interval": "interval resamples",
}


def main() -> None:
    frame = pl.read_csv(ROOT / "docs" / "analysis_coverage.csv", infer_schema_length=0).fill_null("")
    head = "| " + " | ".join(COLUMNS.values()) + " |\n|" + "---|" * len(COLUMNS) + "\n"
    body = "".join(
        "| " + " | ".join(str(r[c]).replace("|", "/") for c in COLUMNS) + " |\n" for r in frame.iter_rows(named=True)
    )
    path = ROOT / "manuscript" / "supplement.md"
    text = path.read_text()
    i, j = text.index(START) + len(START), text.index(END)
    path.write_text(text[:i] + "\n\n" + head + body + "\n" + text[j:])


if __name__ == "__main__":
    main()
