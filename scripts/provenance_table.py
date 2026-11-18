#!/usr/bin/env -S uv run --script
# /// script
# requires-python = ">=3.11"
# dependencies = ["polars"]
# ///
"""Render docs/provenance.csv as Supplementary Table S4 between its markers in manuscript/supplement.md (RM42)."""

from pathlib import Path

import polars as pl

ROOT = Path(__file__).resolve().parents[1]
START, END = "<!-- provenance:start -->", "<!-- provenance:end -->"
COLUMNS = {
    "analysis": "analysis",
    "runs": "runs",
    "claim": "claim",
    "tissue": "tissue / species",
    "methylation": "methylation (study, accession)",
    "platform": "platform",
    "processing": "processing",
    "second_track": "second track",
    "second_source": "its source",
    "batch": "batch",
    "second_samples": "second track from",
    "cell_type": "cell or tissue type of second track",
    "notes": "notes",
}


def render(frame: pl.DataFrame) -> str:
    head = "| " + " | ".join(COLUMNS.values()) + " |\n|" + "---|" * len(COLUMNS) + "\n"
    rows = (
        "| " + " | ".join(str(r[c] or "").replace("|", "/") for c in COLUMNS) + " |\n"
        for r in frame.iter_rows(named=True)
    )
    return head + "".join(rows)


def main() -> None:
    frame = pl.read_csv(ROOT / "docs" / "provenance.csv", infer_schema_length=0).fill_null("")
    path = ROOT / "manuscript" / "supplement.md"
    text = path.read_text()
    i, j = text.index(START) + len(START), text.index(END)
    path.write_text(text[:i] + "\n\n" + render(frame) + "\n" + text[j:])


if __name__ == "__main__":
    main()
