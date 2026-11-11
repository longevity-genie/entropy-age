"""Reference chromatin states for the 450k CpGs: Roadmap Epigenomics chromHMM segmentations.

The 15-state core-marks model (Roadmap Epigenomics Consortium 2015), one segmentation per reference
epigenome, on GRCh37 like the 450k manifest's positions. Each CpG takes the state of the segment it
falls in. The states summarise five histone marks (H3K4me3, H3K4me1, H3K36me3, H3K27me3, H3K9me3), so
matching gene-set CpGs to background CpGs of the same state asks whether chromatin alone accounts for a
difference between them (roadmap RM7).
"""

from __future__ import annotations

import gzip
from dataclasses import dataclass
from pathlib import Path

import numpy as np
import polars as pl

from entropy_age.data.download import download
from entropy_age.paths import ProjectPaths

URL = (
    "https://egg2.wustl.edu/roadmap/data/byFileType/chromhmmSegmentations/ChmmModels/coreMarks/jointModel/final/"
    "{epigenome}_15_coreMarks_mnemonics.bed.gz"
)
#: reference epigenomes used here: peripheral-blood mononuclear cells, neutrophils, dorsolateral prefrontal cortex
EPIGENOMES = {
    "E062": "blood mononuclear cells",
    "E030": "blood neutrophils",
    "E073": "prefrontal cortex",
}
STATES = (
    "1_TssA",
    "2_TssAFlnk",
    "3_TxFlnk",
    "4_Tx",
    "5_TxWk",
    "6_EnhG",
    "7_Enh",
    "8_ZNF/Rpts",
    "9_Het",
    "10_TssBiv",
    "11_BivFlnk",
    "12_EnhBiv",
    "13_ReprPC",
    "14_ReprPCWk",
    "15_Quies",
)


@dataclass(frozen=True, slots=True)
class ChromatinLayout:
    paths: ProjectPaths

    def segmentation(self, epigenome: str) -> Path:
        return self.paths.input / "chromatin" / f"{epigenome}_15_coreMarks_mnemonics.bed.gz"


def fetch_chromatin(paths: ProjectPaths, *, overwrite: bool = False) -> list[Path]:
    layout = ChromatinLayout(paths)
    return [
        download(URL.format(epigenome=e), layout.segmentation(e), overwrite=overwrite, timeout=600) for e in EPIGENOMES
    ]


def read_segmentation(path: Path) -> pl.DataFrame:
    """``chr`` (without the ``chr`` prefix), ``start`` (0-based), ``end`` and ``state`` per segment."""
    with gzip.open(path, "rt") as handle:
        frame = pl.read_csv(
            handle.read().encode(),
            separator="\t",
            has_header=False,
            new_columns=["chr", "start", "end", "state"],
            schema_overrides={"chr": pl.String, "start": pl.Int64, "end": pl.Int64, "state": pl.String},
        )
    unknown = set(frame["state"].unique().to_list()) - set(STATES)
    if unknown:
        raise ValueError(f"{path.name}: states outside the 15-state model: {sorted(unknown)}")
    return frame.with_columns(pl.col("chr").str.strip_prefix("chr"))


def cpg_states(probes: pl.DataFrame, segments: pl.DataFrame) -> pl.Series:
    """The state of the segment holding each probe (1-based ``pos``); null where no segment covers it."""
    out = np.full(probes.height, None, dtype=object)
    chroms = probes["chr"].to_numpy()
    positions = probes["pos"].to_numpy()
    for (chrom,), seg in segments.sort("chr", "start").group_by("chr"):
        rows = np.flatnonzero(chroms == chrom)
        if rows.size == 0:
            continue
        starts, ends, states = seg["start"].to_numpy(), seg["end"].to_numpy(), seg["state"].to_numpy()
        zero_based = positions[rows] - 1
        i = np.searchsorted(starts, zero_based, side="right") - 1
        inside = (i >= 0) & (zero_based < ends[np.clip(i, 0, None)])
        out[rows[inside]] = states[i[inside]]
    return pl.Series("state", out.tolist(), dtype=pl.String)
