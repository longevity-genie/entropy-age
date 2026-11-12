"""The genome-wide CpG index of wgbstools (Loyfer et al. 2023) for hg38, built from a genome FASTA.

wgbstools numbers every CpG of the reference from 1, chromosome after chromosome, and its ``.beta`` and ``.pat``
files address CpGs by that number, not by position. This module rebuilds the numbering from the GRCh38 FASTA so
that array probes and promoters can be found in those files: ``CHROMOSOMES`` in natural order, each CpG at the
position of its C (1-based). ``check`` compares the total with a ``.beta`` file's length (two bytes per CpG).
"""

from __future__ import annotations

import gzip
from pathlib import Path

import numpy as np
import polars as pl

CHROMOSOMES = tuple(f"chr{c}" for c in [*range(1, 23), "X", "Y", "M"])


def cpg_positions(sequence: str) -> np.ndarray:
    """1-based positions of the C of every CG in ``sequence`` (case-insensitive)."""
    s = np.frombuffer(sequence.upper().encode(), dtype=np.uint8)
    hits = np.flatnonzero((s[:-1] == ord("C")) & (s[1:] == ord("G")))
    return hits + 1


def build_index(fasta: Path, chromosomes=CHROMOSOMES) -> pl.DataFrame:
    """``chr``, ``pos`` (1-based C) and ``index`` (1-based, genome-wide, wgbstools order)."""
    wanted = set(chromosomes)
    seqs: dict[str, list[str]] = {}
    name = None
    opener = gzip.open if str(fasta).endswith(".gz") else open
    with opener(fasta, "rt") as handle:
        for line in handle:
            if line.startswith(">"):
                name = line[1:].split()[0]
                if name in wanted:
                    seqs[name] = []
                continue
            if name in seqs:
                seqs[name].append(line.rstrip("\n"))
    frames, start = [], 1
    for chrom in chromosomes:
        if chrom not in seqs:
            continue
        pos = cpg_positions("".join(seqs.pop(chrom)))
        frames.append(pl.DataFrame({"chr": chrom, "pos": pos, "index": np.arange(start, start + pos.size)}))
        start += pos.size
    return pl.concat(frames)


def check(index: pl.DataFrame, beta: Path) -> tuple[int, int]:
    """CpGs in the index, and in a wgbstools ``.beta`` file (uint8 methylated and covered counts per CpG)."""
    return index.height, beta.stat().st_size // 2


def read_beta(beta: Path) -> tuple[np.ndarray, np.ndarray]:
    """Methylated and covered read counts per CpG index (0-based arrays; index i is CpG i + 1)."""
    raw = np.fromfile(beta, dtype=np.uint8).reshape(-1, 2)
    return raw[:, 0].astype(np.int32), raw[:, 1].astype(np.int32)
