"""Transcriptome helpers ported from the BiT age part of the paper (Fig. 2).

The C. elegans count matrix behind Fig. 2 is not distributed upstream (only its
sample sheet, Supplementary Table 1), so these are provided for users who
process their own RNA-seq data the same way.
"""

from __future__ import annotations

import numpy as np
from scipy.special import erf, erfinv


def normalise_counts(counts: np.ndarray) -> np.ndarray:
    """``log10(count + 1)`` then per-sample division by the maximum (samples in rows)."""
    logged = np.log10(np.asarray(counts, dtype=float) + 1.0)
    peak = logged.max(axis=1, keepdims=True)
    return np.divide(logged, peak, out=np.zeros_like(logged), where=peak > 0)


def binarize_by_median(x: np.ndarray, quantile: float = 0.5) -> np.ndarray:
    """BiT age binarisation: per gene, 1 if above its median of non-zero values, else 0.

    Port of upstream ``make_binary``: zeros are masked before taking the
    quantile across samples, and masked or below-median values become 0.
    """
    x = np.asarray(x, dtype=float)
    masked = np.where(x == 0, np.nan, x)
    threshold = np.nanquantile(masked, quantile, axis=0)
    with np.errstate(invalid="ignore"):
        ratio = masked / threshold
    return np.where(np.isnan(ratio) | (ratio <= 1), 0.0, 1.0)


def correct_biological_age(hours: np.ndarray, mean: float = 372.0, spread: float = 192.0 / 3.0) -> np.ndarray:
    """Upstream ``correct_Bio_Age`` (lifespan-distribution correction of worm biological age, hours).

    Written as in the original:
    ``y - erfinv(0.5 - (erf((mean - y) / spread / sqrt 2) / 2 + 0.5) / 2) * 2 sqrt 2 * spread``.
    """
    y = np.asarray(hours, dtype=float)
    cdf = erf((mean - y) / spread / np.sqrt(2)) / 2 + 0.5
    return y - erfinv(0.5 - cdf / 2) * np.sqrt(2) * 2 * spread
