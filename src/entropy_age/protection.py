"""The per-person protection score (roadmap RM6, ``docs/PER_PERSON_SCORE.md``).

Every cohort test compares age slopes of feature groups, which needs many people. The per-person
score asks the same question of one sample: how far its gene-set CpGs have moved from a reference
state, compared with CpGs of the same island status and reference level in *the same* sample.

For a reference with per-CpG level ``m`` and a sample ``x``, the deviation is
``d = |x - m| / sqrt(m (1 - m))`` and the entropy gain is ``H(x) - H(m)``. A set's score ``P`` is the
stratum-weighted mean of ``d`` over its CpGs minus that over all other CpGs of the same stratum,
weights the set's share per stratum (the contrast ``compare_groups`` takes over features); ``B`` is
the matched background deviation and ``F = P / B`` the relative gap. A negative ``P`` means the
sample's set CpGs sit closer to the reference than matched CpGs. Shared sample-wide noise (poor
conversion, a noisy array) raises ``d`` at set and matched CpGs alike and cancels in ``P``.

The sampling error of one sample's ``P`` comes from resampling *genes* within the set (CpGs of one
promoter are correlated, so resampling CpGs would understate it).
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from entropy_age.stochasticity import baseline_strata, binary_entropy

#: reference levels are clipped to this range before the binomial spread is taken
LEVEL_CLIP = 0.01


def deviations(x: np.ndarray, reference: np.ndarray) -> np.ndarray:
    """``|x - m| / sqrt(m (1 - m))`` per sample and CpG, ``m`` clipped to [0.01, 0.99]."""
    m = np.clip(reference, LEVEL_CLIP, 1 - LEVEL_CLIP)
    return np.abs(x - reference) / np.sqrt(m * (1 - m))


def entropy_gain(x: np.ndarray, reference: np.ndarray) -> np.ndarray:
    """``H(x) - H(m)`` per sample and CpG, in bits."""
    return binary_entropy(x) - binary_entropy(reference)


def _codes(strata: np.ndarray) -> tuple[np.ndarray, int]:
    codes = np.unique(strata, return_inverse=True)[1]
    return codes, int(codes.max()) + 1


def contrast_weights(weight_in: np.ndarray, strata: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """The two feature weight vectors whose products with ``d`` give ``P`` and ``B``.

    ``weight_in`` is how many times each feature counts in the set (0 for a feature outside it; a
    resample of genes can count a CpG more than once). Within stratum ``k`` the set's features get
    ``w_k * c_j / n_k`` and the other features ``-w_k / n_out_k``, with ``w_k`` the set's share of its
    counts in ``k``. Strata holding only set or only other features are dropped, as in
    ``compare_groups``. Returns ``(p_weights, b_weights)``: ``P = d @ p_weights``, ``B = d @ b_weights``.
    """
    codes, k = _codes(strata)
    outside = weight_in == 0
    n_in = np.bincount(codes, weights=weight_in.astype(float), minlength=k)
    n_out = np.bincount(codes, weights=outside.astype(float), minlength=k)
    usable = (n_in > 0) & (n_out > 0)
    if not usable.any():
        raise ValueError("no stratum holds both set and other features")
    share = np.where(usable, n_in, 0.0) / n_in[usable].sum()
    with np.errstate(divide="ignore", invalid="ignore"):
        per_in = np.where(usable, share / n_in, 0.0)
        per_out = np.where(usable, share / n_out, 0.0)
    b = np.where(outside, per_out[codes], 0.0)
    p = weight_in * per_in[codes] - b
    return p, b


def matched_scores(d: np.ndarray, in_group: np.ndarray, strata: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """Per sample: the set-minus-matched gap ``P`` and the matched background deviation ``B``."""
    p, b = contrast_weights(np.asarray(in_group, dtype=float), strata)
    return d @ p, d @ b


@dataclass(frozen=True, slots=True)
class Reference:
    """The young state every sample is scored against: per-CpG level, stratum and composition fit.

    ``loadings`` and ``coefficients`` (both ``None`` without composition removal) carry the leading
    components of the reference's background CpGs: a sample's background deviations are projected on
    ``loadings`` to give its component scores, and ``scores @ coefficients`` is removed from every
    CpG before ``d`` is taken.
    """

    level: np.ndarray
    strata: np.ndarray
    samples: int
    background: np.ndarray
    spread: np.ndarray | None = None
    loadings: np.ndarray | None = None
    coefficients: np.ndarray | None = None

    def adjust(self, x: np.ndarray) -> np.ndarray:
        """``x`` with the reference-fitted composition components removed (unchanged without them)."""
        if self.loadings is None:
            return x
        z = (x[:, self.background] - self.level[self.background]) / self.spread
        return x - (z @ self.loadings) @ self.coefficients


def fit_reference(x: np.ndarray, island: np.ndarray, background: np.ndarray, *, components: int = 0) -> Reference:
    """The reference of the samples ``x`` (rows): mean level, island x level-decile stratum, components.

    ``background`` marks the CpGs outside every set; with ``components > 0`` the leading principal
    components of their standardised reference values are fitted, and every CpG's loading on them
    estimated by least squares in the reference, so a sample's composition can be read from its
    background CpGs and removed from all of them (``Reference.adjust``).
    """
    x = np.asarray(x, dtype=float)
    n = x.shape[0]
    level = x.mean(axis=0)
    strata = np.char.add(np.asarray(island).astype(str), baseline_strata(level, 10).astype(str))
    background = np.asarray(background, dtype=bool)
    if components == 0:
        return Reference(level, strata, n, background)
    if not 0 < components < n - 1:
        raise ValueError(f"{components} components need more than {components + 1} reference samples, got {n}")
    centred = x - level
    spread = centred[:, background].std(axis=0)
    spread = np.where(spread > 0, spread, 1.0)
    z = centred[:, background] / spread
    _, _, vt = np.linalg.svd(z, full_matrices=False)
    loadings = vt[:components].T
    scores = z @ loadings
    coefficients, *_ = np.linalg.lstsq(scores, centred, rcond=None)
    return Reference(level, strata, n, background, spread, loadings, coefficients)


def gene_counts(features: list[str], members: dict[str, list[str]]) -> tuple[np.ndarray, list[np.ndarray]]:
    """Which features are in a set, and per gene of it the column indices of its CpGs.

    ``members`` maps gene to its CpGs (a CpG may sit in several genes' promoters). The second value
    lists one index array per gene, the units a resample draws; genes with no CpG among ``features``
    are left out.
    """
    index = {f: j for j, f in enumerate(features)}
    genes = [np.array([index[c] for c in cpgs if c in index], dtype=int) for cpgs in members.values()]
    genes = [g for g in genes if g.size]
    in_set = np.zeros(len(features), dtype=bool)
    for g in genes:
        in_set[g] = True
    return in_set, genes


def resampled_gaps(
    d: np.ndarray, genes: list[np.ndarray], strata: np.ndarray, n_resamples: int, rng: np.random.Generator
) -> np.ndarray:
    """``(n_resamples, samples)`` values of ``P`` with the set's genes drawn with replacement.

    A drawn gene adds its CpGs to the set once per draw; a CpG left undrawn is not counted as
    background either (it stays a set CpG, absent from this resample), so the background is the same
    in every resample.
    """
    in_set = np.zeros(d.shape[1], dtype=bool)
    for g in genes:
        in_set[g] = True
    flat = np.concatenate(genes)
    owner = np.repeat(np.arange(len(genes)), [g.size for g in genes])
    out = np.empty((n_resamples, d.shape[0]))
    for r in range(n_resamples):
        draws = np.bincount(rng.integers(0, len(genes), len(genes)), minlength=len(genes)).astype(float)
        counts = np.bincount(flat, weights=draws[owner], minlength=d.shape[1])
        p, _ = _resample_weights(counts, in_set, strata)
        out[r] = d @ p
    return out


def _resample_weights(counts: np.ndarray, in_set: np.ndarray, strata: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """``contrast_weights`` with undrawn set features left out of both sides."""
    keep = (counts > 0) | ~in_set
    p = np.zeros(counts.size)
    b = np.zeros(counts.size)
    p[keep], b[keep] = contrast_weights(counts[keep], strata[keep])
    return p, b
