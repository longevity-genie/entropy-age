"""Bulk simulations: Gaussian stochastic variation added to a bounded ground state.

This is the model of Fig. 1 and Extended Data Fig. 1 of Meyer & Schumacher (2024).
A sample of *simulated age* ``a`` is the ground state plus a small jitter,
followed by ``a`` independent additions of noise drawn from ``N(0, sigma)``.

Whether an age signal exists depends only on the boundary:

* ``Bound.NONE`` - unbounded random walk; the mean never moves, so no clock can
  be trained (Extended Data Fig. 1c, Pearson ~0).
* ``Bound.CLIP`` - values are clipped to [0, 1] after every step; features near a
  bound can only move inwards, giving *regression to the mean* (Extended Data Fig. 1d).
* ``Bound.LOGIT`` - noise is added in logit space and mapped back with expit at
  the end (Fig. 1). The upstream ``random_epi_logit_everystep`` variant
  (logit -> add -> expit on every step) is mathematically the same model because
  ``logit(expit(x)) == x``; it is kept only as an alias.
* ``Bound.LOGIT_ONLY`` - like ``LOGIT`` but the result stays in logit space.

``sigma`` may be a scalar or a per-feature vector, which is the hook for
feature-specific (gene-specific) stochasticity.
"""

from __future__ import annotations

from collections.abc import Sequence
from enum import StrEnum

import numpy as np
from scipy.special import expit, logit

from entropy_age.simulate.schedule import iterate_age_steps, sort_by_age

_EPS = 1e-6


class Bound(StrEnum):
    NONE = "none"
    CLIP = "clip"
    LOGIT = "logit"
    LOGIT_EVERY_STEP = "logit_every_step"
    LOGIT_ONLY = "logit_only"


class NoiseDistribution(StrEnum):
    """Distribution of the per-step variation, always scaled to standard deviation ``sigma``."""

    NORMAL = "normal"
    LAPLACE = "laplace"
    UNIFORM = "uniform"


def _draw(rng: np.random.Generator, dist: NoiseDistribution, shape: tuple[int, ...]) -> np.ndarray:
    match dist:
        case NoiseDistribution.NORMAL:
            return rng.standard_normal(shape)
        case NoiseDistribution.LAPLACE:
            return rng.laplace(0.0, 1.0 / np.sqrt(2.0), shape)
        case NoiseDistribution.UNIFORM:
            half_width = np.sqrt(3.0)
            return rng.uniform(-half_width, half_width, shape)
    raise ValueError(f"unknown noise distribution {dist!r}")


def random_ground_state(n_features: int, rng: np.random.Generator, resolution: int = 1000) -> np.ndarray:
    """Uniform ground state on a grid of ``1/resolution`` (upstream: ``randint(0, 1000) / 1000``)."""
    return rng.integers(0, resolution, size=n_features) / resolution


def simulate_bulk(
    ground: np.ndarray,
    ages: Sequence[int] | np.ndarray,
    *,
    sigma: float | np.ndarray = 0.05,
    ground_sigma: float = 0.01,
    bound: Bound | str = Bound.CLIP,
    distribution: NoiseDistribution | str = NoiseDistribution.NORMAL,
    rng: np.random.Generator | None = None,
) -> np.ndarray:
    """Simulate one sample per entry of ``ages``; returns an ``(n_samples, n_features)`` array.

    Every sample is an independent trajectory, as in the paper. Instead of
    re-simulating every age from scratch (O(max_age^2) work, as upstream does),
    samples are processed in age order and step ``t`` is applied only to samples
    older than ``t``; the result has the same distribution.
    """
    rng = rng if rng is not None else np.random.default_rng()
    bound = Bound(bound)
    distribution = NoiseDistribution(distribution)
    ground = np.asarray(ground, dtype=float)
    sigma_arr = np.asarray(sigma, dtype=float)
    if sigma_arr.ndim not in (0, 1) or (sigma_arr.ndim == 1 and sigma_arr.shape != ground.shape):
        raise ValueError("sigma must be a scalar or have one value per feature")

    order, sorted_ages = sort_by_age(ages)
    n_samples, n_features = sorted_ages.size, ground.size
    in_logit = bound in (Bound.LOGIT, Bound.LOGIT_EVERY_STEP, Bound.LOGIT_ONLY)

    start = logit(np.clip(ground, _EPS, 1 - _EPS)) if in_logit else ground
    x = np.broadcast_to(start, (n_samples, n_features)).copy()
    x += ground_sigma * _draw(rng, NoiseDistribution.NORMAL, x.shape)
    if bound is Bound.CLIP:
        np.clip(x, 0.0, 1.0, out=x)

    for first in iterate_age_steps(sorted_ages):
        block = x[first:]
        block += sigma_arr * _draw(rng, distribution, block.shape)
        if bound is Bound.CLIP:
            np.clip(block, 0.0, 1.0, out=block)

    if bound in (Bound.LOGIT, Bound.LOGIT_EVERY_STEP):
        x = expit(x)

    out = np.empty_like(x)
    out[order] = x
    return out


def age_grid(max_age: int, per_age: int, *, start: int = 0) -> np.ndarray:
    """``per_age`` samples for every age in ``[start, max_age)`` (upstream: ``range(age_steps)``)."""
    return np.repeat(np.arange(start, max_age), per_age)
