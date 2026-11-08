"""Single-cell DNA methylation as a two-state Markov chain per CpG site.

Each site in each cell is methylated (1) or not (0). Per time step a methylated
site stays methylated with probability ``em`` (maintenance efficiency) and an
unmethylated site becomes methylated with probability ``ed`` (de novo
efficiency). This is the model of Fig. 3-5 of Meyer & Schumacher (2024),
following Pfeifer et al. (1990).

Closed forms used throughout the package (``r = em - ed``):

* equilibrium (fraction methylated):  ``z = ed / (1 - em + ed)``
* expected bulk level after ``t`` steps:  ``z + (m0 - z) * r**t``
* relaxation rate per step:  ``lambda = -ln(r)``  (~ ``(1 - em) + ed`` for small rates)

``lambda`` is the per-site *entropy rate*: how fast a site forgets its ground
state. Gene- or cell-type-specific maintenance and selection act on it.

Simulation detail: the cells at one site are exchangeable and independent, so
the methylated-cell count after one step is exactly
``Binomial(k, em) + Binomial(n_cells - k, ed)``. Simulating counts gives the
same distribution as flipping every cell (the upstream approach) at
``O(sites)`` instead of ``O(sites * cells)`` cost. ``method="cells"`` keeps the
explicit per-cell simulation for validation.
"""

from __future__ import annotations

from collections.abc import Sequence
from dataclasses import dataclass
from enum import StrEnum

import numpy as np

from entropy_age.simulate.schedule import iterate_age_steps, sort_by_age

_EQ_EPS = 1e-4


@dataclass(frozen=True, slots=True)
class MaintenanceRates:
    """Per-site maintenance (``em``) and de novo (``ed``) methylation efficiencies."""

    em: np.ndarray
    ed: np.ndarray

    def __post_init__(self) -> None:
        em = np.atleast_1d(np.asarray(self.em, dtype=float))
        ed = np.atleast_1d(np.asarray(self.ed, dtype=float))
        if em.shape != ed.shape:
            raise ValueError("em and ed must have the same shape")
        if np.any((em < 0) | (em > 1) | (ed < 0) | (ed > 1)):
            raise ValueError("efficiencies are probabilities in [0, 1]")
        object.__setattr__(self, "em", em)
        object.__setattr__(self, "ed", ed)

    def __len__(self) -> int:
        return self.em.size

    # ------------------------------------------------------------------ analytics
    @property
    def equilibrium(self) -> np.ndarray:
        denom = 1.0 - self.em + self.ed
        with np.errstate(invalid="ignore", divide="ignore"):
            return np.where(denom > 0, self.ed / np.where(denom > 0, denom, 1.0), np.nan)

    @property
    def persistence(self) -> np.ndarray:
        """``r = em - ed``: the fraction of the distance to equilibrium kept per step."""
        return self.em - self.ed

    @property
    def relaxation_rate(self) -> np.ndarray:
        """``lambda = -ln(em - ed)`` per step; 0 means perfect maintenance."""
        with np.errstate(divide="ignore"):
            return -np.log(np.clip(self.persistence, 1e-300, None))

    def expected(self, ground: np.ndarray, steps: np.ndarray | int) -> np.ndarray:
        """Expected bulk methylation after ``steps``; broadcasts ``steps`` over a leading axis."""
        steps_arr = np.asarray(steps, dtype=float)
        z = self.equilibrium
        decay = np.power(self.persistence, steps_arr[..., None]) if steps_arr.ndim else self.persistence**steps_arr
        return z + (np.asarray(ground, dtype=float) - z) * decay

    def subset(self, index: np.ndarray) -> MaintenanceRates:
        return MaintenanceRates(self.em[index], self.ed[index])

    # ------------------------------------------------------------- constructors
    @classmethod
    def universal(cls, n_sites: int, em: float, ed: float | None = None) -> MaintenanceRates:
        """Same efficiency for every site; by default ``ed = 1 - em`` (equilibrium 0.5), as upstream."""
        ed_value = 1.0 - em if ed is None else ed
        return cls(np.full(n_sites, em), np.full(n_sites, ed_value))

    @classmethod
    def random(
        cls,
        n_sites: int,
        rng: np.random.Generator,
        *,
        em_limit: float = 0.97,
        ed_limit: float = 0.05,
    ) -> MaintenanceRates:
        """Independent ``em ~ U(em_limit, 1)`` and ``ed ~ U(0, ed_limit)`` (Extended Data Fig. 6)."""
        return cls(rng.uniform(em_limit, 1.0, n_sites), rng.uniform(0.0, ed_limit, n_sites))

    @classmethod
    def from_equilibrium(
        cls,
        equilibrium: np.ndarray,
        rng: np.random.Generator,
        *,
        em_limit: float = 0.97,
        ed_limit: float = 0.05,
    ) -> MaintenanceRates:
        """Empirical rates from a site-specific equilibrium (upstream ``get_noise_Em_all_new``).

        The equilibrium is usually the methylation of the oldest sample. ``em`` is
        fixed to ``em_limit`` and ``ed`` solved from ``z = ed / (1 - em + ed)``;
        where that ``ed`` exceeds ``ed_limit`` a random ``ed`` in ``(0, ed_limit)``
        is drawn and ``em`` solved instead, which keeps ``em >= em_limit``.
        Equilibria are clipped away from exactly 0 and 1 (upstream divides by zero there).
        """
        z = np.clip(np.asarray(equilibrium, dtype=float), _EQ_EPS, 1 - _EQ_EPS)
        ed = (1.0 - em_limit) * z / (1.0 - z)
        too_high = ed > ed_limit
        ed = np.where(too_high, rng.uniform(1e-3, ed_limit, z.size), ed)
        em = np.where(too_high, 1.0 - ed * (1.0 - z) / z, em_limit)
        return cls(np.clip(em, 0.0, 1.0), np.clip(ed, 0.0, 1.0))

    @classmethod
    def from_rate_and_equilibrium(cls, rate: np.ndarray, equilibrium: np.ndarray) -> MaintenanceRates:
        """Invert ``lambda`` and ``z``: ``ed = z (1 - e^-lambda)``, ``1 - em = (1 - z)(1 - e^-lambda)``."""
        loss = 1.0 - np.exp(-np.asarray(rate, dtype=float))
        z = np.clip(np.asarray(equilibrium, dtype=float), 0.0, 1.0)
        return cls(1.0 - (1.0 - z) * loss, z * loss)


class CellMethod(StrEnum):
    BINOMIAL = "binomial"
    CELLS = "cells"


def initial_counts(
    ground: np.ndarray,
    n_samples: int,
    n_cells: int,
    *,
    jitter: float,
    rng: np.random.Generator,
) -> np.ndarray:
    """Methylated-cell counts per sample and site: ``floor(n_cells * clip(ground + N(0, jitter)))``."""
    ground = np.asarray(ground, dtype=float)
    level = np.broadcast_to(ground, (n_samples, ground.size)).copy()
    if jitter > 0:
        level += jitter * rng.standard_normal(level.shape)
    np.clip(level, 0.0, 1.0, out=level)
    return np.floor(n_cells * level).astype(np.int64)


def simulate_methylation(
    ground: np.ndarray,
    ages: Sequence[int] | np.ndarray,
    rates: MaintenanceRates,
    *,
    n_cells: int = 1000,
    jitter: float = 0.01,
    method: CellMethod | str = CellMethod.BINOMIAL,
    rng: np.random.Generator | None = None,
) -> np.ndarray:
    """Bulk methylation (fraction of methylated cells) of independent samples of the given ages.

    Returns an ``(n_samples, n_sites)`` array in the order of ``ages``.
    """
    rng = rng if rng is not None else np.random.default_rng()
    ground = np.asarray(ground, dtype=float)
    if len(rates) != ground.size:
        raise ValueError("one pair of maintenance rates per site is required")
    order, sorted_ages = sort_by_age(ages)
    counts = initial_counts(ground, sorted_ages.size, n_cells, jitter=jitter, rng=rng)

    match CellMethod(method):
        case CellMethod.BINOMIAL:
            for first in iterate_age_steps(sorted_ages):
                k = counts[first:]
                counts[first:] = rng.binomial(k, rates.em) + rng.binomial(n_cells - k, rates.ed)
            levels = counts / n_cells
        case CellMethod.CELLS:
            levels = np.empty(counts.shape, dtype=float)
            for row, (k0, age) in enumerate(zip(counts, sorted_ages, strict=True)):
                cells = np.arange(n_cells)[None, :] < k0[:, None]  # (sites, cells)
                for _ in range(int(age)):
                    u = rng.random(cells.shape)
                    cells ^= np.where(cells, u >= rates.em[:, None], u < rates.ed[:, None])
                levels[row] = cells.mean(axis=1)

    out = np.empty_like(levels)
    out[order] = levels
    return out


def binomial_variance(ground: np.ndarray, steps: int, rates: MaintenanceRates, n_cells: int) -> np.ndarray:
    """Exact variance of the bulk level after ``steps`` given ``floor(n_cells * ground)`` starting cells."""
    k0 = np.floor(n_cells * np.asarray(ground, dtype=float))
    z, r = rates.equilibrium, rates.persistence
    p_from_methylated = z + (1 - z) * r**steps
    p_from_unmethylated = z - z * r**steps
    var = k0 * p_from_methylated * (1 - p_from_methylated)
    var += (n_cells - k0) * p_from_unmethylated * (1 - p_from_unmethylated)
    return var / n_cells**2
