"""Extension: cell turnover and selection acting on stochastic epimutations.

Hypothesis this module is built to explore: the stochastic variation that
aging clocks read is not equally free everywhere. Sites (genes) whose state
matters for cell fitness can be protected from entropy by selection among
cells, and how much protection a site gets depends on the cell type, because
selection needs cell turnover to act.

Model, per tissue sample (one individual):

* ``n_cells`` cells, each a binary vector over sites (1 = methylated / "on").
* Every step, each site of each cell epimutates by the two-state Markov chain of
  :mod:`entropy_age.simulate.methylation` (``em``, ``ed`` per site).
* A fraction ``turnover`` of cells dies at random each step and is replaced by
  daughters of cells drawn with probability proportional to fitness
  ``exp(-strength * sum_j weight_j * deviation_ij)``. With ``reference="founder"``
  a deviation is any site whose state differs from the founder state of that
  cell lineage (selection against change); with ``reference="target"`` it is a
  difference from an explicit functional state.
* Daughters can pick up replication errors: one extra epimutation pass with the
  error probabilities multiplied by ``division_error``.

``turnover = 0`` is a post-mitotic cell type (selection cannot act);
``strength = 0`` is neutral clonal drift. Weighted sites are "constrained"
genes, zero-weight sites are neutral passengers that can still hitchhike with
selected clones.
"""

from __future__ import annotations

from collections.abc import Sequence
from dataclasses import dataclass, field
from enum import StrEnum

import numpy as np

from entropy_age.simulate.methylation import MaintenanceRates
from entropy_age.simulate.schedule import iterate_age_steps, sort_by_age


class Reference(StrEnum):
    FOUNDER = "founder"
    TARGET = "target"


@dataclass(frozen=True, slots=True)
class SelectionModel:
    """Selection among cells. ``weights`` has one non-negative value per site."""

    weights: np.ndarray
    strength: float = 0.0
    turnover: float = 0.0
    division_error: float = 0.0
    reference: Reference = Reference.FOUNDER
    target: np.ndarray | None = None

    def __post_init__(self) -> None:
        weights = np.asarray(self.weights, dtype=float)
        if np.any(weights < 0):
            raise ValueError("selection weights must be non-negative")
        if not 0.0 <= self.turnover <= 1.0:
            raise ValueError("turnover is a per-step fraction in [0, 1]")
        object.__setattr__(self, "weights", weights)
        object.__setattr__(self, "reference", Reference(self.reference))
        if self.reference is Reference.TARGET:
            if self.target is None:
                raise ValueError("reference='target' needs a target state per site")
            object.__setattr__(self, "target", np.asarray(self.target, dtype=bool))

    @classmethod
    def neutral(cls, n_sites: int, turnover: float = 0.0, division_error: float = 0.0) -> SelectionModel:
        return cls(np.zeros(n_sites), 0.0, turnover, division_error)

    @property
    def constrained(self) -> np.ndarray:
        return np.flatnonzero(self.weights > 0)


@dataclass(frozen=True, slots=True)
class CellPopulation:
    """A cell type: its epimutation rates, its selection regime and the cells sampled per tissue."""

    name: str
    rates: MaintenanceRates
    selection: SelectionModel
    n_cells: int = 200
    meta: dict[str, float | str] = field(default_factory=dict)


def _founder_cells(level: np.ndarray, n_cells: int, rng: np.random.Generator) -> np.ndarray:
    """Exactly ``floor(n_cells * level_j)`` methylated cells per site, randomly arranged."""
    counts = np.floor(n_cells * level).astype(np.int64)
    ranks = np.argsort(rng.random((n_cells, level.size)), axis=0).argsort(axis=0)
    return ranks < counts[None, :]


def _epimutate(cells: np.ndarray, rates: MaintenanceRates, rng: np.random.Generator, scale: float = 1.0) -> None:
    loss = np.minimum(1.0, (1.0 - rates.em) * scale)
    gain = np.minimum(1.0, rates.ed * scale)
    u = rng.random(cells.shape)
    cells ^= np.where(cells, u < loss, u < gain)


def _choose_parents(fitness: np.ndarray, rng: np.random.Generator) -> np.ndarray:
    """One fitness-proportional parent per (sample, cell) slot, vectorised over samples."""
    n_samples, n_cells = fitness.shape
    cdf = np.cumsum(fitness, axis=1)
    cdf /= cdf[:, -1:]
    offsets = np.arange(n_samples)[:, None]
    flat = (cdf + offsets).ravel()
    draws = (rng.random((n_samples, n_cells)) + offsets).ravel()
    idx = np.searchsorted(flat, draws, side="right").reshape(n_samples, n_cells)
    return np.minimum(idx - offsets * n_cells, n_cells - 1)


def _turnover_step(
    cells: np.ndarray,
    founders: np.ndarray,
    population: CellPopulation,
    rng: np.random.Generator,
) -> None:
    sel = population.selection
    n_samples, n_cells, _ = cells.shape
    dying = rng.random((n_samples, n_cells)) < sel.turnover
    if not dying.any():
        return
    constrained = sel.constrained
    if sel.strength > 0 and constrained.size:
        reference = founders if sel.reference is Reference.FOUNDER else sel.target[constrained]
        deviation = cells[:, :, constrained] != reference
        fitness = np.exp(-sel.strength * (deviation @ sel.weights[constrained]))
    else:
        fitness = np.ones((n_samples, n_cells))
    parents = _choose_parents(fitness, rng)
    rows = np.arange(n_samples)[:, None]
    newborn_cells = cells[rows, parents]
    cells[dying] = newborn_cells[dying]
    if sel.reference is Reference.FOUNDER and constrained.size:
        founders[dying] = founders[rows, parents][dying]
    if sel.division_error > 0:
        daughters = cells[dying]
        _epimutate(daughters, population.rates, rng, scale=sel.division_error)
        cells[dying] = daughters


def simulate_with_selection(
    ground: np.ndarray,
    ages: Sequence[int] | np.ndarray,
    population: CellPopulation,
    *,
    jitter: float = 0.01,
    rng: np.random.Generator | None = None,
) -> np.ndarray:
    """Bulk levels ``(n_samples, n_sites)`` of independent tissues of the given ages."""
    rng = rng if rng is not None else np.random.default_rng()
    ground = np.asarray(ground, dtype=float)
    n_sites, n_cells = ground.size, population.n_cells
    if len(population.rates) != n_sites or population.selection.weights.size != n_sites:
        raise ValueError("rates and selection weights need one value per site")

    order, sorted_ages = sort_by_age(ages)
    n_samples = sorted_ages.size
    levels = np.clip(ground + jitter * rng.standard_normal((n_samples, n_sites)), 0.0, 1.0)
    cells = np.stack([_founder_cells(level, n_cells, rng) for level in levels])
    constrained = population.selection.constrained
    founders = cells[:, :, constrained].copy()

    for first in iterate_age_steps(sorted_ages):
        # basic slices are views: both helpers update ``cells`` and ``founders`` in place
        _epimutate(cells[first:], population.rates, rng)
        if population.selection.turnover > 0:
            _turnover_step(cells[first:], founders[first:], population, rng)

    bulk = cells.mean(axis=1)
    out = np.empty_like(bulk)
    out[order] = bulk
    return out
