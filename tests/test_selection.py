import numpy as np
import pytest

from entropy_age.simulate.methylation import MaintenanceRates
from entropy_age.simulate.selection import (
    CellPopulation,
    SelectionModel,
    _choose_parents,
    simulate_with_selection,
)


def _setup(n=120, constrained=40):
    rng = np.random.default_rng(11)
    ground = np.where(rng.random(n) < 0.5, 0.03, 0.97)
    weights = np.zeros(n)
    weights[:constrained] = 1.0
    return ground, weights, MaintenanceRates.universal(n, 0.99)


def _drift(model: SelectionModel, seed: int) -> tuple[float, float]:
    ground, weights, rates = _setup()
    ages = np.full(12, 40)
    x = simulate_with_selection(
        ground, ages, CellPopulation("t", rates, model, n_cells=100), rng=np.random.default_rng(seed)
    )
    lost = np.abs(x - ground)
    return lost[:, weights > 0].mean(), lost[:, weights == 0].mean()


def test_post_mitotic_cells_get_no_protection():
    _, weights, _ = _setup()
    constrained, neutral = _drift(SelectionModel(weights, strength=1.0, turnover=0.0), 0)
    assert abs(constrained - neutral) < 0.03


def test_turnover_with_selection_protects_constrained_sites():
    _, weights, _ = _setup()
    constrained, neutral = _drift(SelectionModel(weights, strength=1.0, turnover=0.3), 1)
    assert constrained < 0.6 * neutral


def test_neutral_turnover_does_not_protect():
    _, weights, _ = _setup()
    constrained, neutral = _drift(SelectionModel(weights, strength=0.0, turnover=0.3), 2)
    assert abs(constrained - neutral) < 0.03


def test_parent_choice_follows_fitness():
    rng = np.random.default_rng(0)
    fitness = np.tile(np.array([1.0, 0.0, 3.0, 0.0]), (2000, 1))
    parents = _choose_parents(fitness, rng)
    counts = np.bincount(parents.ravel(), minlength=4) / parents.size
    assert counts[1] == 0 and counts[3] == 0
    assert counts[2] == pytest.approx(0.75, abs=0.01)


def test_validation():
    with pytest.raises(ValueError):
        SelectionModel(np.array([-1.0]))
    with pytest.raises(ValueError):
        SelectionModel(np.zeros(3), turnover=1.5)
    with pytest.raises(ValueError):
        SelectionModel(np.zeros(3), reference="target")
