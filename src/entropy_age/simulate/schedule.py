"""Shared bookkeeping for simulating many independent samples of different ages at once."""

from __future__ import annotations

from collections.abc import Iterator, Sequence

import numpy as np


def sort_by_age(ages: Sequence[int] | np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """Return ``(order, sorted_ages)``; results are written back with ``out[order] = x``."""
    ages_arr = np.asarray(ages)
    if ages_arr.ndim != 1:
        raise ValueError("ages must be one-dimensional")
    if ages_arr.size and (ages_arr.min() < 0 or not np.all(np.equal(np.mod(ages_arr, 1), 0))):
        raise ValueError("ages are numbers of stochastic steps: non-negative integers")
    ages_int = ages_arr.astype(np.int64)
    order = np.argsort(ages_int, kind="stable")
    return order, ages_int[order]


def iterate_age_steps(sorted_ages: np.ndarray) -> Iterator[int]:
    """Yield, for every step ``t = 1..max_age``, the first row index whose age is ``>= t``.

    Rows ``[first:]`` of an age-sorted array are exactly the samples that still
    receive stochastic variation at step ``t``.
    """
    if sorted_ages.size == 0:
        return
    for step in range(1, int(sorted_ages[-1]) + 1):
        yield int(np.searchsorted(sorted_ages, step, side="left"))
