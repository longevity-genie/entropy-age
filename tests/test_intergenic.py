"""RM15: solo-WCGW tagging, directed loss and the age-tracking component."""

import numpy as np
import pytest

from entropy_age.data.annotation import solo_wcgw
from entropy_age.experiments.intergenic import components, directed_loss


def test_solo_wcgw_needs_at_flanks_and_no_cpg_within_35_bp():
    a = "A" * 40
    assert solo_wcgw(a + "[CG]" + a)
    assert not solo_wcgw("A" * 10 + "CG" + "A" * 28 + "[CG]" + a)  # a CpG 28 bp to the left
    assert solo_wcgw("A" * 2 + "CG" + "A" * 36 + "[CG]" + a)  # 36 bp away: outside the window
    assert not solo_wcgw(a[:-1] + "C[CG]" + a)  # C before the CpG
    assert not solo_wcgw(a + "[CG]")  # no context on the right


def _cohort(rng):
    n, p = 120, 600
    ages = rng.uniform(20, 90, n)
    solo = np.arange(p) < 150
    x = 0.7 + rng.normal(0, 0.02, (n, p))
    x[:, solo] -= 0.002 * (ages[:, None] - 20) * rng.uniform(0.5, 1.5, int(solo.sum()))  # coordinated loss
    return x, ages, solo


def test_directed_loss_is_found_at_solo_wcgw_and_not_elsewhere():
    x, ages, solo = _cohort(np.random.default_rng(0))
    rows = {r["cpgs"]: r for r in directed_loss(x, ages, None, solo)}
    assert rows["solo-WCGW"]["per_decade"] == pytest.approx(-0.02, abs=0.004) and rows["solo-WCGW"]["p"] < 1e-10
    assert abs(rows["other open sea"]["per_decade"]) < 0.002


def test_the_leading_component_tracks_age_and_loads_on_solo_wcgw():
    x, ages, solo = _cohort(np.random.default_rng(1))
    rows, _ = components(x, ages, solo, 3)
    first = rows[0]
    assert abs(first["r_age"]) > 0.8
    assert first["solo_weight"] > 2 * first["other_weight"] and first["p_loading"] < 1e-6
    assert first["solo_loading"] < 0  # solo-WCGW CpGs lose methylation as the component rises with age
