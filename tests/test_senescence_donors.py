import numpy as np

from entropy_age.experiments.senescence_donors import arm_effects, squared_share

RS, IRR = "Replicative sensescent", "Irradiation-induced prematurely senescent"


def test_arm_effects_recover_a_treatment_shift_with_donor_and_chip_held():
    rng = np.random.default_rng(0)
    donors, chips, arms = [], [], []
    for d in range(6):
        for k, a in enumerate(["Young", RS, IRR, RS]):
            donors.append(f"d{d}")
            chips.append(f"c{(d + k) % 4}")
            arms.append(a)
    donor, chip, arm = np.array(donors), np.array(chips), np.array(arms)
    donor_shift = np.array([0.05 * int(d[1:]) for d in donor])
    chip_shift = np.array([0.03 * int(c[1:]) for c in chip])
    x = np.column_stack(
        [
            0.3 + donor_shift + chip_shift + 0.2 * (arm == RS) + rng.normal(0, 0.005, arm.size),
            0.5 + donor_shift + chip_shift + rng.normal(0, 0.005, arm.size),
        ]
    )
    effects = arm_effects(x, arm, donor, chip)
    assert set(effects) == {RS, IRR}
    b, se = effects[RS]
    assert abs(b[0] - 0.2) < 0.02 and abs(b[1]) < 0.02
    assert (se > 0).all() and abs(effects[IRR][0]).max() < 0.02


def test_squared_share_is_minus_half_when_the_set_moves_by_a_factor_of_root_two_less():
    in_set = np.array([True, True, False, False])
    strata = np.array(["a", "a", "a", "a"])
    b = np.array([1.0, 1.0, np.sqrt(2.0), np.sqrt(2.0)])
    out = squared_share(b, np.zeros(4), np.ones(4), in_set, strata)
    assert abs(out["F2"] + 0.5) < 1e-9 and abs(out["B2"] - 2.0) < 1e-9
