import numpy as np
import polars as pl

from entropy_age.data.reprocessed import collapse_probes
from entropy_age.experiments.preprocessing_qc import agreement, batch_shares, confounding, eta_squared, pcs


def test_eta_squared_and_batch_share():
    chips = np.array(["a"] * 5 + ["b"] * 5)
    assert eta_squared(np.array([0.0] * 5 + [1.0] * 5), chips) == 1.0
    rng = np.random.default_rng(0)
    x = rng.normal(size=(10, 50)) * 0.01
    x[chips == "b"] += 1.0  # a chip shift dominates
    scores, var = pcs(x, 3)
    assert batch_shares(scores, var, {"chip": chips})["chip"] > 0.9


def test_confounding_flags_age_set_by_chip():
    chips = np.array(["a"] * 6 + ["b"] * 6)
    design = pl.DataFrame({"group": ["young"] * 6 + ["old"] * 6, "sex": ["F", "M"] * 6})
    r = {d["variable"]: d["r2_on_chip"] for d in confounding(design, chips)}
    assert r["group"] == 1.0 and r["sex"] == 0.0


def test_collapse_probes_averages_replicates():
    f = pl.DataFrame({"probe": ["cg1_TC11", "cg1_TC12", "cg2_BC21"], "s": [0.2, 0.4, 0.9]})
    c = collapse_probes(f)
    assert dict(zip(c["probe"], c["s"])) == {"cg1": 0.30000000000000004, "cg2": 0.9}


def test_agreement_matches_by_position():
    rebuilt = pl.DataFrame(
        {
            "probe": [f"cg{i}" for i in range(5)],
            "GSM1": [0.1, 0.2, 0.3, 0.8, 0.9],
            "GSM2": [0.2, 0.2, 0.4, 0.7, 0.9],
            "GSM3": [0.1, 0.3, 0.3, 0.9, 0.8],
        }
    )
    qc = pl.DataFrame(
        {"sample": ["GSM1", "GSM2", "GSM3"], "chip": ["c"] * 3, "position": ["R01C01", "R02C01", "R03C01"]}
    )
    dep = rebuilt.rename({"probe": "cpg", "GSM1": "c_R01C01", "GSM2": "c_R02C01", "GSM3": "c_R03C01"})
    r = agreement(rebuilt, qc, dep)
    assert r["matched_samples"] == 3 and abs(r["sample_r_median"] - 1) < 1e-12
