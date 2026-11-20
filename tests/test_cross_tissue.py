import gzip

import numpy as np

from entropy_age.experiments.cross_tissue import dispersion_gaps, read_people


def test_read_people_takes_donor_from_title(tmp_path):
    path = tmp_path / "m.txt.gz"
    with gzip.open(path, "wt") as handle:
        handle.write('!Sample_title\t"Adult_8_PBMC"\t"Adult_8_Buccal"\t"Adult_9_PBMC"\n')
        handle.write('!Sample_geo_accession\t"GSM1"\t"GSM2"\t"GSM3"\n')
        handle.write('!Sample_characteristics_ch1\t"tissue: PBMC"\t"tissue: Buccal"\t"tissue: PBMC"\n')
        handle.write('!Sample_characteristics_ch1\t"Sex: Female"\t"Sex: Female"\t"Sex: Male"\n')
        handle.write('!Sample_characteristics_ch1\t"age: 30.5"\t"age: 30.5"\t"age: NA"\n')
    people = read_people(path)
    assert people["sample"].to_list() == ["GSM1", "GSM2"]
    assert people["donor"].to_list() == ["Adult_8", "Adult_8"]
    assert people["tissue"].to_list() == ["PBMC", "Buccal"]
    assert people["age"].to_list() == [30.5, 30.5]


def test_dispersion_gap_is_negative_when_the_set_varies_less_between_people():
    rng = np.random.default_rng(0)
    n, m = 60, 400
    in_set = np.arange(m) < 100
    spread = np.where(in_set, 0.01, 0.03)
    x = 0.5 + rng.normal(size=(n, m)) * spread
    island = np.zeros(m, dtype=bool)
    out = dispersion_gaps(x, np.empty((n, 0)), [f"c{j}" for j in range(m)], island, {"s": (in_set, [])})
    gap, sd = out[("dispersion", "s")]
    assert gap < 0 and sd > 0
