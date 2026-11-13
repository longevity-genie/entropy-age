"""RM12: CR's blood effect with CpGs matched on openness too."""

import numpy as np
import pyBigWig

from entropy_age.experiments.cr_chromatin import access_class, openness, with_access


def test_access_class_puts_closed_cpgs_first():
    signal = np.array([0.0, 0.0, 0.0, 1.0, 2.0, 3.0, 4.0, 5.0])
    c = access_class(signal, 4)
    assert c[0] == c[1] == c[2] == 0 and c[-1] == 3 and np.all(np.diff(c) >= 0)


def test_with_access_splits_strata():
    s = with_access(np.array(["True3", "True3"]), np.array([0, 2]))
    assert s.tolist() == ["True3a0", "True3a2"]


def test_openness_averages_tracks(tmp_path):
    paths = []
    for k, level in enumerate((2.0, 4.0)):
        path = tmp_path / f"t{k}.bw"
        bw = pyBigWig.open(str(path), "w")
        bw.addHeader([("chr1", 10_000)])
        bw.addEntries(["chr1"], [1000], ends=[2000], values=[level])
        bw.close()
        paths.append(path)
    v = openness(np.array(["1", "1"]), np.array([1500, 8000]), paths, half=100)
    assert abs(v[0] - 3.0) < 1e-9 and v[1] == 0.0
