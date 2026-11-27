# /// script
# requires-python = ">=3.11,<3.13"
# dependencies = ["pyaging==0.6.3", "twobitreader==3.1.7"]
# [tool.uv]
# extra-index-url = ["https://download.pytorch.org/whl/cpu"]
# index-strategy = "unsafe-best-match"
# ///
"""Write the coefficients of pyaging's mouse clocks as CSV for ``experiments/pseudobulk_clocks.py``.

``uv run --script scripts/export_mouse_clocks.py data/input/clocks/mouse data/input/annotation/mm10/mm10.2bit``
(pyaging 0.6.3, CPU torch; its own environment, so the project's stays without torch) writes, per clock,
``<clock>.csv`` with ``feature, weight, reference`` (the clock's input order; ``reference`` is pyaging's fill value
or NaN) and ``<clock>.json`` with the intercept and any pre- and post-processing constants:

- methylation, RRBS-type CpGs (``chrN:pos``, mm10): Thompson 2018 (582 CpGs, months), Meer 2018 (435, days),
  Petkovich 2017 (90, blood; a power-law conversion to days), Stubbs 2017 (17,992 inputs quantile-normalised to a
  gold standard and scaled, 329 weighted; a quadratic conversion to weeks);
- expression: tAge (Tyshkovskiy 2026, scaled-difference variant; 10,487 mouse Entrez genes, a cohort-relative age
  difference in months) and its gene mapping (``tage_gene_mapping.csv.gz``, source id to mouse Entrez).

The mm10 genome (``mm10.2bit``, UCSC; required, second argument) gives each methylation CpG ``strand``: ``+``
when the feature's 1-based position is the C of a CpG on the top strand, ``-`` when it is the G (the bottom strand's
C), empty otherwise. Coverage files list each strand's C apart, so a CpG's reads are those at ``pos`` and ``pos+1``
(``+``) or ``pos-1`` and ``pos`` (``-``).

Coefficients are those pyaging ships for each publication; the methods of each clock are in its paper.
"""

from __future__ import annotations

import json
import shutil
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import pyaging as pya
from pyaging.utils._hf import download_hf_file

CLOCKS = ("thompson", "meer", "petkovich", "stubbs", "tage")


def _array(value) -> np.ndarray | None:
    if value is None:
        return None
    return np.asarray(value.detach().cpu().numpy() if hasattr(value, "detach") else value, dtype=float)


def strands(features: list[str], genome: Path) -> list[str]:
    """``+``, ``-`` or ``""`` per ``chrN:pos`` feature (1-based), from the mm10 sequence."""
    from twobitreader import TwoBitFile

    seq = TwoBitFile(str(genome))
    out = []
    for f in features:
        chrom, pos = f.split(":")
        p = int(pos)
        tri = seq[chrom][p - 2 : p + 1].upper()  # bases p-1, p, p+1
        out.append("+" if tri[1:] == "CG" else "-" if tri[:2] == "CG" else "")
    return out


def main(out: Path, genome: Path | None = None) -> None:
    out.mkdir(parents=True, exist_ok=True)
    work = out / "_pyaging"
    for name in CLOCKS:
        model = pya.pred.load_clock(name, "cpu", str(work), verbose=False)
        state = model.base_model.state_dict()
        weight = state["linear.weight"].numpy().ravel().astype(float)
        bias = float(state["linear.bias"].numpy().ravel()[0])
        features = list(model.features)
        reference = _array(model.reference_values)
        meta = {
            "clock": name,
            "intercept": bias,
            "citation": model.metadata.get("citation"),
            "doi": model.metadata.get("doi"),
            "unit": model.metadata.get("unit"),
            "preprocess": model.preprocess_name,
            "postprocess": model.postprocess_name,
            "pyaging": pya.__version__ if hasattr(pya, "__version__") else "0.6.3",
        }
        if name == "stubbs":
            stds, index = model.preprocess_dependencies
            meta["used_index"] = [int(i) for i in index]
            frame = pd.DataFrame({"feature": features, "reference": reference, "sd": np.asarray(stds, float)})
            frame["weight"] = 0.0
            frame.loc[meta["used_index"], "weight"] = weight
        else:
            if weight.size != len(features):
                raise ValueError(f"{name}: {weight.size} weights for {len(features)} features")
            frame = pd.DataFrame(
                {"feature": features, "weight": weight, "reference": reference if reference is not None else np.nan}
            )
        if genome is not None and name != "tage":
            frame["strand"] = strands(features, genome)
            meta["strand_counts"] = frame["strand"].value_counts().to_dict()
        frame.to_csv(out / f"{name}.csv", index=False)
        (out / f"{name}.json").write_text(json.dumps(meta, indent=1))
        print(name, len(features), "features,", int((frame["weight"] != 0).sum()), "weighted", file=sys.stderr)
    mapping = download_hf_file("tage_gene_mapping.csv.gz", str(work), repo_id="pyaging/tage")
    (out / "tage_gene_mapping.csv.gz").unlink(missing_ok=True)
    shutil.copyfile(mapping, out / "tage_gene_mapping.csv.gz")
    shutil.rmtree(work, ignore_errors=True)


if __name__ == "__main__":
    if len(sys.argv) != 3:
        sys.exit("usage: export_mouse_clocks.py OUT_DIR mm10.2bit (the scorer needs each CpG's strand)")
    main(Path(sys.argv[1]), Path(sys.argv[2]))
