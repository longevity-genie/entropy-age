# /// script
# requires-python = ">=3.10,<3.13"
# dependencies = [
#   "biolearn @ git+https://github.com/bio-learn/biolearn@0d714f5a0c0ad5743b6315f3edab1f86112d394f",
#   "torch", "seaborn", "matplotlib", "requests", "pyyaml", "appdirs", "scikit-learn", "pyarrow",
# ]
# ///
"""Score clocks with scripts/clocks_biolearn.py with the coverage rule switched off: a labelled check only.

``CLOCKS_ONLY=GrimAge uv run --script scripts/clocks_forced.py IN OUT [META]``. Same scorer, imputation and
measured-share column as clocks_biolearn.py, but a clock is scored however little of its weight is measured, so a
clock that fails coverage on one track can still be set beside the other track to see how far the missing CpGs move
it. Output names must say so (``..._below_coverage_clocks.csv``); such scores are never a clock panel's. In the
sidecar panel a clock that would have failed the coverage rule is ``forced below coverage``, never ``kept``;
``uv run --script scripts/clocks_forced.py --relabel PANEL`` applies that to an existing sidecar.
clock_preprocessing sets these scores beside the other track as its ``forced`` diagnostic.
"""

import sys
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))

FORCED = "forced below coverage"
MIN_ARRAY_SHARE = 0.9  # clocks_biolearn.MIN_ARRAY_SHARE, the rule the forced scores bypass


def relabel(panel: str) -> None:
    """Mark the clocks of a forced sidecar panel that fail the ordinary coverage rule."""
    table = pd.read_csv(panel)
    if "arrays_covered" in table.columns:
        below = (table["status"] == "kept") & (table["arrays_covered"] < MIN_ARRAY_SHARE)
        table.loc[below, "status"] = FORCED
    table.to_csv(panel, index=False)


if __name__ == "__main__":
    if sys.argv[1:2] == ["--relabel"]:
        relabel(sys.argv[2])
        raise SystemExit(0)
    if "below_coverage" not in sys.argv[2]:
        raise SystemExit("clocks_forced: the output name must contain 'below_coverage'")
    import clocks_biolearn

    assert clocks_biolearn.MIN_ARRAY_SHARE == MIN_ARRAY_SHARE
    clocks_biolearn.MIN_ARRAY_SHARE = 0.0
    clocks_biolearn.main(*sys.argv[1:4])
    relabel(sys.argv[2].removesuffix(".tmp").replace(".csv", "_panel.csv"))
