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
it. Output names must say so (``..._below_coverage_clocks.csv``); such scores are never a clock panel's.
"""

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import clocks_biolearn  # noqa: E402

if __name__ == "__main__":
    if "below_coverage" not in sys.argv[2]:
        raise SystemExit("clocks_forced: the output name must contain 'below_coverage'")
    clocks_biolearn.MIN_ARRAY_SHARE = 0.0
    clocks_biolearn.main(*sys.argv[1:4])
