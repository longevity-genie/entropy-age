"""Filesystem layout of the project.

``data/input``   - downloads, exactly as fetched (gitignored; recreate with ``entropy-age fetch``)
``data/interim`` - parsed parquet derived from the inputs (gitignored)
``data/output``  - result tables from the experiments (small CSV summaries are tracked)
``docs/figures`` - rendered figures (tracked)

All locations are :class:`pathlib.Path` objects derived from one root, so a
different checkout or a data directory on another disk is one environment
variable away (``ENTROPY_AGE_HOME`` for the root, ``ENTROPY_AGE_DATA`` for data).
"""

from __future__ import annotations

import os
from dataclasses import dataclass
from pathlib import Path

#: ``src/entropy_age/paths.py`` -> repository root.
_PACKAGE_ROOT = Path(__file__).resolve().parents[2]


@dataclass(frozen=True, slots=True)
class ProjectPaths:
    """Directory layout: downloaded inputs, interim parquet, outputs and figures."""

    root: Path

    @classmethod
    def default(cls) -> ProjectPaths:
        return cls(root=Path(os.environ.get("ENTROPY_AGE_HOME", _PACKAGE_ROOT)).resolve())

    @property
    def data(self) -> Path:
        override = os.environ.get("ENTROPY_AGE_DATA")
        return Path(override).resolve() if override else self.root / "data"

    @property
    def input(self) -> Path:
        return self.data / "input"

    @property
    def interim(self) -> Path:
        return self.data / "interim"

    @property
    def output(self) -> Path:
        return self.data / "output"

    @property
    def figures(self) -> Path:
        return self.root / "docs" / "figures"

    def ensure(self) -> ProjectPaths:
        for path in (self.input, self.interim, self.output, self.figures):
            path.mkdir(parents=True, exist_ok=True)
        return self
