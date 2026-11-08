"""Small, dependency-free downloader with atomic writes."""

from __future__ import annotations

import shutil
import urllib.request
from pathlib import Path

USER_AGENT = "entropy-age (https://github.com/longevity-genie/entropy-age)"


def download(url: str, destination: Path, *, overwrite: bool = False, timeout: float = 120.0) -> Path:
    """Stream ``url`` to ``destination``; an existing file is kept unless ``overwrite``."""
    destination = Path(destination)
    if destination.exists() and not overwrite:
        return destination
    destination.parent.mkdir(parents=True, exist_ok=True)
    partial = destination.with_suffix(destination.suffix + ".part")
    request = urllib.request.Request(url, headers={"User-Agent": USER_AGENT})
    with urllib.request.urlopen(request, timeout=timeout) as response, partial.open("wb") as sink:
        shutil.copyfileobj(response, sink, length=1 << 20)
    partial.replace(destination)
    return destination
