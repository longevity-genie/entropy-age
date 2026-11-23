"""Each main figure with its legend (the caption from manuscript/paper.md) printed underneath, for Drive.

``uv run python scripts/figures_with_legends.py`` writes docs/figures/captioned/figure_<N>.png (and figure_S22.png
for the forest plot, captioned from scripts/forest_series.py's description). Markdown emphasis is stripped.
"""

from __future__ import annotations

import re
import textwrap
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.image as mpimg  # noqa: E402
import matplotlib.pyplot as plt  # noqa: E402

ROOT = Path(__file__).resolve().parents[1]
PAPER = ROOT / "manuscript" / "paper.md"
OUT = ROOT / "docs" / "figures" / "captioned"
EXTRA = {
    "S22": (
        ROOT / "docs" / "figures" / "forest_series.png",
        "Supplementary Figure S22. Every series on one measure. The housekeeping gap (housekeeping promoter CpGs minus "
        "matched CpGs, SD of the relative-entropy slope) in 42 independent cohort-tissue series, 23 human (blue) and "
        "19 mouse (orange), with 95 % intervals that resample donors or animals; chip intercepts held where "
        "identifiable. Large dots: interval excludes zero. Platform and processing are given in docs/forest_series.csv.",
    )
}


def plain(text: str) -> str:
    text = re.sub(r"\*\*(.+?)\*\*", r"\1", text, flags=re.S)
    text = re.sub(r"\*(.+?)\*", r"\1", text, flags=re.S)
    return re.sub(r"\s+", " ", text).strip()


def captions() -> dict[str, tuple[Path, str]]:
    text = PAPER.read_text()
    out = {}
    for m in re.finditer(r"!\[Figure (\d+)\]\(([^)]+)\)\s*\n\s*\n?(\*\*Figure \1\..*?)(?:\n\s*\n)", text, flags=re.S):
        number, rel, caption = m.group(1), m.group(2), m.group(3)
        out[number] = ((PAPER.parent / rel).resolve(), plain(caption))
    return out


def render(image: Path, caption: str, dest: Path) -> None:
    img = mpimg.imread(image)
    h, w = img.shape[:2]
    width_in = 10.0
    lines = textwrap.wrap(caption, 150)
    text_in = 0.17 * len(lines) + 0.3
    img_in = width_in * h / w
    fig = plt.figure(figsize=(width_in, img_in + text_in), dpi=200)
    ax = fig.add_axes((0, text_in / (img_in + text_in), 1, img_in / (img_in + text_in)))
    ax.imshow(img)
    ax.axis("off")
    fig.text(0.02, (text_in - 0.15) / (img_in + text_in), "\n".join(lines), fontsize=7.2, va="top", ha="left",
             family="DejaVu Sans")
    fig.savefig(dest, facecolor="white")
    plt.close(fig)


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    items = {**captions(), **EXTRA}
    for number, (image, caption) in items.items():
        render(image, caption, OUT / f"figure_{number}.png")
    print(f"{len(items)} captioned figures: {', '.join(sorted(items, key=lambda k: (len(k), k)))}")


if __name__ == "__main__":
    main()
