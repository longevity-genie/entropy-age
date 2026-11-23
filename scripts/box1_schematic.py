"""Box 1 schematic: how drift at one CpG is measured and how the gap between a gene set and matched CpGs is formed.

Illustrative values only (simulated), drawn in the house style. ``uv run python scripts/box1_schematic.py`` writes
docs/figures/box1_schematic.png.
"""

from __future__ import annotations

from pathlib import Path

import numpy as np

from entropy_age import plots
from entropy_age.stochasticity import relative_entropy

ROOT = Path(__file__).resolve().parents[1]


def main() -> None:
    plots.apply_style()
    rng = np.random.default_rng(1)
    fig, axes = plots.plt.subplots(1, 3, figsize=(11, 3.4))

    # a: one CpG, young level 5 %, levels in people of three ages
    ax = axes[0]
    young = 0.05
    for i, (age, spread) in enumerate(((25, 0.01), (50, 0.03), (75, 0.06))):
        levels = np.clip(young + rng.normal(0.01 * i, spread, 60), 0, 1)
        ax.scatter(np.full(60, age) + rng.normal(0, 1.2, 60), levels * 100, s=6, color=plots.PALETTE[0], alpha=0.6)
    ax.axhline(young * 100, color=plots.MUTED, ls="--", lw=0.9)
    ax.text(78, young * 100 + 0.4, "young level", fontsize=7, color=plots.MUTED, ha="right")
    ax.set(xlabel="Age (years)", ylabel="Methylation at one CpG (%)", title="a  One CpG drifts with age")

    # b: relative entropy from the young level against age; its slope is the drift rate
    ax = axes[1]
    ages = np.linspace(20, 80, 80)
    for label, rate, colour in (("CpG in an essential promoter", 0.6, plots.PALETTE[0]), ("matched CpG", 1.0, plots.PALETTE[1])):
        level = np.clip(young + rate * 0.0012 * (ages - 20) + rng.normal(0, 0.004, ages.size), 0.001, 0.999)
        kl = relative_entropy(level, np.full(ages.size, young))
        ax.scatter(ages, kl, s=5, color=colour, alpha=0.5)
        coef = np.polyfit(ages, kl, 1)
        ax.plot(ages, np.polyval(coef, ages), color=colour, lw=1.6, label=label)
    ax.set(xlabel="Age (years)", ylabel="Relative entropy from young level (bits)",
           title="b  Drift rate = slope on age")
    ax.legend(fontsize=7, frameon=False, loc="upper left")

    # c: the gap: set drift rates minus matched drift rates, in SD over all CpGs
    ax = axes[2]
    matched = rng.normal(0.0, 1.0, 4000)
    gene_set = rng.normal(-0.4, 1.0, 600)
    bins = np.linspace(-4, 4, 50)
    ax.hist(matched, bins=bins, density=True, color=plots.PALETTE[1], alpha=0.45, label="matched CpGs")
    ax.hist(gene_set, bins=bins, density=True, color=plots.PALETTE[0], alpha=0.55, label="gene-set CpGs")
    ax.axvline(matched.mean(), color=plots.PALETTE[1], lw=1.2)
    ax.axvline(gene_set.mean(), color=plots.PALETTE[0], lw=1.2)
    ax.annotate("", xy=(gene_set.mean(), 0.43), xytext=(matched.mean(), 0.43),
                arrowprops={"arrowstyle": "->", "color": "black", "lw": 1})
    ax.text((gene_set.mean() + matched.mean()) / 2, 0.45, "gap < 0:\ndrifts less", ha="center", fontsize=7)
    ax.set(xlabel="Drift rate (SD over all CpGs)", ylabel="Density", title="c  Gap = set minus matched", ylim=(0, 0.52))
    ax.legend(fontsize=7, frameon=False, loc="upper right")

    fig.tight_layout()
    fig.savefig(ROOT / "docs" / "figures" / "box1_schematic.png")
    print("wrote docs/figures/box1_schematic.png")


if __name__ == "__main__":
    main()
