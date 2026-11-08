"""Matplotlib figures in one consistent, print-friendly style."""

from __future__ import annotations

from collections.abc import Mapping, Sequence
from pathlib import Path

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import polars as pl  # noqa: E402

from entropy_age.clock import ClockMetrics  # noqa: E402

INK = "#1f2933"
MUTED = "#7b8794"
POINT = "#52606d"
ACCENT = "#c2410c"
PALETTE = ("#2563eb", "#c2410c", "#0f766e", "#7c3aed", "#b45309", "#be185d", "#4d7c0f")


def apply_style() -> None:
    plt.rcParams.update(
        {
            "figure.dpi": 110,
            "savefig.dpi": 200,
            "savefig.bbox": "tight",
            "font.size": 9,
            "axes.titlesize": 10,
            "axes.titleweight": "bold",
            "axes.labelsize": 9,
            "axes.edgecolor": MUTED,
            "axes.labelcolor": INK,
            "axes.spines.top": False,
            "axes.spines.right": False,
            "xtick.color": MUTED,
            "ytick.color": MUTED,
            "legend.frameon": False,
            "legend.fontsize": 8,
        }
    )


def panels(n: int, *, ncols: int | None = None, size: float = 3.0) -> tuple[plt.Figure, np.ndarray]:
    ncols = ncols or n
    nrows = int(np.ceil(n / ncols))
    fig, axes = plt.subplots(nrows, ncols, figsize=(size * ncols, size * nrows), squeeze=False)
    for ax in axes.ravel()[n:]:
        ax.set_visible(False)
    return fig, axes.ravel()


def save(fig: plt.Figure, path: Path) -> Path:
    path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(path)
    plt.close(fig)
    return path


def label_panel(ax: plt.Axes, letter: str) -> None:
    ax.text(-0.16, 1.06, letter, transform=ax.transAxes, fontsize=12, fontweight="bold", color=INK, va="bottom")


def prediction_scatter(
    ax: plt.Axes,
    true: np.ndarray,
    pred: np.ndarray,
    *,
    metrics: ClockMetrics | None = None,
    title: str | None = None,
    xlabel: str = "Simulated age",
    ylabel: str = "Predicted age",
    identity: bool = True,
    color: str = POINT,
    size: float = 5.0,
) -> None:
    true = np.asarray(true, dtype=float)
    pred = np.asarray(pred, dtype=float)
    ax.scatter(true, pred, s=size, color=color, alpha=0.55, linewidths=0)
    if np.ptp(true) > 0:
        slope, intercept = np.polyfit(true, pred, 1)
        xs = np.linspace(true.min(), true.max(), 50)
        ax.plot(xs, slope * xs + intercept, color=ACCENT, lw=1.4)
    if identity:
        lo = min(true.min(), pred.min())
        hi = max(true.max(), pred.max())
        ax.plot([lo, hi], [lo, hi], ls=":", color=MUTED, lw=1)
    if metrics is not None:
        summary = (
            f"r = {metrics.pearson_r:.2f}\nR² = {metrics.r2:.2f}\nn = {metrics.n}"
            if np.isfinite(metrics.pearson_r)
            else f"no age signal:\nthe fit shrinks to a constant\nn = {metrics.n}"
        )
        ax.text(
            0.04,
            0.96,
            summary,
            transform=ax.transAxes,
            va="top",
            fontsize=8,
            color=INK,
        )
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    if title:
        ax.set_title(title)


def sweep_box(
    ax: plt.Axes,
    frame: pl.DataFrame,
    *,
    x: str,
    y: str,
    order: Sequence | None = None,
    labels: Sequence[str] | None = None,
    xlabel: str | None = None,
    ylabel: str | None = None,
    title: str | None = None,
    ylim: tuple[float, float] | None = None,
    seed: int = 0,
) -> None:
    rng = np.random.default_rng(seed)
    order = list(order) if order is not None else frame[x].unique(maintain_order=True).to_list()
    data = [frame.filter(pl.col(x) == level)[y].to_numpy() for level in order]
    ax.boxplot(
        data,
        positions=range(len(order)),
        widths=0.55,
        showfliers=False,
        medianprops={"color": ACCENT, "lw": 1.5},
        boxprops={"color": MUTED},
        whiskerprops={"color": MUTED},
        capprops={"color": MUTED},
    )
    for i, values in enumerate(data):
        ax.scatter(i + rng.uniform(-0.15, 0.15, values.size), values, s=8, color=INK, alpha=0.7, linewidths=0)
    ax.set_xticks(range(len(order)), labels if labels is not None else [str(v) for v in order], rotation=45, ha="right")
    ax.set_xlabel(x if xlabel is None else xlabel)
    ax.set_ylabel(ylabel or y)
    if ylim:
        ax.set_ylim(*ylim)
    if title:
        ax.set_title(title)


def series_scatter(
    ax: plt.Axes,
    x: np.ndarray,
    series: Mapping[str, np.ndarray],
    *,
    xlabel: str = "Simulated age",
    ylabel: str = "Predicted age",
    title: str | None = None,
    size: float = 4.0,
) -> None:
    for (name, y), color in zip(series.items(), PALETTE, strict=False):
        ax.scatter(x, y, s=size, color=color, alpha=0.5, linewidths=0, label=name)
    ax.legend(markerscale=2.5, loc="upper left")
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    if title:
        ax.set_title(title)


def mean_band(ax: plt.Axes, x: np.ndarray, y: np.ndarray, *, label: str, color: str) -> None:
    """Mean of ``y`` per unique ``x`` with an interquartile band."""
    levels = np.unique(x)
    mean = np.array([y[x == v].mean() for v in levels])
    lo = np.array([np.quantile(y[x == v], 0.25) for v in levels])
    hi = np.array([np.quantile(y[x == v], 0.75) for v in levels])
    ax.plot(levels, mean, color=color, lw=1.6, label=label)
    ax.fill_between(levels, lo, hi, color=color, alpha=0.15, linewidth=0)
