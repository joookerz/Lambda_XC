"""Plotting utilities."""

from __future__ import annotations

from pathlib import Path
from typing import Iterable

import matplotlib.pyplot as plt
import numpy as np
import scipy.stats as ss

from .analysis import PanelResult


def qq_panels(
    results: Iterable[PanelResult],
    *,
    dpi: int = 150,
    output: Path | None = None,
) -> plt.Figure:
    results = list(results)
    fig, axes = plt.subplots(1, len(results), figsize=(4.5 * len(results), 4), dpi=dpi)
    if len(results) == 1:
        axes = [axes]

    for ax, res in zip(axes, results):
        chi2 = np.sort(res.chi2_values)
        expected = ss.chi2.ppf((np.arange(1, chi2.size + 1) - 0.5) / chi2.size, 1)
        ax.scatter(expected, chi2, s=4, alpha=0.7)
        maxv = max(expected.max(), chi2.max())
        ax.plot([0, maxv], [0, maxv], lw=1, color="black")
        ax.set_title(res.label)
        ax.set_xlabel("Expected χ²₁")
        ax.set_ylabel("Observed χ²")
        ax.grid(alpha=0.3)

    fig.tight_layout()
    if output:
        path = Path(output).expanduser().resolve()
        path.parent.mkdir(parents=True, exist_ok=True)
        fig.savefig(path, bbox_inches="tight")
    return fig
