"""Plot heatmaps of detection rates from power_summary.csv.

Generates two heatmaps per summary:
  - λ_XC z-test detection rate (lambda_detect_rate)
  - PCA (TW) detection rate (pca_detect_rate)

Rows: sample size N; Columns: effect (F_ST for mixture; tau for cline).
"""

from __future__ import annotations

import argparse
from pathlib import Path
from typing import Optional

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


def _heatmap(ax, data: np.ndarray, row_labels, col_labels, cmap="viridis", vmin=0.0, vmax=1.0, annotate=True, fmt=".2f"):
    im = ax.imshow(data, aspect="auto", cmap=cmap, vmin=vmin, vmax=vmax)
    ax.set_xticks(range(len(col_labels)))
    ax.set_yticks(range(len(row_labels)))
    ax.set_xticklabels(col_labels)
    ax.set_yticklabels([str(int(x)) for x in row_labels])
    ax.set_xlabel("effect")
    ax.set_ylabel("N")
    if annotate:
        for i in range(data.shape[0]):
            for j in range(data.shape[1]):
                ax.text(j, i, format(data[i, j], fmt), ha="center", va="center", color="white" if data[i, j] > 0.5 else "black", fontsize=8)
    return im


def main() -> None:
    ap = argparse.ArgumentParser(description="Plot detection-rate heatmaps from power_summary.csv")
    ap.add_argument("--summary", type=Path, required=True, help="Path to power_summary.csv")
    ap.add_argument("--out-dir", type=Path, help="Directory to write PNGs (defaults to summary's dir)")
    ap.add_argument("--title-prefix", type=str, default="", help="Optional title prefix for figures")
    args = ap.parse_args()

    df = pd.read_csv(args.summary)
    if "model" in df.columns and df["model"].nunique() == 1:
        model = df["model"].iloc[0]
    else:
        model = ""

    # Ensure numeric
    df["N"] = df["N"].astype(int)
    df["effect"] = df["effect"].astype(float)

    # Pivot matrices
    piv_lam = df.pivot(index="N", columns="effect", values="lambda_detect_rate").sort_index().sort_index(axis=1)
    piv_pca = df.pivot(index="N", columns="effect", values="pca_detect_rate").sort_index().sort_index(axis=1)

    out_dir = (args.out_dir or args.summary.parent).expanduser().resolve()
    out_dir.mkdir(parents=True, exist_ok=True)

    # λ_XC heatmap
    fig, ax = plt.subplots(figsize=(1.8 * len(piv_lam.columns), 0.8 * len(piv_lam.index) + 1.5), dpi=160)
    im = _heatmap(ax, piv_lam.values, piv_lam.index.values, [f"{c:.3f}" for c in piv_lam.columns])
    cbar = fig.colorbar(im, ax=ax)
    cbar.set_label("λ_XC detect rate")
    title = f"{args.title_prefix} λ_XC z-test detect rate"
    if model:
        title = f"{model} - {title}"
    ax.set_title(title)
    fig.tight_layout()
    fig.savefig(out_dir / "lambda_detect_rate_heatmap.png", bbox_inches="tight")
    plt.close(fig)

    # PCA heatmap
    fig, ax = plt.subplots(figsize=(1.8 * len(piv_pca.columns), 0.8 * len(piv_pca.index) + 1.5), dpi=160)
    im = _heatmap(ax, piv_pca.values, piv_pca.index.values, [f"{c:.3f}" for c in piv_pca.columns], cmap="magma")
    cbar = fig.colorbar(im, ax=ax)
    cbar.set_label("PCA(TW) detect rate")
    title = f"{args.title_prefix} PCA(TW) detect rate"
    if model:
        title = f"{model} - {title}"
    ax.set_title(title)
    fig.tight_layout()
    fig.savefig(out_dir / "pca_detect_rate_heatmap.png", bbox_inches="tight")
    plt.close(fig)

    print(f"Saved: {(out_dir / 'lambda_detect_rate_heatmap.png').as_posix()} and {(out_dir / 'pca_detect_rate_heatmap.png').as_posix()}")


if __name__ == "__main__":
    main()

