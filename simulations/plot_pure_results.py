"""Visualise pure λ_XC simulation results."""

from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


def plot_heatmap(df: pd.DataFrame, path: Path | None) -> plt.Figure:
    pivot = df.pivot(index="sample_size", columns="fst_target", values="lambda_ratio")
    fst_vals = pivot.columns.to_numpy()
    sample_vals = pivot.index.to_numpy()

    vmin = float(np.min(pivot.values))
    vmax = float(np.max(pivot.values))
    center = 1.0
    spread = max(abs(vmin - center), abs(vmax - center))
    if spread == 0:
        spread = 0.01

    fig, ax = plt.subplots(figsize=(1.8 * len(fst_vals), 1.05 * len(sample_vals)))
    im = ax.imshow(
        pivot.values,
        aspect="auto",
        cmap="coolwarm",
        vmin=center - spread,
        vmax=center + spread,
    )
    cbar = fig.colorbar(im, ax=ax)
    cbar.set_label("λ_emp / λ_theory")

    ax.set_xticks(range(len(fst_vals)))
    ax.set_xticklabels([f"{x:.3f}" for x in fst_vals], rotation=45, ha="right")
    ax.set_yticks(range(len(sample_vals)))
    ax.set_yticklabels([f"{int(x)}" for x in sample_vals])
    ax.set_xlabel("F_ST target")
    ax.set_ylabel("Sample size (diploids)")
    ax.set_title("Ratio: λ_emp / λ_theory")

    for i, sn in enumerate(sample_vals):
        for j, fst in enumerate(fst_vals):
            value = pivot.iloc[i, j]
            ax.text(
                j,
                i,
                f"{value:.3f}",
                ha="center",
                va="center",
                color="black",
                fontsize=9,
            )

    if path:
        path.parent.mkdir(parents=True, exist_ok=True)
        fig.savefig(path, dpi=200, bbox_inches="tight")
    return fig


def plot_lambda_panels(df: pd.DataFrame, path: Path | None) -> plt.Figure:
    sample_sizes = sorted(df["sample_size"].unique())

    n_cols = 2
    n_rows = int(np.ceil(len(sample_sizes) / n_cols))
    fig, axes = plt.subplots(
        n_rows,
        n_cols,
        figsize=(5.2 * n_cols, 3.6 * n_rows),
        sharex=True,
        sharey=False,
    )
    axes = np.array(axes).reshape(n_rows, n_cols)

    for idx, sample in enumerate(sample_sizes):
        ax = axes.flat[idx]
        sub = df[df["sample_size"] == sample].sort_values("fst_target")
        fst = sub["fst_target"].to_numpy()
        lam_emp = sub["lambda_emp"].to_numpy()
        lam_theory = sub["lambda_theory"].to_numpy()
        ci_lo = sub["lambda_emp_ci_lo"].to_numpy()
        ci_hi = sub["lambda_emp_ci_hi"].to_numpy()
        err = np.vstack((lam_emp - ci_lo, ci_hi - lam_emp))

        ax.plot(
            fst,
            lam_theory,
            linestyle="--",
            color="black",
            linewidth=1.1,
            label="λ_theory",
        )
        ax.errorbar(
            fst,
            lam_emp,
            yerr=err,
            fmt="o",
            markersize=5,
            capsize=4,
            elinewidth=1.0,
            linewidth=1.0,
            color="#1f77b4",
            label="λ_emp",
        )

        for idx_pt, (x, y) in enumerate(zip(fst, lam_emp)):
            theory_val = lam_theory[idx_pt]
            ax.text(
                x,
                y + 0.035,
                f"{y:.3f}\n({theory_val:.3f})",
                ha="center",
                va="bottom",
                fontsize=7,
                color="#1f77b4",
            )

        ax.set_title(f"Sample size N = {int(sample)}")
        ax.set_xlabel("F_ST target")
        ax.set_ylabel("λ")
        ax.grid(alpha=0.3)
        if idx == 0:
            ax.legend(fontsize=8)

    for ax in axes.flat[len(sample_sizes) :]:
        ax.axis("off")

    fig.suptitle("Empirical vs theoretical λ", y=0.98)
    fig.tight_layout()

    if path:
        path.parent.mkdir(parents=True, exist_ok=True)
        fig.savefig(path, dpi=200, bbox_inches="tight")
    return fig


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Plot pure λ_XC simulation results.")
    parser.add_argument(
        "--summary",
        type=Path,
        required=True,
        help="Path to pure_simulation_summary.csv.",
    )
    parser.add_argument(
        "--heatmap",
        type=Path,
        help="Optional output path for the λ difference heatmap PNG.",
    )
    parser.add_argument(
        "--lines",
        type=Path,
        help="Optional output path for the λ vs F_ST comparison PNG.",
    )
    parser.add_argument(
        "--show",
        action="store_true",
        help="Display figures interactively (may require GUI support).",
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    df = pd.read_csv(args.summary)
    df["sample_size"] = df["sample_size"].astype(int)
    df["fst_target"] = df["fst_target"].astype(float)
    if "lambda_diff" not in df.columns:
        df["lambda_diff"] = df["lambda_emp"] - df["lambda_theory"]

    heatmap_path = args.heatmap.expanduser().resolve() if args.heatmap else None
    lines_path = args.lines.expanduser().resolve() if args.lines else None

    heatmap_fig = plot_heatmap(df, heatmap_path)
    lines_fig = plot_lambda_panels(df, lines_path)

    if args.show:
        plt.show()
    else:
        plt.close(heatmap_fig)
        plt.close(lines_fig)


if __name__ == "__main__":
    main()
