"""Command-line helpers for lambda_xc."""

from __future__ import annotations

import argparse
import json
from pathlib import Path
from typing import List

import pandas as pd

from .pipeline import PanelConfig, SimulationSettings, run_panel_analysis
from .plots import qq_panels


def parse_args(argv: List[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Run cross-chromosome inflation analyses for HapMap3 panels."
    )
    parser.add_argument(
        "--data-root",
        type=Path,
        required=True,
        help="Directory containing genotype files for the requested panels.",
    )
    parser.add_argument(
        "--genotype-format",
        choices=["auto", "hapmap", "plink", "vcf"],
        default="auto",
        help=(
            "File format for the genotype panels. "
            "Use 'auto' to infer from the data directory."
        ),
    )
    parser.add_argument(
        "--panels",
        nargs="+",
        default=["CHB", "CEU", "CEU+TSI", "CEU+YRI"],
        help="Panels to analyse (comma-separated populations when mixing).",
    )
    parser.add_argument("--nsnp-per-chr", type=int, default=2_000)
    parser.add_argument("--prune-kb", type=int, default=100)
    parser.add_argument("--pairs", type=int, default=1_000_000)
    parser.add_argument("--permutations", type=int, default=80)
    parser.add_argument("--bootstraps", type=int, default=100)
    parser.add_argument("--seed", type=int, default=511)
    parser.add_argument(
        "--qq-path",
        type=Path,
        help="Optional path to save QQ plots (PNG).",
    )
    parser.add_argument(
        "--results-json",
        type=Path,
        help="Optional path to store the raw results as JSON.",
    )
    parser.add_argument(
        "--results-csv",
        type=Path,
        help="Optional path to store a tabular summary as CSV.",
    )
    return parser.parse_args(argv)


def parse_panel_names(names: List[str]) -> List[PanelConfig]:
    configs: List[PanelConfig] = []
    for label in names:
        pops = tuple(part.strip() for part in label.split("+"))
        configs.append(PanelConfig(label=label, populations=pops))
    return configs


def format_summary(results):
    rows = []
    for res in results:
        rows.append(
            {
                "panel": res.label,
                "N": res.sample_size,
                "lambda_theory": res.lambda_theory,
                "lambda_emp": res.lambda_emp,
                "lambda_emp_sd": res.lambda_emp_sd,
                "lambda_emp_ci_lo": res.lambda_emp_ci[0],
                "lambda_emp_ci_hi": res.lambda_emp_ci[1],
                "baseline": res.baseline,
                "baseline_sd": res.baseline_sd,
                "baseline_ci_lo": res.baseline_ci[0],
                "baseline_ci_hi": res.baseline_ci[1],
                "lambda_corr": res.lambda_corr,
                "lambda_corr_sd": res.lambda_corr_sd,
                "lambda_corr_ci_lo": res.lambda_corr_ci[0],
                "lambda_corr_ci_hi": res.lambda_corr_ci[1],
            }
        )
    df = pd.DataFrame(rows)
    return df


def main(argv: List[str] | None = None) -> None:
    args = parse_args(argv)

    panels = parse_panel_names(args.panels)
    settings = SimulationSettings(
        nsnp_per_chr=args.nsnp_per_chr,
        prune_kb=args.prune_kb,
        n_pairs=args.pairs,
        n_perm=args.permutations,
        n_boot=args.bootstraps,
        seed=args.seed,
    )

    results = run_panel_analysis(
        args.data_root,
        panels,
        settings,
        genotype_format=args.genotype_format,
    )
    df = format_summary(results)
    print(
        df.to_string(
            index=False,
            float_format=lambda x: f"{x:.3f}",
        )
    )

    if args.qq_path:
        qq_panels(results, output=args.qq_path)
        print(f"QQ plot saved to {args.qq_path.resolve()}")

    if args.results_json:
        payload = [
            {
                "label": res.label,
                "sample_size": res.sample_size,
                "lambda_theory": res.lambda_theory,
                "lambda_emp": res.lambda_emp,
                "lambda_emp_sd": res.lambda_emp_sd,
                "lambda_emp_ci": res.lambda_emp_ci,
                "baseline": res.baseline,
                "baseline_sd": res.baseline_sd,
                "baseline_ci": res.baseline_ci,
                "lambda_corr": res.lambda_corr,
                "lambda_corr_sd": res.lambda_corr_sd,
                "lambda_corr_ci": res.lambda_corr_ci,
            }
            for res in results
        ]
        args.results_json.parent.mkdir(parents=True, exist_ok=True)
        args.results_json.write_text(json.dumps(payload, indent=2))
        print(f"JSON results written to {args.results_json.resolve()}")

    if args.results_csv:
        args.results_csv.parent.mkdir(parents=True, exist_ok=True)
        df.to_csv(args.results_csv, index=False)
        print(f"CSV summary written to {args.results_csv.resolve()}")


if __name__ == "__main__":
    main()
