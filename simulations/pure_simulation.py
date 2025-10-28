"""Pure simulation of cross-chromosome inflation without baseline correction.

This module generates synthetic two-population mixtures with configurable
sample sizes and F_ST targets, computes the empirical λ_XC (median χ² scaled
by the χ²_1 median), and compares it against the theoretical prediction

    λ_theory = 1 + 1.096 · N · \bar{F}^2 ,

where \bar{F} is the realised mean per-locus F estimate from the simulated
allele frequencies.  For each (sample_size, F_ST) combination one synthetic
dataset is generated and bootstrap resampling of the χ² values provides
standard errors and percentile confidence intervals.
"""

from __future__ import annotations

import argparse
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable, Iterator, Sequence, Tuple

import numpy as np
import pandas as pd

CHI2_MEDIAN = 0.4559364


@dataclass
class SimulationConfig:
    sample_sizes: Tuple[int, ...]
    fst_values: Tuple[float, ...]
    n_chrom: int = 22
    snps_per_chr: int = 800
    n_pairs: int = 80_000
    n_boot: int = 400
    min_maf: float = 0.05
    seed: int = 2025
    store_bootstrap: bool = False


def _random_pairs(chr_labels: np.ndarray, n_pairs: int, rng: np.random.Generator) -> Iterator[Tuple[int, int]]:
    """Yield indices for random cross-chromosome SNP pairs."""
    bins: dict[int, list[int]] = {}
    for idx, chrom in enumerate(chr_labels):
        bins.setdefault(int(chrom), []).append(int(idx))
    chromosomes = np.array(list(bins.keys()), dtype=int)
    for _ in range(n_pairs):
        c1, c2 = rng.choice(chromosomes, 2, replace=False)
        yield rng.choice(bins[int(c1)]), rng.choice(bins[int(c2)])


def _chi2_pairs(
    standardized: np.ndarray,
    valid: np.ndarray,
    pairs: Iterable[Tuple[int, int]],
    n_samples: int,
) -> np.ndarray:
    """Return χ² = n · r² values for the provided SNP pairs."""
    stats: list[float] = []
    denom = float(n_samples - 1)
    for i, j in pairs:
        if not (valid[i] and valid[j]):
            continue
        r = float(np.dot(standardized[i], standardized[j]) / denom)
        stats.append((r * r) * n_samples)
    return np.asarray(stats, dtype=float)


def _bootstrap_lambda(
    chi2_values: np.ndarray,
    lambda_theory: float,
    n_boot: int,
    rng: np.random.Generator,
) -> Tuple[np.ndarray, np.ndarray]:
    """Generate bootstrap replicates of λ_emp and λ_emp / λ_theory."""
    samples = []
    ratios = []
    size = chi2_values.size
    for _ in range(n_boot):
        resampled = rng.choice(chi2_values, size=size, replace=True)
        lam = float(np.median(resampled) / CHI2_MEDIAN)
        samples.append(lam)
        ratios.append(lam / lambda_theory if lambda_theory != 0 else np.nan)
    return np.asarray(samples, dtype=float), np.asarray(ratios, dtype=float)


def simulate_once(
    sample_size: int,
    fst: float,
    cfg: SimulationConfig,
    rng: np.random.Generator,
) -> Tuple[dict[str, float], np.ndarray, np.ndarray]:
    """Simulate a single dataset and return summary stats plus bootstrap draws."""
    if sample_size % 2 != 0:
        raise ValueError("Sample size must be even to split into two populations.")
    half = sample_size // 2

    n_snps = cfg.n_chrom * cfg.snps_per_chr
    chr_labels = np.repeat(np.arange(1, cfg.n_chrom + 1), cfg.snps_per_chr)

    base_freq = rng.uniform(cfg.min_maf, 1.0 - cfg.min_maf, size=n_snps)
    delta = np.sqrt(2.0 * fst * base_freq * (1.0 - base_freq))
    p1 = np.clip(base_freq + 0.5 * delta, 1e-4, 1 - 1e-4)
    p2 = np.clip(base_freq - 0.5 * delta, 1e-4, 1 - 1e-4)

    geno1 = rng.binomial(2, p1[:, None], size=(n_snps, half))
    geno2 = rng.binomial(2, p2[:, None], size=(n_snps, half))
    geno = np.concatenate([geno1, geno2], axis=1).astype(float)

    means = geno.mean(axis=1)
    centered = geno - means[:, None]
    std = centered.std(axis=1, ddof=1)
    valid = std > 0
    standardized = np.zeros_like(centered)
    standardized[valid] = centered[valid] / std[valid, None]

    pairs = list(_random_pairs(chr_labels, cfg.n_pairs, rng))
    chi2_values = _chi2_pairs(standardized, valid, pairs, sample_size)
    if chi2_values.size == 0:
        raise RuntimeError("No χ² statistics were produced. Consider increasing SNP density or pair count.")

    lambda_emp = float(np.median(chi2_values) / CHI2_MEDIAN)

    mix_freq = 0.5 * (p1 + p2)
    denom = 2.0 * mix_freq * (1.0 - mix_freq)
    with np.errstate(divide="ignore", invalid="ignore"):
        F_vals = np.divide(
            (p1 - p2) ** 2,
            denom,
            out=np.zeros_like(mix_freq),
            where=denom > 0,
        )
    F_effective = float(F_vals.mean())
    lambda_theory = 1.0 + 1.096 * sample_size * (F_effective ** 2)

    boot_emp, boot_ratio = _bootstrap_lambda(chi2_values, lambda_theory, cfg.n_boot, rng)

    lambda_ratio = lambda_emp / lambda_theory if lambda_theory != 0 else np.nan

    summary = {
        "sample_size": float(sample_size),
        "fst_target": float(fst),
        "fst_effective": F_effective,
        "snps_total": float(n_snps),
        "snps_retained": float(int(valid.sum())),
        "chi2_pairs": float(chi2_values.size),
        "lambda_emp": lambda_emp,
        "lambda_theory": lambda_theory,
        "lambda_diff": lambda_emp - lambda_theory,
        "lambda_ratio": lambda_ratio,
        "lambda_emp_sd": float(boot_emp.std(ddof=1)),
        "lambda_emp_ci_lo": float(np.percentile(boot_emp, 2.5)),
        "lambda_emp_ci_hi": float(np.percentile(boot_emp, 97.5)),
        "lambda_ratio_sd": float(boot_ratio.std(ddof=1)),
        "lambda_ratio_ci_lo": float(np.percentile(boot_ratio, 2.5)),
        "lambda_ratio_ci_hi": float(np.percentile(boot_ratio, 97.5)),
    }
    return summary, boot_emp, boot_ratio


def run_simulation(cfg: SimulationConfig, output_dir: Path) -> None:
    summaries = []
    bootstrap_rows = []

    base_ss = np.random.SeedSequence(cfg.seed)
    total = len(cfg.sample_sizes) * len(cfg.fst_values)
    scenario_seeds = base_ss.spawn(total)
    index = 0

    for sample_size in cfg.sample_sizes:
        for fst in cfg.fst_values:
            rng = np.random.default_rng(scenario_seeds[index])
            index += 1
            summary, boot_emp, boot_ratio = simulate_once(sample_size, fst, cfg, rng)
            summaries.append(summary)
            if cfg.store_bootstrap:
                for idx, (val, ratio) in enumerate(zip(boot_emp, boot_ratio)):
                    bootstrap_rows.append(
                        {
                            "sample_size": sample_size,
                            "fst_target": fst,
                            "bootstrap_index": idx,
                            "lambda_emp_boot": val,
                            "lambda_ratio_boot": ratio,
                        }
                    )

    summary_df = pd.DataFrame(summaries)
    summary_path = output_dir / "pure_simulation_summary.csv"
    summary_df.sort_values(["sample_size", "fst_target"]).to_csv(summary_path, index=False)
    print(f"Summary written to {summary_path}")

    if cfg.store_bootstrap and bootstrap_rows:
        boot_df = pd.DataFrame(bootstrap_rows)
        boot_path = output_dir / "pure_simulation_bootstrap.csv"
        boot_df.to_csv(boot_path, index=False)
        print(f"Bootstrap replicates written to {boot_path}")


def parse_args(argv: Sequence[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Pure λ_XC simulation without baseline correction.")
    parser.add_argument(
        "--sample-grid",
        type=int,
        nargs="+",
        required=True,
        help="List of total diploid sample sizes (even numbers).",
    )
    parser.add_argument(
        "--fst-grid",
        type=float,
        nargs="+",
        required=True,
        help="List of target F_ST values.",
    )
    parser.add_argument(
        "--snps-per-chr",
        type=int,
        default=800,
        help="Number of SNPs simulated per chromosome (default: 800).",
    )
    parser.add_argument(
        "--pairs",
        type=int,
        default=80_000,
        help="Number of random cross-chromosome SNP pairs per scenario (default: 80000).",
    )
    parser.add_argument(
        "--bootstraps",
        type=int,
        default=400,
        help="Bootstrap replicates for λ_emp standard errors (default: 400).",
    )
    parser.add_argument(
        "--seed",
        type=int,
        default=2025,
        help="Seed for the RNG (default: 2025).",
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=Path("simulations/output"),
        help="Directory where CSV outputs will be written.",
    )
    parser.add_argument(
        "--store-bootstrap",
        action="store_true",
        help="Persist per-bootstrap λ samples to CSV in addition to the summary.",
    )
    parser.add_argument(
        "--min-maf",
        type=float,
        default=0.05,
        help="Lower bound on baseline allele frequencies (default: 0.05).",
    )
    return parser.parse_args(argv)


def main(argv: Sequence[str] | None = None) -> None:
    args = parse_args(argv)
    cfg = SimulationConfig(
        sample_sizes=tuple(int(x) for x in args.sample_grid),
        fst_values=tuple(float(x) for x in args.fst_grid),
        snps_per_chr=int(args.snps_per_chr),
        n_pairs=int(args.pairs),
        n_boot=int(args.bootstraps),
        seed=int(args.seed),
        min_maf=float(args.min_maf),
        store_bootstrap=bool(args.store_bootstrap),
    )

    output_dir = args.output_dir.expanduser().resolve()
    output_dir.mkdir(parents=True, exist_ok=True)

    run_simulation(cfg, output_dir)


if __name__ == "__main__":
    main()
