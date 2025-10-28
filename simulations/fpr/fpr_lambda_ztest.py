"""Null (no-structure) FPR for λ_XC using a low-cost z-test with permutations.

Procedure per replicate (single panmictic population):
  1) Compute λ_emp from median of χ² = n r² across random cross-chromosome pairs.
     Also bootstrap χ² to estimate sd_emp (K resamples, default 100).
  2) Do B column permutations (default 20), compute λ_perm^(b) with the same
     pair indices, and record their mean μ_perm and sample SD s_perm.
  3) Form z = (λ_emp − μ_perm) / sqrt(sd_emp² + s_perm² / B).
     One-sided p = 1 − Φ(z). Declare detection if z > 1.645 (α=0.05).

Outputs per-replicate details and per-N summary (FPR under null).
"""

from __future__ import annotations

import argparse
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable, Iterator, Sequence, Tuple, List

import numpy as np
import pandas as pd
import scipy.stats as ss

CHI2_MEDIAN = 0.4559364


def _random_pairs(chr_labels: np.ndarray, n_pairs: int, rng: np.random.Generator) -> List[Tuple[int, int]]:
    bins: dict[int, List[int]] = {}
    for idx, chrom in enumerate(chr_labels):
        bins.setdefault(int(chrom), []).append(int(idx))
    chromosomes = np.array(list(bins.keys()), dtype=int)
    pairs: List[Tuple[int, int]] = []
    for _ in range(n_pairs):
        c1, c2 = rng.choice(chromosomes, 2, replace=False)
        pairs.append((rng.choice(bins[int(c1)]), rng.choice(bins[int(c2)])))
    return pairs


def _chi2_pairs(G: np.ndarray, pairs: Iterable[Tuple[int, int]], n_samples: int) -> np.ndarray:
    stats: List[float] = []
    denom = float(n_samples - 1)
    for i, j in pairs:
        r = float(np.dot(G[i], G[j]) / denom)
        stats.append((r * r) * n_samples)
    return np.asarray(stats, dtype=float)


def _standardize(geno: np.ndarray) -> np.ndarray:
    means = geno.mean(axis=1)
    centered = geno - means[:, None]
    std = centered.std(axis=1, ddof=1)
    std[std == 0] = 1.0
    return centered / std[:, None]


def simulate_null(n: int, snps_per_chr: int, rng: np.random.Generator, min_maf: float = 0.05) -> tuple[np.ndarray, np.ndarray]:
    n_chrom = 22
    m = n_chrom * snps_per_chr
    chr_labels = np.repeat(np.arange(1, n_chrom + 1), snps_per_chr)
    p = rng.uniform(min_maf, 1 - min_maf, size=m)
    geno = rng.binomial(2, p[:, None], size=(m, n)).astype(float)
    return geno, chr_labels


@dataclass
class Config:
    N_grid: tuple[int, ...]
    snps_per_chr: int
    pairs: int
    permutations: int
    bootstraps: int
    replicates: int
    lambda_alpha_z: float
    seed: int
    output_dir: Path


def run(cfg: Config) -> None:
    ssq = np.random.SeedSequence(cfg.seed)
    seeds = ssq.spawn(len(cfg.N_grid) * cfg.replicates)
    rows = []
    idx = 0
    for N in cfg.N_grid:
        for rep in range(cfg.replicates):
            rng = np.random.default_rng(seeds[idx]); idx += 1
            geno, chr_labels = simulate_null(N, cfg.snps_per_chr, rng)
            Z = _standardize(geno)
            pairs = _random_pairs(chr_labels, cfg.pairs, rng)

            # λ_emp and bootstrap SD
            chi2 = _chi2_pairs(Z, pairs, N)
            lam_emp = float(np.median(chi2) / CHI2_MEDIAN)
            # Bootstrap over chi2 samples
            if chi2.size < 2:
                sd_emp = np.nan
            else:
                boots = []
                for _ in range(cfg.bootstraps):
                    res = rng.choice(chi2, size=chi2.size, replace=True)
                    boots.append(float(np.median(res) / CHI2_MEDIAN))
                boots = np.asarray(boots, dtype=float)
                sd_emp = float(boots.std(ddof=1))

            # Permutation baseline (reuse pairs; just permute columns of Z)
            lam_perm_vals = []
            for _ in range(cfg.permutations):
                # Row-wise independent column permutations to fully destroy cross-locus correlation
                Zp = Z.copy()
                for r in range(Zp.shape[0]):
                    idx_cols = rng.permutation(Zp.shape[1])
                    Zp[r] = Zp[r, idx_cols]
                chi2p = _chi2_pairs(Zp, pairs, N)
                lam_perm_vals.append(float(np.median(chi2p) / CHI2_MEDIAN))
            lam_perm_vals = np.asarray(lam_perm_vals, dtype=float)
            mu_perm = float(lam_perm_vals.mean())
            s_perm = float(lam_perm_vals.std(ddof=1))

            # z-test combining uncertainties
            var_mu_perm = (s_perm ** 2) / max(cfg.permutations, 1)
            var_emp = 0.0 if np.isnan(sd_emp) else (sd_emp ** 2)
            se = float(np.sqrt(var_emp + var_mu_perm))
            if se > 0:
                z = (lam_emp - mu_perm) / se
                p_one = 1.0 - float(ss.norm.cdf(z))
                det = int(z > cfg.lambda_alpha_z)
            else:
                z = np.nan
                p_one = np.nan
                det = 0

            rows.append({
                "N": N,
                "replicate": rep,
                "lambda_emp": lam_emp,
                "mu_perm": mu_perm,
                "sd_emp": sd_emp,
                "s_perm": s_perm,
                "z": z,
                "p_one_sided": p_one,
                "detect": det,
            })

    out = pd.DataFrame(rows)
    cfg.output_dir.mkdir(parents=True, exist_ok=True)
    out.to_csv(cfg.output_dir / "fpr_ztest_details.csv", index=False)
    summ = (
        out.groupby(["N"], as_index=False)
        .agg(
            fpr_ztest=("detect", "mean"),
            lambda_median=("lambda_emp", "median"),
            mu_perm_mean=("mu_perm", "mean"),
            z_median=("z", "median"),
        )
        .sort_values(["N"])
    )
    summ.to_csv(cfg.output_dir / "fpr_ztest_summary.csv", index=False)
    print(f"Wrote {cfg.output_dir/'fpr_ztest_details.csv'} and {cfg.output_dir/'fpr_ztest_summary.csv'}")


def parse_args(argv: Sequence[str] | None = None) -> argparse.Namespace:
    p = argparse.ArgumentParser(description="Null FPR for λ_XC via low-cost z-test with permutations")
    p.add_argument("--N-grid", type=int, nargs="+", required=True)
    p.add_argument("--snps-per-chr", type=int, default=400)
    p.add_argument("--pairs", type=int, default=10000)
    p.add_argument("--permutations", type=int, default=20)
    p.add_argument("--bootstraps", type=int, default=100)
    p.add_argument("--replicates", type=int, default=50)
    p.add_argument("--z-alpha-threshold", type=float, default=1.645, help="one-sided z threshold (default 1.645 for α=0.05)")
    p.add_argument("--seed", type=int, default=2025)
    p.add_argument("--output-dir", type=Path, default=Path("simulations/fpr/output_null_ztest"))
    return p.parse_args(argv)


def main(argv: Sequence[str] | None = None) -> None:
    args = parse_args(argv)
    cfg = Config(
        N_grid=tuple(args.N_grid),
        snps_per_chr=int(args.snps_per_chr),
        pairs=int(args.pairs),
        permutations=int(args.permutations),
        bootstraps=int(args.bootstraps),
        replicates=int(args.replicates),
        lambda_alpha_z=float(args.z_alpha_threshold),
        seed=int(args.seed),
        output_dir=args.output_dir.expanduser().resolve(),
    )
    run(cfg)


if __name__ == "__main__":
    main()
