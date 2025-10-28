"""Null (no-structure) false positive rate for λ_XC vs PCA (optional).

Simulates a single panmictic population (no mixture, fst=0),
computes λ_emp from cross-chromosome pairs, and estimates:
  - FPR_lam_vs_one: Pr(λ_emp >= 1.0)
  - FPR_lam_perm: Pr(λ_emp >= q_{0.975}(λ_perm)) using column permutations
  - (optional) PCA FPR via TW threshold

Outputs per-replicate details and a summary by sample size.
"""

from __future__ import annotations

import argparse
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable, Iterator, Sequence, Tuple

import numpy as np
import pandas as pd

CHI2_MEDIAN = 0.4559364


def _random_pairs(chr_labels: np.ndarray, n_pairs: int, rng: np.random.Generator) -> Iterator[Tuple[int, int]]:
    bins: dict[int, list[int]] = {}
    for idx, chrom in enumerate(chr_labels):
        bins.setdefault(int(chrom), []).append(int(idx))
    chromosomes = np.array(list(bins.keys()), dtype=int)
    for _ in range(n_pairs):
        c1, c2 = rng.choice(chromosomes, 2, replace=False)
        yield rng.choice(bins[int(c1)]), rng.choice(bins[int(c2)])


def _chi2_pairs(G: np.ndarray, valid: np.ndarray, pairs: Iterable[Tuple[int, int]], n_samples: int) -> np.ndarray:
    stats: list[float] = []
    denom = float(n_samples - 1)
    for i, j in pairs:
        if not (valid[i] and valid[j]):
            continue
        r = float(np.dot(G[i], G[j]) / denom)
        stats.append((r * r) * n_samples)
    return np.asarray(stats, dtype=float)


def _standardize(geno: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    means = geno.mean(axis=1)
    centered = geno - means[:, None]
    std = centered.std(axis=1, ddof=1)
    valid = std > 0
    G = np.zeros_like(centered)
    G[valid] = centered[valid] / std[valid, None]
    return G, valid


def simulate_null(n: int, snps_per_chr: int, rng: np.random.Generator, min_maf: float = 0.05) -> tuple[np.ndarray, np.ndarray]:
    """Single-population null: independent SNPs with constant allele freq across individuals."""
    n_chrom = 22
    m = n_chrom * snps_per_chr
    chr_labels = np.repeat(np.arange(1, n_chrom + 1), snps_per_chr)
    p = rng.uniform(min_maf, 1 - min_maf, size=m)
    geno = rng.binomial(2, p[:, None], size=(m, n)).astype(float)
    return geno, chr_labels


def lambda_emp(geno: np.ndarray, chr_labels: np.ndarray, n_pairs: int, rng: np.random.Generator) -> float:
    G, valid = _standardize(geno)
    pairs = _random_pairs(chr_labels, n_pairs, rng)
    chi2 = _chi2_pairs(G, valid, pairs, geno.shape[1])
    if chi2.size == 0:
        return np.nan
    return float(np.median(chi2) / CHI2_MEDIAN)


def lambda_perm_quantile(geno: np.ndarray, chr_labels: np.ndarray, n_pairs: int, n_perm: int, q: float, rng: np.random.Generator) -> float:
    vals = []
    for _ in range(n_perm):
        Gperm = geno.copy()
        for r in range(Gperm.shape[0]):
            idx = rng.permutation(Gperm.shape[1])
            Gperm[r] = Gperm[r, idx]
        vals.append(lambda_emp(Gperm, chr_labels, n_pairs, rng))
    return float(np.nanpercentile(np.asarray(vals, float), q * 100.0))


@dataclass
class Config:
    N_grid: tuple[int, ...]
    snps_per_chr: int
    pairs: int
    permutations: int
    replicates: int
    seed: int
    output_dir: Path


def run(cfg: Config) -> None:
    ss = np.random.SeedSequence(cfg.seed)
    seeds = ss.spawn(len(cfg.N_grid) * cfg.replicates)
    rows = []
    k = 0
    for N in cfg.N_grid:
        for rep in range(cfg.replicates):
            rng = np.random.default_rng(seeds[k]); k += 1
            geno, chr_labels = simulate_null(N, cfg.snps_per_chr, rng)
            lam = lambda_emp(geno, chr_labels, cfg.pairs, rng)
            thr = lambda_perm_quantile(geno, chr_labels, cfg.pairs, cfg.permutations, 0.975, rng)
            det_vs_one = int(lam >= 1.0)
            det_vs_perm = int(lam >= thr)
            rows.append({
                "N": N,
                "replicate": rep,
                "lambda_emp": lam,
                "lambda_thr_perm97": thr,
                "detect_vs_one": det_vs_one,
                "detect_vs_perm": det_vs_perm,
            })

    out_dir = cfg.output_dir
    out_dir.mkdir(parents=True, exist_ok=True)
    details = pd.DataFrame(rows)
    details.to_csv(out_dir / "fpr_details.csv", index=False)

    summary = (
        details.groupby(["N"], as_index=False)
        .agg(
            fpr_vs_one=("detect_vs_one", "mean"),
            fpr_vs_perm=("detect_vs_perm", "mean"),
            lambda_median=("lambda_emp", "median"),
        )
        .sort_values(["N"])
    )
    summary.to_csv(out_dir / "fpr_summary.csv", index=False)
    print(f"Wrote {out_dir/'fpr_details.csv'} and {out_dir/'fpr_summary.csv'}")


def parse_args(argv: Sequence[str] | None = None) -> argparse.Namespace:
    p = argparse.ArgumentParser(description="Null FPR for λ_XC detection")
    p.add_argument("--N-grid", type=int, nargs="+", required=True)
    p.add_argument("--snps-per-chr", type=int, default=400)
    p.add_argument("--pairs", type=int, default=20000)
    p.add_argument("--permutations", type=int, default=50)
    p.add_argument("--replicates", type=int, default=100)
    p.add_argument("--seed", type=int, default=2025)
    p.add_argument("--output-dir", type=Path, default=Path("simulations/fpr/output_null"))
    return p.parse_args(argv)


def main(argv: Sequence[str] | None = None) -> None:
    args = parse_args(argv)
    cfg = Config(
        N_grid=tuple(args.N_grid),
        snps_per_chr=int(args.snps_per_chr),
        pairs=int(args.pairs),
        permutations=int(args.permutations),
        replicates=int(args.replicates),
        seed=int(args.seed),
        output_dir=args.output_dir.expanduser().resolve(),
    )
    run(cfg)


if __name__ == "__main__":
    main()
