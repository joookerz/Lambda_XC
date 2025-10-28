"""Compare detection power of λ_XC vs PCA under weak structure.

Two models:
  - mixture: two-population mixture with target F_ST
  - cline: continuous (IBD-like) ancestry gradient across individuals

Detection rules:
  - λ_XC: z-test vs permutation baseline mean using
           z = (λ_emp − μ_perm) / sqrt(sd_emp² + s_perm²/B)
           one-sided at α via z > z_alpha
  - PCA : Tracy–Widom (TW) threshold for the leading eigenvalue of GRM
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


# -------------------------- utilities --------------------------

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


def _lambda_emp(G: np.ndarray, valid: np.ndarray, chr_labels: np.ndarray, n_pairs: int, n: int, rng: np.random.Generator) -> float:
    pairs = _random_pairs(chr_labels, n_pairs, rng)
    chi2 = _chi2_pairs(G, valid, pairs, n)
    if chi2.size == 0:
        return np.nan
    return float(np.median(chi2) / CHI2_MEDIAN)


def _lambda_perm_dist(geno: np.ndarray, chr_labels: np.ndarray, n_pairs: int, n_perm: int, rng: np.random.Generator) -> np.ndarray:
    vals = []
    for _ in range(n_perm):
        # row-wise independent permutations across individuals
        Gperm = geno.copy()
        for r in range(Gperm.shape[0]):
            idx = rng.permutation(Gperm.shape[1])
            Gperm[r] = Gperm[r, idx]
        Gstd, valid = _standardize(Gperm)
        val = _lambda_emp(Gstd, valid, chr_labels, n_pairs, Gperm.shape[1], rng)
        vals.append(val)
    return np.asarray(vals, dtype=float)


def _grm_leading_eig(G: np.ndarray) -> float:
    # GRM: K = Z^T Z / M ; leading eig of K equals (s1^2) / M where s1 is leading singular value of Z
    # To avoid full SVD, compute eigenvalues of N x N matrix via np.linalg.eigvalsh on K
    M = G.shape[0]
    K = (G.T @ G) / float(M)
    w = np.linalg.eigvalsh(K)
    return float(w[-1])


# PCA detection uses Tracy–Widom threshold (no permutation baseline).


def _tw_quantile(p: float) -> float:
    """Approximate upper-tail quantiles for Tracy–Widom (β=1)."""
    table = {
        0.90: 0.4501,
        0.95: 0.9793,
        0.975: 1.2670,
        0.99: 2.0234,
    }
    if p in table:
        return table[p]
    keys = np.array(sorted(table.keys()))
    return float(table[float(keys[np.argmin(np.abs(keys - p))])])


def _pca_tw_threshold(m_valid: int, n: int, quantile: float = 0.975) -> float:
    """TW-based upper threshold for top eigenvalue of K=(Z^T Z)/m_valid."""
    gamma = n / float(m_valid)
    mu = (1.0 + np.sqrt(gamma)) ** 2
    sigma = (1.0 + np.sqrt(gamma)) * ((1.0 / np.sqrt(m_valid)) + (1.0 / np.sqrt(n))) ** (2.0 / 3.0)
    tw_q = _tw_quantile(quantile)
    return float(mu + sigma * tw_q)


# -------------------------- models --------------------------

def simulate_mixture(n: int, snps_per_chr: int, fst: float, rng: np.random.Generator, min_maf: float = 0.05) -> tuple[np.ndarray, np.ndarray]:
    if n % 2 != 0:
        raise ValueError("N must be even for mixture model.")
    half = n // 2
    n_chrom = 22
    m = n_chrom * snps_per_chr
    chr_labels = np.repeat(np.arange(1, n_chrom + 1), snps_per_chr)

    base = rng.uniform(min_maf, 1 - min_maf, size=m)
    delta = np.sqrt(2.0 * fst * base * (1.0 - base))
    p1 = np.clip(base + 0.5 * delta, 1e-4, 1 - 1e-4)
    p2 = np.clip(base - 0.5 * delta, 1e-4, 1 - 1e-4)

    g1 = rng.binomial(2, p1[:, None], size=(m, half))
    g2 = rng.binomial(2, p2[:, None], size=(m, half))
    geno = np.concatenate([g1, g2], axis=1).astype(float)
    return geno, chr_labels


def simulate_cline(n: int, snps_per_chr: int, tau: float, rng: np.random.Generator, min_maf: float = 0.05) -> tuple[np.ndarray, np.ndarray]:
    # tau controls cline strength (analog of sqrt(mean Δ^2 / [2p(1-p)])).
    n_chrom = 22
    m = n_chrom * snps_per_chr
    chr_labels = np.repeat(np.arange(1, n_chrom + 1), snps_per_chr)

    base = rng.uniform(min_maf, 1 - min_maf, size=m)
    sign = rng.choice([-1.0, 1.0], size=m)
    Delta = sign * tau * np.sqrt(2.0 * base * (1.0 - base))

    A = rng.uniform(0.0, 1.0, size=n)
    A = A - A.mean()
    A = A / (A.std(ddof=1) + 1e-12)
    A = A / 4.0  # small gradient (tunable via tau)

    geno = np.empty((m, n), dtype=float)
    for i in range(m):
        p_i = np.clip(base[i] + Delta[i] * A, 1e-4, 1 - 1e-4)
        geno[i] = rng.binomial(2, p_i)
    return geno, chr_labels


# -------------------------- experiment --------------------------

@dataclass
class Config:
    model: str
    N_grid: tuple[int, ...]
    effect_grid: tuple[float, ...]
    snps_per_chr: int
    n_pairs: int
    n_perm: int
    replicates: int
    seed: int
    output_dir: Path
    lambda_boot: int = 100         # bootstrap resamples for λ_emp SD
    lambda_z_alpha: float = 1.645  # one-sided z threshold (α=0.05)
    pca_tw_quantile: float = 0.975 # upper-tail quantile for TW threshold


def run_compare(cfg: Config) -> None:
    rng_master = np.random.default_rng(cfg.seed)
    seed_seq = np.random.SeedSequence(cfg.seed)
    seeds = seed_seq.spawn(len(cfg.N_grid) * len(cfg.effect_grid) * cfg.replicates)

    rows = []
    idx = 0
    for N in cfg.N_grid:
        for eff in cfg.effect_grid:
            for rep in range(cfg.replicates):
                rng = np.random.default_rng(seeds[idx])
                idx += 1

                if cfg.model == "mixture":
                    geno, chr_labels = simulate_mixture(N, cfg.snps_per_chr, eff, rng)
                elif cfg.model == "cline":
                    geno, chr_labels = simulate_cline(N, cfg.snps_per_chr, eff, rng)
                else:
                    raise ValueError("model must be 'mixture' or 'cline'")

                G, valid = _standardize(geno)

                # λ_XC z-test vs permutation baseline mean
                # 1) Pairs (reuse for permutations)
                pairs_list: List[Tuple[int, int]] = list(_random_pairs(chr_labels, cfg.n_pairs, rng))
                chi2 = _chi2_pairs(G, valid, pairs_list, N)
                lam_emp = float(np.median(chi2) / CHI2_MEDIAN)
                # 2) Bootstrap SD for λ_emp
                if chi2.size > 1 and cfg.lambda_boot > 0:
                    boots = []
                    for _ in range(cfg.lambda_boot):
                        res = rng.choice(chi2, size=chi2.size, replace=True)
                        boots.append(float(np.median(res) / CHI2_MEDIAN))
                    sd_emp = float(np.asarray(boots, dtype=float).std(ddof=1))
                else:
                    sd_emp = np.nan
                # 3) Permutation baseline mean and SD
                lam_perm_vals = []
                for _ in range(max(cfg.n_perm, 1)):
                    # independent column permutation per row to destroy cross-locus correlation
                    Zp = G.copy()
                    for r in range(Zp.shape[0]):
                        idx_cols = rng.permutation(Zp.shape[1])
                        Zp[r] = Zp[r, idx_cols]
                    chi2p = _chi2_pairs(Zp, valid, pairs_list, N)
                    lam_pv = float(np.median(chi2p) / CHI2_MEDIAN)
                    lam_perm_vals.append(lam_pv)
                lam_perm_vals = np.asarray(lam_perm_vals, dtype=float)
                mu_perm = float(lam_perm_vals.mean())
                s_perm = float(lam_perm_vals.std(ddof=1))
                # 4) z-test
                var_emp = 0.0 if np.isnan(sd_emp) else (sd_emp ** 2)
                var_mu = (s_perm ** 2) / max(cfg.n_perm, 1)
                se = float(np.sqrt(var_emp + var_mu))
                if se > 0:
                    z_val = (lam_emp - mu_perm) / se
                    lam_p = 1.0 - float(ss.norm.cdf(z_val))
                    det_xc = bool(z_val > cfg.lambda_z_alpha)
                else:
                    z_val = np.nan
                    lam_p = np.nan
                    det_xc = False

                # PCA detection: TW threshold
                G_ = G[valid]
                eig1_obs = _grm_leading_eig(G_)
                m_valid = int(G_.shape[0])
                eig1_thr = _pca_tw_threshold(m_valid, N, quantile=cfg.pca_tw_quantile)
                eig1_p = np.nan
                det_pca = bool(eig1_obs >= eig1_thr)

                rows.append(
                    {
                        "model": cfg.model,
                        "N": N,
                        "effect": eff,
                        "replicate": rep,
                        "lambda_emp": lam_emp,
                        "lambda_sd_emp": sd_emp,
                        "lambda_mu_perm": mu_perm,
                        "lambda_s_perm": s_perm,
                        "lambda_z": z_val,
                        "lambda_p_one_sided": lam_p,
                        "detect_lambda": int(det_xc),
                        "eig1_obs": eig1_obs,
                        "eig1_thr97": eig1_thr,
                        "eig1_p": eig1_p,
                        "detect_pca": int(det_pca),
                    }
                )

    out_dir = cfg.output_dir
    out_dir.mkdir(parents=True, exist_ok=True)
    details = pd.DataFrame(rows)
    details.to_csv(out_dir / "power_details.csv", index=False)

    summary = (
        details.groupby(["model", "N", "effect"], as_index=False)
        .agg(
            lambda_detect_rate=("detect_lambda", "mean"),
            pca_detect_rate=("detect_pca", "mean"),
            lambda_median=("lambda_emp", "median"),
            eig1_median=("eig1_obs", "median"),
        )
        .sort_values(["N", "effect"])
    )
    summary.to_csv(out_dir / "power_summary.csv", index=False)
    print(f"Wrote {out_dir/'power_details.csv'} and {out_dir/'power_summary.csv'}")


def parse_args(argv: Sequence[str] | None = None) -> argparse.Namespace:
    p = argparse.ArgumentParser(description="Compare λ_XC vs PCA detection power.")
    p.add_argument("--model", choices=["mixture", "cline"], required=True)
    p.add_argument("--N-grid", type=int, nargs="+", required=True)
    p.add_argument("--effect-grid", type=float, nargs="+", required=True, help="mixture: F_ST; cline: tau")
    p.add_argument("--snps-per-chr", type=int, default=800)
    p.add_argument("--pairs", type=int, default=100000)
    p.add_argument("--permutations", type=int, default=20, help="Permutations for λ_XC baseline mean (default 20)")
    p.add_argument("--replicates", type=int, default=30)
    p.add_argument("--seed", type=int, default=2026)
    p.add_argument("--output-dir", type=Path, default=Path("simulations/compare/output"))
    p.add_argument("--lambda-bootstraps", type=int, default=100, help="Bootstrap resamples for λ_emp SD (default 100)")
    p.add_argument("--lambda-z-alpha", type=float, default=1.645, help="One-sided z threshold for λ_XC (default 1.645)")
    p.add_argument("--pca-tw-quantile", type=float, default=0.975, help="Upper-tail quantile for TW threshold (default 0.975)")
    return p.parse_args(argv)


def main(argv: Sequence[str] | None = None) -> None:
    args = parse_args(argv)
    cfg = Config(
        model=args.model,
        N_grid=tuple(args.N_grid),
        effect_grid=tuple(args.effect_grid),
        snps_per_chr=int(args.snps_per_chr),
        n_pairs=int(args.pairs),
        n_perm=int(args.permutations),
        replicates=int(args.replicates),
        seed=int(args.seed),
        output_dir=args.output_dir.expanduser().resolve(),
        lambda_boot=int(args.lambda_bootstraps),
        lambda_z_alpha=float(args.lambda_z_alpha),
        pca_tw_quantile=float(args.pca_tw_quantile),
    )
    run_compare(cfg)


if __name__ == "__main__":
    main()
