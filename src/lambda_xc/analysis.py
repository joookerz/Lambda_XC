"""Core simulation routines for the cross-chromosome inflation project."""

from __future__ import annotations

from dataclasses import dataclass
from typing import Dict, Iterable, Iterator, List, Sequence, Tuple

import numpy as np
import scipy.stats as ss

from .data import HapMapDataset

CHI2_MEDIAN = 0.4559364


def _median_lambda(values: np.ndarray) -> float:
    return float(np.median(values) / CHI2_MEDIAN)


def _random_pairs(chr_labels: np.ndarray, n_pairs: int, rng: np.random.Generator) -> Iterator[Tuple[int, int]]:
    bins: Dict[int, List[int]] = {}
    for idx, chrom in enumerate(chr_labels):
        bins.setdefault(int(chrom), []).append(int(idx))
    chromosomes = np.array(list(bins.keys()), dtype=int)

    for _ in range(n_pairs):
        c1, c2 = rng.choice(chromosomes, 2, replace=False)
        yield rng.choice(bins[int(c1)]), rng.choice(bins[int(c2)])


def _chi2_pairs(
    geno_matrix: np.ma.MaskedArray,
    pairs: Iterable[Tuple[int, int]],
    min_samples: int = 4,
) -> np.ndarray:
    stats: List[float] = []
    for i, j in pairs:
        g1 = geno_matrix[int(i)]
        g2 = geno_matrix[int(j)]
        ok = (~g1.mask) & (~g2.mask)
        n = int(ok.sum())
        if n < min_samples:
            continue
        r, _ = ss.pearsonr(g1[ok], g2[ok])
        stats.append((r * r) * n)
    return np.asarray(stats, dtype=float)


def permutation_baseline(
    geno_matrix: np.ma.MaskedArray,
    chr_labels: np.ndarray,
    n_pairs: int,
    n_perm: int,
    rng: np.random.Generator,
) -> Tuple[float, float, float, float]:
    """Estimate the finite-sample baseline through column permutations."""
    medians = []
    for _ in range(n_perm):
        permuted = geno_matrix[:, rng.permutation(geno_matrix.shape[1])]
        stats = _chi2_pairs(permuted, _random_pairs(chr_labels, n_pairs, rng))
        if stats.size == 0:
            raise RuntimeError("No χ² statistics produced during permutation baseline.")
        medians.append(_median_lambda(stats))
    values = np.asarray(medians, dtype=float)
    mean = float(values.mean())
    sd = float(values.std(ddof=1))
    lo, hi = np.percentile(values, [2.5, 97.5])
    return mean, sd, float(lo), float(hi)


def bootstrap_inflation(
    chi2_values: np.ndarray,
    baseline_mean: float,
    n_boot: int,
    rng: np.random.Generator,
) -> Tuple[float, float, float, float, float, float, float]:
    """Return λ_emp, λ_corr and their bootstrap SD / percentile CI."""
    lam_emp = _median_lambda(chi2_values)
    lam_corr = float(lam_emp / baseline_mean)

    emp_samples = []
    cor_samples = []
    for _ in range(n_boot):
        resampled = rng.choice(chi2_values, size=chi2_values.size, replace=True)
        lam = _median_lambda(resampled)
        emp_samples.append(lam)
        cor_samples.append(lam / baseline_mean)

    def sd_ci(samples: Sequence[float]) -> Tuple[float, float, float]:
        arr = np.asarray(samples, dtype=float)
        sd = float(arr.std(ddof=1))
        lo, hi = np.percentile(arr, [2.5, 97.5])
        return sd, float(lo), float(hi)

    sd_emp, lo_emp, hi_emp = sd_ci(emp_samples)
    sd_cor, lo_cor, hi_cor = sd_ci(cor_samples)
    return lam_emp, sd_emp, lo_emp, hi_emp, lam_corr, sd_cor, lo_cor, hi_cor


def fst_average(dataset: HapMapDataset, populations: Sequence[str], snp_indices: np.ndarray) -> float:
    """Average-of-ratios FST estimate for a pair of populations."""
    if len(populations) != 2:
        raise ValueError("FST requires exactly two populations.")
    p1, p2 = populations
    g1 = dataset.load_population(p1)[snp_indices]
    g2 = dataset.load_population(p2)[snp_indices]
    n1 = (~g1.mask).sum(axis=1).astype(float)
    n2 = (~g2.mask).sum(axis=1).astype(float)
    p1_hat = np.divide(g1.filled(0).sum(axis=1), 2 * n1, where=n1 > 0, out=np.zeros_like(n1))
    p2_hat = np.divide(g2.filled(0).sum(axis=1), 2 * n2, where=n2 > 0, out=np.zeros_like(n2))
    pq = 0.5 * (p1_hat + p2_hat)
    pq = pq * (1 - pq)
    delta2 = (p1_hat - p2_hat) ** 2
    inv1 = np.divide(1.0, n1, where=n1 > 0, out=np.zeros_like(n1))
    inv2 = np.divide(1.0, n2, where=n2 > 0, out=np.zeros_like(n2))
    num = delta2 - (pq / 2.0) * (inv1 + inv2)
    den = 2.0 * pq
    with np.errstate(divide="ignore", invalid="ignore"):
        ratio = np.divide(num, den, where=(den > 0), out=np.full_like(num, np.nan))
    ratio = ratio[~np.isnan(ratio)]
    return float(np.nan if ratio.size == 0 else ratio.mean())


def theoretical_lambda(populations: Sequence[str], fst: float, sample_size: int) -> float:
    if len(populations) == 1:
        return 1.0
    return float(1.0 + 1.096 * sample_size * (fst ** 2))


@dataclass
class PanelConfig:
    label: str
    populations: Tuple[str, ...]


@dataclass
class PanelResult:
    label: str
    sample_size: int
    lambda_theory: float
    lambda_emp: float
    lambda_emp_sd: float
    lambda_emp_ci: Tuple[float, float]
    baseline: float
    baseline_sd: float
    baseline_ci: Tuple[float, float]
    lambda_corr: float
    lambda_corr_sd: float
    lambda_corr_ci: Tuple[float, float]
    chi2_values: np.ndarray


def analyse_panel(
    dataset: HapMapDataset,
    panel: PanelConfig,
    snp_indices: np.ndarray,
    chr_labels: np.ndarray,
    n_pairs: int,
    n_perm: int,
    n_boot: int,
    rng: np.random.Generator,
) -> PanelResult:
    mats = [dataset.load_population(pop)[snp_indices] for pop in panel.populations]
    geno_matrix = mats[0] if len(mats) == 1 else np.ma.concatenate(mats, axis=1)
    sample_size = int(geno_matrix.shape[1])

    baseline_mean, baseline_sd, base_lo, base_hi = permutation_baseline(
        geno_matrix, chr_labels, n_pairs, n_perm, rng
    )
    chi2_values = _chi2_pairs(geno_matrix, _random_pairs(chr_labels, n_pairs, rng))
    if chi2_values.size == 0:
        raise RuntimeError(f"No χ² statistics for panel {panel.label}")
    (
        lambda_emp,
        lambda_emp_sd,
        lam_lo,
        lam_hi,
        lambda_corr,
        lambda_corr_sd,
        cor_lo,
        cor_hi,
    ) = bootstrap_inflation(chi2_values, baseline_mean, n_boot, rng)

    if len(panel.populations) == 2:
        fst = fst_average(dataset, panel.populations, snp_indices)
    else:
        fst = np.nan
    lambda_th = theoretical_lambda(panel.populations, fst, sample_size)

    return PanelResult(
        label=panel.label,
        sample_size=sample_size,
        lambda_theory=lambda_th,
        lambda_emp=lambda_emp,
        lambda_emp_sd=lambda_emp_sd,
        lambda_emp_ci=(lam_lo, lam_hi),
        baseline=baseline_mean,
        baseline_sd=baseline_sd,
        baseline_ci=(base_lo, base_hi),
        lambda_corr=lambda_corr,
        lambda_corr_sd=lambda_corr_sd,
        lambda_corr_ci=(cor_lo, cor_hi),
        chi2_values=chi2_values,
    )
