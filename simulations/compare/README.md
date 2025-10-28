# λ_XC vs PCA: Detection Power Comparison

This folder contains a simulation-based comparison of cross-chromosome inflation (λ_XC)
vs. PCA in detecting weak population structure. We include two models:

- Two-population mixtures with target F_ST
- Continuous cline (IBD-like gradient) with tunable strength

Detection criteria:
- λ_XC: z-test vs permutation baseline mean using B permutations and K bootstraps
  - z = (λ_emp − μ_perm) / sqrt(sd_emp² + s_perm² / B)
  - one-sided α via z > z_alpha (default 1.645)
- PCA (standard): Tracy–Widom threshold for the leading eigenvalue of GRM K = Z^T Z / M.

A replicate is called positive if the observed statistic exceeds the 97.5th percentile of its permutation baseline (one-sided, α=0.05). Detection rate is the fraction of positive replicates.

## Quick start

Run the cline (continuous structure) comparison with modest N and weak effects:

```bash
python simulations/compare/compare_power.py \
  --model cline \
  --N-grid 100 150 200 \
  --effect-grid 0.005 0.010 0.015 \
  --snps-per-chr 800 \
  --pairs 100000 \
  --permutations 20 \
  --lambda-bootstraps 100 \
  --replicates 30 \
  --seed 2026 \
  --output-dir simulations/compare/output \
  --lambda-z-alpha 1.645 \
  --pca-tw-quantile 0.975
```

Run the two-pop mixture comparison:

```bash
python simulations/compare/compare_power.py \
  --model mixture \
  --N-grid 100 150 200 \
  --effect-grid 0.002 0.004 0.006 0.010 \
  --snps-per-chr 800 \
  --pairs 100000 \
  --permutations 20 \
  --lambda-bootstraps 100 \
  --replicates 30 \
  --seed 2026 \
  --output-dir simulations/compare/output \
  --lambda-z-alpha 1.645 \
  --pca-tw-quantile 0.975
```

Outputs:
- `power_details.csv`: per-replicate stats and p-values
- `power_summary.csv`: detection rates for λ_XC and PCA across the grid

## Rationale

- In the weak-structure regime, λ_XC - 1 ∝ N, and aggregation over many cross-chromosome pairs allows precise baselines via permutations.
- PCA’s largest eigenvalue detection threshold follows Tracy–Widom/Marchenko–Pastur behavior, which becomes conservative with large M/N and when structure is diffuse/continuous. Under modest N and continuous clines, λ_XC can remain sensitive while PCA fails to achieve significance.

Adjust `--effect-grid` to locate regions where λ_XC outperforms PCA.
