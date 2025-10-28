#!/usr/bin/env bash
set -euo pipefail

# 100-replicate power test (mixture) on three representative grid points.
# Uses λ_XC threshold = 1.0 and PCA Tracy–Widom at 0.975.

python simulations/compare/compare_power.py \
  --model mixture \
  --N-grid 100 300 1000 \
  --effect-grid 0.002 0.006 0.010 \
  --snps-per-chr 400 \
  --pairs 20000 \
  --permutations 20 \
  --replicates 100 \
  --lambda-threshold 1.0 \
  --pca-tw-quantile 0.975 \
  --seed 7002 \
  --output-dir simulations/compare/output/mixture_power100

