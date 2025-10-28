#!/usr/bin/env bash
set -euo pipefail

# Grid over N and cline strength (tau), using TW for PCA and λ_XC threshold 1.0

python simulations/compare/compare_power.py \
  --model cline \
  --N-grid 100 300 600 1000 \
  --effect-grid 0.005 0.010 0.015 \
  --snps-per-chr 400 \
  --pairs 50000 \
  --permutations 50 \
  --replicates 10 \
  --lambda-threshold 1.0 \
  --pca-tw-quantile 0.975 \
  --seed 3031 \
  --output-dir simulations/compare/output/cline_grid
