#!/usr/bin/env bash
set -euo pipefail

# High-precision cline grid with λ_XC z-test (B=20, K=100), replicates=50

python simulations/compare/compare_power.py \
  --model cline \
  --N-grid 100 150 200 300 \
  --effect-grid 0.005 0.010 0.015 \
  --snps-per-chr 800 \
  --pairs 50000 \
  --permutations 20 \
  --lambda-bootstraps 100 \
  --replicates 50 \
  --lambda-z-alpha 1.645 \
  --pca-tw-quantile 0.975 \
  --seed 8123 \
  --output-dir simulations/compare/output/cline_z_grid_rep50

