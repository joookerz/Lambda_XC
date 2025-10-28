#!/usr/bin/env bash
set -euo pipefail

# High-precision mixture grid with λ_XC z-test (B=20, K=100), replicates=50

python simulations/compare/compare_power.py \
  --model mixture \
  --N-grid 100 300 600 1000 \
  --effect-grid 0.002 0.004 0.006 0.008 0.010 \
  --snps-per-chr 600 \
  --pairs 50000 \
  --permutations 20 \
  --lambda-bootstraps 100 \
  --replicates 50 \
  --lambda-z-alpha 1.645 \
  --pca-tw-quantile 0.975 \
  --seed 8124 \
  --output-dir simulations/compare/output/mixture_z_grid_rep50

