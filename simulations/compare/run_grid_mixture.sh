#!/usr/bin/env bash
set -euo pipefail

# Finer grid over N and F_ST for two-population mixture, using TW for PCA

python simulations/compare/compare_power.py \
  --model mixture \
  --N-grid 100 200 300 500 800 1000 \
  --effect-grid 0.001 0.002 0.004 0.006 0.008 0.010 0.015 \
  --snps-per-chr 600 \
  --pairs 50000 \
  --permutations 50 \
  --replicates 10 \
  --lambda-threshold 1.0 \
  --pca-tw-quantile 0.975 \
  --seed 4042 \
  --output-dir simulations/compare/output/mixture_grid

