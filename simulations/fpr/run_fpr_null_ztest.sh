#!/usr/bin/env bash
set -euo pipefail

# Null FPR for λ_XC using low-cost z-test with B=20 permutations, pairs=10k

python simulations/fpr/fpr_lambda_ztest.py \
  --N-grid 100 300 1000 \
  --snps-per-chr 400 \
  --pairs 10000 \
  --permutations 20 \
  --bootstraps 100 \
  --replicates 50 \
  --z-alpha-threshold 1.645 \
  --seed 2025 \
  --output-dir simulations/fpr/output_null_ztest

