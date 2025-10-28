#!/usr/bin/env bash
set -euo pipefail

# Null FPR for λ_XC with and without permutation baseline.

python simulations/fpr/fpr_lambda.py \
  --N-grid 100 300 1000 \
  --snps-per-chr 400 \
  --pairs 20000 \
  --permutations 50 \
  --replicates 100 \
  --seed 2025 \
  --output-dir simulations/fpr/output_null
