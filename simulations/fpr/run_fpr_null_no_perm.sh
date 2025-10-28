#!/usr/bin/env bash
set -euo pipefail

# Null FPR for λ_XC without permutations (threshold vs 1.0 by default)

python simulations/fpr/fpr_lambda_no_perm.py \
  --N-grid 100 300 1000 \
  --snps-per-chr 400 \
  --pairs 20000 \
  --replicates 100 \
  --lambda-threshold 1.0 \
  --seed 2025 \
  --output-dir simulations/fpr/output_null_no_perm
