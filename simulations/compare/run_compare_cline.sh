#!/usr/bin/env bash
set -euo pipefail

python simulations/compare/compare_power.py \
  --model cline \
  --N-grid 100 150 200 \
  --effect-grid 0.005 0.010 0.015 \
  --snps-per-chr 800 \
  --pairs 100000 \
  --permutations 50 \
  --replicates 30 \
  --seed 2026 \
  --output-dir simulations/compare/output

