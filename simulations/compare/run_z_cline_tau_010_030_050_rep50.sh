#!/usr/bin/env bash
set -euo pipefail

# Cline power test (rep=50) with stronger tau grid: 0.010, 0.030, 0.050
# λ_XC uses z-test (B=20 perms, K=100 boot), PCA uses TW (q=0.975)

OUT="simulations/compare/output/cline_z_tau_010_030_050_rep50"
mkdir -p "$OUT"

nohup python simulations/compare/compare_power.py \
  --model cline \
  --N-grid 100 300 600 1000 \
  --effect-grid 0.010 0.030 0.050 \
  --snps-per-chr 800 \
  --pairs 50000 \
  --permutations 20 \
  --lambda-bootstraps 100 \
  --replicates 50 \
  --lambda-z-alpha 1.645 \
  --pca-tw-quantile 0.975 \
  --seed 9201 \
  --output-dir "$OUT" \
  > "$OUT/run.log" 2>&1 &

echo "Started cline (tau=0.010,0.030,0.050; rep=50). Tail log with:"
echo "  tail -f $OUT/run.log"

