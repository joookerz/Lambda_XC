#!/usr/bin/env bash
set -euo pipefail

# Cline power test split by N (rep=50) with tau grid 0.010, 0.030, 0.050

BASE="simulations/compare/output/cline_z_tau_010_030_050_rep50_split"
mkdir -p "$BASE"

for N in 100 300 600 1000; do
  OUT="$BASE/N$N"
  mkdir -p "$OUT"
  nohup python simulations/compare/compare_power.py \
    --model cline \
    --N-grid $N \
    --effect-grid 0.010 0.030 0.050 \
    --snps-per-chr 800 \
    --pairs 50000 \
    --permutations 20 \
    --lambda-bootstraps 100 \
    --replicates 50 \
    --lambda-z-alpha 1.645 \
    --pca-tw-quantile 0.975 \
    --seed $((9201 + N)) \
    --output-dir "$OUT" \
    > "$OUT/run.log" 2>&1 &
  echo "Launched N=$N -> tail -f $OUT/run.log"
done

