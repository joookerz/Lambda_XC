#!/usr/bin/env bash
# Wrapper for the pure λ_XC simulation with tuned parameters.

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(cd "${SCRIPT_DIR}/.." && pwd)"

python "${SCRIPT_DIR}/pure_simulation.py" \
  --sample-grid 200 500 1000 2000 \
  --fst-grid 0.001 0.01 0.03 0.05 0.07 0.10 \
  --snps-per-chr 800 \
  --pairs 80000 \
  --bootstraps 400 \
  --seed 2025 \
  --output-dir "${PROJECT_ROOT}/simulations/output" \
  --store-bootstrap

