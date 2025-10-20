# λ<sub>XC</sub> Simulation Toolkit

Refactored codebase for the simulations in *“A Cross-Chromosome Inflation
Factor for Detecting Subtle Population Structure.”*

## Features

- Distance-based SNP pruning with fixed density per chromosome.
- Cross-chromosome χ² calculations, permutation baseline, and bootstrap CIs.
- CLI wrapper for batch runs and QQ-plot generation.

## Layout

```
src/lambda_xc/      core package (data loading, analysis, plots, CLI)
scripts/            helper entry points (e.g. run_analysis.py)
data/               place HapMap3 genotype files here
pyproject.toml      dependencies and console script
```

## Usage

```bash
cd lambda_xc_project
python -m venv .venv
source .venv/bin/activate
pip install -e .

lambda-xc-run \
  --data-root /path/to/hapmap3 \
  --qq-path results/qq.png \
  --results-csv results/lambda.csv
```

Adjust flags such as `--panels`, `--pairs`, or `--permutations` to explore
other scenarios.
