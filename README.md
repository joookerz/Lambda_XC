# λ<sub>XC</sub> Simulation Toolkit

Refactored codebase for the simulations in *“A Cross-Chromosome Inflation
Factor for Detecting Subtle Population Structure.”*

## Features

- Distance-based SNP pruning with fixed density per chromosome.
- Cross-chromosome χ² calculations, permutation baseline, and bootstrap CIs.
- CLI wrapper for batch runs and QQ-plot generation.

## Requirements

- Python ≥ 3.9
- NumPy ≥ 1.22, SciPy ≥ 1.9, Pandas ≥ 1.5, Matplotlib ≥ 3.5
- HapMap3 `.snp/.geno/.ind` files for the populations of interest (e.g. CHB, CEU, TSI, YRI)

## Repository layout

```
src/lambda_xc/      core package (data loading, analysis, plots, CLI)
scripts/            helper entry points (e.g. run_analysis.py)
data/               place HapMap3 genotype files here
pyproject.toml      dependencies and console script
```

## Installation

```bash
git clone https://github.com/joookerz/Lambda_XC.git
cd Lambda_XC/lambda_xc_project
python -m venv .venv
source .venv/bin/activate  # Windows: .venv\Scripts\activate
pip install -e .
```

## Preparing the data

1. Obtain HapMap3 genotype bundles (or the course-supplied dataset).
2. Extract the `.geno`, `.ind`, and `.snp` files into a directory, e.g. `lambda_xc_project/data/hapmap3`.
3. Ensure `HapMap3.snp` is present in the same directory (used as the SNP reference).

The CLI will accept any path via `--data-root`, so you can keep the data outside the repository if preferred.

## Quick start

```bash
lambda-xc-run \
  --data-root data/hapmap3 \
  --qq-path results/qq_panels.png \
  --results-csv results/lambda_summary.csv
```

This command reproduces the paper’s four panels (CHB, CEU, CEU+TSI, CEU+YRI) with default parameters:

- 2,000 SNPs per chromosome after pruning (100 kb window)
- 1,000,000 random cross-chromosome SNP pairs
- 80 column permutations for the baseline
- 100 bootstrap resamples

Results:

- Tabulated summary printed to stdout (theoretical λ, empirical λ, baseline, corrected λ, and CIs)
- Optional CSV/JSON files mirroring the summary (`--results-csv` / `--results-json`)
- Optional QQ plot PNG (`--qq-path`)

### Command-line options

| Flag | Description | Default |
| ---- | ----------- | ------- |
| `--data-root PATH` | Directory containing `.geno/.ind/.snp` files | *required* |
| `--panels` | Panels to analyse (e.g. `CEU`, `CEU+TSI`) | `CHB CEU CEU+TSI CEU+YRI` |
| `--nsnp-per-chr` | Target SNPs per chromosome after pruning | `2000` |
| `--prune-kb` | Distance threshold (kb) for pruning | `100` |
| `--pairs` | Number of random cross-chromosome pairs | `1000000` |
| `--permutations` | Number of column permutations | `80` |
| `--bootstraps` | Bootstrap replicates for CIs | `100` |
| `--seed` | RNG seed | `511` |
| `--qq-path PATH` | Save QQ plot to file | `None` |
| `--results-json PATH` | Write detailed JSON summary | `None` |
| `--results-csv PATH` | Write CSV summary | `None` |

To run without installing the console script you can call the module directly:

```bash
python -m lambda_xc.cli --data-root data/hapmap3 [options...]
```

## Interpreting the output

- **`lambda_emp`**: raw cross-chromosome inflation factor (median χ² scaled by the χ²₁ median).
- **`baseline`**: permutation-derived finite sample baseline (`λ̄_perm`) with SD and percentile CI.
- **`lambda_corr`**: bias-corrected inflation (`λ_emp / λ̄_perm`).
- **Confidence intervals**: 95% percentile intervals from bootstrap (for λ values) or permutations (for baseline).

The CSV/JSON helpers make it easy to incorporate the metrics into tables or figures for the accompanying manuscript.

## Advanced usage

### Custom panel definitions

Create a small Python script that imports `run_panel_analysis` and defines custom `PanelConfig` objects for arbitrary population mixes. See `scripts/run_analysis.py` for an example wrapper.

### Reproducible environments

Use the provided `pyproject.toml` with `pip install -e .` or export requirements via:

```bash
pip install pip-tools
pip-compile pyproject.toml
```

Then store the resulting `requirements.txt` in the repository if needed.
