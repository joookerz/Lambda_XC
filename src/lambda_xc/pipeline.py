"""High-level orchestration for the cross-chromosome inflation simulations."""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Iterable, List, Sequence

import numpy as np

from .analysis import PanelConfig, PanelResult, analyse_panel
from .data import HapMapDataset
from .pruning import distance_prune


@dataclass
class SimulationSettings:
    nsnp_per_chr: int = 2_000
    prune_kb: int = 100
    n_pairs: int = 1_000_000
    n_perm: int = 80
    n_boot: int = 100
    seed: int = 511


def run_panel_analysis(
    data_root: Path | str,
    panels: Sequence[PanelConfig],
    settings: SimulationSettings | None = None,
    genotype_format: str | None = None,
) -> List[PanelResult]:
    """Execute the full pipeline for the requested panels."""
    cfg = settings or SimulationSettings()
    dataset = HapMapDataset(
        Path(data_root),
        genotype_format=genotype_format or "auto",
    )

    snp_indices = distance_prune(dataset, cfg.nsnp_per_chr, cfg.prune_kb)
    chr_labels = dataset.snp_table.iloc[snp_indices]["chromosome"].to_numpy()
    rng = np.random.default_rng(cfg.seed)

    results: List[PanelResult] = []
    for panel in panels:
        result = analyse_panel(
            dataset=dataset,
            panel=panel,
            snp_indices=snp_indices,
            chr_labels=chr_labels,
            n_pairs=cfg.n_pairs,
            n_perm=cfg.n_perm,
            n_boot=cfg.n_boot,
            rng=rng,
        )
        results.append(result)
    return results
