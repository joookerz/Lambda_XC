"""Marker pruning helpers."""

from __future__ import annotations

from typing import List

import numpy as np

from .data import HapMapDataset


def distance_prune(
    dataset: HapMapDataset,
    nsnp_per_chr: int,
    prune_kb: int,
) -> np.ndarray:
    """Return SNP indices after simple distance-based pruning."""
    keep: List[int] = []
    snp_table = dataset.snp_table
    chromosomes = sorted(snp_table["chromosome"].unique())

    for chrom in chromosomes:
        start, stop = dataset.chromosome_range(int(chrom))
        rows = np.arange(start, stop, dtype=int)
        pos = dataset.chromosome_positions(int(chrom))

        last = -np.inf
        picked = 0
        for row, position in zip(rows, pos):
            if position - last >= prune_kb * 1_000:
                keep.append(int(row))
                last = position
                picked += 1
                if picked == nsnp_per_chr:
                    break

    keep = np.array(sorted(set(keep)), dtype=int)
    if keep.size == 0:
        raise RuntimeError("Distance pruning produced no SNPs; check parameters.")
    return keep
