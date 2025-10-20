"""Data-loading utilities for the HapMap3 panels used in the simulations."""

from __future__ import annotations

from dataclasses import dataclass
from functools import lru_cache
from pathlib import Path
from typing import Dict, Iterable, Tuple

import numpy as np
import pandas as pd


MISSING = 9


def _read_snp_table(path: Path) -> pd.DataFrame:
    """Read a PLINK .snp file into a DataFrame."""
    return pd.read_table(
        path,
        sep=r"\s+",
        names=["rsid", "chromosome", "morgans", "position", "ref", "alt"],
        index_col=0,
    )


def _read_geno(path: Path) -> np.ma.MaskedArray:
    """Read a PLINK .geno file as a masked array."""
    return np.genfromtxt(
        path,
        dtype="uint8",
        delimiter=1,
        missing_values=MISSING,
        usemask=True,
    )


@dataclass(frozen=True)
class HapMapDataset:
    """Convenience wrapper around the HapMap genotype files."""

    data_root: Path
    snp_reference: str = "HapMap3.snp"

    def __post_init__(self) -> None:
        root = self.data_root.expanduser().resolve()
        object.__setattr__(self, "data_root", root)
        snp_path = root / self.snp_reference
        if not snp_path.exists():
            raise FileNotFoundError(
                f"SNP reference file '{self.snp_reference}' not found under {root}"
            )
        snp_table = _read_snp_table(snp_path)
        object.__setattr__(self, "_snp_table", snp_table)

        chrom_ranges: Dict[int, Tuple[int, int]] = {}
        chrom_positions: Dict[int, np.ndarray] = {}
        chromosomes = snp_table["chromosome"].to_numpy()
        for chrom in np.unique(chromosomes):
            idx = np.where(chromosomes == chrom)[0]
            chrom_ranges[int(chrom)] = (int(idx[0]), int(idx[-1] + 1))
            chrom_positions[int(chrom)] = snp_table.iloc[idx]["position"].to_numpy()

        object.__setattr__(self, "_chrom_ranges", chrom_ranges)
        object.__setattr__(self, "_chrom_positions", chrom_positions)

    @property
    def snp_table(self) -> pd.DataFrame:
        return self._snp_table.copy()

    def chromosome_range(self, chromosome: int) -> Tuple[int, int]:
        try:
            return self._chrom_ranges[int(chromosome)]
        except KeyError as exc:
            raise ValueError(f"Chromosome {chromosome} not present in SNP table") from exc

    def chromosome_positions(self, chromosome: int) -> np.ndarray:
        try:
            return self._chrom_positions[int(chromosome)].copy()
        except KeyError as exc:
            raise ValueError(f"Chromosome {chromosome} not present in SNP table") from exc

    @lru_cache(maxsize=None)
    def load_population(self, population: str) -> np.ma.MaskedArray:
        geno_path = self.data_root / f"{population}.geno"
        if not geno_path.exists():
            raise FileNotFoundError(f"Genotype file not found for population '{population}'")
        return _read_geno(geno_path)

    def load_population_slice(
        self, population: str, chromosome: int
    ) -> np.ma.MaskedArray:
        start, stop = self.chromosome_range(chromosome)
        geno = self.load_population(population)
        return geno[start:stop]

    def load_population_subset(
        self, population: str, snp_indices: Iterable[int]
    ) -> np.ma.MaskedArray:
        geno = self.load_population(population)
        return geno[np.asarray(list(snp_indices))]
