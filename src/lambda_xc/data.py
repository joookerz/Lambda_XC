"""Data-loading utilities for the genotype panels used in the simulations."""

from __future__ import annotations

import gzip
from dataclasses import dataclass
from functools import lru_cache
from io import TextIOWrapper
from pathlib import Path
from typing import Dict, Iterable, List, Tuple

import numpy as np
import pandas as pd


MISSING = 9
SUPPORTED_FORMATS = {"hapmap", "plink", "vcf"}


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


def _read_bim(path: Path) -> pd.DataFrame:
    """Read a PLINK .bim table."""
    table = pd.read_table(
        path,
        sep=r"\s+",
        names=["chromosome", "rsid", "morgans", "position", "allele1", "allele2"],
    )
    table = table.set_index("rsid")
    return table


def _read_fam(path: Path) -> List[str]:
    """Read the sample IDs from a PLINK .fam file."""
    fam = pd.read_table(
        path,
        sep=r"\s+",
        header=None,
        names=["family", "sample", "father", "mother", "sex", "phenotype"],
    )
    return fam["sample"].tolist()


def _read_plink_bed(prefix: Path, n_snps: int, n_samples: int) -> np.ma.MaskedArray:
    """Read a PLINK binary bed file (SNP-major)."""

    bed_path = prefix.with_suffix(".bed")
    with bed_path.open("rb") as handle:
        header = handle.read(3)
        if len(header) != 3 or header[0:2] != b"\x6C\x1B" or header[2] != 0x01:
            raise ValueError(
                "Only SNP-major PLINK .bed files are supported. "
                "Ensure the file is in the modern binary format."
            )
        bytes_per_snp = (n_samples + 3) // 4
        data = handle.read()

    expected = bytes_per_snp * n_snps
    if len(data) < expected:
        raise ValueError(
            "Incomplete PLINK .bed file; expected at least "
            f"{expected} bytes but found {len(data)}."
        )

    geno = np.zeros((n_snps, n_samples), dtype=np.uint8)
    mask = np.zeros_like(geno, dtype=bool)

    for snp in range(n_snps):
        offset = snp * bytes_per_snp
        chunk = data[offset : offset + bytes_per_snp]
        sample_idx = 0
        for byte in chunk:
            for shift in range(0, 8, 2):
                if sample_idx >= n_samples:
                    break
                code = (byte >> shift) & 0b11
                if code == 0b00:
                    geno[snp, sample_idx] = 0
                elif code == 0b11:
                    geno[snp, sample_idx] = 2
                elif code == 0b10:
                    geno[snp, sample_idx] = 1
                else:  # 0b01 is missing
                    mask[snp, sample_idx] = True
                    geno[snp, sample_idx] = 0
                sample_idx += 1

    return np.ma.MaskedArray(geno, mask=mask)


def _open_text_file(path: Path) -> TextIOWrapper:
    """Return a text handle that transparently reads gzipped files."""

    if path.suffix == ".gz":
        return TextIOWrapper(gzip.open(path, "rb"))
    return TextIOWrapper(path.open("rb"))


def _parse_vcf_genotypes(path: Path) -> Tuple[np.ma.MaskedArray, pd.DataFrame]:
    """Parse genotypes and metadata from a VCF file."""

    genotypes: List[List[int]] = []
    masks: List[List[bool]] = []
    records: List[Tuple[str, int, str, str, str]] = []

    with _open_text_file(path) as handle:
        samples: List[str] = []
        for raw in handle:
            line = raw.strip()
            if not line:
                continue
            if line.startswith("##"):
                continue
            if line.startswith("#CHROM"):
                header = line.split("\t")
                samples = header[9:]
                continue

            fields = line.split("\t")
            chrom, pos, rsid, ref, alt = fields[:5]
            fmt = fields[8].split(":") if len(fields) > 8 else []
            try:
                gt_index = fmt.index("GT")
            except ValueError:
                gt_index = -1

            row: List[int] = []
            row_mask: List[bool] = []
            for entry in fields[9:]:
                if gt_index == -1:
                    row.append(0)
                    row_mask.append(True)
                    continue
                parts = entry.split(":")
                gt = parts[gt_index] if gt_index < len(parts) else "./."
                if gt in {".", "./.", ".|."}:
                    row.append(0)
                    row_mask.append(True)
                    continue
                if "|" in gt:
                    alleles = gt.split("|")
                else:
                    alleles = gt.split("/")
                valid = True
                count = 0
                for allele in alleles:
                    if allele == ".":
                        valid = False
                        break
                    try:
                        count += int(allele)
                    except ValueError:
                        valid = False
                        break
                if not valid:
                    row.append(0)
                    row_mask.append(True)
                else:
                    row.append(count)
                    row_mask.append(False)

            genotypes.append(row)
            masks.append(row_mask)

            if rsid == ".":
                rsid = f"{chrom}:{pos}"
            records.append((rsid, chrom, int(pos), ref, alt))

    if not genotypes:
        raise ValueError(f"VCF file '{path}' does not contain any variants.")

    geno_matrix = np.asarray(genotypes, dtype=np.uint8)
    mask_matrix = np.asarray(masks, dtype=bool)
    metadata = pd.DataFrame.from_records(
        records,
        columns=["rsid", "chromosome", "position", "ref", "alt"],
    ).set_index("rsid")
    metadata.insert(1, "morgans", np.nan)

    return np.ma.MaskedArray(geno_matrix, mask=mask_matrix), metadata


def _vcf_metadata(path: Path) -> pd.DataFrame:
    """Return the variant table without materialising the genotype matrix."""

    _, metadata = _parse_vcf_genotypes(path)
    return metadata


@dataclass(frozen=True)
class HapMapDataset:
    """Convenience wrapper around genotype files."""

    data_root: Path
    snp_reference: str | None = "HapMap3.snp"
    genotype_format: str = "auto"

    def __post_init__(self) -> None:
        root = self.data_root.expanduser().resolve()
        object.__setattr__(self, "data_root", root)

        fmt = self.genotype_format.lower()
        if fmt not in SUPPORTED_FORMATS | {"auto"}:
            raise ValueError(
                f"Unsupported genotype format '{self.genotype_format}'. "
                f"Supported formats: {', '.join(sorted(SUPPORTED_FORMATS))}."
            )
        if fmt == "auto":
            fmt = self._auto_detect_format()
        object.__setattr__(self, "_format", fmt)

        snp_table = self._load_variant_table(fmt)
        object.__setattr__(self, "_snp_table", snp_table)

        chrom_ranges: Dict[int, Tuple[int, int]] = {}
        chrom_positions: Dict[int, np.ndarray] = {}
        chromosomes = snp_table["chromosome"].to_numpy()
        try:
            chrom_values = chromosomes.astype(int)
        except ValueError as exc:
            raise ValueError(
                "Chromosome labels must be numeric for pruning operations."
            ) from exc
        for chrom in np.unique(chrom_values):
            idx = np.where(chrom_values == chrom)[0]
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
        fmt = self._format
        if fmt == "hapmap":
            path = self.data_root / f"{population}.geno"
            if not path.exists():
                raise FileNotFoundError(
                    f"Genotype file not found for population '{population}'"
                )
            geno = _read_geno(path)
            if geno.shape[0] != self.snp_table.shape[0]:
                raise ValueError(
                    "Genotype matrix does not match SNP reference dimensions."
                )
            return geno

        if fmt == "plink":
            prefix = self.data_root / population
            bim_path = prefix.with_suffix(".bim")
            fam_path = prefix.with_suffix(".fam")
            if not (prefix.with_suffix(".bed").exists() and bim_path.exists() and fam_path.exists()):
                raise FileNotFoundError(
                    f"PLINK files (.bed/.bim/.fam) not found for population '{population}'"
                )
            samples = _read_fam(fam_path)
            bim = _read_bim(bim_path)
            if bim.shape[0] != self.snp_table.shape[0] or not bim[["chromosome", "position"]].equals(
                self._snp_table[["chromosome", "position"]]
            ):
                raise ValueError(
                    "Variant order in PLINK files does not match the reference variant table."
                )
            geno = _read_plink_bed(prefix, self.snp_table.shape[0], len(samples))
            return geno

        if fmt == "vcf":
            candidates = [
                self.data_root / f"{population}.vcf",
                self.data_root / f"{population}.vcf.gz",
            ]
            path = next((p for p in candidates if p.exists()), None)
            if path is None:
                raise FileNotFoundError(
                    f"VCF file (.vcf or .vcf.gz) not found for population '{population}'"
                )
            geno, metadata = _parse_vcf_genotypes(path)
            if geno.shape[0] != self.snp_table.shape[0]:
                raise ValueError(
                    "VCF genotype matrix does not match reference variant count."
                )
            # Ensure metadata aligns for sanity.
            if not metadata[["chromosome", "position"]].equals(
                self._snp_table[["chromosome", "position"]]
            ):
                raise ValueError(
                    "Variant order in VCF does not match the reference variant table."
                )
            return geno

        raise RuntimeError(f"Unsupported genotype format '{fmt}'")

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

    # Internal helpers -------------------------------------------------

    def _auto_detect_format(self) -> str:
        geno_files = list(self.data_root.glob("*.geno"))
        if geno_files:
            return "hapmap"
        if list(self.data_root.glob("*.bed")) and list(self.data_root.glob("*.bim")):
            return "plink"
        if list(self.data_root.glob("*.vcf")) or list(self.data_root.glob("*.vcf.gz")):
            return "vcf"
        raise ValueError(
            "Unable to infer genotype format from data directory. "
            "Provide --genotype-format explicitly."
        )

    def _load_variant_table(self, fmt: str) -> pd.DataFrame:
        if fmt == "hapmap":
            if self.snp_reference is None:
                raise ValueError("A SNP reference file is required for HapMap data.")
            snp_path = self.data_root / self.snp_reference
            if not snp_path.exists():
                raise FileNotFoundError(
                    f"SNP reference file '{self.snp_reference}' not found under {self.data_root}"
                )
            return _read_snp_table(snp_path)

        if fmt == "plink":
            if self.snp_reference is not None:
                snp_path = self.data_root / self.snp_reference
                if snp_path.exists():
                    return _read_snp_table(snp_path)
            bim_files = sorted(self.data_root.glob("*.bim"))
            if not bim_files:
                raise FileNotFoundError(
                    "Could not locate a PLINK .bim file to use as variant reference."
                )
            return _read_bim(bim_files[0])

        if fmt == "vcf":
            if self.snp_reference is not None:
                candidate = self.data_root / self.snp_reference
                if candidate.exists():
                    # Accept either .snp-style tables or VCF files.
                    if candidate.suffix in {".vcf", ".gz"}:
                        return _vcf_metadata(candidate)
                    return _read_snp_table(candidate)
            vcf_files = sorted(
                list(self.data_root.glob("*.vcf"))
                + list(self.data_root.glob("*.vcf.gz"))
            )
            if not vcf_files:
                raise FileNotFoundError(
                    "Could not locate a VCF file to derive variant metadata."
                )
            return _vcf_metadata(vcf_files[0])

        raise RuntimeError(f"Unsupported genotype format '{fmt}'")
