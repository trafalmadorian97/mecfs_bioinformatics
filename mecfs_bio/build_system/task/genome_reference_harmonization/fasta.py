"""
Memory-mapped access to an uncompressed, faidx-indexed genome FASTA.

The uncompressed file is memory-mapped, and (chromosome, position) is translated to a
byte offset with the .fai index, so nothing is loaded up front. The operating system
page cache serves the touched pages, and those pages are reclaimable rather than
anonymous memory. Alleles are compared as prefixes of the reference starting at the
variant position, ignoring case (soft-masked bases are lowercase).
"""

import re
from collections.abc import Mapping
from pathlib import Path

import numpy as np
import polars as pl
from attrs import frozen

from mecfs_bio.constants.gwaslab_constants import GWASLAB_CHROM_CODE_FOR_NAME

FASTA_FILENAME = "genome.fa"
FAI_SUFFIX = ".fai"
# Upper bound on the byte size of one gather's offset matrix, so that very long alleles
# are processed in row slices instead of materializing a huge matrix.
DEFAULT_MAX_GATHER_BYTES = 64 * 1024 * 1024
_OFFSET_BYTES = 8
_ASCII_UPPERCASE_MASK = 0xDF
_MAIN_CONTIG = re.compile(r"^(?:chr)?([0-9]+|X|Y|M|MT)$")


@frozen
class FaiEntry:
    length: int # Number of actual nucleotide bases
    offset: int # Location in fasta file where first base starts
    line_bases: int # Each normal line contains this many bases
    line_bytes: int # Each normal line contains this many bytes


def contig_to_gwaslab_code(name: str) -> int | None:
    """gwaslab numeric chromosome code for a main contig name (chr1, 1, chrX, chrM, MT), else None."""
    match = _MAIN_CONTIG.match(name)
    if match is None:
        return None
    token = match.group(1)
    if token == "M":
        token = "MT"
    if token in GWASLAB_CHROM_CODE_FOR_NAME:
        return GWASLAB_CHROM_CODE_FOR_NAME[token]
    return int(token)


@frozen
class IndexedFasta:
    """An uncompressed FASTA and its .fai entries, keyed by gwaslab chromosome code."""

    fasta_path: Path
    entries: Mapping[int, FaiEntry]

    @classmethod
    def open(cls, directory: Path) -> "IndexedFasta":
        fasta_path = directory / FASTA_FILENAME
        fai_path = directory / (FASTA_FILENAME + FAI_SUFFIX)
        assert fasta_path.is_file(), f"missing {fasta_path}"
        assert fai_path.is_file(), f"missing {fai_path}"
        entries: dict[int, FaiEntry] = {}
        for line in fai_path.read_text().splitlines():
            name, length, offset, line_bases, line_bytes = line.split("\t")[:5]
            code = contig_to_gwaslab_code(name)
            if code is None:
                continue
            assert code not in entries, f"two FASTA contigs map to chromosome {code}"
            entries[code] = FaiEntry(
                length=int(length),
                offset=int(offset),
                line_bases=int(line_bases),
                line_bytes=int(line_bytes),
            )
        return cls(fasta_path=fasta_path, entries=entries)


def reference_matches(
    fasta: IndexedFasta,
    chrom: int,
    positions: np.ndarray,
    alleles: pl.Series,
    max_gather_bytes: int = DEFAULT_MAX_GATHER_BYTES,
) -> np.ndarray:
    """Boolean array: allele i equals the reference sequence starting at 1-based positions[i]."""
    assert positions.ndim == 1 and len(positions) == len(alleles), (
        "positions and alleles must be one-dimensional and the same length"
    )
    assert chrom in fasta.entries, f"chromosome {chrom} is not in {fasta.fasta_path}"
    assert max_gather_bytes > 0
    result = np.zeros(len(positions), dtype=bool)
    if len(positions) == 0:
        return result
    entry = fasta.entries[chrom]
    genome = np.memmap(fasta.fasta_path, dtype=np.uint8, mode="r")
    lengths = alleles.str.len_bytes().to_numpy()
    assert (lengths > 0).all(), "alleles must be non-empty"
    order = np.argsort(lengths, kind="stable")
    boundaries = np.flatnonzero(np.diff(lengths[order])) + 1
    for group in np.split(order, boundaries):
        length = int(lengths[group[0]])
        rows_per_slice = max(1, max_gather_bytes // (length * _OFFSET_BYTES))
        for start in range(0, len(group), rows_per_slice):
            rows = group[start : start + rows_per_slice]
            result[rows] = _match_equal_length(
                genome=genome,
                entry=entry,
                positions=positions[rows],
                alleles=alleles.gather(rows),
                length=length,
            )
    return result


def _match_equal_length(
    genome: np.ndarray,
    entry: FaiEntry,
    positions: np.ndarray,
    alleles: pl.Series,
    length: int,
) -> np.ndarray:
    start = positions.astype(np.int64) - 1
    in_bounds = (start >= 0) & (start + length <= entry.length)
    clipped = np.clip(start, 0, max(entry.length - length, 0))
    base_index = clipped[:, None] + np.arange(length, dtype=np.int64)[None, :]
    file_offsets = (
        entry.offset
        + base_index
        + (base_index // entry.line_bases) * (entry.line_bytes - entry.line_bases)
    )
    reference = genome[file_offsets] & _ASCII_UPPERCASE_MASK
    observed = (
        np.frombuffer("".join(alleles.to_list()).encode("ascii"), dtype=np.uint8).reshape(
            len(positions), length
        )
        & _ASCII_UPPERCASE_MASK
    )
    return in_bounds & (reference == observed).all(axis=1)
