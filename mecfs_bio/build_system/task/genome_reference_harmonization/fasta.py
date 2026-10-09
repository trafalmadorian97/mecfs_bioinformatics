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

from mecfs_bio.constants.gwaslab_constants import (
    GWASLAB_CHROM_CODE_FOR_NAME,
    GWASLAB_CHROM_NAME_FOR_CODE,
)

FASTA_FILENAME = "genome.fa"
FAI_SUFFIX = ".fai"
# Upper bound on the byte size of one gather's offset matrix, so that very long alleles
# are processed in row slices instead of materializing a huge matrix.
DEFAULT_MAX_GATHER_BYTES = 64 * 1024 * 1024
_OFFSET_BYTES = 8
_ASCII_UPPERCASE_MASK = 0xDF
_ACGT_CODES = np.frombuffer(b"ACGT", dtype=np.uint8)
_MAIN_CONTIG = re.compile(r"^(?:chr)?([0-9]+|X|Y|M|MT)$")


@frozen(slots=True)
class FaiEntry:
    length: int  # Number of actual nucleotide bases
    offset: int  # Location in fasta file where first base starts
    line_bases: int  # Each normal line contains this many bases
    line_bytes: int  # Each normal line contains this many bytes


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


def gwaslab_code_to_contig_name(code: int) -> str:
    """Bare contig name (1-22, X, Y, MT) for a gwaslab numeric chromosome code."""
    return GWASLAB_CHROM_NAME_FOR_CODE.get(code, str(code))


@frozen(slots=True)
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


@frozen(slots=True)
class _LengthSlice:
    """Row indices whose spans share one length, few enough for one gather."""

    rows: np.ndarray
    length: int


@frozen(slots=True)
class _ReferenceBases:
    """Uppercased reference bases of equal-length spans, and which spans are in bounds."""

    in_bounds: np.ndarray
    bases: np.ndarray

    def __attrs_post_init__(self):
        assert self.in_bounds.ndim == 1 and self.in_bounds.dtype == np.bool_
        assert self.bases.ndim == 2 and self.bases.dtype == np.uint8
        assert self.bases.shape[0] == len(self.in_bounds)


def _length_slices(lengths: np.ndarray, max_gather_bytes: int) -> list[_LengthSlice]:
    """Group rows by span length, split so no gather's offset matrix exceeds the budget."""
    slices: list[_LengthSlice] = []
    order = np.argsort(lengths, kind="stable")
    boundaries = np.flatnonzero(np.diff(lengths[order])) + 1
    for group in np.split(order, boundaries):
        length = int(lengths[group[0]])
        rows_per_slice = max(1, max_gather_bytes // (length * _OFFSET_BYTES))
        for start in range(0, len(group), rows_per_slice):
            slices.append(
                _LengthSlice(rows=group[start : start + rows_per_slice], length=length)
            )
    return slices


def _gather_reference(
    genome: np.ndarray, entry: FaiEntry, positions: np.ndarray, length: int
) -> _ReferenceBases:
    start = positions.astype(np.int64) - 1
    in_bounds = (start >= 0) & (start + length <= entry.length)
    clipped = np.clip(start, 0, max(entry.length - length, 0))
    base_index = clipped[:, None] + np.arange(length, dtype=np.int64)[None, :]
    file_offsets = (
        entry.offset
        + base_index
        + (base_index // entry.line_bases) * (entry.line_bytes - entry.line_bases)
    )
    return _ReferenceBases(
        in_bounds=in_bounds, bases=genome[file_offsets] & _ASCII_UPPERCASE_MASK
    )


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
    for piece in _length_slices(lengths, max_gather_bytes):
        reference = _gather_reference(
            genome, entry, positions[piece.rows], piece.length
        )
        observed = (
            np.frombuffer(
                "".join(alleles.gather(piece.rows).to_list()).encode("ascii"),
                dtype=np.uint8,
            ).reshape(len(piece.rows), piece.length)
            & _ASCII_UPPERCASE_MASK
        )
        result[piece.rows] = reference.in_bounds & (reference.bases == observed).all(
            axis=1
        )
    return result


def reference_is_acgt(
    fasta: IndexedFasta,
    chrom: int,
    positions: np.ndarray,
    lengths: np.ndarray,
    max_gather_bytes: int = DEFAULT_MAX_GATHER_BYTES,
) -> np.ndarray:
    """Boolean array: the lengths[i] reference bases from 1-based positions[i] are in bounds
    and each is A, C, G or T in either case. False where the span holds N or another IUPAC
    ambiguity code."""
    assert (
        positions.ndim == 1 and lengths.ndim == 1 and len(positions) == len(lengths)
    ), "positions and lengths must be one-dimensional and the same length"
    assert chrom in fasta.entries, f"chromosome {chrom} is not in {fasta.fasta_path}"
    assert max_gather_bytes > 0
    result = np.zeros(len(positions), dtype=bool)
    if len(positions) == 0:
        return result
    assert (lengths > 0).all(), "spans must be non-empty"
    entry = fasta.entries[chrom]
    genome = np.memmap(fasta.fasta_path, dtype=np.uint8, mode="r")
    for piece in _length_slices(lengths, max_gather_bytes):
        reference = _gather_reference(
            genome, entry, positions[piece.rows], piece.length
        )
        result[piece.rows] = reference.in_bounds & np.isin(
            reference.bases, _ACGT_CODES
        ).all(axis=1)
    return result
