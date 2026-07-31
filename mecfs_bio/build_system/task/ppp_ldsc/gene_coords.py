"""
Gene coordinates per Olink protein, read from the Sun et al. 2023 supplementary table ST3.

The cis-excluded variant set needs to know where each protein's gene sits, so that variants near
it can be dropped and the remaining heritability or genetic correlation reflects trans effects
only. ST3 supplies those hg38 coordinates.

This module takes ST3 already extracted into a dataframe (ExtractSheetFromExelFileTask) rather
than reading the spreadsheet itself, which keeps the spreadsheet's quirks -- the two title rows
above the header, the mixed numeric/text columns -- at the point where ST3 is wired in, and leaves
the parsing here as a plain dataframe-to-dictionary function.
"""

from __future__ import annotations

import polars as pl
import structlog

from mecfs_bio.constants.ppp_database_constants import Oid
from mecfs_bio.constants.vocabulary_classes.genomic_interval import GenomicInterval

logger = structlog.get_logger()

# ST3 column labels, as they appear in the extracted sheet.
ST3_OLINK_ID_COL = "Olink ID"
ST3_GENE_CHROM_COL = "Gene CHROM"
ST3_GENE_START_COL = "Gene start"
ST3_GENE_END_COL = "Gene end"

_X_CHROMOSOME = 23


def parse_chrom(value: str | None) -> int | None:
    """Chromosome label -> int (X -> 23); None for unparseable/non-standard labels.

    NOTE:
        - There are a few "proteins" in PPP that are actually multi-protein complexes.
        - These multi-protein complexes do not correspond to a unique gene on a single chromosome.
        - In the Excel sheet, these proteins have entries in the chromosome column like chrom1;chrom2
        - We return None for these proteins.

    """
    if value is None:
        return None
    text = str(value).strip()
    if text.isdigit():
        return int(text)
    if text.upper() == "X":
        return _X_CHROMOSOME
    return None


def _parse_position(value: str | None) -> int | None:
    """A gene start/end -> int; None when absent or non-numeric (the multi-gene complexes
    described in parse_chrom carry several positions joined by semicolons)."""
    if value is None:
        return None
    text = str(value).strip()
    return int(text) if text.isdigit() else None


def read_gene_coords(st3: pl.DataFrame) -> dict[Oid, GenomicInterval]:
    """Map Olink ID (OID) -> hg38 gene coordinates from the extracted ST3 sheet.

    Proteins whose chromosome or coordinates are non-standard are skipped: they simply have no
    cis window, and callers treat them accordingly."""
    for column in (
        ST3_OLINK_ID_COL,
        ST3_GENE_CHROM_COL,
        ST3_GENE_START_COL,
        ST3_GENE_END_COL,
    ):
        assert column in st3.columns, (
            f"extracted ST3 sheet has no '{column}' column; found {st3.columns}"
        )

    coords: dict[Oid, GenomicInterval] = {}
    for row in st3.select(
        pl.col(ST3_OLINK_ID_COL).cast(pl.String),
        pl.col(ST3_GENE_CHROM_COL).cast(pl.String),
        pl.col(ST3_GENE_START_COL).cast(pl.String),
        pl.col(ST3_GENE_END_COL).cast(pl.String),
    ).iter_rows(named=True):
        oid = row[ST3_OLINK_ID_COL]
        chrom = parse_chrom(row[ST3_GENE_CHROM_COL])
        start = _parse_position(row[ST3_GENE_START_COL])
        end = _parse_position(row[ST3_GENE_END_COL])
        if oid is None or chrom is None or start is None or end is None:
            continue
        coords[Oid(oid)] = GenomicInterval(chrom=chrom, start=start, end=end)

    logger.debug(f"Read gene coordinates for {len(coords)} of {st3.height} ST3 rows")
    return coords
