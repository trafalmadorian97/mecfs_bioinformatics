"""Unit tests for reading protein gene coordinates from the extracted Sun et al. ST3 sheet."""

import polars as pl
import pytest

from mecfs_bio.build_system.task.ppp_ldsc.gene_coords import (
    ST3_GENE_CHROM_COL,
    ST3_GENE_END_COL,
    ST3_GENE_START_COL,
    ST3_OLINK_ID_COL,
    parse_chrom,
    read_gene_coords,
)
from mecfs_bio.constants.ppp_database_constants import Oid


def _st3(rows: list[tuple[str, str, str, str]]) -> pl.DataFrame:
    return pl.DataFrame(
        {
            ST3_OLINK_ID_COL: [row[0] for row in rows],
            ST3_GENE_CHROM_COL: [row[1] for row in rows],
            ST3_GENE_START_COL: [row[2] for row in rows],
            ST3_GENE_END_COL: [row[3] for row in rows],
        }
    )


@pytest.mark.parametrize(
    ("label", "expected"),
    [
        ("19", 19),
        (" 7 ", 7),
        ("X", 23),
        ("x", 23),
        ("1;1;1", None),  # a multi-gene protein complex has no single locus
        ("Y", None),
        ("", None),
        (None, None),
    ],
)
def test_parse_chrom(label: str | None, expected: int | None):
    assert parse_chrom(label) == expected


def test_read_gene_coords():
    coords = read_gene_coords(
        _st3(
            [
                ("OID30771", "19", "58345178", "58353492"),
                ("OID20921", "X", "51968510", "51970000"),
            ]
        )
    )
    assert set(coords) == {Oid("OID30771"), Oid("OID20921")}
    assert coords[Oid("OID30771")].chrom == 19
    assert coords[Oid("OID30771")].start == 58345178
    assert coords[Oid("OID30771")].end == 58353492
    assert coords[Oid("OID20921")].chrom == 23


def test_read_gene_coords_skips_proteins_without_a_single_locus():
    """Multi-protein complexes carry semicolon-joined values in every coordinate column; they
    get no cis window rather than a nonsensical one."""
    coords = read_gene_coords(
        _st3(
            [
                ("OID30707", "1;1;1", "103655760;103687415", "103664554;103696454"),
                ("OID30771", "19", "58345178", "58353492"),
            ]
        )
    )
    assert set(coords) == {Oid("OID30771")}


def test_read_gene_coords_requires_expected_columns():
    with pytest.raises(AssertionError):
        read_gene_coords(pl.DataFrame({ST3_OLINK_ID_COL: ["OID30771"]}))
