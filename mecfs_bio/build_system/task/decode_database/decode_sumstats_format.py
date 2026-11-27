"""
Reading the text of deCODE SomaScan per-protein sumstats files.

The files are tab-separated with the header DECODE_HEADER and "." for missing values.
Chromosomes are written chr1..chr22, chrX and coded 1..22, 23.
"""

import io

import polars as pl

from mecfs_bio.constants.decode_somascan_constants import (
    DECODE_BETA_COL,
    DECODE_CHR_X_CODE,
    DECODE_CHROM_COL,
    DECODE_EFFECT_ALLELE_COL,
    DECODE_HEADER,
    DECODE_IMP_MAF_COL,
    DECODE_MLOG10P_COL,
    DECODE_N_COL,
    DECODE_NAME_COL,
    DECODE_OTHER_ALLELE_COL,
    DECODE_POS_COL,
    DECODE_PVAL_COL,
    DECODE_RSIDS_COL,
    DECODE_SE_COL,
)

DECODE_SCHEMA = pl.Schema(
    {
        DECODE_CHROM_COL: pl.String,
        DECODE_POS_COL: pl.Int64,
        DECODE_NAME_COL: pl.String,
        DECODE_RSIDS_COL: pl.String,
        DECODE_EFFECT_ALLELE_COL: pl.String,
        DECODE_OTHER_ALLELE_COL: pl.String,
        DECODE_BETA_COL: pl.Float64,
        DECODE_PVAL_COL: pl.Float64,
        DECODE_MLOG10P_COL: pl.Float64,
        DECODE_SE_COL: pl.Float64,
        DECODE_N_COL: pl.Float64,
        DECODE_IMP_MAF_COL: pl.Float64,
    }
)
assert tuple(DECODE_SCHEMA.names()) == DECODE_HEADER

_NULL_VALUE = "."
_CHROM_CODES: dict[str, int] = {f"chr{i}": i for i in range(1, 23)} | {
    "chrX": DECODE_CHR_X_CODE
}


def chrom_code(chrom: str) -> int:
    """chr1..chr22 -> 1..22, chrX -> 23."""
    assert chrom in _CHROM_CODES, f"unexpected deCODE chromosome {chrom!r}"
    return _CHROM_CODES[chrom]


def chrom_code_expr(col: str) -> pl.Expr:
    """chrom_code over a string column, as Int8; raises on any other value."""
    return pl.col(col).replace_strict(_CHROM_CODES, return_dtype=pl.Int8)


def check_header(line: bytes) -> None:
    """Assert that line (without its newline) is the deCODE header."""
    fields = tuple(line.rstrip(b"\r\n").decode().split("\t"))
    assert fields == DECODE_HEADER, f"unexpected deCODE header: {fields}"


def parse_decode_lines(text: bytes) -> pl.DataFrame:
    """Parse header-less deCODE lines into DECODE_SCHEMA."""
    if not text:
        return DECODE_SCHEMA.to_frame()
    return pl.read_csv(
        io.BytesIO(text),
        separator="\t",
        has_header=False,
        schema=DECODE_SCHEMA,
        null_values=[_NULL_VALUE],
        quote_char=None,
    )
