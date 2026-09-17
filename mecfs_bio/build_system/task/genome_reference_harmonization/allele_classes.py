"""
Per-row allele classes relative to the genome reference.

Equal-length variants (SNVs, MNPs) are compared on the plus strand. The reverse
complement is consulted only when neither plus-strand allele matches, as gwaslab does.
Different-length variants (indels) are compared on the plus strand only: a
reverse-complemented left-anchored indel loses its anchor base, so a minus-strand match
would be spurious.

AlleleClass is what the type checker sees; the CLASS_* constants are what code uses. An
import-time assertion keeps the two in step (see outcomes.py for why this is not a
StrEnum).
"""

from typing import Final, Literal, get_args

import numpy as np
import polars as pl

from mecfs_bio.build_system.task.genome_reference_harmonization.fasta import (
    IndexedFasta,
    reference_matches,
)
from mecfs_bio.constants.gwaslab_constants import (
    GWASLAB_EFFECT_ALLELE_COL,
    GWASLAB_NON_EFFECT_ALLELE_COL,
    GWASLAB_POS_COL,
)

AlleleClass = Literal[
    "nea_ref",
    "ea_ref",
    "nea_ref_rc",
    "ea_ref_rc",
    "not_on_reference",
    "indel_both",
    "indel_nea_only",
    "indel_ea_only",
    "indel_not_on_reference",
]
CLASS_NEA_REF: Final[AlleleClass] = "nea_ref"
CLASS_EA_REF: Final[AlleleClass] = "ea_ref"
CLASS_NEA_REF_RC: Final[AlleleClass] = "nea_ref_rc"
CLASS_EA_REF_RC: Final[AlleleClass] = "ea_ref_rc"
CLASS_NOT_ON_REFERENCE: Final[AlleleClass] = "not_on_reference"
CLASS_INDEL_BOTH: Final[AlleleClass] = "indel_both"
CLASS_INDEL_NEA_ONLY: Final[AlleleClass] = "indel_nea_only"
CLASS_INDEL_EA_ONLY: Final[AlleleClass] = "indel_ea_only"
CLASS_INDEL_NOT_ON_REFERENCE: Final[AlleleClass] = "indel_not_on_reference"
assert set(get_args(AlleleClass)) == {
    CLASS_NEA_REF,
    CLASS_EA_REF,
    CLASS_NEA_REF_RC,
    CLASS_EA_REF_RC,
    CLASS_NOT_ON_REFERENCE,
    CLASS_INDEL_BOTH,
    CLASS_INDEL_NEA_ONLY,
    CLASS_INDEL_EA_ONLY,
    CLASS_INDEL_NOT_ON_REFERENCE,
}

ALLELE_CLASS_COL = "_allele_class"
IS_PALINDROMIC_SNV_COL = "_is_palindromic_snv"
IS_PALINDROMIC_MNP_COL = "_is_palindromic_mnp"
_VALID_ALLELE_PATTERN = "^[ACGT]+$"
_NEA_RC_COL = "_nea_rc"
_EA_RC_COL = "_ea_rc"
_BASES = ["A", "C", "G", "T"]
_COMPLEMENTS = ["T", "G", "C", "A"]


def reverse_complement_expr(column: str) -> pl.Expr:
    return pl.col(column).str.reverse().str.replace_many(_BASES, _COMPLEMENTS)


def prepare_alleles(frame: pl.DataFrame) -> pl.DataFrame:
    """Alleles as uppercase strings and positions as Int64."""
    return frame.with_columns(
        pl.col(GWASLAB_EFFECT_ALLELE_COL).cast(pl.String).str.to_uppercase(),
        pl.col(GWASLAB_NON_EFFECT_ALLELE_COL).cast(pl.String).str.to_uppercase(),
        pl.col(GWASLAB_POS_COL).cast(pl.Int64),
    )


def valid_alleles_expr() -> pl.Expr:
    ea = pl.col(GWASLAB_EFFECT_ALLELE_COL)
    nea = pl.col(GWASLAB_NON_EFFECT_ALLELE_COL)
    return (
        ea.str.contains(_VALID_ALLELE_PATTERN)
        & nea.str.contains(_VALID_ALLELE_PATTERN)
        & (ea != nea)
    )


def classify_alleles(
    frame: pl.DataFrame, fasta: IndexedFasta, chrom: int, max_gather_bytes: int
) -> pl.DataFrame:
    """Add ALLELE_CLASS_COL and the palindrome flags to rows with valid, uppercase alleles."""
    pos = frame[GWASLAB_POS_COL].to_numpy()
    nea_match = reference_matches(
        fasta,
        chrom=chrom,
        positions=pos,
        alleles=frame[GWASLAB_NON_EFFECT_ALLELE_COL],
        max_gather_bytes=max_gather_bytes,
    )
    ea_match = reference_matches(
        fasta,
        chrom=chrom,
        positions=pos,
        alleles=frame[GWASLAB_EFFECT_ALLELE_COL],
        max_gather_bytes=max_gather_bytes,
    )
    equal_length = (
        frame[GWASLAB_EFFECT_ALLELE_COL].str.len_bytes()
        == frame[GWASLAB_NON_EFFECT_ALLELE_COL].str.len_bytes()
    ).to_numpy()
    needs_rc = equal_length & ~nea_match & ~ea_match
    rc_nea_match = np.zeros(frame.height, dtype=bool)
    rc_ea_match = np.zeros(frame.height, dtype=bool)
    rc_rows = np.flatnonzero(needs_rc)
    if len(rc_rows):
        rc = frame.select(
            reverse_complement_expr(GWASLAB_NON_EFFECT_ALLELE_COL).alias(_NEA_RC_COL),
            reverse_complement_expr(GWASLAB_EFFECT_ALLELE_COL).alias(_EA_RC_COL),
        )
        rc_nea_match[rc_rows] = reference_matches(
            fasta,
            chrom=chrom,
            positions=pos[rc_rows],
            alleles=rc[_NEA_RC_COL].gather(rc_rows),
            max_gather_bytes=max_gather_bytes,
        )
        rc_ea_match[rc_rows] = reference_matches(
            fasta,
            chrom=chrom,
            positions=pos[rc_rows],
            alleles=rc[_EA_RC_COL].gather(rc_rows),
            max_gather_bytes=max_gather_bytes,
        )
    classes = np.select(
        [
            equal_length & nea_match,
            equal_length & ea_match,
            needs_rc & rc_nea_match,
            needs_rc & rc_ea_match,
            equal_length,
            nea_match & ea_match,
            nea_match,
            ea_match,
        ],
        [
            CLASS_NEA_REF,
            CLASS_EA_REF,
            CLASS_NEA_REF_RC,
            CLASS_EA_REF_RC,
            CLASS_NOT_ON_REFERENCE,
            CLASS_INDEL_BOTH,
            CLASS_INDEL_NEA_ONLY,
            CLASS_INDEL_EA_ONLY,
        ],
        default=CLASS_INDEL_NOT_ON_REFERENCE,
    )
    ea_length = pl.col(GWASLAB_EFFECT_ALLELE_COL).str.len_bytes()
    same_length = ea_length == pl.col(GWASLAB_NON_EFFECT_ALLELE_COL).str.len_bytes()
    palindromic = same_length & (
        reverse_complement_expr(GWASLAB_NON_EFFECT_ALLELE_COL)
        == pl.col(GWASLAB_EFFECT_ALLELE_COL)
    )
    return frame.with_columns(
        pl.Series(ALLELE_CLASS_COL, classes, dtype=pl.String),
        (palindromic & (ea_length == 1)).alias(IS_PALINDROMIC_SNV_COL),
        (palindromic & (ea_length > 1)).alias(IS_PALINDROMIC_MNP_COL),
    )
