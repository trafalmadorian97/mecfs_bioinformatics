"""
Column names and status values of the trait x protein colocalization outputs.

The dense per-locus window frames use gwaslab column names (CHR, POS, EA, NEA, BETA, SE, EAF,
N); the constants below name everything else.
"""

from typing import Final

import polars as pl

from mecfs_bio.constants.gwaslab_constants import (
    GWASLAB_BETA_COL,
    GWASLAB_CHROM_COL,
    GWASLAB_EFFECT_ALLELE_COL,
    GWASLAB_EFFECT_ALLELE_FREQ_COL,
    GWASLAB_NON_EFFECT_ALLELE_COL,
    GWASLAB_POS_COL,
    GWASLAB_SAMPLE_SIZE_COLUMN,
    GWASLAB_SE_COL,
)

LOCUS_ID_COL: Final = "locus_id"
PROTEIN_ID_COL: Final = "protein_id"

# Schema of a dense summary-statistics window (trait or protein), positions in GRCh38.
WINDOW_SCHEMA: Final = pl.Schema(
    {
        GWASLAB_CHROM_COL: pl.Int64,
        GWASLAB_POS_COL: pl.Int64,
        GWASLAB_EFFECT_ALLELE_COL: pl.String,
        GWASLAB_NON_EFFECT_ALLELE_COL: pl.String,
        GWASLAB_BETA_COL: pl.Float64,
        GWASLAB_SE_COL: pl.Float64,
        GWASLAB_EFFECT_ALLELE_FREQ_COL: pl.Float64,
        GWASLAB_SAMPLE_SIZE_COLUMN: pl.Float64,
    }
)
