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

# Columns of an aligned trait x protein window: the variant key (CHR, POS, EA, NEA) followed by
# each side's statistics.
BETA_TRAIT_COL: Final = "BETA_TRAIT"
SE_TRAIT_COL: Final = "SE_TRAIT"
EAF_TRAIT_COL: Final = "EAF_TRAIT"
N_TRAIT_COL: Final = "N_TRAIT"
BETA_PROTEIN_COL: Final = "BETA_PROTEIN"
SE_PROTEIN_COL: Final = "SE_PROTEIN"
EAF_PROTEIN_COL: Final = "EAF_PROTEIN"
N_PROTEIN_COL: Final = "N_PROTEIN"

# Protein hit-segment table: one row per (protein, chromosome, segment of significant variants).
HIT_PROTEIN_ID_COL: Final = "PROTEIN_ID"
HIT_GENE_COL: Final = "GENE"
HIT_CHR_COL: Final = "CHR"
HIT_SEG_START_COL: Final = "SEG_START"
HIT_SEG_END_COL: Final = "SEG_END"
HIT_LEAD_POS_COL: Final = "LEAD_POS"
HIT_LEAD_MLOG10P_COL: Final = "LEAD_MLOG10P"
HIT_N_SIG_VARIANTS_COL: Final = "N_SIG_VARIANTS"
