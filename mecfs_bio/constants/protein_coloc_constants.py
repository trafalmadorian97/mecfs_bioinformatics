"""
Column names and status values of the trait x protein colocalization outputs.

The dense per-locus window frames use gwaslab column names (CHR, POS, EA, NEA, BETA, SE, EAF,
N); the constants below name everything else.
"""

from typing import Final, Literal

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

LOCUS_ID_COL: Final = "LOCUS_ID"
PROTEIN_ID_COL: Final = "PROTEIN_ID"
GENE_COL: Final = "GENE"

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
HIT_PROTEIN_ID_COL: Final = PROTEIN_ID_COL
HIT_GENE_COL: Final = GENE_COL
HIT_CHR_COL: Final = "CHR"
HIT_SEG_START_COL: Final = "SEG_START"
HIT_SEG_END_COL: Final = "SEG_END"
HIT_LEAD_POS_COL: Final = "LEAD_POS"
HIT_LEAD_MLOG10P_COL: Final = "LEAD_MLOG10P"
HIT_N_SIG_VARIANTS_COL: Final = "N_SIG_VARIANTS"

# Trait loci: one row per (merged) window around a genome-wide significant trait lead variant.
LOCUS_CHR_COL: Final = "CHR"
LOCUS_START_COL: Final = "START"
LOCUS_END_COL: Final = "END"
LOCUS_LEAD_POS_COL: Final = "LEAD_POS"
LOCUS_LEAD_EA_COL: Final = "LEAD_EA"
LOCUS_LEAD_NEA_COL: Final = "LEAD_NEA"
LOCUS_LEAD_MLOG10P_COL: Final = "LEAD_MLOG10P"
STATUS_COL: Final = "STATUS"

# Candidate (locus, protein) pairs and their selection scores.
PROTEIN_MLOG10P_COL: Final = "PROTEIN_MLOG10P"
CIS_STATUS_COL: Final = "CIS_STATUS"
N_OTHER_PROTEINS_COL: Final = "N_OTHER_PROTEINS"
# Distance from the trait lead to the protein's gene (0 inside it); null when not on its chromosome.
GENE_LEAD_DISTANCE_COL: Final = "GENE_LEAD_DISTANCE_BP"
SCORE_COL: Final = "SCORE"
SCORE_TRAIT_COL: Final = "SCORE_TRAIT"
SCORE_PROTEIN_COL: Final = "SCORE_PROTEIN"
SCORE_PENALTY_COL: Final = "SCORE_PENALTY"

CisStatus = Literal["cis", "trans", "unknown"]
CIS_STATUS_CIS: Final[CisStatus] = "cis"
CIS_STATUS_TRANS: Final[CisStatus] = "trans"
CIS_STATUS_UNKNOWN: Final[CisStatus] = "unknown"

ColocStatus = Literal[
    "candidate",
    "excluded_mhc",
    "no_candidates",
    "selected",
    "not_selected_budget",
    "tested",
    "failed",
]
STATUS_CANDIDATE: Final[ColocStatus] = "candidate"
STATUS_EXCLUDED_MHC: Final[ColocStatus] = "excluded_mhc"
STATUS_NO_CANDIDATES: Final[ColocStatus] = "no_candidates"
STATUS_SELECTED: Final[ColocStatus] = "selected"
STATUS_NOT_SELECTED_BUDGET: Final[ColocStatus] = "not_selected_budget"
STATUS_TESTED: Final[ColocStatus] = "tested"
STATUS_FAILED: Final[ColocStatus] = "failed"

# Raw colocalization results: candidate columns plus these per tested pair.
N_TRAIT_VARIANTS_COL: Final = "N_TRAIT_VARIANTS"
N_PROTEIN_VARIANTS_COL: Final = "N_PROTEIN_VARIANTS"
N_ALIGNED_COL: Final = "N_ALIGNED"
N_DUPLICATE_KEYS_DROPPED_COL: Final = "N_DUPLICATE_KEYS_DROPPED"
N_DROPPED_PREFIX: Final = "N_DROPPED_"
PROTEIN_LEAD_POS_COL: Final = "PROTEIN_LEAD_POS"
PROTEIN_LEAD_MLOG10P_COL: Final = "PROTEIN_LEAD_MLOG10P"
PP_H0_COL: Final = "PP_H0"
PP_H1_COL: Final = "PP_H1"
PP_H2_COL: Final = "PP_H2"
PP_H3_COL: Final = "PP_H3"
PP_H4_COL: Final = "PP_H4"
PRIOR_P1_COL: Final = "P1"
PRIOR_P2_COL: Final = "P2"
PRIOR_P12_COL: Final = "P12"
METHOD_COL: Final = "METHOD"
FAILURE_REASON_COL: Final = "FAILURE_REASON"

# Files inside a colocalization result directory.
RAW_RESULTS_FILENAME: Final = "raw_results.parquet"
WINDOWS_DIRNAME: Final = "windows"


def pair_file_stem(locus_id: str, protein_id: str) -> str:
    """File stem for one (locus, protein) pair, free of characters Windows filesystems reject."""
    return f"{locus_id}__{protein_id}".replace(":", "_")


# Docs-facing display table and plot index.
DISPLAY_LOCUS_COL: Final = "LOCUS"
DISPLAY_LEAD_VARIANT_COL: Final = "LEAD_VARIANT"
PLOT_PNG_FILENAME_COL: Final = "PNG_FILENAME"
PLOT_SVG_FILENAME_COL: Final = "SVG_FILENAME"
