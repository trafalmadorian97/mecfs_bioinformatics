"""
Constants for the deCODE Iceland SomaScan proteomics GWAS (Eldjarn et al. 2023,
doi 10.1038/s41586-023-06563-x) and the database built from it.

Two groups of per-protein files exist: "smp" (SMP-normalized SomaScan levels) and "raw"
(non-normalized, PC0). Every file has the same 33.4M rows in the same order (chromosomes in
text order: chr1, chr10, ..., chr9, chrX), GRCh38.
"""

import math
from typing import Literal

DecodeSomascanGroup = Literal["smp", "raw"]

# Columns of a deCODE per-protein sumstats file, in file order.
DECODE_CHROM_COL = "Chrom"
DECODE_POS_COL = "Pos"
DECODE_NAME_COL = "Name"
DECODE_RSIDS_COL = "rsids"
DECODE_EFFECT_ALLELE_COL = "effectAllele"
DECODE_OTHER_ALLELE_COL = "otherAllele"
DECODE_BETA_COL = "Beta"
DECODE_PVAL_COL = "Pval"
DECODE_MLOG10P_COL = "minus_log10_pval"
DECODE_SE_COL = "SE"
DECODE_N_COL = "N"
DECODE_IMP_MAF_COL = "ImpMAF"
DECODE_HEADER: tuple[str, ...] = (
    DECODE_CHROM_COL,
    DECODE_POS_COL,
    DECODE_NAME_COL,
    DECODE_RSIDS_COL,
    DECODE_EFFECT_ALLELE_COL,
    DECODE_OTHER_ALLELE_COL,
    DECODE_BETA_COL,
    DECODE_PVAL_COL,
    DECODE_MLOG10P_COL,
    DECODE_SE_COL,
    DECODE_N_COL,
    DECODE_IMP_MAF_COL,
)
# otherAllele value deCODE uses for a variant with no usable second allele.
DECODE_MISSING_ALLELE = "!"
# Chromosome code for chrX (chr1..chr22 map to 1..22).
DECODE_CHR_X_CODE = 23

# Parquet key-value metadata key of the variant index: data rows in every protein file.
DECODE_TOTAL_ROWS_METADATA_KEY = "total_rows"

# Columns of the deduplicated variant_info_iceland table (allele frequencies per allele).
VI_CHR_COL = "CHR"
VI_POS_COL = "POS"
VI_A1_COL = "A1"
VI_FREQ_COL = "ImpFreqA1"
VI_INFO_COL = "ImpInfo"

# Two frequencies agree when they differ by at most this (deCODE prints 6 decimals).
EAF_TOLERANCE = 1e-5

# The per-protein task keeps every row with p <= 5e-8 (any chromosome) in a supplement file.
DECODE_SIGNIFICANT_ROWS_MIN_MLOG10P = -math.log10(5e-8)

# Variant index columns beyond the gwaslab ones: the row's data row number in every protein
# file (0 = first line after the header), and deCODE's variant Name, checked on every gather.
DECODE_INDEX_SOURCE_ROW_COL = "SOURCE_ROW"
DECODE_INDEX_NAME_COL = "NAME"
