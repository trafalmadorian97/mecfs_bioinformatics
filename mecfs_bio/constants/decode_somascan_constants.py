"""
Constants for the deCODE Iceland SomaScan proteomics GWAS (Eldjarn et al. 2023,
doi 10.1038/s41586-023-06563-x) and the database built from it.

Two groups of per-protein files exist: "smp" (SMP-normalized SomaScan levels) and "raw"
(non-normalized, PC0). Each file has about 33.4M rows, GRCh38, with chromosomes in text order
(chr1, chr10, ..., chr9, chrX). Files mostly share one row list and order, but some lack a few
rows, so rows are matched across files by Name, never by row number.
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

# Variant index column beyond the gwaslab ones: deCODE's variant Name
# (chrom:pos:effectAllele:otherAllele), by which protein files are matched to the index.
DECODE_INDEX_NAME_COL = "NAME"

# Per-protein heritability table columns beyond the PPP ones: the protein id (seq id), and the
# 5th and 95th percentiles of the per-variant N over the regression variants.
DECODE_H2_PROTEIN_ID_COL = "protein_id"
DECODE_H2_N_P05_COL = "n_p05"
DECODE_H2_N_P95_COL = "n_p95"

# Columns of the typed Eldjarn et al. 2023 ST36 table (highly pleiotropic pQTLs; see
# St36PleiotropicPqtlsPipe): chromosome code, positions (hg38) and pleiotropy flags.
ST36_CHROM_COL = "chrom"
ST36_OLINK_POS_COL = "olink_pos"
ST36_PLEIOTROPIC_OLINK_COL = "highly_pleiotropic_olink"
ST36_PLEIOTROPIC_SMP_COL = "highly_pleiotropic_soma_smp"
ST36_PLEIOTROPIC_PC0_COL = "highly_pleiotropic_soma_pc0"
ST36_SMP_POS_COL = "soma_smp_pos"
ST36_SMP_NAME_COL = "soma_smp_name"
ST36_SMP_RSNAME_COL = "soma_smp_rsname"
ST36_SMP_NPROT_COL = "soma_smp_nprot"
ST36_PC0_POS_COL = "soma_pc0_pos"
