"""
Column-name constants for the Western et al. 2024 CSF pQTL LDSC analysis: the
per-aptamer heritability table.

One row per aptamer. Only all-variants heritability is produced today, but the table
still carries a variant_set column (constant "all_variants") so a future cis-excluded set
can be added as extra rows without changing the schema. The aptamer is identified by its
SomaScan analyte id; uniprot and gene_symbol name the protein target for readability but
are not unique keys.
"""

from typing import Literal

CSF_H2_ANALYTE_COL = "analyte"
CSF_H2_UNIPROT_COL = "uniprot"
CSF_H2_GENE_SYMBOL_COL = "gene_symbol"
CSF_H2_VARIANT_SET_COL = "variant_set"
CSF_H2_H2_COL = "h2"
CSF_H2_H2_SE_COL = "h2_se"
CSF_H2_P_COL = "p"
CSF_H2_INTERCEPT_COL = "intercept"
CSF_H2_MEAN_CHI2_COL = "mean_chi2"
CSF_H2_LAMBDA_GC_COL = "lambda_gc"
CSF_H2_N_SNPS_COL = "n_snps"
CSF_H2_N_BAR_COL = "n_bar"

# Values of CSF_H2_VARIANT_SET_COL. Only the full regression SNP set exists today; the
# literal is kept open for a future cis-excluded set (matching the PPP convention's
# string values so the two databases' tables stay comparable).
CsfVariantSet = Literal["all_variants", "cis_excluded"]
CSF_VARIANT_SET_ALL: CsfVariantSet = "all_variants"
