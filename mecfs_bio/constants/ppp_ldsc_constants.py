"""
Column-name constants for the UKB-PPP LDSC analysis tasks: the per-protein sample-size
table and the per-protein heritability table.

The extended-MHC exclusion region used by these tasks is NOT defined here; it reuses the
shared genomic constant (extended_mhc_interval / EXTENDED_MHC_BUILD_38 in
genomic_coordinate_constants), applied on the index's primary hg38 POS.
"""

from typing import Literal

# --- Sample-size table (one row per protein) ---
PPP_N_OID_COL = "oid"
PPP_N_GENE_COL = "gene"
PPP_N_SYNID_COL = "synapse_id"
PPP_N_SAMPLE_SIZE_COL = "n"

# --- Heritability table (two rows per protein: one per variant set) ---
PPP_H2_OID_COL = "oid"
PPP_H2_GENE_COL = "gene"
PPP_H2_VARIANT_SET_COL = "variant_set"
PPP_H2_H2_COL = "h2"
PPP_H2_H2_SE_COL = "h2_se"
PPP_H2_INTERCEPT_COL = "intercept"
PPP_H2_MEAN_CHI2_COL = "mean_chi2"
PPP_H2_LAMBDA_GC_COL = "lambda_gc"
PPP_H2_N_SNPS_COL = "n_snps"
PPP_H2_N_BAR_COL = "n_bar"

# Values of PPP_H2_VARIANT_SET_COL: the full regression SNP set, versus the set with
# variants within the cis window (+/- CIS_WINDOW_BP of the protein's gene) removed.
PppVariantSet = Literal["all_variants", "cis_excluded"]
PPP_VARIANT_SET_ALL: PppVariantSet = "all_variants"
PPP_VARIANT_SET_CIS_EXCLUDED: PppVariantSet = "cis_excluded"

# --- Cross-trait genetic-correlation table (one row per protein) ---
# The trait-vs-protein genetic correlation (rg) with its jackknife SE and the two
# heritabilities it is built from. gcov is the genetic covariance (rg numerator) and
# gcov_intercept its cross-trait LDSC intercept (a sample-overlap diagnostic).
PPP_RG_OID_COL = "oid"
PPP_RG_GENE_COL = "gene"
PPP_RG_VARIANT_SET_COL = "variant_set"
PPP_RG_RG_COL = "rg"
PPP_RG_RG_SE_COL = "rg_se"
PPP_RG_RG_Z_COL = "rg_z"
PPP_RG_RG_P_COL = "rg_p"
PPP_RG_GCOV_COL = "gcov"
PPP_RG_GCOV_INTERCEPT_COL = "gcov_intercept"
PPP_RG_H2_TRAIT_COL = "h2_trait"
PPP_RG_H2_PROTEIN_COL = "h2_protein"
PPP_RG_N_SNPS_COL = "n_snps"
PPP_RG_N_BAR_TRAIT_COL = "n_bar_trait"
PPP_RG_N_BAR_PROTEIN_COL = "n_bar_protein"
