"""
Column and file names for the fast all-pairs CT-LDSC outputs.

The pairs table is read by GeneticCorrelationClustermapTask (p1, p2, rg, p) and by LAVA's
sample-overlap reader (p1, p2, gcov_int, h2_int); the heritability table is read by the same LAVA
reader (p, Intercept). Those names are therefore fixed by the readers.
"""

# Pairs table: one row per unordered trait pair.
FAST_CT_LDSC_P1_COL = "p1"
FAST_CT_LDSC_P2_COL = "p2"
FAST_CT_LDSC_RG_COL = "rg"
FAST_CT_LDSC_RG_SE_COL = "se"
FAST_CT_LDSC_RG_Z_COL = "z"
FAST_CT_LDSC_RG_P_COL = "p"
FAST_CT_LDSC_GCOV_COL = "gcov"
FAST_CT_LDSC_GCOV_SE_COL = "gcov_se"
FAST_CT_LDSC_GCOV_INT_COL = "gcov_int"
FAST_CT_LDSC_GCOV_INT_SE_COL = "gcov_int_se"
FAST_CT_LDSC_H2_OBS_P1_COL = "h2_obs_p1"
FAST_CT_LDSC_H2_OBS_P2_COL = "h2_obs_p2"
FAST_CT_LDSC_H2_INT_COL = "h2_int"
FAST_CT_LDSC_N_SNPS_COL = "n_snps"

# Heritability table: one row per trait.
FAST_CT_LDSC_TRAIT_COL = "p"
FAST_CT_LDSC_H2_OBS_COL = "h2_obs"
FAST_CT_LDSC_H2_OBS_SE_COL = "h2_obs_se"
FAST_CT_LDSC_H2_LIAB_COL = "h2_liab"
FAST_CT_LDSC_H2_LIAB_SE_COL = "h2_liab_se"
FAST_CT_LDSC_INTERCEPT_COL = "Intercept"
FAST_CT_LDSC_INTERCEPT_SE_COL = "Intercept_se"
FAST_CT_LDSC_LAMBDA_GC_COL = "lambda_gc"
FAST_CT_LDSC_MEAN_CHI2_COL = "mean_chi2"
FAST_CT_LDSC_N_BAR_COL = "n_bar"
FAST_CT_LDSC_COVERAGE_COL = "coverage"
FAST_CT_LDSC_ALLELE_MISMATCH_FRACTION_COL = "allele_mismatch_fraction"
FAST_CT_LDSC_POSITION_CONCORDANCE_COL = "position_concordance"

# Aligned-trait directory: z and n in context row order, plus the alignment diagnostics.
FAST_CT_LDSC_Z_COL = "z"
FAST_CT_LDSC_N_COL = "n"
FAST_CT_LDSC_ALIGNED_FILENAME = "aligned.parquet"
FAST_CT_LDSC_DIAGNOSTICS_FILENAME = "diagnostics.json"

# Columns of the HapMap3 w_hm3.snplist.
HAPMAP3_SNPLIST_RSID_COL = "SNP"
HAPMAP3_SNPLIST_A1_COL = "A1"
HAPMAP3_SNPLIST_A2_COL = "A2"
