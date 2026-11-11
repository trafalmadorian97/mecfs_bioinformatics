"""
Column names for the GWAS-SSF v1.0 standard (GWAS Catalog Summary Statistics
Format), as emitted in the harmonized .tsv.gz files served from the EBI GWAS
Catalog FTP.

These names are a published public standard, not specific to any one deposit, so
they live here rather than alongside the CSF-database tasks that happen to be the
first consumer.
"""

GWAS_SSF_CHROM_COL = "chromosome"
GWAS_SSF_POS_COL = "base_pair_location"  # 1-based, on the study's genome build
GWAS_SSF_EFFECT_ALLELE_COL = "effect_allele"
GWAS_SSF_OTHER_ALLELE_COL = "other_allele"
GWAS_SSF_BETA_COL = "beta"
GWAS_SSF_SE_COL = "standard_error"
GWAS_SSF_EFFECT_ALLELE_FREQ_COL = "effect_allele_frequency"
GWAS_SSF_NEG_LOG10_P_COL = "neg_log_10_p_value"
GWAS_SSF_VARIANT_ID_COL = "variant_id"
GWAS_SSF_RSID_COL = "rs_id"
GWAS_SSF_N_COL = "n"
