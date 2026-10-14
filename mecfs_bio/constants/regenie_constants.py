"""
Column names emitted by regenie, as they appear in raw regenie summary-statistics
files (e.g. the per-chromosome UKB-PPP files).
"""

REGENIE_CHROM_COL = "CHROM"
REGENIE_GENPOS_COL = "GENPOS"  # base-pair position on the analysis genome build
REGENIE_ID_COL = "ID"  # variant marker name, e.g. CHROM:POS:ALLELE0:ALLELE1:...
REGENIE_ALLELE0_COL = "ALLELE0"  # reference / non-effect allele
REGENIE_ALLELE1_COL = "ALLELE1"  # alternative / effect (tested) allele
REGENIE_A1FREQ_COL = "A1FREQ"  # frequency of ALLELE1
REGENIE_INFO_COL = "INFO"
REGENIE_N_COL = "N"
REGENIE_BETA_COL = "BETA"
REGENIE_SE_COL = "SE"
REGENIE_CHISQ_COL = "CHISQ"
REGENIE_LOG10P_COL = "LOG10P"

# Binary-trait output columns. gwaslab passes them through without renaming, so they
# survive into gwaslab-format tables (DecodeME carries all of them).
REGENIE_A1FREQ_CASES_COL = "A1FREQ_CASES"  # frequency of ALLELE1 in cases
REGENIE_A1FREQ_CONTROLS_COL = "A1FREQ_CONTROLS"  # frequency of ALLELE1 in controls
REGENIE_N_CASES_COL = "N_CASES"
REGENIE_N_CONTROLS_COL = "N_CONTROLS"
REGENIE_TEST_COL = "TEST"
REGENIE_EXTRA_COL = "EXTRA"
