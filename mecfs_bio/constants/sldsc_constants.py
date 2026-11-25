"""
Constants for cell-type-specific stratified LD score regression (S-LDSC).
"""

from typing import Literal

# exact: per-cell-type regression weights, reproducing ldsc's cell_type_specific.
# shared_controls: one weight vector per control group (target column excluded from the
# weight formula), enabling a Frisch-Waugh shortcut; faster but approximate.
SLDSCWeighting = Literal["exact", "shared_controls"]

# What to do with sumstats rsIDs that occur more than once among the S-LDSC regression SNPs
# (after the HapMap3 allele filter). raise: fail fast. drop: remove every row of each such rsID,
# since which row is correct cannot be told. gwaslab instead keeps them all, double-counting.
SLDSCDuplicateRsidPolicy = Literal["raise", "drop"]

# Result-table columns written by ldsc's cell_type_specific and by the vectorized Task.
SLDSC_CTS_NAME_COL = "Name"
SLDSC_CTS_COEF_COL = "Coefficient"
SLDSC_CTS_COEF_SE_COL = "Coefficient_std_error"
SLDSC_CTS_P_COL = "Coefficient_P_value"

# Minimum number of sumstats SNPs in the HapMap3 regression-SNP set for S-LDSC to run. Fewer
# usually indicates an upstream processing problem (for example, lost rsIDs).
SLDSC_DEFAULT_MIN_HAPMAP3_SNPS = 100_000
