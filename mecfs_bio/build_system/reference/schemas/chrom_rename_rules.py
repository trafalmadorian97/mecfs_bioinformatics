import pandas as pd

CHROM_RENAME_RULES = {
    "chrX": 23,
    "chrY": 24,
    "chrM": 25,
} | {f"chr{i}": i for i in range(1, 25)}

CHROM_RENAME_RULES_NUMERIC_X_Y_M = {
    "X": 23,
    "Y": 24,
    "M": 25,
} | {str(i): i for i in range(1, 25)}

CHROM_RENAME_DF_NUMERIC_X_Y_M = pd.DataFrame(
    {
        "CHROM": CHROM_RENAME_RULES_NUMERIC_X_Y_M.keys(),
        "int_chrom": CHROM_RENAME_RULES_NUMERIC_X_Y_M.values(),
    },
)
