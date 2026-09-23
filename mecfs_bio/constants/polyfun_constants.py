"""Column names used by PolyFun-provided dataframes.

The baseline-LF annotation members, LD-score members, derived annotation matrix,
and precomputed prior all key variants by CHR, BP, SNP, A1, A2. A1 and A2 are
REF-oriented: A1 is the hg19 reference allele (the gwaslab non-effect allele) and
A2 is the alternate allele (the gwaslab effect allele).
"""

from mecfs_bio.constants.gwaslab_constants import (
    GWASLAB_CHROM_COL,
    GWASLAB_EFFECT_ALLELE_COL,
    GWASLAB_NON_EFFECT_ALLELE_COL,
    GWASLAB_POS_COL,
)

POLYFUN_CHR_COL = "CHR"
POLYFUN_BP_COL = "BP"
POLYFUN_SNP_COL = "SNP"
POLYFUN_A1_COL = "A1"
POLYFUN_A2_COL = "A2"

# Renames a PolyFun variant key to the gwaslab four-part join key.
POLYFUN_TO_GWASLAB_KEY_RENAME: dict[str, str] = {
    POLYFUN_CHR_COL: GWASLAB_CHROM_COL,
    POLYFUN_BP_COL: GWASLAB_POS_COL,
    POLYFUN_A1_COL: GWASLAB_NON_EFFECT_ALLELE_COL,
    POLYFUN_A2_COL: GWASLAB_EFFECT_ALLELE_COL,
}
