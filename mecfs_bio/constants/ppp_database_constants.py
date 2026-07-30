"""
Column names specific to the shared UKB-PPP variant index (the parquet table whose
row order is the canonical alignment slot for the per-protein beta/se files).

The index otherwise uses gwaslab-standard column names (see gwaslab_constants):
CHR, POS, EA, NEA, rsID, EAF. The names below are the index-specific additions.
"""

from typing import NewType

# Secondary hg19 position (primary POS is hg38). Complete when the index is
# templated off a PPP protein, since the hg19 position is carried in the regenie ID.
PPP_INDEX_POS_HG19_COL = "POS_HG19"

# Strand-ambiguous (palindromic A/T or C/G) flag, for downstream cross-dataset use.
PPP_INDEX_IS_STRAND_AMBIGUOUS_COL = "is_strand_ambiguous"

# Internal, order-agnostic allele-set key used only during index construction.
PPP_INDEX_ALLELE_KEY_COL = "__allele_key__"


Oid = NewType("Oid", str)
SynID = NewType("SynID", str)
