"""
Identifiers and index-specific column names for the Western et al. 2024 CSF pQTL
database (the parquet table whose row order is the canonical alignment slot for the
per-aptamer beta/se/N files).

The index otherwise uses gwaslab-standard column names (see gwaslab_constants):
CHR, POS, EA, NEA, rsID, EAF. The names below are the index-specific additions, and
the NewTypes give the aptamer/accession identifiers distinct static types.
"""

from typing import NewType

# Strand-ambiguous (palindromic A/T or C/G) flag, for downstream cross-dataset use.
CSF_INDEX_IS_STRAND_AMBIGUOUS_COL = "is_strand_ambiguous"

# Internal, order-agnostic allele-set key used only during index construction and
# alignment.
CSF_INDEX_ALLELE_KEY_COL = "__allele_key__"

# SomaScan analyte id, the aptamer primary key, e.g. "X13681.173". Gene symbol is not
# unique across aptamers, so this (not the gene) namespaces every per-aptamer asset.
Analyte = NewType("Analyte", str)

# SomaScan sequence id, e.g. "13681-173" (the analyte with a different separator).
SeqId = NewType("SeqId", str)

# GWAS Catalog study accession, e.g. "GCST90421540".
GcstAccession = NewType("GcstAccession", str)

# GWAS Catalog trait label for a study, e.g. "Beta-crystallin B2 levels". Unique per
# study; it is the only key the study metadata and the aptamer table share, so the
# manifest resolver maps it to an Analyte (see regenerate_csf_manifest.py).
Trait = NewType("Trait", str)
