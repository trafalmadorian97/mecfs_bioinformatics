"""
Column layout shared by all FinnGen R13 summary statistics files.

Missing rsids are encoded as NA in the raw files, and read as null.
"""

from mecfs_bio.build_system.task.gwaslab.gwaslab_create_sumstats_task import (
    GWASLabColumnSpecifiers,
)

FINNGEN_R13_RSID_COL = "rsids"

FINNGEN_R13_COLUMN_SPECIFIERS = GWASLabColumnSpecifiers(
    rsid=FINNGEN_R13_RSID_COL,
    chrom="#chrom",
    pos="pos",
    ea="alt",
    nea="ref",
    beta="beta",
    se="sebeta",
    p="pval",
    eaf="af_alt",
)
