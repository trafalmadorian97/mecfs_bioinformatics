"""
Note:
    According to the paper Horton, Roger, et al.
    "Gene map of the extended human MHC." Nature Reviews Genetics 5.12 (2004): 889-899.
    The extended MHC region is defined as being between the HIST1H2AA gene and the RPL12P1 pseudo-gene.

    According to https://www.genecards.org/cgi-bin/carddisp.pl?gene=H2AC1,
    HIST1H2AA starts at position 25726291 in genome build 37

    According to https://www.genecards.org/cgi-bin/carddisp.pl?gene=RPL12P1

    RPL12P1 ends at genome position 33368421 in genome build 37
"""

from typing import Literal

from mecfs_bio.constants.vocabulary_classes.genomic_interval import GenomicInterval

EXTENDED_MHC_START_BUILD_37 = 25726291
EXTENDED_MHC_END_BUILD_37 = 33368421


MHCRegion = Literal["classical", "extended"]


EXTENDED_MHC_BUILD_37 = GenomicInterval(
    start=EXTENDED_MHC_START_BUILD_37,
    end=EXTENDED_MHC_END_BUILD_37,
    chrom=6,
)
GenomeBuild = Literal["19", "38"]
