"""
Note:
    According to the paper Horton, Roger, et al.
    "Gene map of the extended human MHC." Nature Reviews Genetics 5.12 (2004): 889-899.
    The extended MHC region is defined as being between the HIST1H2AA gene and the RPL12P1 pseudo-gene.

    According to https://www.genecards.org/cgi-bin/carddisp.pl?gene=H2AC1,
    HIST1H2AA starts at position 25726291 in genome build 37

    According to https://www.genecards.org/cgi-bin/carddisp.pl?gene=RPL12P1

    RPL12P1 ends at genome position 33368421 in genome build 37

    The build-38 bounds use the same two boundary genes in GRCh38 (GeneCards):
    H2AC1 (HIST1H2AA) spans chr6:25723397-25726562 and RPL12P1 spans
    chr6:33400015-33400644, so the extended MHC brackets [25723397, 33400644].
    Excluding this region is standard for LD-score regression (its long-range,
    atypical LD violates the LD-score model; cf. Bulik-Sullivan et al. 2015).
"""

from typing import Literal

from mecfs_bio.constants.vocabulary_classes.genomic_interval import GenomicInterval

EXTENDED_MHC_START_BUILD_37 = 25726291
EXTENDED_MHC_END_BUILD_37 = 33368421

EXTENDED_MHC_START_BUILD_38 = 25723397
EXTENDED_MHC_END_BUILD_38 = 33400644


MHCRegion = Literal["classical", "extended"]


EXTENDED_MHC_BUILD_37 = GenomicInterval(
    start=EXTENDED_MHC_START_BUILD_37,
    end=EXTENDED_MHC_END_BUILD_37,
    chrom=6,
)
EXTENDED_MHC_BUILD_38 = GenomicInterval(
    start=EXTENDED_MHC_START_BUILD_38,
    end=EXTENDED_MHC_END_BUILD_38,
    chrom=6,
)
GenomeBuild = Literal["19", "38"]


def extended_mhc_interval(build: GenomeBuild) -> GenomicInterval:
    """Return the extended MHC region for the given genome build (19/GRCh37 or
    38/GRCh38), bracketed by the HIST1H2AA and RPL12P1 boundary genes."""
    if build == "19":
        return EXTENDED_MHC_BUILD_37
    if build == "38":
        return EXTENDED_MHC_BUILD_38
    raise ValueError(f"Extended MHC region not implemented for genome build {build}")
