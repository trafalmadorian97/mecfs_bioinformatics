"""Record of genome-reference harmonization: which build an asset is oriented to,
and which columns carry that orientation."""

from attrs import frozen

from mecfs_bio.constants.genomic_coordinate_constants import GenomeBuild


@frozen
class HarmonizationInfo:
    """Records that, at the positions in pos_col, ref_allele_col equals the reference
    base of the given genome build."""

    build: GenomeBuild
    ref_allele_col: str
    pos_col: str
