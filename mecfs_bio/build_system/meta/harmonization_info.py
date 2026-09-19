"""Provenance of genome-reference harmonization: which build an asset is oriented to,
and which columns carry that orientation."""

from attrs import frozen

from mecfs_bio.constants.genomic_coordinate_constants import GenomeBuild


@frozen
class HarmonizationInfo:
    """Records that, at the positions in pos_col, ref_allele_col equals the reference
    base of the given genome build. build is the load-bearing field for construction-time
    build-match checks; the column names document the claim and support an optional
    data-level verification."""

    build: GenomeBuild
    ref_allele_col: str
    pos_col: str
