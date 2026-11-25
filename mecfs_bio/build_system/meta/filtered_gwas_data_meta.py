"""
Metadata describing an asset resulting from the processing of GWAS data.
"""

from pathlib import PurePath

from attrs import frozen

from mecfs_bio.build_system.meta.asset_id import AssetId
from mecfs_bio.build_system.meta.base_meta import FileMeta
from mecfs_bio.build_system.meta.harmonization_info import HarmonizationInfo
from mecfs_bio.build_system.meta.read_spec.read_spec import ReadSpec
from mecfs_bio.constants.genomic_coordinate_constants import GenomeBuild


@frozen(slots=True)
class FilteredGWASDataMeta(FileMeta):
    """
    genome_build records the genome build of the position column, without any claim
    about allele orientation. harmonization_info makes the stronger claim that the
    non-effect allele column equals the reference base of its build. When
    harmonization_info is set, genome_build must be set to the same build.
    """

    @property
    def asset_id(self) -> AssetId:
        return self.id

    id: AssetId
    trait: str
    project: str
    sub_dir: str | PurePath
    read_spec: ReadSpec | None = None
    extension: str = ".parquet"
    harmonization_info: HarmonizationInfo | None = None
    genome_build: GenomeBuild | None = None

    def __attrs_post_init__(self) -> None:
        if self.harmonization_info is not None:
            assert self.genome_build == self.harmonization_info.build, (
                f"genome_build={self.genome_build!r} must equal "
                f"harmonization_info.build={self.harmonization_info.build!r}"
            )
