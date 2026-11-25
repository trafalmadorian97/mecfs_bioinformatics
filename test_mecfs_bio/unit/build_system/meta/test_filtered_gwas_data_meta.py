import pytest

from mecfs_bio.build_system.meta.asset_id import AssetId
from mecfs_bio.build_system.meta.filtered_gwas_data_meta import FilteredGWASDataMeta
from mecfs_bio.build_system.meta.harmonization_info import HarmonizationInfo
from mecfs_bio.constants.genomic_coordinate_constants import GenomeBuild

_HARMONIZED_TO_19 = HarmonizationInfo(build="19", ref_allele_col="NEA", pos_col="POS")


def _meta(
    genome_build: GenomeBuild | None, harmonization_info: HarmonizationInfo | None
) -> FilteredGWASDataMeta:
    return FilteredGWASDataMeta(
        id=AssetId("a"),
        trait="trait",
        project="project",
        sub_dir="processed",
        harmonization_info=harmonization_info,
        genome_build=genome_build,
    )


def test_harmonization_info_without_genome_build_is_rejected() -> None:
    with pytest.raises(AssertionError):
        _meta(genome_build=None, harmonization_info=_HARMONIZED_TO_19)


def test_genome_build_contradicting_harmonization_info_is_rejected() -> None:
    with pytest.raises(AssertionError):
        _meta(genome_build="38", harmonization_info=_HARMONIZED_TO_19)


def test_genome_build_agreeing_with_harmonization_info_is_accepted() -> None:
    _meta(genome_build="19", harmonization_info=_HARMONIZED_TO_19)


def test_genome_build_without_harmonization_info_is_accepted() -> None:
    _meta(genome_build="38", harmonization_info=None)
