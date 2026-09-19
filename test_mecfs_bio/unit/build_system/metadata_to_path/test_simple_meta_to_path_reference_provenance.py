from pathlib import PurePath

from mecfs_bio.build_system.meta.asset_id import AssetId
from mecfs_bio.build_system.meta.harmonization_info import HarmonizationInfo
from mecfs_bio.build_system.meta.reference_meta.fasta_meta import FASTAMeta
from mecfs_bio.build_system.meta.reference_meta.harmonizable_reference_table_meta import (
    HarmonizableReferenceTableMeta,
)
from mecfs_bio.build_system.meta.reference_meta.reference_data_directory_meta import (
    ReferenceDataDirectoryMeta,
)
from mecfs_bio.build_system.meta.reference_meta.reference_file_meta import (
    ReferenceFileMeta,
)
from mecfs_bio.build_system.rebuilder.metadata_to_path.simple_meta_to_path import (
    simple_meta_to_relative_path,
)


def test_fasta_meta_resolves_like_reference_directory():
    fasta = FASTAMeta(
        group="genome_sequence",
        sub_group="ucsc_hg19",
        sub_folder=PurePath("processed"),
        id=AssetId("idx"),
        build="19",
    )
    twin = ReferenceDataDirectoryMeta(
        group="genome_sequence",
        sub_group="ucsc_hg19",
        sub_folder=PurePath("processed"),
        id=AssetId("idx"),
    )
    assert simple_meta_to_relative_path(fasta) == simple_meta_to_relative_path(twin)


def test_harmonizable_table_resolves_like_reference_file():
    table = HarmonizableReferenceTableMeta(
        group="g",
        sub_group="s",
        sub_folder=PurePath("processed"),
        extension=".parquet",
        id=AssetId("t"),
        filename="f",
        harmonization_info=HarmonizationInfo(
            build="19", ref_allele_col="REF", pos_col="POS"
        ),
    )
    twin = ReferenceFileMeta(
        group="g",
        sub_group="s",
        sub_folder=PurePath("processed"),
        extension=".parquet",
        id=AssetId("t"),
        filename="f",
    )
    assert simple_meta_to_relative_path(table) == simple_meta_to_relative_path(twin)
