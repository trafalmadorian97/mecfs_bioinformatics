from pathlib import PurePath

import pytest

from mecfs_bio.build_system.meta.asset_id import AssetId
from mecfs_bio.build_system.meta.filtered_gwas_data_meta import FilteredGWASDataMeta
from mecfs_bio.build_system.meta.harmonization_info import HarmonizationInfo
from mecfs_bio.build_system.meta.reference_meta.harmonizable_reference_table_meta import (
    HarmonizableReferenceTableMeta,
)
from mecfs_bio.build_system.task.fake_task import FakeTask
from mecfs_bio.build_system.task.r_tasks.susie_r_finemap_task import (
    BroadInstituteFormatLDMatrix,
    SusieRFinemapTask,
)
from mecfs_bio.constants.genomic_coordinate_constants import GenomeBuild


def _gwas(build: GenomeBuild) -> FakeTask:
    return FakeTask(
        meta=FilteredGWASDataMeta(
            id=AssetId("g"),
            trait="decode_me",
            project="me_cfs",
            sub_dir=PurePath("processed"),
            harmonization_info=HarmonizationInfo(
                build=build, ref_allele_col="NEA", pos_col="POS"
            ),
        )
    )


def _ld_labels(build: GenomeBuild) -> FakeTask:
    return FakeTask(
        meta=HarmonizableReferenceTableMeta(
            group="ukbb_reference_ld",
            sub_group="chr1",
            sub_folder=PurePath("processed"),
            extension=".parquet",
            id=AssetId("ld"),
            harmonization_info=HarmonizationInfo(
                build=build, ref_allele_col="NEA", pos_col="POS"
            ),
        )
    )


def test_build_mismatch_raises():
    with pytest.raises(AssertionError):
        SusieRFinemapTask.create(
            asset_id="s",
            gwas_data_task=_gwas("19"),
            ld_labels_task=_ld_labels("38"),
            ld_matrix_source=BroadInstituteFormatLDMatrix(
                FakeTask(meta=_ld_labels("38").meta)
            ),
            effective_sample_size=1000,
        )


def test_matching_builds_construct():
    task = SusieRFinemapTask.create(
        asset_id="s",
        gwas_data_task=_gwas("19"),
        ld_labels_task=_ld_labels("19"),
        ld_matrix_source=BroadInstituteFormatLDMatrix(
            FakeTask(meta=_ld_labels("19").meta)
        ),
        effective_sample_size=1000,
    )
    assert isinstance(task, SusieRFinemapTask)
