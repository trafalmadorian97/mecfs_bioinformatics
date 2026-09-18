from pathlib import PurePath

import pytest

from mecfs_bio.build_system.meta.asset_id import AssetId
from mecfs_bio.build_system.meta.harmonization_info import HarmonizationInfo
from mecfs_bio.build_system.meta.read_spec.dataframe_read_spec import (
    DataFrameParquetFormat,
    DataFrameReadSpec,
)
from mecfs_bio.build_system.meta.reference_meta.harmonizable_reference_table_meta import (
    HarmonizableReferenceTableMeta,
)
from mecfs_bio.build_system.task.fake_task import FakeTask
from mecfs_bio.build_system.task.rename_cols_task import RenameColsTask


def _harmonizable_source(task_id: str) -> FakeTask:
    return FakeTask(
        meta=HarmonizableReferenceTableMeta(
            group="ukbb_reference_ld",
            sub_group="chr1_1_2",
            sub_folder=PurePath("raw"),
            extension=".gz",
            id=AssetId(task_id),
            filename="chr1_1_2",
            read_spec=DataFrameReadSpec(DataFrameParquetFormat()),
            harmonization_info=HarmonizationInfo(
                build="19", ref_allele_col="allele1", pos_col="position"
            ),
        )
    )


def test_create_remaps_harmonization_info():
    task = RenameColsTask.create(
        source_task=_harmonizable_source("src"),
        asset_id="renamed",
        renames={"allele1": "NEA", "allele2": "EA", "position": "POS"},
    )
    assert isinstance(task.meta, HarmonizableReferenceTableMeta)
    info = task.meta.harmonization_info
    assert info is not None
    assert (info.build, info.ref_allele_col, info.pos_col) == ("19", "NEA", "POS")


def test_create_rejects_colliding_provenance_columns():
    with pytest.raises(AssertionError):
        RenameColsTask.create(
            source_task=_harmonizable_source("src"),
            asset_id="renamed",
            renames={
                "allele1": "position"
            },  # ref-allele col collides onto pos col name
        )
