"""Rename dataframe columns. A deliberately narrow alternative to PipeDataFrameTask: because the
rename map is the whole operation, create() can update HarmonizationInfo (ref/pos columns)
accurately across the rename, which a general pipe cannot."""

from collections.abc import Mapping
from pathlib import Path, PurePath

import narwhals
from attrs import frozen

from mecfs_bio.build_system.asset.base_asset import Asset
from mecfs_bio.build_system.asset.file_asset import FileAsset
from mecfs_bio.build_system.meta.asset_id import AssetId
from mecfs_bio.build_system.meta.filtered_gwas_data_meta import FilteredGWASDataMeta
from mecfs_bio.build_system.meta.harmonization_info import HarmonizationInfo
from mecfs_bio.build_system.meta.meta import Meta
from mecfs_bio.build_system.meta.read_spec.read_dataframe import (
    ValidBackend,
    scan_dataframe_asset,
)
from mecfs_bio.build_system.meta.reference_meta.harmonizable_reference_table_meta import (
    HarmonizableReferenceTableMeta,
)
from mecfs_bio.build_system.meta.reference_meta.reference_file_meta import (
    ReferenceFileMeta,
)
from mecfs_bio.build_system.rebuilder.fetch.base_fetch import Fetch
from mecfs_bio.build_system.task.base_task import Task
from mecfs_bio.build_system.task.dataframe_output import (
    OutFormat,
    ParquetOutFormat,
    get_extension_and_read_spec_from_format,
    write_df_according_to_format,
)
from mecfs_bio.build_system.wf.base_wf import WF

_OUT_FILENAME = "renamed"


def remap_harmonization_info(
    info: HarmonizationInfo, renames: Mapping[str, str]
) -> HarmonizationInfo:
    new_ref = renames.get(info.ref_allele_col, info.ref_allele_col)
    new_pos = renames.get(info.pos_col, info.pos_col)
    assert new_ref != new_pos, (
        f"rename collides the ref-allele column and position column onto {new_ref!r}"
    )
    return HarmonizationInfo(build=info.build, ref_allele_col=new_ref, pos_col=new_pos)


@frozen
class RenameColsTask(Task):
    meta: Meta
    source_task: Task
    renames: Mapping[str, str]
    out_format: OutFormat = ParquetOutFormat()
    backend: ValidBackend = "polars"

    @property
    def deps(self) -> list[Task]:
        return [self.source_task]

    def execute(self, scratch_dir: Path, fetch: Fetch, wf: WF) -> Asset:
        source = fetch(self.source_task.asset_id)
        frame: narwhals.LazyFrame = scan_dataframe_asset(
            asset=source, meta=self.source_task.meta, parquet_backend=self.backend
        )
        frame = frame.rename(dict(self.renames))
        out_path = scratch_dir / _OUT_FILENAME
        write_df_according_to_format(
            df=frame, out_path=out_path, out_format=self.out_format
        )
        return FileAsset(out_path)

    @classmethod
    def create(
        cls,
        source_task: Task,
        asset_id: str,
        renames: Mapping[str, str],
        out_format: OutFormat = ParquetOutFormat(),
        backend: ValidBackend = "polars",
    ) -> "RenameColsTask":
        source_meta = source_task.meta
        extension, read_spec = get_extension_and_read_spec_from_format(
            out_format=out_format
        )
        meta: Meta
        if isinstance(source_meta, HarmonizableReferenceTableMeta):
            meta = HarmonizableReferenceTableMeta(
                group=source_meta.group,
                sub_group=source_meta.sub_group,
                sub_folder=PurePath("processed"),
                extension=extension,
                id=AssetId(asset_id),
                filename=None,
                read_spec=read_spec,
                harmonization_info=(
                    remap_harmonization_info(source_meta.harmonization_info, renames)
                    if source_meta.harmonization_info is not None
                    else None
                ),
            )
        elif isinstance(source_meta, FilteredGWASDataMeta):
            meta = FilteredGWASDataMeta(
                id=AssetId(asset_id),
                trait=source_meta.trait,
                project=source_meta.project,
                sub_dir=source_meta.sub_dir,
                read_spec=read_spec,
                harmonization_info=(
                    remap_harmonization_info(source_meta.harmonization_info, renames)
                    if source_meta.harmonization_info is not None
                    else None
                ),
            )
        elif isinstance(source_meta, ReferenceFileMeta):
            meta = ReferenceFileMeta(
                group=source_meta.group,
                sub_group=source_meta.sub_group,
                sub_folder=PurePath("processed"),
                extension=extension,
                id=AssetId(asset_id),
                read_spec=read_spec,
            )
        else:
            raise ValueError(
                f"RenameColsTask: unsupported source meta {type(source_meta).__name__}"
            )
        return cls(
            meta=meta,
            source_task=source_task,
            renames=dict(renames),
            out_format=out_format,
            backend=backend,
        )
