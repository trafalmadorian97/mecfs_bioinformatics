"""
Task to filter SNPs in a GWAS according to another table of SNPs passing quality control.
Implemented as a DataFrame join.
"""

from pathlib import Path

from attrs import frozen

from mecfs_bio.build_system.asset.file_asset import FileAsset
from mecfs_bio.build_system.meta.filtered_gwas_data_meta import FilteredGWASDataMeta
from mecfs_bio.build_system.meta.meta import Meta
from mecfs_bio.build_system.meta.read_spec.read_dataframe import (
    scan_dataframe_asset,
)
from mecfs_bio.build_system.rebuilder.fetch.base_fetch import Fetch
from mecfs_bio.build_system.task.base_task import Task
from mecfs_bio.build_system.wf.base_wf import WF


@frozen
class FilterSNPsTask(Task):
    _raw_gwas_task: Task
    _snp_list_task: Task
    _meta: FilteredGWASDataMeta
    _col_in_raw_data: str = "ID"
    _col_in_filter_data: str = "ID"

    @property
    def meta(self) -> Meta:
        return self._meta

    @property
    def deps(self) -> list["Task"]:
        return [self._raw_gwas_task, self._snp_list_task]

    def execute(self, scratch_dir: Path, fetch: Fetch, wf: WF) -> FileAsset:
        df_1 = scan_dataframe_asset(
            asset=fetch(self._raw_gwas_task.asset_id), meta=self._raw_gwas_task.meta
        )
        df_2 = scan_dataframe_asset(
            asset=fetch(self._snp_list_task.asset_id), meta=self._snp_list_task.meta
        )
        result = df_1.join(
            df_2, left_on=self._col_in_raw_data, right_on=self._col_in_filter_data
        )
        target_path = scratch_dir / "tmp.parqet"
        result.sink_parquet(target_path)
        return FileAsset(target_path)
