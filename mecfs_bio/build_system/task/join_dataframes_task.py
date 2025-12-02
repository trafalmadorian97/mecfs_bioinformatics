from pathlib import Path
from typing import Sequence

import structlog
from attrs import frozen
from narwhals.typing import JoinStrategy

from mecfs_bio.build_system.asset.base_asset import Asset
from mecfs_bio.build_system.asset.file_asset import FileAsset
from mecfs_bio.build_system.meta.asset_id import AssetId
from mecfs_bio.build_system.meta.meta import Meta
from mecfs_bio.build_system.meta.read_spec.dataframe_read_spec import (
    DataFrameReadSpec,
    DataFrameTextFormat,
)
from mecfs_bio.build_system.meta.read_spec.read_dataframe import (
    scan_dataframe_asset,
)
from mecfs_bio.build_system.meta.result_table_meta import ResultTableMeta
from mecfs_bio.build_system.rebuilder.fetch.base_fetch import Fetch
from mecfs_bio.build_system.task.base_task import Task
from mecfs_bio.build_system.wf.base_wf import WF

logger = structlog.get_logger()


@frozen
class JoinDataFramesTask(Task):
    """
    Task to load to dataframes, join them, and then write out the resulting dataframe.

    By default, writes as csv.
    """

    _df_1_task: Task
    _df_2_task: Task
    how: JoinStrategy
    left_on: Sequence[str]
    right_on: Sequence[str]
    _meta: Meta

    @property
    def _df_1_id(self) -> AssetId:
        return self._df_1_task.asset_id

    @property
    def _df_12_id(self) -> AssetId:
        return self._df_2_task.asset_id

    @property
    def _df_1_meta(self) -> Meta:
        return self._df_1_task.meta

    @property
    def _df_2_meta(self) -> Meta:
        return self._df_2_task.meta

    @property
    def meta(self) -> Meta:
        return self._meta

    @property
    def deps(self) -> list["Task"]:
        return [self._df_1_task, self._df_2_task]

    def execute(self, scratch_dir: Path, fetch: Fetch, wf: WF) -> Asset:
        result_path = scratch_dir / "result.csv"
        asset_1 = fetch(self._df_1_id)
        df_1 = scan_dataframe_asset(
            asset_1,
            meta=self._df_1_meta,
        )
        asset_2 = fetch(self._df_12_id)
        df_2 = scan_dataframe_asset(
            asset_2,
            meta=self._df_2_meta,
        )
        joined = df_1.join(
            df_2, how=self.how, left_on=list(self.left_on), right_on=list(self.right_on)
        )
        result = joined.collect().to_pandas()
        # logger.debug(f"\nresult of join:\n {result.to_markdown(index=False)}")
        result.to_csv(result_path, index=False)
        return FileAsset(result_path)

    @classmethod
    def create_from_result_df(
        cls,
        asset_id: str,
        result_df_task: Task,
        reference_df_task,
        how: JoinStrategy,
        left_on: Sequence[str],
        right_on: Sequence[str],
    ):
        """
        Join a result dataframe to a reference dataframe.
        """

        source_meta = result_df_task.meta
        assert isinstance(source_meta, ResultTableMeta)
        meta = ResultTableMeta(
            asset_id=AssetId(asset_id),
            trait=source_meta.trait,
            project=source_meta.project,
            extension=".csv",
            read_spec=DataFrameReadSpec(DataFrameTextFormat(separator=",")),
        )
        return cls(
            df_1_task=result_df_task,
            df_2_task=reference_df_task,
            how=how,
            left_on=left_on,
            right_on=right_on,
            meta=meta,
        )
