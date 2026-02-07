from pathlib import Path
from typing import Sequence

import structlog
from attrs import frozen
from narwhals.typing import JoinStrategy

from mecfs_bio.build_system.asset.base_asset import Asset
from mecfs_bio.build_system.asset.file_asset import FileAsset
from mecfs_bio.build_system.meta.asset_id import AssetId
from mecfs_bio.build_system.meta.filtered_gwas_data_meta import FilteredGWASDataMeta
from mecfs_bio.build_system.meta.gwaslab_meta.gwaslab_lead_variants_meta import (
    GWASLabLeadVariantsMeta,
)
from mecfs_bio.build_system.meta.meta import Meta
from mecfs_bio.build_system.meta.read_spec.read_dataframe import (
    ValidBackend,
    scan_dataframe_asset,
)
from mecfs_bio.build_system.meta.result_table_meta import ResultTableMeta
from mecfs_bio.build_system.rebuilder.fetch.base_fetch import Fetch
from mecfs_bio.build_system.task.base_task import Task
from mecfs_bio.build_system.task.pipe_dataframe_task import (
    CSVOutFormat,
    OutFormat,
    ParquetOutFormat,
    get_extension_and_read_spec_from_format,
)
from mecfs_bio.build_system.task.pipes.data_processing_pipe import DataProcessingPipe
from mecfs_bio.build_system.task.pipes.identity_pipe import IdentityPipe
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

    out_format: OutFormat = CSVOutFormat(sep=",")
    df_1_pipe: DataProcessingPipe = IdentityPipe()
    df_2_pipe: DataProcessingPipe = IdentityPipe()
    out_pipe: DataProcessingPipe = IdentityPipe()
    backend: ValidBackend = "polars"

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
            asset_1, meta=self._df_1_meta, parquet_backend=self.backend
        )
        df_1 = self.df_1_pipe.process(df_1)
        asset_2 = fetch(self._df_12_id)
        df_2 = scan_dataframe_asset(
            asset_2, meta=self._df_2_meta, parquet_backend=self.backend
        )
        df_2 = self.df_2_pipe.process(df_2)
        joined = df_1.join(
            df_2, how=self.how, left_on=list(self.left_on), right_on=list(self.right_on)
        )
        joined = self.out_pipe.process(joined)
        if isinstance(self.out_format, CSVOutFormat):
            result = joined.collect().to_pandas()
            result.to_csv(result_path, index=False, sep=self.out_format.sep)
        elif isinstance(self.out_format, ParquetOutFormat):
            joined.sink_parquet(result_path)
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
        out_format: OutFormat = CSVOutFormat(sep=","),
        df_1_pipe: DataProcessingPipe = IdentityPipe(),
        df_2_pipe: DataProcessingPipe = IdentityPipe(),
        out_pipe: DataProcessingPipe = IdentityPipe(),
        backend: ValidBackend = "polars",
    ):
        """
        Join a result dataframe to a reference dataframe.
        """
        extension, read_spec = get_extension_and_read_spec_from_format(
            out_format=out_format
        )
        source_meta = result_df_task.meta
        meta: Meta
        if isinstance(source_meta, ResultTableMeta):
            meta = ResultTableMeta(
                asset_id=AssetId(asset_id),
                trait=source_meta.trait,
                project=source_meta.project,
                extension=extension,
                read_spec=read_spec,
                # read_spec=DataFrameReadSpec(DataFrameTextFormat(separator=",")),
            )
        elif isinstance(source_meta, FilteredGWASDataMeta):
            meta = FilteredGWASDataMeta(
                id=AssetId(asset_id),
                trait=source_meta.trait,
                project=source_meta.project,
                extension=extension,
                read_spec=read_spec,
                sub_dir=source_meta.sub_dir,
            )
        elif isinstance(source_meta, GWASLabLeadVariantsMeta):
            meta = ResultTableMeta(
                asset_id=AssetId(asset_id),
                trait=source_meta.trait,
                project=source_meta.project,
                extension=extension,
                read_spec=read_spec,
                sub_dir=source_meta.sub_dir,
            )
        else:
            raise ValueError(f"Source meta has unknown type: {source_meta}")
        return cls(
            df_1_task=result_df_task,
            df_2_task=reference_df_task,
            how=how,
            left_on=left_on,
            right_on=right_on,
            meta=meta,
            df_1_pipe=df_1_pipe,
            df_2_pipe=df_2_pipe,
            out_format=out_format,
            backend=backend,
            out_pipe=out_pipe,
        )
