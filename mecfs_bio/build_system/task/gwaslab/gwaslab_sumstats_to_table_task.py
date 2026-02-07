from pathlib import Path

import structlog

logger = structlog.get_logger()

import narwhals
import pandas as pd
from attrs import frozen

from mecfs_bio.build_system.asset.base_asset import Asset
from mecfs_bio.build_system.asset.file_asset import FileAsset
from mecfs_bio.build_system.meta.asset_id import AssetId
from mecfs_bio.build_system.meta.filtered_gwas_data_meta import FilteredGWASDataMeta
from mecfs_bio.build_system.meta.gwaslab_meta.gwaslab_sumstats_meta import (
    GWASLabSumStatsMeta,
)
from mecfs_bio.build_system.meta.meta import Meta
from mecfs_bio.build_system.meta.read_spec.dataframe_read_spec import (
    DataFrameParquetFormat,
    DataFrameReadSpec,
)
from mecfs_bio.build_system.meta.read_spec.read_sumstats import read_sumstats
from mecfs_bio.build_system.rebuilder.fetch.base_fetch import Fetch
from mecfs_bio.build_system.task.base_task import Task
from mecfs_bio.build_system.task.pipes.data_processing_pipe import DataProcessingPipe
from mecfs_bio.build_system.task.pipes.identity_pipe import IdentityPipe
from mecfs_bio.build_system.wf.base_wf import WF


@frozen
class GwasLabSumstatsToTableTask(Task):
    """
    Task to write a sumstats object to a plain table for further processing.
    """

    _meta: Meta
    source_sumstats_task: Task
    pipe: DataProcessingPipe = IdentityPipe()

    @property
    def meta(self) -> Meta:
        return self._meta

    @property
    def source_meta(self) -> Meta:
        return self.source_sumstats_task.meta

    @property
    def source_asset_id(self) -> AssetId:
        return self.source_meta.asset_id

    @property
    def deps(self) -> list["Task"]:
        return [self.source_sumstats_task]

    def execute(self, scratch_dir: Path, fetch: Fetch, wf: WF) -> Asset:
        asset = fetch(self.source_asset_id)
        sumstats = read_sumstats(asset)
        df: pd.DataFrame = sumstats.data
        logger.debug(f"Post sumstats df has shape {df.shape}")
        n_df = narwhals.from_native(df).lazy()
        df = self.pipe.process(n_df).collect().to_pandas()
        logger.debug(f"Post pipe df has shape {df.shape}")
        out_loc = scratch_dir / "data.parquet"
        df.to_parquet(out_loc, index=False)
        return FileAsset(path=out_loc)

    @classmethod
    def create_from_source_task(
        cls,
        source_tsk: Task,
        asset_id: str,
        sub_dir: str,
        pipe: DataProcessingPipe = IdentityPipe(),
    ):
        source_meta = source_tsk.meta
        assert isinstance(source_meta, GWASLabSumStatsMeta)
        meta = FilteredGWASDataMeta(
            id=AssetId(asset_id),
            trait=source_meta.trait,
            project=source_meta.project,
            sub_dir=sub_dir,
            read_spec=DataFrameReadSpec(format=DataFrameParquetFormat()),
        )
        return cls(meta=meta, source_sumstats_task=source_tsk, pipe=pipe)
