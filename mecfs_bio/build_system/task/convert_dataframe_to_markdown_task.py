from pathlib import Path

import numpy as np
import pandas as pd
from attrs import frozen

from mecfs_bio.build_system.asset.base_asset import Asset
from mecfs_bio.build_system.asset.file_asset import FileAsset
from mecfs_bio.build_system.meta.asset_id import AssetId
from mecfs_bio.build_system.meta.filtered_gwas_data_meta import FilteredGWASDataMeta
from mecfs_bio.build_system.meta.meta import Meta
from mecfs_bio.build_system.meta.read_spec.read_dataframe import scan_dataframe_asset
from mecfs_bio.build_system.meta.result_table_meta import ResultTableMeta
from mecfs_bio.build_system.rebuilder.fetch.base_fetch import Fetch
from mecfs_bio.build_system.task.base_task import Task
from mecfs_bio.build_system.task.pipes.data_processing_pipe import DataProcessingPipe
from mecfs_bio.build_system.task.pipes.identity_pipe import IdentityPipe
from mecfs_bio.build_system.wf.base_wf import WF


@frozen
class ConvertDataFrameToMarkdownTask(Task):
    """
    Task to convert a dataframe to markdown format.
    Useful for writing up results
    """

    _meta: Meta
    _df_task: Task
    pipe: DataProcessingPipe = IdentityPipe()

    def __attrs_post_init__(self):
        assert self._source_meta.read_spec() is not None

    @property
    def _source_meta(self) -> Meta:
        return self._df_task.meta

    @property
    def _source_id(self) -> AssetId:
        return self._source_meta.asset_id

    @property
    def meta(self) -> Meta:
        return self._meta

    @property
    def deps(self) -> list["Task"]:
        return [self._df_task]

    def execute(self, scratch_dir: Path, fetch: Fetch, wf: WF) -> Asset:
        source_asset = fetch(self._source_id)
        nw_df = self.pipe.process(
            scan_dataframe_asset(asset=source_asset, meta=self._source_meta)
        )

        df = nw_df.collect().to_pandas()
        df = _array_to_list(df)
        out_path = scratch_dir / "output.md"
        markdown_str = df.to_markdown(index=False)
        with open(out_path, "w") as f:
            f.write(markdown_str)
        return FileAsset(out_path)

    @classmethod
    def create_from_result_table_task(
        cls, source_task: Task, asset_id: str, pipe: DataProcessingPipe = IdentityPipe()
    ):
        source_meta = source_task.meta
        assert isinstance(source_meta, ResultTableMeta)
        meta = FilteredGWASDataMeta(
            short_id=AssetId(asset_id),
            trait=source_meta.trait,
            project=source_meta.project,
            sub_dir=source_meta.sub_dir,
            extension=".md",
        )
        return cls(meta=meta, df_task=source_task, pipe=pipe)


def _array_to_list(df: pd.DataFrame) -> pd.DataFrame:
    df = df.copy()
    for col in df.columns:
        df[col] = df[col].apply(lambda x: list(x) if isinstance(x, np.ndarray) else x)
    return df
