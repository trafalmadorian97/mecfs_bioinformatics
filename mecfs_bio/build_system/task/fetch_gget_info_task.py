from pathlib import Path

import gget
import pandas as pd
import structlog
from attrs import frozen

from mecfs_bio.build_system.asset.base_asset import Asset
from mecfs_bio.build_system.asset.file_asset import FileAsset
from mecfs_bio.build_system.meta.asset_id import AssetId
from mecfs_bio.build_system.meta.meta import Meta
from mecfs_bio.build_system.meta.read_spec.dataframe_read_spec import (
    DataFrameReadSpec,
    DataFrameTextFormat,
)
from mecfs_bio.build_system.meta.read_spec.read_dataframe import scan_dataframe_asset
from mecfs_bio.build_system.meta.reference_meta.reference_file_meta import (
    ReferenceFileMeta,
)
from mecfs_bio.build_system.meta.result_table_meta import ResultTableMeta
from mecfs_bio.build_system.rebuilder.fetch.base_fetch import Fetch
from mecfs_bio.build_system.task.base_task import Task
from mecfs_bio.build_system.wf.base_wf import WF

logger = structlog.getLogger()

_dummy_gget_result = pd.DataFrame(
    {
        "ensembl_description": [None],
        "uniprot_description": [None],
        "ncbi_description": [None],
        "subcellular_localisation": [None],
    }
)


@frozen
class FetchGGetInfoTask(Task):
    source_df_task: Task
    ensembl_id_col: str
    _meta: Meta
    genes_to_use: int | None = None

    @property
    def source_meta(self) -> Meta:
        return self.source_df_task.meta

    @property
    def source_id(self) -> AssetId:
        return self.source_df_task.asset_id

    @property
    def meta(self) -> Meta:
        return self._meta

    @property
    def deps(self) -> list["Task"]:
        return [self.source_df_task]

    def execute(self, scratch_dir: Path, fetch: Fetch, wf: WF) -> Asset:
        source_asset = fetch(self.source_id)
        df = (
            scan_dataframe_asset(source_asset, meta=self.source_meta)
            .collect()
            .to_pandas()
        )
        genes = df[self.ensembl_id_col].tolist()
        if self.genes_to_use is not None:
            genes = genes[: self.genes_to_use]
        logger.debug(f"Using gget to retrieve info on {len(genes)} genes")
        logger.debug(f"Genes are: {genes}")
        gget_result = gget.info(genes)
        if isinstance(gget_result, pd.DataFrame):
            gene_info = gget_result
        else:
            gene_info = _dummy_gget_result
        result_df = pd.merge(
            df, gene_info, left_on=self.ensembl_id_col, right_index=True, how="left"
        )
        out_path = scratch_dir / f"{self.source_id}.csv"
        result_df.to_csv(out_path, index=False)
        return FileAsset(out_path)

    @classmethod
    def create(
        cls,
        asset_id: str,
        source_df_task: Task,
        ensembl_id_col: str,
        genes_to_use: int | None = None,
    ):
        source_meta = source_df_task.meta
        meta: Meta
        if isinstance(source_meta, ResultTableMeta):
            meta = ResultTableMeta(
                asset_id=asset_id,
                trait=source_meta.trait,
                project=source_meta.project,
                extension=".csv",
                read_spec=DataFrameReadSpec(DataFrameTextFormat(separator=",")),
            )
        elif isinstance(source_meta, ReferenceFileMeta):
            meta = ReferenceFileMeta(
                asset_id=AssetId(asset_id),
                group=source_meta.group,
                sub_group=source_meta.sub_group,
                extension=".csv",
                sub_folder=source_meta.sub_folder,
            )
        else:
            raise ValueError("unknown source meta type")

        return cls(
            source_df_task=source_df_task,
            ensembl_id_col=ensembl_id_col,
            meta=meta,
            genes_to_use=genes_to_use,
        )
