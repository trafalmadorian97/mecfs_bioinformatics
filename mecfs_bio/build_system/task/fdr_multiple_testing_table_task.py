from pathlib import Path
from typing import Literal

import pandas as pd
import structlog
from attrs import frozen
from statsmodels.stats.multitest import multipletests

from mecfs_bio.build_system.asset.base_asset import Asset
from mecfs_bio.build_system.asset.directory_asset import DirectoryAsset
from mecfs_bio.build_system.asset.file_asset import FileAsset
from mecfs_bio.build_system.meta.asset_id import AssetId
from mecfs_bio.build_system.meta.base_meta import DirMeta
from mecfs_bio.build_system.meta.meta import Meta
from mecfs_bio.build_system.meta.processed_gwas_data_directory_meta import (
    ProcessedGwasDataDirectoryMeta,
)
from mecfs_bio.build_system.meta.read_spec.dataframe_read_spec import (
    DataFrameReadSpec,
    DataFrameTextFormat,
    DataFrameWhiteSpaceSepTextFormat,
)
from mecfs_bio.build_system.meta.read_spec.read_dataframe import (
    scan_dataframe,
    scan_dataframe_asset,
)
from mecfs_bio.build_system.meta.result_table_meta import ResultTableMeta
from mecfs_bio.build_system.rebuilder.fetch.base_fetch import Fetch
from mecfs_bio.build_system.task.base_task import Task
from mecfs_bio.build_system.task.magma.magma_gene_analysis_task import (
    GENE_ANALYSIS_OUTPUT_STEM_NAME,
    MagmaGeneAnalysisTask,
)
from mecfs_bio.build_system.wf.base_wf import WF

logger = structlog.get_logger()

Method = Literal["fdr_bh", "fdr_by", "bonferroni"]
REJECT_NULL_LABEL = "Reject Null"


@frozen
class MultipleTestingTableTask(Task):
    """
    Task to read a table of statistical test results with p values,
    filter it according to some multiple-test-correction procedure, and then write the filtered table.

    table_in_dir_filename/ table_in_dir_readspec:
    - set these when the source task provides a directory, and we need to find a results table within that directory

    For information on multiple testing procedures, see: https://www.statsmodels.org/dev/generated/statsmodels.stats.multitest.multipletests.html
    """

    _meta: Meta
    _table_source_task: Task
    p_value_column: str
    alpha: float  # uncorrected error rate
    method: Method
    apply_filter: bool = True
    table_in_dir_filename: str | None = None
    table_in_dir_readspec: DataFrameReadSpec | None = None

    def __attrs_post_init__(self):
        if (
            self.table_in_dir_filename is not None
            or self.table_in_dir_readspec is not None
        ):
            assert isinstance(self._table_source_meta, DirMeta)

    @property
    def _table_source_id(self) -> AssetId:
        return self._table_source_task.asset_id

    @property
    def _table_source_meta(self) -> Meta:
        return self._table_source_task.meta

    @property
    def meta(self) -> Meta:
        return self._meta

    @property
    def deps(self) -> list["Task"]:
        return [self._table_source_task]

    def _get_frame_from_source_asset(self, asset: Asset) -> pd.DataFrame:
        if isinstance(asset, FileAsset):
            return (
                scan_dataframe_asset(asset, meta=self._table_source_meta)
                .collect()
                .to_pandas()
            )
        if isinstance(asset, DirectoryAsset):
            assert self.table_in_dir_readspec is not None
            assert self.table_in_dir_filename is not None
            return (
                scan_dataframe(
                    asset.path / self.table_in_dir_filename,
                    spec=self.table_in_dir_readspec,
                )
                .collect()
                .to_pandas()
            )
        raise ValueError("unknown asset type")

    def execute(self, scratch_dir: Path, fetch: Fetch, wf: WF) -> Asset:
        asset = fetch(self._table_source_id)
        frame = self._get_frame_from_source_asset(asset)
        p_vals = frame[self.p_value_column]
        logger.debug(
            f"\nApplying multiple testing correction method {self.method} to dataframe with alpha {self.alpha}"
        )
        logger.debug(f"Initial dataframe contains {len(frame)} rows")
        reject, corrected, _, _ = multipletests(
            pvals=p_vals, alpha=self.alpha, method=self.method
        )
        # logger.debug(f"Corrected p values {corrected}")
        frame[CORRECTED_COL_LABEL] = corrected
        if self.apply_filter:
            frame = frame.loc[reject]
            assert (frame[CORRECTED_COL_LABEL] <= self.alpha).all()
        else:
            frame[REJECT_NULL_LABEL] = reject
        frame: pd.DataFrame = frame.sort_values(by=[self.p_value_column])
        out_path = scratch_dir / "filtered_frame.csv"
        frame.to_csv(out_path, index=False)
        return FileAsset(out_path)

    @classmethod
    def create_from_magma_gene_analysis_task(
        cls,
        alpha: float,
        method: Method,
        asset_id: str,
        source_task: MagmaGeneAnalysisTask,
    ):
        source_meta = source_task.meta
        assert isinstance(source_meta, ProcessedGwasDataDirectoryMeta)
        meta = ResultTableMeta(
            asset_id=asset_id,
            trait=source_meta.trait,
            project=source_meta.project,
            extension=".csv",
            read_spec=DataFrameReadSpec(DataFrameTextFormat(separator=",")),
        )
        return cls(
            meta=meta,
            table_source_task=source_task,
            p_value_column="P",
            alpha=alpha,
            method=method,
            table_in_dir_filename=GENE_ANALYSIS_OUTPUT_STEM_NAME + ".genes.out",
            table_in_dir_readspec=DataFrameReadSpec(
                DataFrameWhiteSpaceSepTextFormat(comment_code="#")
            ),
        )

    @classmethod
    def create_from_result_table_task(
        cls,
        alpha: float,
        method: Method,
        asset_id: str,
        p_value_column: str,
        source_task: Task,
        apply_filter: bool = True,
    ):
        source_meta = source_task.meta
        assert isinstance(source_meta, ResultTableMeta)
        meta = ResultTableMeta(
            asset_id=asset_id,
            trait=source_meta.trait,
            project=source_meta.project,
            extension=".csv",
            read_spec=DataFrameReadSpec(DataFrameTextFormat(separator=",")),
        )
        return cls(
            meta=meta,
            table_source_task=source_task,
            p_value_column=p_value_column,
            alpha=alpha,
            method=method,
            apply_filter=apply_filter,
        )


CORRECTED_COL_LABEL = "_Corrected P Value_"
