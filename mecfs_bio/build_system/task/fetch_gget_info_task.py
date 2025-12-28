"""
Task to use gget to annotate gene lists with annotations from genetics databases.
"""

import math
from pathlib import Path

import gget
import narwhals
import numpy as np
import pandas as pd
import structlog
from attrs import frozen

from mecfs_bio.build_system.asset.base_asset import Asset
from mecfs_bio.build_system.asset.file_asset import FileAsset
from mecfs_bio.build_system.meta.asset_id import AssetId
from mecfs_bio.build_system.meta.meta import Meta
from mecfs_bio.build_system.meta.read_spec.read_dataframe import scan_dataframe_asset
from mecfs_bio.build_system.meta.reference_meta.reference_file_meta import (
    ReferenceFileMeta,
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

logger = structlog.getLogger()

_dummy_gget_result = pd.DataFrame(
    {
        "ensembl_description": [None],
        "uniprot_description": [None],
        "ncbi_description": [None],
        "subcellular_localisation": [None],
    }
)

UNIPROT_ID_COL = "uniprot_id"
UNIPROT_DESCRIPTION = "uniprot_description"
PROTEIN_NAMES_COL = "protein_names"


@frozen
class FetchGGetInfoTask(Task):
    """
    Task to use gget (https://github.com/pachterlab/gget) to retrieve database information about a list of genes from a dataframe
    Useful for analyzing GWAS results.

    Listen to an interview with the primary developer of gget here:
    https://podcasts.apple.com/nz/podcast/99-laura-luebbert-gget-hunting-viruses-and/id1534473511?i=1000664104787

    Sometimes gget returns dataframes with inconsistent formating.
    e.g.: some columns are partly lists, and partly singleton values.
    Thus this file also contains functionality to munge the output of gget into a more consistent format.

    """

    source_df_task: Task
    ensembl_id_col: str
    _meta: Meta
    genes_to_use: int | None = None
    out_format: OutFormat = CSVOutFormat(",")
    post_pipe: DataProcessingPipe = IdentityPipe()

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
        result_df = (
            self.post_pipe.process(narwhals.from_native(result_df).lazy())
            .collect()
            .to_pandas()
        )
        result_df = _preprocess_columns(result_df)
        if isinstance(self.out_format, CSVOutFormat):
            result_df.to_csv(out_path, index=False, sep=self.out_format.sep)
        elif isinstance(self.out_format, ParquetOutFormat):
            result_df.to_parquet(
                out_path,
            )
        else:
            raise ValueError(f"Unsupported output format: {self.out_format}")
        return FileAsset(out_path)

    @classmethod
    def create(
        cls,
        asset_id: str,
        source_df_task: Task,
        ensembl_id_col: str,
        genes_to_use: int | None = None,
        post_pipe: DataProcessingPipe = IdentityPipe(),
        out_format: OutFormat = CSVOutFormat(","),
    ):
        source_meta = source_df_task.meta
        meta: Meta
        extension, read_spec = get_extension_and_read_spec_from_format(out_format)
        if isinstance(source_meta, ResultTableMeta):
            meta = ResultTableMeta(
                asset_id=asset_id,
                trait=source_meta.trait,
                project=source_meta.project,
                extension=extension,
                read_spec=read_spec,
            )
        elif isinstance(source_meta, ReferenceFileMeta):
            meta = ReferenceFileMeta(
                asset_id=AssetId(asset_id),
                group=source_meta.group,
                sub_group=source_meta.sub_group,
                extension=extension,
                read_spec=read_spec,
                sub_folder=source_meta.sub_folder,
            )
        else:
            raise ValueError("unknown source meta type")

        return cls(
            source_df_task=source_df_task,
            ensembl_id_col=ensembl_id_col,
            meta=meta,
            genes_to_use=genes_to_use,
            post_pipe=post_pipe,
            out_format=out_format,
        )


def _preprocess_columns(df: pd.DataFrame) -> pd.DataFrame:
    df = df.copy()
    if UNIPROT_ID_COL in df.columns:
        df[UNIPROT_ID_COL] = _wrap_col_in_list(df[UNIPROT_ID_COL])
    if PROTEIN_NAMES_COL in df.columns:
        df[PROTEIN_NAMES_COL] = _wrap_col_in_list(df[PROTEIN_NAMES_COL])
    if UNIPROT_DESCRIPTION in df.columns:
        df[UNIPROT_DESCRIPTION] = _wrap_col_in_list(df[UNIPROT_DESCRIPTION])
    for col in df.columns:
        df[col] = _array_to_list(df[col])
        df[col] = _unnest_list(df[col])
        df[col] = _clear_lists(df[col])

    return df


def _wrap_col_in_list(ser: pd.Series) -> pd.Series:
    return ser.apply(
        lambda x: [x] if ((not isinstance(x, list)) and (x is not None)) else x
    )


def _array_to_list(ser: pd.Series) -> pd.Series:
    return ser.apply(lambda x: list(x) if isinstance(x, np.ndarray) else x)


def _clear_lists(ser: pd.Series) -> pd.Series:
    def _clean(r: list) -> list:
        return [item for item in r if (item != "" and not _is_all_nan(item))]

    return ser.apply(lambda x: _clean(x) if isinstance(x, list) else x)


def _is_all_nan(l) -> bool:
    if not isinstance(l, list):
        return False
    return all(math.isnan(x) for x in l)


def _unnest_list(s: pd.Series) -> pd.Series:
    def _unest(x):
        if (
            isinstance(x, list)
            and len(x) == 1
            and (isinstance(x[0], list) or isinstance(x[0], np.ndarray))
        ):
            return list(x[0])
        return x

    return s.apply(_unest)
