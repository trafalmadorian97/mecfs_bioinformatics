"""
Task to combine gene lists from multiple sources
"""

from pathlib import Path
from typing import Sequence

import narwhals
import pandas as pd
from attrs import frozen

from mecfs_bio.build_system.asset.base_asset import Asset
from mecfs_bio.build_system.asset.file_asset import FileAsset
from mecfs_bio.build_system.meta.asset_id import AssetId
from mecfs_bio.build_system.meta.meta import Meta
from mecfs_bio.build_system.meta.read_spec.read_dataframe import scan_dataframe_asset
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

ENSEMBL_ID_LABEL = "Ensembl ID"


@frozen
class SrcGeneList:
    task: Task
    name: str
    ensemble_id_column: str
    pipe: DataProcessingPipe = IdentityPipe()


@frozen
class CombineGeneListsTask(Task):
    """
    Task to aggregate gene lists from multiple sources

    Example use case: I have one gene list from MAGMA and another from Gwaslab, and I want to combine them
    to create master gene list for the trait of interest.
    """

    _meta: Meta
    src_gene_lists: Sequence[SrcGeneList]
    out_format: OutFormat = CSVOutFormat(sep=",")
    out_pipe: DataProcessingPipe = IdentityPipe()

    def __attrs_post_init__(self):
        assert len(self.src_gene_lists) > 0
        names = set([gl.name for gl in self.src_gene_lists])
        assert len(names) == len(self.src_gene_lists)

    @property
    def meta(self) -> Meta:
        return self._meta

    @property
    def deps(self) -> list["Task"]:
        return [gl.task for gl in self.src_gene_lists]

    def execute(self, scratch_dir: Path, fetch: Fetch, wf: WF) -> Asset:
        gene_dict: dict[str, list[str]] = {}
        for gl in self.src_gene_lists:
            asset = fetch(gl.task.asset_id)
            df = (
                scan_dataframe_asset(
                    asset,
                    gl.task.meta,
                )
                .collect()
                .to_pandas()
            )
            gene_dict = _combine_gene_dicts(
                gene_dict,
                _get_gene_dict_from_df(
                    df, name_col=gl.ensemble_id_column, method_name=gl.name
                ),
            )
        result_df = pd.DataFrame(
            {ENSEMBL_ID_LABEL: gene_dict.keys(), "sources": gene_dict.values()}
        )
        result_df = (
            self.out_pipe.process(narwhals.from_native(result_df).lazy())
            .collect()
            .to_pandas()
        )
        out_path = scratch_dir / (self._meta.asset_id + ".csv")
        if isinstance(self.out_format, CSVOutFormat):
            result_df.to_csv(out_path, index=False, sep=self.out_format.sep)
        elif isinstance(self.out_format, ParquetOutFormat):
            result_df.to_parquet(out_path)
        else:
            raise ValueError(f"Unsupported output format: {self.out_format}")
        return FileAsset(out_path)

    @classmethod
    def create(
        cls,
        asset_id: str,
        src_gene_lists: Sequence[SrcGeneList],
        out_format: OutFormat,
        out_pipe: DataProcessingPipe = IdentityPipe(),
    ):
        assert len(src_gene_lists) > 0
        first_gene_list = src_gene_lists[0]
        src_meta = first_gene_list.task.meta
        extension, read_spec = get_extension_and_read_spec_from_format(out_format)
        assert isinstance(src_meta, ResultTableMeta)
        meta = ResultTableMeta(
            asset_id=AssetId(asset_id),
            trait=src_meta.trait,
            project=src_meta.project,
            extension=extension,
            read_spec=read_spec,
        )
        return cls(
            src_gene_lists=src_gene_lists,
            meta=meta,
            out_format=out_format,
            out_pipe=out_pipe,
        )


def _get_gene_dict_from_df(
    df: pd.DataFrame, name_col: str, method_name: str
) -> dict[str, list[str]]:
    return {item: [method_name] for item in df[name_col]}


def _combine_gene_dicts(gd1: dict[str, list[str]], gd2: dict[str, list[str]]):
    gd1 = dict(gd1)
    for name, src_list in gd2.items():
        if name not in gd1:
            gd1[name] = src_list
        else:
            gd1[name].extend(src_list)
    return gd1
