from pathlib import Path
from typing import Sequence

import pandas as pd
from attrs import frozen

from mecfs_bio.build_system.asset.base_asset import Asset
from mecfs_bio.build_system.asset.file_asset import FileAsset
from mecfs_bio.build_system.meta.meta import Meta
from mecfs_bio.build_system.meta.read_spec.read_dataframe import scan_dataframe_asset
from mecfs_bio.build_system.rebuilder.fetch.base_fetch import Fetch
from mecfs_bio.build_system.task.base_task import Task
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
    """

    _meta: Meta
    src_gene_lists: Sequence[SrcGeneList]

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
        out_path = scratch_dir / (self._meta.asset_id + ".csv")
        result_df.to_csv(out_path, index=False)
        return FileAsset(out_path)


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
