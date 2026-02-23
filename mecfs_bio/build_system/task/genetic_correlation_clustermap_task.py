import xarray as xr
from pathlib import Path

from attrs import frozen
from fontTools.misc.transform import Identity

from mecfs_bio.build_system.asset.base_asset import Asset
from mecfs_bio.build_system.meta.meta import Meta
from mecfs_bio.build_system.meta.read_spec.read_dataframe import scan_dataframe_asset
from mecfs_bio.build_system.rebuilder.fetch.base_fetch import Fetch
from mecfs_bio.build_system.task.base_task import Task
from mecfs_bio.build_system.task.pipes.data_processing_pipe import DataProcessingPipe
from mecfs_bio.build_system.task.pipes.identity_pipe import IdentityPipe
from mecfs_bio.build_system.wf.base_wf import WF

XR_GENETIC_CORR_ARRAY ="genetic_correlation"
XR_GENETIC_CORR_P_VALUE_ARRAY= "genetic_correlation_p_value"
XR_TRAIT_1_DIM = "trait_1"
XR_TRAIT_2_DIM = "trait_2"

@frozen
class GeneticCorrSource:
    task: Task
    trait_1_col: str
    trait_2_col: str
    rg_col: str
    p_col: str
    df_pipe: DataProcessingPipe=IdentityPipe()


def load_xr_corr_dataset(
    src: GeneticCorrSource,
    fetch: Fetch,
)-> xr.Dataset:
    asset = fetch(src.task.asset_id)
    df = (
        src.df_pipe.process(
            scan_dataframe_asset(
                asset,
                meta=src.task.meta,
            )
        )
        .collect()
        .to_pandas()
    )
    pivoted_rg = df.pivot(
        index=src.trait_1_col,
        columns=src.trait_2_col,
        values=src.rg_col,
    )
    pivoted_p = df



@frozen
class GeneticCorrelationClustermapTask(Task):
    _meta: Meta
    genetic_corr_table_task: Task
    @property
    def meta(self) -> Meta:
        return self._meta

    @property
    def deps(self) -> list["Task"]:
        return [self.genetic_corr_table_task]

    def execute(self, scratch_dir: Path, fetch: Fetch, wf: WF) -> Asset:
        pass