from pathlib import Path
from typing import Any, Sequence

import numpy as np
import pandas as pd
import structlog
from attrs import frozen

from mecfs_bio.build_system.asset.base_asset import Asset
from mecfs_bio.build_system.asset.directory_asset import DirectoryAsset
from mecfs_bio.build_system.asset.file_asset import FileAsset
from mecfs_bio.build_system.meta.asset_id import AssetId
from mecfs_bio.build_system.meta.meta import Meta
from mecfs_bio.build_system.meta.processed_gwas_data_directory_meta import (
    ProcessedGwasDataDirectoryMeta,
)
from mecfs_bio.build_system.meta.read_spec.dataframe_read_spec import (
    DataFrameReadSpec,
    DataFrameTextFormat,
)
from mecfs_bio.build_system.meta.result_table_meta import ResultTableMeta
from mecfs_bio.build_system.rebuilder.fetch.base_fetch import Fetch
from mecfs_bio.build_system.task.base_task import Task
from mecfs_bio.build_system.task.magma.magma_gene_set_analysis_task import (
    GENE_SET_ANALYSIS_OUTPUT_STEM_NAME,
)
from mecfs_bio.build_system.wf.base_wf import WF
from mecfs_bio.constants.magma_constants import MAGMA_VARIABLE_COLUMN, MAGMA_MODEL_COLUMN, MAGMA_P_COLUMN

VARCODE_COLUMN = "VARCODE"
P_MARG_COLUMN = "P_MARG"
P_COND_COLUMN = "P_COND"

RETAINED_CLUSTERS_COLUMN = "Retained_clusters"

logger = structlog.get_logger()


@frozen
class MagmaForwardStepwiseSelectTask(Task):
    """
    - Since the gene expression patterns of many cell types are highly correlated, it can be difficult
     to distinguish the independent significant signals on a typical plot of MAGMA results.
    - To help resolve this, this function attempts to identify a set of relatively independent significant cell types, in
    the sense that the p value of one does  not decline too much when conditioning on another.
    This is loosely based on: https://github.com/Integrative-Mental-Health-Lab/linking_cell_types_to_brain_phenotypes/blob/675b5c9b58b8762934183a3ca61ae49ad587934a/MAGMA/5.forward_selection_condition_results.md
    """

    _meta: Meta
    magma_marginal_output_task: Task
    magma_conditional_output_task: Task
    min_prop_sig: float = 0.5
    significance_threshold: float = 0.05

    @property
    def marginal_asset_id(self) -> AssetId:
        return self.magma_marginal_output_task.asset_id

    @property
    def conditional_asset_id(self) -> AssetId:
        return self.magma_conditional_output_task.asset_id

    @property
    def meta(self) -> Meta:
        return self._meta

    @property
    def deps(self) -> list["Task"]:
        return [self.magma_marginal_output_task, self.magma_conditional_output_task]

    def execute(self, scratch_dir: Path, fetch: Fetch, wf: WF) -> Asset:
        marginal_asset = fetch(self.marginal_asset_id)
        assert isinstance(marginal_asset, DirectoryAsset)
        conditional_asset = fetch(self.conditional_asset_id)
        assert isinstance(conditional_asset, DirectoryAsset)
        marginal_path = marginal_asset.path / (
            GENE_SET_ANALYSIS_OUTPUT_STEM_NAME + ".gsa.out"
        )
        conditional_path = conditional_asset.path / (
            GENE_SET_ANALYSIS_OUTPUT_STEM_NAME + ".gsa.out"
        )
        df_marg = pd.read_csv(marginal_path, comment="#", sep=r"\s+")
        df_cond = pd.read_csv(conditional_path, comment="#", sep=r"\s+")
        df_wide = generate_wide_dataframe(df_cond=df_cond, df_marg=df_marg)
        p_marg_dict, prop_sig_dict = generate_mappers_from_wide_dataframe(
            df_wide=df_wide
        )

        all_cluster_list = df_marg.sort_values(by=MAGMA_P_COLUMN)[
            MAGMA_VARIABLE_COLUMN
        ].tolist()
        sig_cluster_set = set(
            df_marg.loc[
                df_marg[MAGMA_P_COLUMN] <= (self.significance_threshold / len(df_marg))
                ][MAGMA_VARIABLE_COLUMN].tolist()
        )
        sig_cluster_list = [
            item for item in all_cluster_list if item in sig_cluster_set
        ]

        retained_cluster_list = get_retained_clusters(
            all_cluster_list=sig_cluster_list,
            prop_sig_dict=prop_sig_dict,
            min_prop_sig=self.min_prop_sig,
        )
        logger.debug(f"Number of retained clusters: {len(retained_cluster_list)}")
        retained_cluster_df = pd.DataFrame(
            {
                RETAINED_CLUSTERS_COLUMN: retained_cluster_list,
            }
        )
        out_path = scratch_dir / self.asset_id
        retained_cluster_df.to_csv(out_path, index=False)
        return FileAsset(out_path)

    @classmethod
    def create(
        cls,
        asset_id: str,
        magma_marginal_output_task: Task,
        magma_conditional_output_task: Task,
        min_prop_sig: float = 0.5,
        significance_threshold=0.05,
    ):
        source_meta = magma_marginal_output_task.meta
        assert isinstance(source_meta, ProcessedGwasDataDirectoryMeta)
        meta = ResultTableMeta(
            asset_id=AssetId(asset_id),
            trait=source_meta.trait,
            project=source_meta.project,
            extension=".csv",
            read_spec=DataFrameReadSpec(DataFrameTextFormat(separator=",")),
        )
        return cls(
            meta=meta,
            magma_marginal_output_task=magma_marginal_output_task,
            magma_conditional_output_task=magma_conditional_output_task,
            min_prop_sig=min_prop_sig,
            significance_threshold=significance_threshold,
        )


def get_retained_clusters(
    all_cluster_list: Sequence[str],
    prop_sig_dict: dict[tuple[str, str], float],
    min_prop_sig: float,
) -> list[str]:
    to_retain = []
    to_retain.append(all_cluster_list[0])
    for cluster in all_cluster_list[1:]:
        if pairwise_independent(
            cluster1=cluster,
            clusters_to_retain=to_retain,
            min_prop_sig=min_prop_sig,
            prop_sig_dict=prop_sig_dict,
        ):
            to_retain.append(cluster)
    return to_retain


def pairwise_independent(
    cluster1: str,
    clusters_to_retain: list[str],
    prop_sig_dict: dict[tuple[str, str], float],
    min_prop_sig: float,
) -> bool:
    for retained_cluster in clusters_to_retain:
        if not independent(
            cluster1,
            retained_cluster,
            prop_sig_dict=prop_sig_dict,
            min_prop_sig=min_prop_sig,
        ):
            return False
    return True


def independent(
    cluster1: str,
    cluster2: str,
    prop_sig_dict: dict[tuple[str, str], float],
    min_prop_sig: float,
) -> bool:
    if (
        prop_sig_dict[(cluster1, cluster2)] >= min_prop_sig
        and prop_sig_dict[(cluster2, cluster1)] >= min_prop_sig
    ):
        return True
    return False


def generate_wide_dataframe(
    df_cond: pd.DataFrame, df_marg: pd.DataFrame
) -> pd.DataFrame:
    """
    Based on https://github.com/Integrative-Mental-Health-Lab/linking_cell_types_to_brain_phenotypes/blob/675b5c9b58b8762934183a3ca61ae49ad587934a/MAGMA/5.forward_selection_condition_results.md
    """
    df_cond = (
        df_cond.copy()
        .loc[:, [MAGMA_VARIABLE_COLUMN, MAGMA_MODEL_COLUMN, MAGMA_P_COLUMN]]
        .rename(columns={MAGMA_P_COLUMN: P_COND_COLUMN})
    )
    df_marg = (
        df_marg.copy()
        .loc[:, [MAGMA_VARIABLE_COLUMN, MAGMA_P_COLUMN]]
        .rename(columns={MAGMA_P_COLUMN: P_MARG_COLUMN})
        .sort_values(by=P_MARG_COLUMN)
    )
    df_comb = df_cond.merge(df_marg, on=MAGMA_VARIABLE_COLUMN, how="left")
    df_comb[VARCODE_COLUMN] = (
        df_comb.groupby(MAGMA_MODEL_COLUMN)[MAGMA_VARIABLE_COLUMN]
        .cumcount()
        .map({0: "a", 1: "b"})
    )
    df_wide: Any = df_comb.pivot(
        index=MAGMA_MODEL_COLUMN,
        columns=VARCODE_COLUMN,
        values=[MAGMA_VARIABLE_COLUMN, P_COND_COLUMN, P_MARG_COLUMN],
    )
    df_wide.columns = [f"{v}_{c}" for v, c in df_wide.columns]
    for col in ["P_COND_a", "P_COND_b", "P_MARG_a", "P_MARG_b"]:
        df_wide[col] = df_wide[col].astype(np.float64)
    df_wide["PS_a"] = np.log10(df_wide["P_COND_a"]) / np.log10(df_wide["P_MARG_a"])
    df_wide["PS_b"] = np.log10(df_wide["P_COND_b"]) / np.log10(df_wide["P_MARG_b"])

    assert (df_wide["P_MARG_a"] <= df_wide["P_MARG_b"]).all()
    return df_wide


def generate_mappers_from_wide_dataframe(
    df_wide: pd.DataFrame,
) -> tuple[dict[str, float], dict[tuple[str, str], float]]:
    p_marg_dict = {}
    prop_sig_dict = {}
    for i in range(len(df_wide)):
        p_marg_dict[df_wide["VARIABLE_a"].iloc[i]] = df_wide["P_MARG_a"].iloc[i]
        p_marg_dict[df_wide["VARIABLE_b"].iloc[i]] = df_wide["P_MARG_b"].iloc[i]
        prop_sig_dict[
            (df_wide["VARIABLE_a"].iloc[i], df_wide["VARIABLE_b"].iloc[i])
        ] = df_wide["PS_a"].iloc[i]
        prop_sig_dict[
            (df_wide["VARIABLE_b"].iloc[i], df_wide["VARIABLE_a"].iloc[i])
        ] = df_wide["PS_b"].iloc[i]
    return p_marg_dict, prop_sig_dict
