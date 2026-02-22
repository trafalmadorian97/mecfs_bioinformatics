from pathlib import Path

import numpy as np
import pandas as pd
import pytest

from mecfs_bio.build_system.asset.base_asset import Asset
from mecfs_bio.build_system.asset.directory_asset import DirectoryAsset
from mecfs_bio.build_system.asset.file_asset import FileAsset
from mecfs_bio.build_system.meta.asset_id import AssetId
from mecfs_bio.build_system.meta.simple_directory_meta import SimpleDirectoryMeta
from mecfs_bio.build_system.meta.simple_file_meta import SimpleFileMeta
from mecfs_bio.build_system.task.fake_task import FakeTask
from mecfs_bio.build_system.task.magma.magma_forward_stepwise_select_task import (
    MagmaForwardStepwiseSelectTask,
    generate_mappers_from_wide_dataframe,
    generate_wide_dataframe,
)
from mecfs_bio.build_system.task.magma.magma_gene_set_analysis_task import (
    GENE_SET_ANALYSIS_OUTPUT_STEM_NAME,
)
from mecfs_bio.build_system.wf.base_wf import SimpleWF
from mecfs_bio.constants.magma_constants import (
    MAGMA_MODEL_COLUMN,
    MAGMA_P_COLUMN,
    MAGMA_VARIABLE_COLUMN,
)


def test_proportional_significance():
    """
    Test that a proportional significance value is consistent with manual calculation
    """
    df_marg = pd.read_csv(
        "test_mecfs_bio/unit/build_system/task/magma/dummy_marg_magma_output.txt",
        sep=r"\s+",
    )
    df_cond = pd.read_csv(
        "test_mecfs_bio/unit/build_system/task/magma/dummy_cond_magma_output.txt",
        sep=r"\s+",
    )
    df_wide = generate_wide_dataframe(df_cond, df_marg=df_marg)
    marg_dict, prop_sig_dict = generate_mappers_from_wide_dataframe(df_wide)
    assert prop_sig_dict[("Cluster242", "Cluster239")] == pytest.approx(
        np.log10(0.00048536) / np.log10(1.4544e-15)
    )


def test_forward_stepwise_when_missing_data(tmp_path: Path):
    empty_cond_df = pd.DataFrame(
        {MAGMA_VARIABLE_COLUMN: [], MAGMA_MODEL_COLUMN: [], MAGMA_P_COLUMN: []}
    )
    empty_marg_df = pd.DataFrame({MAGMA_VARIABLE_COLUMN: [], MAGMA_P_COLUMN: []})
    conf_dir_path = tmp_path / "cond_df.parquet"
    marg_dir_path = tmp_path / "marg_df.parquet"
    marg_dir_path.mkdir()
    conf_dir_path.mkdir()
    empty_marg_df.to_csv(
        marg_dir_path / (GENE_SET_ANALYSIS_OUTPUT_STEM_NAME + ".gsa.out"), sep="\t"
    )
    empty_cond_df.to_csv(
        conf_dir_path / (GENE_SET_ANALYSIS_OUTPUT_STEM_NAME + ".gsa.out"), sep="\t"
    )
    scratch = tmp_path / "scratch"
    scratch.mkdir()
    tsk = MagmaForwardStepwiseSelectTask(
        meta=SimpleFileMeta("stepwise"),
        magma_marginal_output_task=FakeTask(
            SimpleDirectoryMeta(
                "marginal",
            )
        ),
        magma_conditional_output_task=FakeTask(SimpleDirectoryMeta("conditional")),
    )

    def fetch(asset_id: AssetId) -> Asset:
        if asset_id == "marginal":
            return DirectoryAsset(marg_dir_path)
        if asset_id == "conditional":
            return DirectoryAsset(conf_dir_path)
        raise ValueError("Unknown asset id")

    result = tsk.execute(fetch=fetch, scratch_dir=scratch, wf=SimpleWF())
    assert isinstance(result, FileAsset)
    result_df = pd.read_csv(result.path)
    assert len(result_df) == 0
